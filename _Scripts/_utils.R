library(dplyr)
library(eyelinker)


# NOTE: The custom 'sync_timestamps' function is only necessary since the
# sample rate for this example file is < 1000 Hz. Basically, since at 250 Hz the
# tracker checks for messages every millisecond but gaze samples every 4 ms, you
# often get messages that don't line up exactly with the sample timestamps. The
# function below fixes this by rounding each message timestamp to the nearest
# sample in the gaze data.
#
# If your data was recorded at 1000 Hz or 2000 Hz you can skip this bit (you can
# check by checking dat$info$sample.rate).

sync_timestamps <- function(df, raw) {

  # Get tracker sample rate and block onset timestamps from raw pupil data
  ms_per_sample <- median(raw$time[2:11] - raw$time[1:10])
  block_onsets <- raw %>%
    group_by(block) %>%
    summarize(onset = time[1])

  # Adjust df so that its timestamps map onto the nearest sample timestamp
  df <- df %>%
    left_join(block_onsets, by = "block") %>%
    group_by(block) %>%
    mutate(
      time = round((time - onset) / ms_per_sample) * ms_per_sample + onset
    ) %>%
    select(-c(onset))

  df
}


pupil_lowpass <- function(x, high_freq, srate) {

  # Create lowpass butterworth filter for pupil data
  lowpass_w <- high_freq / (srate / 2)
  bfilt_lo <- gsignal::butter(3, w = lowpass_w, type = "low", output = "Sos")

  # Apply the filter to the pupil samples
  gsignal::filtfilt(bfilt_lo, x)
}


interpolate_blinks <- function(time, ps, in_blink, smooth = FALSE) {
  if (!any(in_blink)) {
    return(ps)
  }
  if (smooth) {
    ps_fun <- splinefun(time[!in_blink], ps[!in_blink], method = "natural")
  } else {
    ps_fun <- approxfun(time[!in_blink], ps[!in_blink])
  }
  ifelse(in_blink, ps_fun(time), ps)
}


which_phase <- function(time, phases) {
  phase_idx <- which_interval(time, cbind(phases$start, phases$end))
  phase <- ifelse(is.na(phase_idx), NA, as.character(phases$phase[phase_idx]))
  phase
}


in_blink <- function(time, blinks, lpad = 50, rpad = 60) {
  time %within% cbind(blinks$stime - lpad, blinks$etime + rpad)
}


get_slope <- function(ps_dat) {
  fit <- lm(ps ~ time, data = ps_dat)$coefficients
  tibble(intercept = fit[1], slope = fit[2])
}


get_epochs <- function(msg, start, end, pad = c(0, 0)) {
  epochs <- msg %>%
    group_by(block) %>%
    filter(grepl(start, text) | grepl(end, text)) %>%
    filter(n() >= 2) %>%
    summarize(
      start = time[grepl(start, text)] - pad[1],
      end = time[grepl(end, text)] + pad[2],
    )
  epochs
}


within_epoch <- function(time, epochs) {
  time %within% cbind(epochs$start, epochs$end)
}


get_pupil_data <- function(
  asc, start = "cue_on", end = "probe_on", zero = NA, pad = c(100, 100)
) {

  # First, get map of ASC 'blocks' to actual blocks/trials of task

  events <- asc$msg %>%
    filter(grepl(" b\\d ", text)) %>%
    separate(text, c("phase", "block_num", "trial_num"), sep = " ") %>%
    mutate(
      phase = as.factor(phase),
      block_num = as.numeric(gsub("b", "", block_num)),
      trial_num = as.numeric(gsub("t", "", trial_num))
    )

  asc_trial_map <- events %>%
    subset(phase == "trial_start") %>%
    select(c(block, block_num, trial_num))


  # Get the start/end times for the window of interest for each trial

  trial_windows <- get_epochs(asc$msg, start, end, pad)


  # Find error trials & trials with blinks within 100 ms of the windows of
  # interest, and flag them for exclusion

  blink_epochs <- asc$blinks %>%
    rename(start = stime, end = etime) %>%
    mutate(
      start = start - 100,
      end = end + 100
    )

  blink_during_cue <- events %>%
    subset(phase == "cue_on") %>%
    filter(within_epoch(time, blink_epochs))

  blink_during_probe <- events %>%
    subset(phase == "probe_on") %>%
    filter(within_epoch(time, blink_epochs))

  incomplete_trials <- events %>%
    group_by(block) %>%
    summarize(
      incomplete = !("probe_on" %in% phase)
    ) %>%
    filter(incomplete)

  immediate_cue_trials <- events %>%
    # Due to an occasional EyeLink init delay that klibs didn't anticipate
    group_by(block) %>%
    summarize(
      cue_onset = ifelse(
        any(phase == "cue_on"),
        time[phase == "cue_on"] - time[phase == "trial_start"], NA
      )
    ) %>%
    filter(cue_onset < 700)

  bad_trials <- bind_rows(
    select(incomplete_trials, c(block)),
    select(immediate_cue_trials, c(block)),
    select(blink_during_cue, c(block)),
    select(blink_during_probe, c(block))
  )$block


  # Append event markers to raw pupil data

  ps_dat <- asc$raw %>%
    select(c(block, time, ps)) %>%
    filter(!(block %in% bad_trials)) %>%
    left_join(select(events, c(time, phase)), by = "time") %>%
    group_by(block) %>%
    fill(phase, .direction = "down") %>%
    filter(phase %in% c("trial_start", "cue_on", "probe_on"))


  # Subset raw pupil data to main window of interest for each trial

  ps_dat <- ps_dat %>%
    left_join(trial_windows, by = "block") %>%
    group_by(block) %>%
    filter(time >= start & time < end) %>%
    select(-c(start, end))


  # Identify and repair blinks using linear interpolation, discarding any blinks
  # regions at start/end of data

  blink_epochs <- asc$blinks %>%
    rename(start = stime, end = etime) %>%
    mutate(
      start = start - 75,
      end = end + 75
    )

  ps_dat <- ps_dat %>%
    mutate(
      in_blink = within_epoch(time, blink_epochs)
    ) %>%
    group_by(block) %>%
    mutate(
      ps = interpolate_blinks(time, ps, in_blink)
    ) %>%
    filter(!is.na(ps))


  # Apply 10 Hz lowpass filter to pupil signal to remove high-frequency noise

  ps_dat <- ps_dat %>%
    group_by(block) %>%
    mutate(
      ps = pupil_lowpass(ps, 10, 1000)
    )


  # If requested, adjust time = 0 to be relative to a given event

  if (!is.na(zero)) {
    ps_dat <- ps_dat %>%
      group_by(block) %>%
      mutate(
        time = time - time[phase == zero][1]
      )
  }

  # Reinterpolate blink regions to be smooth based on surrounding filtered data
  # and downsample from 1000 Hz to 100 Hz

  ps_dat <- ps_dat %>%
    group_by(block) %>%
      mutate(
      ps = interpolate_blinks(time, ps, in_blink, smooth = TRUE)
    ) %>%
    filter(time %% 10 == 0) %>%
    select(-in_blink)


  # Map ASC 'blocks' to actual task blocks/trials

  ps_dat <- ps_dat %>%
    ungroup() %>%
    left_join(asc_trial_map, by = "block") %>%
    select(-block)

  ps_dat
}


get_pupil_slopes <- function(asc, window = NA) {

  # Get downsampled/filtered windows of pupil size for each trial

  ps_dat <- get_pupil_data(asc) %>%
    filter(phase == "cue_on") %>%
    group_by(block_num, trial_num) %>%
    mutate(
      time = time - time[1]
    )

  if (!is.na(window)) {
    ps_dat <- ps_dat %>%
      group_by(block_num, trial_num) %>%
      filter(time < window)
  }


  # Standardize pupil size for the participant

  ps_dat <- ps_dat %>%
    ungroup() %>%
    mutate(ps = (ps - mean(ps)) / sd(ps))


  # Fit a regession line to the pupil size changes during the 'cue on' epoch
  # for each trial, getting the slope and intercept of the model

  ps_slopes <- ps_dat %>%
    group_by(block_num, trial_num) %>%
    group_modify(~ get_slope(.x)) %>%
    select(c(block_num, trial_num, intercept, slope))

  ps_slopes
}
