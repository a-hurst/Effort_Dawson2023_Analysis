############################################
### ColourWheelEffort data import script ###
############################################

# NOTE: For this to work, you need to set your current R directory to the
# folder containing the _Data folder, either via RStudio or with setwd()


### Import required packages & load project settings ###

library(readr)
library(purrr)
library(tibble)
library(dplyr)
library(tidyr)
library(eyelinker)

source("./_Scripts/_utils.R")
options(readr.show_progress = FALSE)



### Import task data ###

# Get file list

datapath <- "./_Data/"

datafiles <- list.files(
  datapath, pattern = "*.txt",
  recursive = TRUE, full.names = TRUE
)

ascfiles <- list.files(
  datapath, pattern = "*.asc.gz",
  recursive = TRUE, full.names = TRUE
)


# Import trial-by-trial task data

col_overrides <- cols(
  .default = col_double(),
  participant = col_character(),
  gender = col_factor(levels = c("m", "f", "n")),
  handedness = col_factor(levels = c("l", "r", "a")),
  random_seed = col_skip(),
  practice = col_logical(),
  difficulty = col_factor(levels = c("easy", "difficult")),
  catch_trial = col_logical(),
  probe_loc = col_factor(levels = c("L", "R")),
  trial_err = col_character(),
  probe_col = col_character(),
  response_col = col_character(),
  cue_validity = col_factor(levels = c("valid", "invalid", "neutral"))
)

a <- map_df(datafiles, function(f) {
  df <- read_tsv(f, comment = "#", col_types = col_overrides)
  id_num <- as.numeric(gsub("^p(\\d+).*", "\\1", basename(f)))
  add_column(df, asc_id = id_num, .before = 2)
})


# Import raw pupil data from .asc files

# NOTE: This will take a few minutes

pupildat_raw <- map_df(ascfiles, function(f) {

  # Get ID number from ASC file
  id_num <- as.numeric(gsub("^p(\\d+).*", "\\1", basename(f)))

  # If we have task data for that ID, import and preprocess their per-trial
  # pre-target pupil data
  if (id_num %in% unique(a$asc_id)) {
    cat(paste0("Processing ", basename(f), "...\n"))
    asc <- read_asc(f)
    # Preprocess data provided they have > 10 trials
    if (max(asc$msg$block) > 10) {
      df <- get_pupil_data(
        asc, start = "trial_start", end = "probe_on",
        zero = "cue_on", pad = c(0, 100)
      )
      add_column(df, asc_id = id_num, .before = 1)
    }
  }
})


# Summarize pupil changes

weird_range <- pupildat_raw %>%
  group_by(asc_id, block_num, trial_num) %>%
  summarize(
    range = max(ps) - min(ps)
  ) %>%
  group_by(asc_id) %>%
  mutate(
    range_madm = (range - median(range)) / mad(range)
  ) %>%
  subset(abs(range_madm) > 3)

pupildat <- pupildat_raw %>%
  filter(phase != "probe_on" & time > -200) %>%
  anti_join(weird_range, by = c("asc_id", "block_num", "trial_num")) %>%
  group_by(asc_id) %>%
  mutate(
    ps_z = (ps - mean(ps)) / sd(ps)
  ) %>%
  group_by(asc_id, block_num, trial_num) %>%
  summarize(
    ps_avg_0 = mean(ps_z[time >= -20 & time <= 20]),
    ps_avg_1000 = mean(ps_z[time > 950 & time <= 1000]),
    ps_avg_pretarget = mean(ps_z[time >= (max(time) - 40)]),
    minmax = max(ps_z) - min(ps_z),
    minmax_abs = max(ps) - min(ps)
  )


# Join pupil change data to task data

a2 <- a %>%
  left_join(pupildat, by = c("asc_id", "block_num", "trial_num"))
