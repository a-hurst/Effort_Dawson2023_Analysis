###############################################
### ColourWheelEffort data modelling script ###
###############################################


### Import required packages ###

library(mgcv)
library(brms)
library(bayestestR)
library(stringr)



### Set sampling settings for all models ###

n_warmup <- 3000
n_samples <- 8000
n_chains <- 6
n_cores <- 6

set.seed(530453080)



### Prepare and run the reaction time model ###

# Calculate rates of trial errors, anticipatory responses, & blinks per id

bad_trials <- a2 %>%
  subset(!practice & !catch_trial) %>%
  group_by(asc_id) %>%
  summarize(
    gaze_err = mean(!is.na(trial_err) & trial_err == "gaze_err") * 100,
    no_resp = mean(!is.na(trial_err) & trial_err == "no_resp") * 100,
    too_soon = mean(!is.na(trial_err) & trial_err == "too_soon") * 100,
    anticipatory = mean(is.na(trial_err) & probe_rt < 200) * 100,
    no_pupil = (
      mean(is.na(trial_err) & probe_rt > 200 & is.na(ps_avg_0)) * 100
    ),
    total = (
      mean(!is.na(trial_err) | probe_rt < 200 | is.na(ps_avg_0)) * 100
    ),
    usable = (
      sum(is.na(trial_err) & probe_rt > 200 & !is.na(ps_avg_0))
    ),
    n_trials = n()
  )


# Flag and exclude problem participants

problem_ids <- bad_trials %>%
  subset(total > 40 | n_trials < 96)


# Filter practice trials, catch trials, anticipatory responses, etc.

dat_ps1 <- a2 %>%
  anti_join(problem_ids, by = "asc_id") %>%
  subset(!practice & !catch_trial) %>%
  subset(is.na(trial_err) & probe_rt > 200 & !is.na(ps_avg_0))


# Standardize/log-transform RTs for model and calculate pupil change metrics

# NOTE: We use z-scored pupil change here becaues it makes it easier to
# calculate effect sizes for our effects of interest (e.g. one SD increase in
# pupil change = 0.2 SDs increase in RT).

dat_ps <- dat_ps1 %>%
  ungroup() %>%
  mutate(
    lrt = log(probe_rt),
    rt_z = (probe_rt - mean(probe_rt)) / sd(probe_rt),
    lrt_z = (lrt - mean(lrt)) / sd(lrt),
    ps_change = ps_avg_pretarget - ps_avg_0,
    trial_n = (block_num - 3) * 120 + trial_num
  ) %>%
  mutate(
    difficulty = as.factor(difficulty),
    block = as.factor(block_num),
    lrt_zz = (lrt - mean(lrt)) / sd(lrt),
    ps_change_z = (ps_change - mean(ps_change)) / sd(ps_change),
  )

# Set equal-prior contrasts on categorical predictors (helps with Bayes models)

contrasts(dat_ps$difficulty) <- contr.equalprior_pairs(2)
contrasts(dat_ps$cue_validity) <- contr.equalprior_pairs(3)


# Define the reaction time model and pre-compile it (takes a minute)

bf_ps <- bf(
  lrt_z ~ difficulty * cue_validity * ps_change_z  +
    (difficulty * cue_validity * ps_change_z | asc_id)
)

mod_ps_tmp <- brm(
  bf_ps, data = dat_ps,
  prior = c(
    prior(normal(0, 2), class = b),
    prior(exponential(1), class = sd),
    prior(exponential(1), class = sigma)
  ),
  chains = 0
)


# Actually run the reaction time model (will take a few minutes)

mod_ps <- update(
  mod_ps_tmp, newdata = dat_ps,
  recompile = FALSE,
  iter = n_samples, warmup = n_warmup,
  chains = n_chains, cores = n_cores
)



#### GAM for pupil changes between trial start & probe onset ####

# Gather and prepare raw pupil size data

trialinfo <- a %>%
  select(c(asc_id, block_num, trial_num, difficulty, cue_validity))

bad_trials <- a %>%
  select(c(asc_id, block_num, trial_num, trial_err)) %>%
  filter(trial_err %in% c("too_soon", "gaze_err"))

pupildat_cleaned <- pupildat_raw %>%
  filter(phase != "probe_on") %>%
  anti_join(weird_range, by = c("asc_id", "block_num", "trial_num")) %>%
  anti_join(problem_ids, by = "asc_id") %>%
  subset(block_num >= 3) %>%
  group_by(asc_id) %>%
  mutate(
    ps_z = (ps - mean(ps)) / sd(ps)
  )


# Check for and flag trials with abnormal baseline

psz_onset <- pupildat_cleaned %>%
  group_by(asc_id, block_num, trial_num) %>%
  summarize(
    psz_baseline = mean(ps_z[time <= -950])
  )

extreme_onset <- psz_onset %>%
  subset(!is.na(psz_baseline)) %>%
  group_by(asc_id) %>%
  mutate(
    ps_z_blz = scale(psz_baseline)
  ) %>%
  subset(abs(psz_baseline) > 2)


# Prepare and summarize data for the model

pupildat_time <- pupildat_cleaned %>%
  group_by(asc_id, block_num, trial_num) %>%
  mutate(
    ps_z = ps_z - mean(ps_z[time <= -950])
  ) %>%
  subset(!is.na(ps_z)) %>%
  anti_join(extreme_onset, by = c("asc_id", "block_num", "trial_num")) %>%
  anti_join(bad_trials, by = c("asc_id", "block_num", "trial_num")) %>%
  left_join(trialinfo, by = c("asc_id", "block_num", "trial_num"))

pupildat_id <- pupildat_time %>%
  subset(block_num >= 3) %>%
  group_by(asc_id, difficulty, time) %>%
  summarize(
    ps_z = mean(ps_z)
  ) %>%
  mutate(
    difficulty = as.factor(str_to_title(difficulty))
  )

contrasts(pupildat_id$difficulty) <- contr.sum(2) / 2


# Actually run the GAM

ps_time_mod <- gam(
  ps_z ~ difficulty + s(time, by = difficulty), data = pupildat_id
)



#### GAM for pupil changes between trial start & probe onset ####

# Check for and flag trials with abnormal baseline

psz_onset_cue <- pupildat_cleaned %>%
  group_by(asc_id, block_num, trial_num) %>%
  summarize(
    psz_baseline = mean(ps_z[time >= -20 & time <= 20])
  )

extreme_onset_cue <- psz_onset_cue %>%
  subset(!is.na(psz_baseline)) %>%
  group_by(asc_id) %>%
  mutate(
    ps_z_blz = scale(psz_baseline)
  ) %>%
  subset(abs(psz_baseline) > 2)


# Prepare and summarize data for the model

pupildat_postcue <- pupildat_cleaned %>%
  group_by(asc_id, block_num, trial_num) %>%
  mutate(
    ps_z = ps_z - mean(ps_z[time >= -20 & time <= 20])
  ) %>%
  subset(!is.na(ps_z)) %>%
  subset(time >= 0 & time < 1000) %>%
  anti_join(extreme_onset_cue, by = c("asc_id", "block_num", "trial_num")) %>%
  anti_join(bad_trials, by = c("asc_id", "block_num", "trial_num")) %>%
  left_join(trialinfo, by = c("asc_id", "block_num", "trial_num"))

pupildat_id_cue <- pupildat_postcue %>%
  mutate(
    spatial_cue = cue_validity != "neutral"
  ) %>%
  group_by(asc_id, difficulty, spatial_cue, time) %>%
  summarize(
    ps_z = mean(ps_z)
  ) %>%
  mutate(
    difficulty = as.factor(str_to_title(difficulty)),
    spatial_cue = as.factor(ifelse(spatial_cue, "True", "False"))
  )

contrasts(pupildat_id_cue$difficulty) <- contr.sum(2) / 2
contrasts(pupildat_id_cue$spatial_cue) <- contr.sum(2) / 2


# Actually run and visualze the GAM

ps_cue_mod_easy <- gam(
  ps_z ~ spatial_cue + s(time, by = spatial_cue),
  data = subset(pupildat_id_cue, difficulty == "Easy")
)

ps_cue_mod_hard <- gam(
  ps_z ~ spatial_cue + s(time, by = spatial_cue),
  data = subset(pupildat_id_cue, difficulty == "Difficult")
)


### Prepare and run colour wheel fidelity model ###

# Gather & prepare data from colour wheel trials

colour_dat <- a2 %>%
  anti_join(problem_ids, by = "asc_id") %>%
  subset(!practice & probe_rt > 200) %>%
  subset(!is.na(response_angle) & abs(angle_err) < 90) %>%
  mutate(
    angle_err_rad = angle_err * (pi / 180)
  )

contrasts(colour_dat$cue_validity) <- contr.equalprior_pairs(3)


# Define the Von Mises model and pre-compile it (takes a minute)

bf_col <- bf(
  angle_err_rad ~ 0,
  kappa ~ cue_validity + (cue_validity | asc_id)
) + von_mises()

mod_col_tmp <- brm(
  bf_col, data = colour_dat,
  family = von_mises(),
  prior = c(
    # Since kappa is on log scale, prior of 2.3 ~= kappa of 10
    prior(normal(2.3, 0.8), class = "Intercept", dpar = "kappa"),
    prior(normal(0, 2), class = "b", dpar = "kappa"),
    prior(exponential(1), class = "sd", dpar = "kappa"),
    prior(lkj(2), class = "cor")
  ),
  chains = 0
)


# Actually run the colour model (will take a minute or two)

mod_col <- update(
  mod_col_tmp, newdata = colour_dat,
  recompile = FALSE,
  iter = n_samples, warmup = n_warmup,
  chains = n_chains, cores = n_cores
)
