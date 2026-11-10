################################################
### ColourWheelEffort model reporting script ###
################################################

### Import required packages ###

library(emmeans)
library(tidybayes)
library(bayestestR)
library(marginaleffects)
library(modelr)



### Get demographics for reporting ###

pre_exclusions <- a %>%
  group_by(asc_id) %>%
  summarize(
    gender = gender[1],
    age = age[1],
    handedness = handedness[1]
  ) %>%
  subset(age < 100) # Removes testing id

post_exclusions <- a %>%
  anti_join(problem_ids, by = "asc_id") %>%
  group_by(asc_id) %>%
  summarize(
    gender = gender[1],
    age = age[1],
    handedness = handedness[1]
  ) %>%
  subset(age < 100) # Removes testing id



### Summarize reaction time model ###

# Allow back-transforming model estimates to raw RT values

tran_custom <- make.tran("scale", y = dat_ps$probe_rt, inner = "log")

mod_pst <- update(ref_grid(mod_ps), tran = tran_custom)


# Calculate overall estimated marginal means across factors

emm_ps <- emmeans(
  mod_pst, ~ cue_validity * difficulty, level = 0.9, type = "response"
)


# Get main effects & interactions on transformed scale (z-scored log-rt)

emmeans(
  mod_ps, ~ cue_validity, level = 0.9
) %>%
  pairs() %>%
  describe_posterior(ci = 0.9)

emmeans(
  mod_ps, ~ difficulty, level = 0.9
) %>%
  pairs() %>%
  describe_posterior(ci = 0.9)

emmeans(
  mod_ps, ~ cue_validity * difficulty, level = 0.9
) %>%
  pairs(interaction = "consec") %>%
  describe_posterior(ci = 0.9)


# Get main effects & interactions on response scale (msec)

emmeans(
  mod_pst, ~ cue_validity, level = 0.9, type = "response"
) %>%
  regrid() %>%
  pairs() %>%
  describe_posterior(ci = 0.9)

emmeans(
  mod_pst, ~ difficulty, level = 0.9, type = "response"
) %>%
  regrid() %>%
  pairs() %>%
  describe_posterior(ci = 0.9)

emmeans(
  mod_pst, ~ cue_validity * difficulty, level = 0.9, type = "response"
) %>%
  regrid() %>%
  pairs(interaction = "consec") %>%
  describe_posterior(ci = 0.9)

emmeans(
  mod_pst, ~ cue_validity | difficulty, level = 0.9, type = "response"
) %>%
  regrid() %>%
  pairs() %>%
  describe_posterior(ci = 0.9)


# Get slopes for reaction time ~ pupil size effect

emtrends(
  mod_ps, ~ difficulty, level = 0.9, var = "ps_change_z"
) %>%
  pairs() %>%
  describe_posterior(ci = 0.9)

emtrends(
  mod_ps, ~ cue_validity, level = 0.9, var = "ps_change_z"
) %>%
  pairs() %>%
  describe_posterior(ci = 0.9)

emtrends(
  mod_ps, ~ cue_validity * difficulty, level = 0.9, var = "ps_change_z"
) %>%
  pairs(interaction = "trt.vs.ctrl") %>%
  describe_posterior(ci = 0.9)

emtrends(
  mod_ps, ~ difficulty | cue_validity, level = 0.9, var = "ps_change_z"
) %>%
  pairs() %>%
  describe_posterior(ci = 0.9)

emtrends(
  mod_ps, ~ cue_validity | difficulty, level = 0.9, var = "ps_change_z"
) %>%
  describe_posterior(ci = 0.9)



### Summarize GAM outputs ###

# Get actual results from pupil-across-task model

summary(ps_time_mod)


# Get predictions from pupil-across-task model

ps_time_grid <- pupildat_id %>%
  group_by(difficulty) %>%
  data_grid(time = seq(-1000, 1500, by = 10))

ps_time_pred <- predictions(ps_time_mod, newdata = ps_time_grid)


# Get estimates at specific timepoints for pupil-on-task model

as_tibble(ps_time_pred) %>%
  subset(time == 0)

as_tibble(ps_time_pred) %>%
  group_by(difficulty) %>%
  summarize(
    diff = estimate[time == 1000] - estimate[time == 0]
  )


# Get actual results from post-cue pupil model

summary(ps_cue_mod_easy)
summary(ps_cue_mod_hard)


# Get predictions from post-cue pupil model

ps_cue_grid <- pupildat_id_cue %>%
  group_by(spatial_cue) %>%
  data_grid(time = seq(0, 1000, by = 10))

ps_cue_pred_easy <- predictions(ps_cue_mod_easy, newdata = ps_cue_grid) %>%
  as_tibble() %>%
  mutate(difficulty = "Easy")

ps_cue_pred_hard <- predictions(ps_cue_mod_hard, newdata = ps_cue_grid) %>%
  as_tibble() %>%
  mutate(difficulty = "Difficult")

ps_cue_pred <- bind_rows(ps_cue_pred_easy, ps_cue_pred_hard) %>%
  mutate(
    difficulty = factor(difficulty, levels = c("Easy", "Difficult")),
    cue_type = factor(
      ifelse(spatial_cue == "True", "Spatial", "Neutral"),
      levels = c("Spatial", "Neutral")
    )
  )


# Get estimates at specific timepoints for post-cue pupil model

as_tibble(ps_cue_pred) %>%
  group_by(difficulty, cue_type) %>%
  summarize(
    diff = estimate[time == 1000] - estimate[time == 600]
  )



### Summarize colour wheel model results ###

# Get estimated marginal means of kappa across cue conditions

emm_col <- emmeans(
  mod_col, ~ cue_validity, dpar = "kappa", level = 0.9, regrid = "response"
)


# Get pairwise comparisons between cueing conditions

emm_col_pairs <- emm_col %>% pairs()
emm_col_pairs_log <- emmeans(mod_col, ~ cue_validity, dpar = "kappa") %>%
  pairs()

describe_posterior(emm_col_pairs, ci = 0.9)

# NOTE: ROPE range based on visual inpection of differences in Von Mises
#       distributions at different log(kappa) values. Across the values
#       tested, a difference in log(kappa) of 0.03 appears to correspond
#       to a negligible effect (and 0.06 appears to be a small effect)

describe_posterior(
  emm_col_pairs_log, ci = 0.9, rope_range = c(-0.03, 0.03)
)



### Summarize colour wheel vs pupil change model results ###

# Get slopes for fidelity ~ pupil size effect

col_ps_trends <- emtrends(
  mod_col_ps, ~ cue_validity, var = "ps_change_z", dpar = "kappa", level = 0.9
)

describe_posterior(col_ps_trends, ci = 0.9)

col_ps_trends %>%
  pairs() %>%
  describe_posterior(ci = 0.9)
