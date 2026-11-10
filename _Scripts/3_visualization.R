################################################
### ColourWheelEffort plots & visualizations ###
################################################

### Import required packages ###

library(stringr)
library(emmeans)
library(tidybayes)
library(marginaleffects)
library(modelr)
library(ggplot2)
library(ggdist)
library(cowplot)



### Plot main reaction time effects ###

plotdat_rt <- emm_ps %>%
  describe_posterior(ci = c(0.9, 0.6), test = c()) %>%
  as_tibble() %>%
  mutate(
    difficulty = factor(
      str_to_title(difficulty),
      levels = c("Easy", "Difficult")
    ),
    cue_validity = factor(
      str_to_title(cue_validity),
      levels = c("Valid", "Neutral", "Invalid")
    ),
    .width = CI
  )

pd <- position_dodge(0.4)
ggplot(plotdat_rt, aes(x = cue_validity, y = Median, color = difficulty)) +
  geom_line(aes(group = difficulty), lty = 6, position = pd, alpha = 0.5) +
  geom_pointinterval(
    aes(ymin = CI_low, ymax = CI_high), position = pd
  ) +
  labs(color = "Trial Type") +
  ylab("Reaction Time (ms)") +
  xlab("Cue Type") +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  theme_classic() +
  theme(
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7))
  )



### Plot pupil change across task between trial difficulties ###

ggplot(
  ps_time_pred,
  aes(x = time, y = estimate, group = difficulty, color = difficulty)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.5, color = NA, fill = "grey70"
  ) +
  geom_line() +
  labs(color = "Trial Type") +
  ylab("Pupil Size Change (Standardized)") +
  xlab("Time Relative To Cue Onset (ms)") +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  theme(
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7))
  )



### Plot post-cue pupil changes across cue type & difficulty ###

ggplot(
  ps_cue_pred,
  aes(x = time, y = estimate, group = cue_type, color = cue_type)
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.5, color = NA, fill = "grey70"
  ) +
  geom_line() +
  labs(color = "Cue Type") +
  ylab("Pupil Size Change (Standardized)") +
  xlab("Time Relative To Cue Onset (ms)") +
  scale_color_manual(values = c("#004488", "#BB5566")) +
  theme(
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7))
  ) +
  facet_wrap(~ difficulty)



### Plot pupil-change ~ reaction time slopes across factors ###

mod_ps_draws <- dat_ps %>%
  group_by(cue_validity, difficulty) %>%
  data_grid(ps_change_z = seq(-2, 2, by = 0.2)) %>%
  add_epred_draws(mod_ps, re_formula = NA) %>%
  mutate(
    lrt = (.epred * sd(dat_ps$lrt)) + mean(dat_ps$lrt),
    rt = exp(lrt),
    ps_change = (
      (ps_change_z * sd(dat_ps$ps_change)) + mean(dat_ps$ps_change)
    ),
    difficulty = as.factor(str_to_title(difficulty)),
    cue_validity = as.factor(str_to_title(cue_validity))
  )

ggplot(mod_ps_draws, aes(x = ps_change, y = rt, color = cue_validity)) +
  stat_lineribbon(
    aes(y = rt), .width = 0, alpha = 1, fill = "grey70"
  ) +
  facet_wrap(~ difficulty) +
  labs(color = "Cue Type") +
  ylab("Reaction Time (ms)") +
  xlab("Post-Cue Pupil Change (Standardized Units)") +
  scale_color_manual(values = c("#009E73", "#D55E00", "#0072B2")) +
  theme(
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7))
  ) +
  facet_wrap(~ difficulty)



### Plot cueing effects at +/- 1 SD of pupil change ###

mod_ps_draws2 <- dat_ps %>%
  group_by(cue_validity, difficulty) %>%
  data_grid(ps_change_z = c(-1, 1)) %>%
  add_epred_draws(mod_ps, re_formula = NA) %>%
  ungroup() %>%
  select(-c(.chain, .iteration, .row)) %>%
  mutate(
    .epred = exp((.epred * sd(dat_ps$lrt)) + mean(dat_ps$lrt)),
    ps_change_z = as.factor(as.character(ps_change_z)),
    difficulty = as.factor(str_to_title(difficulty))
  ) %>%
  spread(cue_validity, .epred) %>%
  mutate(
    cue_eff = invalid - valid,
    ps_change_z = factor(
      ps_change_z, levels = c("-1", "1"), labels = c("-1 SD", "+1 SD")
    ),
    `Trial Type` = factor(difficulty, levels = c("Easy", "Difficult")),
  )

mod_ps_hdis <- mod_ps_draws2 %>%
  group_by(ps_change_z, `Trial Type`) %>%
  point_interval(cue_eff, .width = c(0.6, 0.9))

pd <- position_dodge(0.3)
ggplot(mod_ps_hdis, aes(x = ps_change_z, y = cue_eff, color = `Trial Type`)) +
  geom_line(aes(group = `Trial Type`), lty = 6, position = pd, alpha = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_pointinterval(
    aes(ymin = .lower, ymax = .upper), position = pd
  ) +
  ylab("Cueing Effect (ms)") +
  xlab("Post-Cue Pupil Change (Standardized)") +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  theme_classic() +
  theme(
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7))
  )



### Plot kappa estimates & pairwise effects from colour wheel model ###

# Visualize kappa estimates across cue type

cue_colors <- c("#009E73", "#D55E00", "#0072B2")

col_draws <- emm_col %>%
  gather_emmeans_draws() %>%
  mutate(
    cue_validity = as.factor(str_to_title(cue_validity))
  )

colplot1 <- ggplot(
  col_draws,
  aes(x = cue_validity, y = .value, fill = cue_validity)
) +
  geom_violin(width = 0.7, color = NA) +
  geom_violin(width = 0.7, fill = "white", alpha = 0.6, color = NA) +
  stat_pointinterval(.width = c(0.6, 0.9), aes(color = cue_validity)) +
  ylab("Fidelity of Encoding (k)") +
  xlab("Cue Validity") +
  scale_x_discrete(limits = c("Invalid", "Neutral", "Valid")) +
  scale_fill_manual(values = cue_colors) +
  scale_color_manual(values = cue_colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7))
  )


# Visualize pairwise comparisons

col_draws_pairs <- emm_col_pairs %>%
  gather_emmeans_draws() %>%
  mutate(
    contrast = factor(case_when(
      contrast == "invalid - neutral" ~ "I - N",
      contrast == "valid - neutral" ~ "V - N",
      contrast == "valid - invalid" ~ "V - I",
    ),
      levels = c("I - N", "V - N", "V - I")
    )
  )

contrast_colors <- c("#DC267F", "#FFB000", "#648FFF")

colplot2 <- ggplot(
  col_draws_pairs, aes(x = contrast, y = .value, fill = contrast)
) +
  geom_hline(yintercept = 0, lty = 2) +
  stat_pointinterval(.width = c(0.6, 0.9), aes(color = contrast)) +
  ylab("Difference In Fidelity (k)") +
  xlab("Cue Comparison") +
  scale_fill_manual(values = contrast_colors) +
  scale_color_manual(values = contrast_colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7))
  )

plot_grid(
  colplot1, colplot2,
  labels = c("A", "B"), label_size = 12
)


### Plot pupil change ~ cue validity slopes across cue types ###

mod_col_draws <- colour_ps_dat %>%
  group_by(cue_validity) %>%
  data_grid(ps_change_z = seq(-2, 2, by = 0.5)) %>%
  add_epred_draws(mod_col_ps, re_formula = NA, dpar = "kappa") %>%
  mutate(
    cue_validity = as.factor(str_to_title(cue_validity))
  )

colplot3 <- ggplot(
  mod_col_draws, aes(x = ps_change_z, y = kappa, color = cue_validity)
) +
  stat_lineribbon(
    aes(y = kappa), .width = 0, alpha = 1, fill = "grey70"
  ) +
  scale_color_manual(values = c("#009E73", "#D55E00", "#0072B2")) +
  theme(
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7))
  ) +
  labs(color = "Cue Type") +
  ylab("Fidelity of Encoding (k)") +
  xlab("Post-Cue Pupil Change (Standardized Units)")



### Plot fidelity cueing effects at +/- 1 SD of pupil change ###

mod_ps_draws_col <- colour_ps_dat %>%
  group_by(cue_validity) %>%
  data_grid(ps_change_z = c(-1, 1)) %>%
  add_epred_draws(mod_col_ps, re_formula = NA, dpar = "kappa") %>%
  ungroup() %>%
  select(-c(.chain, .iteration, .row)) %>%
  mutate(
    ps_change_z = as.factor(as.character(ps_change_z))
  ) %>%
  spread(cue_validity, kappa) %>%
  mutate(
    cue_eff = valid - invalid,
    ps_change_z = factor(
      ps_change_z, levels = c("-1", "1"), labels = c("-1 SD", "+1 SD")
    )
  )

mod_col_ps_hdis <- mod_ps_draws_col %>%
  group_by(ps_change_z) %>%
  point_interval(cue_eff, .width = c(0.6, 0.9))

colplot4 <- ggplot(
  mod_col_ps_hdis, aes(x = ps_change_z, y = cue_eff, color = "#F8766D")
) +
  geom_line(aes(group = 1), lty = 6, alpha = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_pointinterval(
    aes(ymin = .lower, ymax = .upper)
  ) +
  ylab("Cueing Effect (k)") +
  xlab("Post-Cue Pupil Change") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7))
  )

plot_grid(
  colplot3, colplot4,
  labels = c("A", "B"), label_size = 12,
  rel_widths = c(0.6, 0.4)
)
