# Figure 3 - Comparison of scenario predictions against descendant empirical
# data

# Load libraries
library(tidyverse)

# Read in all Blue Genes experimental data. This is the derived dataset that has
# already been cleaned and formatted.
bg_full <- read_csv(url("https://raw.githubusercontent.com/mlvahsen/BlueGenes/main/derived_data/All_Trait_Data.csv"))

# Create 'not in' operator
`%notin%` <- Negate(`%in%`)

# Subset data for levels 1-4, no competition, and remove pots that had no agb
# Also do this for ONLY descendant genotypes
bg_full %>% 
  filter(level < 5 & comp == 0 & agb_scam > 0 & age == "modern") %>% 
  mutate(root_shoot = total_bg/agb_scam)-> bg_sub

# Filter out pots we don't want because we had to harvest rhizomes or root to
# shoots don't make sense (r:s > 5 observations have very little agb growth and
# reflect large rhizomes at set-up so an experimental artifact).
bg_sub %>% 
  filter(pot_no %notin% c(165, 176) &
           root_shoot < 5 & complete.cases(root_shoot)) -> bg_rs

# Reading in NOAA tide gauge data and setting marsh elevation for 2019 Annapolis
msl_2019 <- mean(c(6.2, 1.8, 8.5, 14.2, 25.7, 22, 23.1, 26.7, 32.1, 34.3, 14.3, 5.8))
mhw_2019 <- mean(c(22.6, 16, 22.5, 30.2, 40.1, 36.2, 38.2, 41.7, 48.1, 50.9, 31.8, 20.3))

# Create z* column in bg_rs data
bg_rs %>% 
  mutate(z_star = (elevation*100 - msl_2019) / (mhw_2019 - msl_2019)) -> bg_rs

# Read in simulation predictions for 2020
preds_2020 <- read_rds("outputs/preds_2020.rds")

bg_rs %>% 
  ggplot(aes(x = z_star, y = log(root_shoot))) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = F, color = "black") +
  # Add predictions from each simulation
  annotate(
    geom = "point",
    x = preds_2020 %>% filter(scenario == "constant") %>% pull(zstar_2020),
    y = preds_2020 %>% filter(scenario == "constant") %>% pull(lnrs_2020),
    size = 3.5, shape = 4, stroke = 2, color = "#502688") +
  annotate(
    geom = "point",
    x = preds_2020 %>% filter(scenario == "plastic") %>% pull(zstar_2020),
    y = preds_2020 %>% filter(scenario == "plastic") %>% pull(lnrs_2020),
    size = 3.5, shape = 4, stroke = 2, color = "#8F61DB") +
  annotate(
    geom = "point",
    x = preds_2020 %>% filter(scenario == "evo_plastic") %>% pull(zstar_2020),
    y = preds_2020 %>% filter(scenario == "evo_plastic") %>% pull(lnrs_2020),
    size = 3.5, shape = 4, stroke = 2, color = "#FC90AF") +
  # Add curved arrows
  annotate(geom = "curve", x = 0.9, y = 1, xend = 0.6, yend = 0.80, 
    curvature = .3, arrow = arrow(length = unit(2, "mm")), color = "#502688") +
  annotate(geom = "curve", x = 0.25, y = 0.85, xend = 0.34, yend = 0.43, 
    curvature = .3, arrow = arrow(length = unit(2, "mm")), color = "#FC90AF") +
  annotate(geom = "curve", x = 0.50, y = -0.20, xend = 0.32, yend = 0.05, 
    curvature = -.3, arrow = arrow(length = unit(2, "mm")), color = "#8F61DB") +
  annotate(geom = "text", x = 1.16, y = 1.01, label = "constant", color = "#502688", size = 3.5) +
  annotate(geom = "text", x = 0.12, y = 0.96, label = "evo + plastic", color = "#FC90AF", size = 3.5) +
  annotate(geom = "text", x = 0.7, y = -0.18, label = "plastic", color = "#8F61DB", size = 3.5) +
  theme_bw(base_size = 10) +
  labs(y = "ln(root-to-shoot ratio)",
       x = "E* (relative tidal elevation)") -> validation_plot

png("figs/Fig3_validation.png", width = 9, height = 6, res = 300, units = "cm")
validation_plot
dev.off()

# Get mean observed root-to-shoot ratio at same zstar for plastic scenario
plastic_obs_mod <- lm(log(root_shoot) ~ z_star, data = bg_rs)
exp(predict(plastic_obs_mod, newdata = tibble(z_star = preds_2020$zstar_2020[2])))
