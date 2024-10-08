# Figure S2 = plot organic and sediment accretion rates over time for each scenario

# Load libraries
library(tidyverse)

# Read in simulation output data
simulation_preds <- read_rds("outputs/simulation_preds.rds")

# Set basic, plastic, and evo_plastic scenarios
basic <- simulation_preds[[1]]
plastic <- simulation_preds[[2]]
evo_plastic <- simulation

set_colors <- c("#502688", "#FC90AF", "#8F61DB")

png("Figs/FigS2_pools.png", height = 4, width = 12, res = 300, units = "in")
tibble(
  sediment = c(basic$dsdt_out, plastic$dsdt_out, evo_plastic$dsdt_out),
  organic = c(basic$dodt_out, plastic$dodt_out, evo_plastic$dodt_out),
  run = rep(c("constant", "plastic", "evo + plastic"), each = length(basic$dodt_out)),
  time = rep(1920:2020, 3)
) %>% 
  gather(key = pool, value = rate, sediment:organic) %>% 
  ggplot(aes(x = time, y = rate, color = run, group = pool, linetype = pool)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~run) +
  labs(y = expression(paste("accretion rate (mm ", yr^-1, ")")),
       x = "year",
       linetype = "mass pool") +
  scale_color_manual(values = set_colors, guide = "none") +
  theme_bw(base_size = 16) +
  scale_x_continuous(breaks = c(1920, 1950, 1980, 2010))
dev.off()
# Warning is related to NAs at t = 1