# Figure S1 - Plot for depicting how plasticity and environmental sensitivity of
# selection parameters were derived

# Read in all Blue Genes experimental data. This is the derived dataset that has
# already been cleaned and formatted.
bg_full <- read_csv(url("https://raw.githubusercontent.com/mlvahsen/BlueGenes/main/derived_data/All_Trait_Data.csv"))

# Create 'not in' operator
`%notin%` <- Negate(`%in%`)

# Subset data for levels 1-4, no competition, and remove pots that had no agb
# Also do this for ONLY ancestral genotypes
bg_full %>% 
  filter(level < 5 & comp == 0 & agb_scam > 0 & age == "ancestral") %>% 
  mutate(root_shoot = total_bg/agb_scam)-> bg_sub

# Filter out pots we don't want because we had to harvest rhizomes or root to
# shoots don't make sense (r:s > 5 observations have very little agb growth and
# reflect large rhizomes at set-up so an experimental artifact).
bg_sub %>% 
  filter(pot_no %notin% c(165, 176) &
           root_shoot < 5 & complete.cases(root_shoot)) -> bg_rs

# Create z* column in bg_rs data
bg_rs %>% 
  mutate(z_star = (elevation*100 - msl_2019) / (mhw_2019 - msl_2019)) -> bg_rs

# Create ln(root:shoot) column
bg_rs$ln_rs <- log(bg_rs$root_shoot)

# Plot up plasticity
bg_rs %>% 
  ggplot(aes(x = z_star, y = ln_rs)) +
  geom_point(alpha = 0.3) +
  geom_textsmooth(method = "lm", color = "black", se = F, label = "mean", size = 5, linewidth = 1.2) +
  ylab("ln(root-to-shoot ratio)") +
  xlab("E* (relative tidal elevation)") +
  theme_bw(base_size = 14) +
  ylim(-1, 1.6) +
  xlim(-0.1473627,2.3844481)-> a

# Split data into four different groups and fit quadratic regressions for each
quad_mod1 <- lm(agb_scam ~ ln_rs + I(ln_rs^2), data = bg_rs %>% filter(level == 1))
quad_mod2 <- lm(agb_scam ~ ln_rs + I(ln_rs^2), data = bg_rs %>% filter(level == 2))
quad_mod3 <- lm(agb_scam ~ ln_rs + I(ln_rs^2), data = bg_rs %>% filter(level == 3))
quad_mod4 <- lm(agb_scam ~ ln_rs + I(ln_rs^2), data = bg_rs %>% filter(level == 4))

# Plot up each quadratic regression
bg_rs %>% 
  mutate(flooding_level = case_when(level == 1 ~ "flooding level 1",
                                    level == 2 ~ "flooding level 2",
                                    level == 3 ~ "flooding level 3",
                                    level == 4 ~ "flooding level 4")) -> bg_rs

# Collect linear and quadratic coefficients for each
tibble(level = 1:4,
       beta = c(coef(quad_mod1)[2],coef(quad_mod2)[2],
                coef(quad_mod3)[2],coef(quad_mod4)[2]),
       gamma = c(coef(quad_mod1)[3],coef(quad_mod2)[3],
                 coef(quad_mod3)[3],coef(quad_mod4)[3])) -> quad_coefs

# Create a function to calculate optimal phenotype
calc_optimal_phen <- function(beta, gamma){
  z = -beta / (2*gamma)
  return(z)
}

# Apply function to all levels
optimal_phenotypes <- calc_optimal_phen(quad_coefs$beta, quad_coefs$gamma)

lnrs_opt1 <- predict(quad_mod1, newdata = data.frame(ln_rs = optimal_phenotypes[1]))
lnrs_opt2 <- predict(quad_mod2, newdata = data.frame(ln_rs = optimal_phenotypes[2]))
lnrs_opt3 <- predict(quad_mod3, newdata = data.frame(ln_rs = optimal_phenotypes[3]))
lnrs_opt4 <- predict(quad_mod4, newdata = data.frame(ln_rs = optimal_phenotypes[4]))

bg_rs %>% 
  ggplot(aes(x = ln_rs, y = agb_scam, color = flooding_level)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_smooth(aes(color = flooding_level), method = "lm", formula = y ~ x + I(x^2), se = F, linewidth = 1.5) +
  xlab("ln(root-to-shoot ratio)") +
  ylab(paste("aboveground biomass (fitness)")) +
  annotate(geom = "point", x = optimal_phenotypes[1], y = lnrs_opt1, size = 5, color = "gold") +
  annotate(geom = "point", x = optimal_phenotypes[2], y = lnrs_opt2, size = 5, color = "gold") +
  annotate(geom = "point", x = optimal_phenotypes[3], y = lnrs_opt3, size = 5, color = "gold") +
  annotate(geom = "point", x = optimal_phenotypes[4], y = lnrs_opt4, size = 5, color = "gold") +
  theme_bw(base_size = 14) + 
  scale_color_manual(values = c("#bdc9e1","#74a9cf", "#2b8cbe", "#045a8d")) +
  theme(legend.position = "none") +
  annotate("label", x= 1.20, y = 6.5, label = "least flooded", size = 3.5, 
           fill = "#bdc9e1", fontface = "bold", color = "white") +
  annotate("label", x= -0.5, y = 7.5, label = "most flooded", size = 3.5,
           fill = "#045a8d", fontface = "bold", color = "white") -> b
  

# Plot optimal phenotype as a function of elevation
tibble(z_star = bg_rs %>% group_by(level) %>%
         dplyr::summarize(mean = mean(z_star)) %>% pull(mean),
       ln_rs = optimal_phenotypes) -> optimal_df 

optimal_df %>% 
  ggplot(aes(x = z_star, y = ln_rs)) + 
  geom_point(size = 5, color = "gold") +
  geom_textsmooth(method = "lm", se = F, color = "gray45",
                  label = "optimal", size = 5, linewidth = 1.2) +
  theme_bw(base_size = 14) + 
  ylab("ln(root-to-shoot ratio)") +
  xlab("E* (relative tidal elevation)") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  ylim(-1, 1.6) +
  xlim(-0.1473627,2.3844481)-> c

# Optimal and mean reaction norm plot
tibble(lnrs = c(predict(lnrs_mod_optimal, newdata = data.frame(z_star = seq(-0.1473627,2.3844481,0.01))),
                predict(lnrs_mod, newdata = data.frame(z_star = seq(-0.1473627,2.3844481,0.01)))),
       z_star = rep(seq(-0.1473627,2.3844481,0.01),2),
       phenotype = rep(c("optimal phenotype", "mean phenotype"), each = length(seq(-0.1473627,2.3844481,0.01)))) %>% 
  ggplot(aes(x = z_star, y =lnrs, color = phenotype)) +
  geom_textline(aes(label = phenotype), hjust = 0.25, linewidth = 1.2) +
  ylab("ln(root-to-shoot ratio)") +
  xlab("E* (relative tidal elevation)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  geom_textvline(aes(xintercept = zstar_init, label = "starting elevation"),
                 hjust = 0.5, color = "black", linetype = "dotted") +
  scale_color_manual(values = c("black", "gray45")) +
  xlim(-0.1473627,2.3844481)-> d

png("Figs/FigS1_parameterization_plot.png", height = 8.5, width = 9, units = "in", res = 300)
(d + a + plot_layout(widths = c(2,1))) / (b + c + plot_layout(widths = c(2,1))) +
  plot_layout(heights = c(3,2)) +
  plot_annotation(tag_levels = "a")
dev.off()
