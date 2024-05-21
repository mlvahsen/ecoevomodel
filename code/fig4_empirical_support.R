# Figure 4 - empirical support of differences in slopes across cohorts 

# Read in all Blue Genes experimental data. This is the derived dataset that has
# already been cleaned and formatted.
bg_full <- read_csv(url("https://raw.githubusercontent.com/mlvahsen/BlueGenes/main/derived_data/All_Trait_Data.csv"))

# Create 'not in' operator
`%notin%` <- Negate(`%in%`)

# Subset data for levels 1-4, no competition, and remove pots that had no agb
bg_full %>% 
  filter(level < 5 & comp == 0 & agb_scam > 0) %>% 
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

png("Figs/Fig4_empirical.png", height = 11, width = 11, units = "cm", res = 300)
bg_rs %>% 
  mutate(`age cohort` = case_when(age == "ancestral" ~ "ancestral (1895-1947)",
                                  T ~ "descendant (2003-2016)")) %>% 
  ggplot(aes(x = z_star, y = ln_rs, color = `age cohort`)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "lm", linewidth = 1.2, se = F) +
  theme_bw(base_size = 14) +
  ylab("ln(root-to-shoot ratio)") +
  xlab("E* (relative tidal elevation)") +
  theme(legend.position = c(0.70, 0.13),
        legend.box.background = element_rect(colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.margin = margin(1,1,1,1)) +
  scale_color_manual(values = c("#E69F00","#009E73")) 
dev.off()

