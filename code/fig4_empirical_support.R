png("Figs/empirical.png", height = 5.9, width = 5.9, units = "in", res = 300)
bg_rs %>% 
  mutate(`age cohort` = case_when(age == "ancestral" ~ "ancestral (1895-1947)",
                                  T ~ "descendant (2003-2016)")) %>% 
  ggplot(aes(x = z_star, y = ln_rs, color = `age cohort`)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "lm", linewidth = 1.2, se = F) +
  theme_bw(base_size = 14) +
  ylab("ln(root-to-shoot ratio)") +
  xlab("E* (relative tidal elevation)") +
  theme(legend.position = c(0.25, 0.875),
        legend.box.background = element_rect(colour = "black")) +
  scale_color_manual(values = c("#E69F00","#009E73")) 
dev.off()

