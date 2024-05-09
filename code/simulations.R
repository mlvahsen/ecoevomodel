source("code/set_parameters.R")
source("code/model.R")

set_colors <- c("#502688", "#FC90AF", "#8F61DB")

## Run simulations ####
# Run model with no plasticity and no evolution
predict_marsh(years = years, z_init = z_init, rs_int = rs_int, rs_int_opt = rs_int_opt,
              rs_slope = rs_slope, rs_slope_opt = rs_slope_opt,
              h2 = 0, # Set model to not allow evolution
              sigma2p = sigma2p,
              strength_selec = strength_selec,
              msl_vec = msl_vec, mhw_vec = mhw_vec, mlw_vec = mlw_vec,
              q = q, ssc = ssc, n_tides = n_tides, rho_m = rho_m, 
              zmaxstar = zmaxstar, zminstar = zminstar,
              bmax = bmax, kr = kr, bg_tr = bg_tr, rho_o = rho_o,
              # Set model to not allow plasticity
              no_plastic = TRUE) -> basic

# Run model with plasticity but no evolution
predict_marsh(years = years, z_init = z_init, rs_int = rs_int, rs_int_opt = rs_int_opt,
              rs_slope = rs_slope, rs_slope_opt = rs_slope_opt,
              h2 = 0, # Set model to not allow evolution
              sigma2p = sigma2p,
              strength_selec = strength_selec,
              msl_vec = msl_vec, mhw_vec = mhw_vec, mlw_vec = mlw_vec,
              q = q, ssc = ssc, n_tides = n_tides, rho_m = rho_m, 
              zmaxstar = zmaxstar, zminstar = zminstar,
              bmax = bmax, kr = kr, bg_tr = bg_tr, rho_o = rho_o) -> plastic

# Run model with plasticity AND evolution 
predict_marsh(years = years, z_init = z_init, rs_int = rs_int, rs_int_opt = rs_int_opt,
              rs_slope = rs_slope, rs_slope_opt = rs_slope_opt,
              h2 = h2,
              sigma2p = sigma2p,
              strength_selec = strength_selec,
              msl_vec = msl_vec, mhw_vec = mhw_vec, mlw_vec = mlw_vec,
              q = q, ssc = ssc, n_tides = n_tides, rho_m = rho_m, 
              zmaxstar = zmaxstar, zminstar = zminstar,
              bmax = bmax, kr = kr, bg_tr = bg_tr, rho_o = rho_o) -> evo_plastic

## Make plots ####
# Plot of trait change over time
tibble(constant = basic$lnrs_store[2:101],
       plastic = plastic$lnrs_store[2:101],
       `evo + plastic` = evo_plastic$lnrs_store[2:101],
       time = 1921:2020) %>%
  gather(key = type, value = rs, constant:`evo + plastic`) %>% 
  ggplot(aes(x = time, y = rs, group = type, color = type, label = type)) +
  geom_textline(size = 4, hjust = 0.9, linewidth = 1.2) +
  theme_bw(base_size = 14) +
  ylab("ln(root-to-shoot ratio)") +
  xlab("year") +
  theme(legend.position = "none") +
  scale_color_manual(values = set_colors) +
  scale_x_continuous(breaks = seq(1920,2020,length.out = 6))-> a

# Plot organic versus inorganic accretion rates over time
tibble(basic = (basic$dsdt_out[2:101] + basic$dodt_out[2:101])*10,
       plastic = (plastic$dsdt_out[2:101] + plastic$dodt_out[2:101])*10,
       evo_plastic = (evo_plastic$dsdt_out[2:101] + evo_plastic$dodt_out[2:101])*10,
       time = 1921:2020) %>%
  gather(key = type,
         value = value,
         basic:evo_plastic) %>%
  ggplot(aes(x = time, y = value, color = type)) +
  geom_line(linewidth = 1.2) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  ylab(expression(paste("accretion rate (mm ", yr^-1, ")"))) +
  xlab("year") +
  theme(legend.position = "none") +
  scale_color_manual(values = set_colors)+
  scale_x_continuous(breaks = seq(1920,2020,length.out = 6)) -> b

# Plot carbon sequestration over time
tibble(basic = basic$carbon[2:years]*1e4,
       plastic = plastic$carbon[2:years]*1e4,
       evo_plastic = evo_plastic$carbon[2:years]*1e4,
       time = 1921:2019) %>%
  gather(key = type, value = carbon, basic:evo_plastic) %>% 
  ggplot(aes(x = time, y = carbon, group = type, color = type)) +
  geom_line(linewidth = 1.2) +
  ylab(expression(paste("carbon accum. rate (g C ", m^-2, yr^-1,")"))) +
  theme_bw(base_size = 14) +
  xlab("year") +
  theme(legend.position = "none") +
  scale_color_manual(values = set_colors)+
  scale_x_continuous(breaks = seq(1920,2020,length.out = 6))-> c

png("figs/Fig2_simulations.png", height = 3.5, width = 11.5, units = "in", res = 300)
(a + b + c) + plot_annotation(tag_levels = "a") & theme(plot.margin = margin(3,3,3,3)) 
dev.off()

## Calculate effect sizes ####
# Calculate percent difference in root-to-shoot ratios at end of simulation
rs100_basic <- exp(basic$lnrs_store[100])
rs100_plastic <- exp(plastic$lnrs_store[100])
rs100_evo <- exp(evo_plastic$lnrs_store[100])

round((rs100_basic - rs100_plastic)*100/ rs100_basic,1) # 42.8
round((rs100_basic - rs100_evo)*100/ rs100_basic,1) # 21.6
round((rs100_evo - rs100_plastic)*100/ rs100_plastic,1) # 37.2

# Calculate percent difference in average vertical accretion rates
accV_basic <- mean(basic$dzdt_out, na.rm = T)
accV_plastic <- mean(plastic$dzdt_out, na.rm = T)
accV_evo <- mean(evo_plastic$dzdt_out, na.rm = T)

round((accV_basic - accV_plastic)*100/accV_basic,1) # 20.8
round((accV_basic - accV_evo)*100/accV_basic,1) # 12.4

# Calculate percent difference in average carbon accumulation rates
accC_basic <- mean(basic$carbon, na.rm = T)
accC_plastic <- mean(plastic$carbon, na.rm = T)
accC_evo <- mean(evo_plastic$carbon, na.rm = T)

round((accC_basic - accC_plastic)*100/accC_basic,1) # 24.1
round((accC_basic - accC_evo)*100/accC_basic,1) # 14.2

## Save predictions at t = 2020 for Figure 3 validation plot ####
tibble(scenario = c("constant", "plastic", "evo_plastic"),
       lnrs_2020 = c(basic$lnrs_store[101],
                     plastic$lnrs_store[101],
                     evo_plastic$lnrs_store[101]),
       zstar_2020 = c(basic$zstar_vec[100],
                      plastic$zstar_vec[100],
                      evo_plastic$zstar_vec[100])) -> preds_2020

write_rds(preds_2020, "outputs/preds_2020.rds")
