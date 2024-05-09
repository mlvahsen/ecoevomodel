## Derive and set parameters for eco-evolutionary model. Use data from Vahsen et
## al. (2023) https://doi.org/10.1111/nph.19117 marsh organ experiment

## Load libraries ####
library(tidyverse); library(geomtextpath); library(patchwork); library(readr);
library(lme4)

## Set initial tidal drivers (cm) and elevation ####
# Reading in NOAA tide gauge data and setting marsh elevation for 2019 Annapolis
msl_2019 <- mean(c(6.2, 1.8, 8.5, 14.2, 25.7, 22, 23.1, 26.7, 32.1, 34.3, 14.3, 5.8))
mhw_2019 <- mean(c(22.6, 16, 22.5, 30.2, 40.1, 36.2, 38.2, 41.7, 48.1, 50.9, 31.8, 20.3))
mlw_2019 <- mean(c(-8.5, -13.7, -6.3, 0.2, 10, 6.2, 6.4, 10.4, 16.1, 20.3, 0.7, -9.4))
# Set elevation of the marsh in 2019
z_init_2019 <- 25

# Get initial water levels 100 years ago by back-calculating given 2019 values
# and historic average rate of sea-level rise (3.75mm/yr * 100 yrs = 37.5cm)
msl <- msl_2019 - 37.5
mhw <- mhw_2019 - 37.5
mlw <- mlw_2019 - 37.5

# Get initial elevation using Kirkpatrick soil core data used in Vahsen et al.
# 2021 to calculate how much accretion has occurred over 100 years
acc_data <- read_csv(url("https://raw.githubusercontent.com/mlvahsen/SchoenoplectusGermination/main/supp_data/Vahsen_etal_Pb210.csv"))

# Find depth that is closest to 100 years
acc_data %>% 
  mutate(age = core_year - year,
         # Set to segment top depths instead of median
         depth_cm = depth_cm - 1) %>% 
  group_by(depth_cm) %>% 
  summarize(mean_age = mean(age)) %>% 
  mutate(diff = 100 - mean_age) %>% 
  arrange(abs(diff)) 
# 18 cm depth is the closest to 100 years of accretion 

# Given this we can assume that elevation 100 years ago is roughly equal the
# current elevation minus 18 cm (this is a 0.18 mm/yr average accretion rate
# which seems reasonable)
z_init <- z_init_2019 - 18

## Derive parameter values for plasticity and environmental sensitivity of selection ####

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

# Fit a linear model of the relationship between z* and ln(root-to-shoot ratio)
lnrs_mod <- lm(ln_rs ~ z_star, data = bg_rs)

# Calculate optimal phenotypes at each flooding level for environmental
# sensitivity of selection
quad_mod1 <- lm(agb_scam ~ ln_rs + I(ln_rs^2),
                data = bg_rs %>% filter(level == 1))
quad_mod2 <- lm(agb_scam ~ ln_rs + I(ln_rs^2),
                data = bg_rs %>% filter(level == 2))
quad_mod3 <- lm(agb_scam ~ ln_rs + I(ln_rs^2),
                data = bg_rs %>% filter(level == 3))
quad_mod4 <- lm(agb_scam ~ ln_rs + I(ln_rs^2),
                data = bg_rs %>% filter(level == 4))

# Collect linear and quadratic coefficients for each flooding level
tibble(level = 1:4,
       beta = c(coef(quad_mod1)[2],coef(quad_mod2)[2],
                coef(quad_mod3)[2],coef(quad_mod4)[2]),
       gamma = c(coef(quad_mod1)[3],coef(quad_mod2)[3],
                 coef(quad_mod3)[3],coef(quad_mod4)[3])) -> quad_coefs

# Create a function to calculate optimal phenotype (i.e., value of r:s that
# maximizes agb)
calc_optimal_phen <- function(beta, gamma){
  z = -beta / (2*gamma)
  return(z)
}

# Apply function to all levels
optimal_phenotypes <- calc_optimal_phen(quad_coefs$beta, quad_coefs$gamma)

# Oragnize mean z_star and ln_rs values for optimum at each level
tibble(z_star = bg_rs %>% group_by(level) %>%
         summarize(mean = mean((elevation*100 - msl_2019) / (mhw_2019 - msl_2019))) %>% pull(mean),
       ln_rs = optimal_phenotypes) -> optimal_df 

# Fit regression to these four points to capture environmental sensitivity of
# selection
lnrs_mod_optimal <- lm(ln_rs ~ z_star, data = optimal_df)

# Create a function to calculate Z* (relative tidal elevation) from Z
z_to_zstar <- function(z, msl, mhw){
  zstar <- (z - msl) / (mhw - msl)
  return(zstar)
}

# Convert initial elevation in 1920 z to zstar
zstar_init <- z_to_zstar(z_init, msl, mhw)

## Set parameters ####

# Get x-intercepts from biomass parabola
biomass_elevation <- lm(agb_scam ~ elevation + I(elevation^2), data = bg_rs)
quad <- function(a, b, c){
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  if(all(Im(answer) == 0)) answer <- Re(answer)
  if(answer[1] == answer[2]) return(answer[1])
  answer}
x_ints <- quad(coef(biomass_elevation)[3], coef(biomass_elevation)[2], coef(biomass_elevation)[1])
# Convert x-intercepts to cm
x_ints_cm <- as.numeric(x_ints) * 100

# Get maximum aboveground biomass (peak of the parabola)
z_bmax <- calc_optimal_phen(coef(biomass_elevation)[2], coef(biomass_elevation)[3])
bmax <- predict(biomass_elevation, data.frame(elevation = z_bmax))
# Convert bmax to g/cm2
pot_area_cm2 <- pi * 5.08^2
bmax_cm2 <- bmax/pot_area_cm2

## CMEM PARAMETERS ##

# Minimum vegetation elevation (cm)
zmin <- x_ints_cm[1]
# Maximum vegetation elevation (cm)
zmax <- x_ints_cm[2]
# Length of simulation
years <- 100
# Total cm of SLR over simulation
total_slr <- 37.5
# Initial rate of sea level rise (cm/yr)
initial_rslr <- 0.375
# Capture rate
q <- 2.8
# Suspended sediment concentration (g/cm3)
ssc <- 3e-5
# Number of tides per year
n_tides <- 704
# Mineral self-packing density (g/cm3)
rho_m <- 1.99
# Organic self-packing density (g/cm3)
rho_o <- 0.085
# Peak aboveground biomass (g/cm2)
bmax <- bmax_cm2
# Recalcitrant fraction
kr <- 0.2
# Belowground turnover rate
bg_tr <- 0.55

## EVOLUTION PARAMETERS ##

# Intercept of plasticity model
rs_int <- as.numeric(coef(lnrs_mod)[1])
# Slope of plasticity model
rs_slope <- as.numeric(coef(lnrs_mod)[2])
# Intercept of optimal phenotype model
rs_int_opt <- as.numeric(coef(lnrs_mod_optimal)[1])
# Slope of optimal phenotype model
rs_slope_opt <- as.numeric(coef(lnrs_mod_optimal)[2])
# Strength of selection
strength_selec_mod <- lm(agb_scam ~ ln_rs + I(ln_rs^2), data = bg_rs)
strength_selec <- as.numeric(coef(strength_selec_mod)[3])

# Phenotypic variance - fit a model that takes into account all environmental
# treatments (e.g., flooding, co2) and exogenous variation due to experimental
# set-up (e.g. cloning group, initial weight). Residual variance from this model
# is phenotypic variance.
mod_env <- lm(ln_rs ~ co2*elevation*salinity + I(elevation^2) + weight_init +
                date_cloned_grp + origin_lab + location, data = bg_rs)
sigma2p <- var(residuals(mod_env))

# Heritability - fit LMM to the residuals of the environmental model and
# calculate ICC using genotype as a random intercept to get heritability.
bg_rs$resid <- residuals(mod_env)
mod_gen <- lmer(resid ~ (1|genotype), data = bg_rs)
h2 <- performance::icc(mod_gen)$ICC_adjusted

## Set sea level rise scenario ####
beta <- (total_slr/years - initial_rslr) / (total_slr - 1)
alpha <- initial_rslr - beta
msl_vec <- msl + alpha*(0:(years-1)) + beta*(0:(years-1))^2
mhw_vec <- msl_vec + (mhw - msl)
mlw_vec <- msl_vec - (msl - mlw)

## Convert zmin and zmax to relative tidal elevations ####
zmaxstar <- z_to_zstar(zmax, msl_2019, mhw_2019)
zminstar <- z_to_zstar(zmin, msl_2019, mhw_2019)
