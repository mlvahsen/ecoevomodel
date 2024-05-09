## Function for eco-evolutionary marsh model ####

predict_marsh <- function(years, z_init, rs_int, rs_slope, rs_int_opt, rs_slope_opt,
                          h2, sigma2p, strength_selec, msl_vec, mhw_vec, mlw_vec,
                          q, ssc, n_tides, rho_m, zmaxstar, zminstar, bmax, kr,
                          bg_tr, rho_o, no_plastic = FALSE){
  
  ## Create storage vectors to hold outputs ##
  
  # ln(root-to-shoot ratio)
  lnrs_store <- rep(NA, years)
  # intercept of ln(root-to-shoot ratio) reaction norm
  lnrs_int_store <- rep(NA, years)
  # intercept of optimal reaction norm
  lnrs_int_opt_store <- rep(NA, years)
  # inorganic accretion rate
  dsdt_out <- rep(NA, years)
  # organic accretion rate
  dodt_out <- rep(NA, years)
  # total vertical accretion rate
  dzdt_out <- rep(NA, years)
  # carbon accumulation rate
  carbon <- rep(NA, years)
  # elevation
  z_vec <- rep(NA, years)
  # relative tidal elevation
  zstar_vec <- rep(NA, years)
  # set initial elevation
  z_vec[1] <- z_init
  # aboveground biomass
  agb_vec <- rep(NA, years)
  # zstar
  zstar <- rep(NA, years)
  
  # Set first value of lnrs_int_store to be rs_int (the intercept of the
  # reaction norm)
  lnrs_int_store[1] <- rs_int
  # Do the same for the intercept for optimal phenotype
  lnrs_int_opt_store[1] <- rs_int_opt
  
  # Loop through model accounting for plasticity and evolution
  for (t in 1:years){
    
    # Convert current elevation to z* (relative tidal elevation)
    zstar <- z_to_zstar(z_vec[t], msl_vec[t], mhw_vec[t])
    zstar_vec[t] <- zstar
    
    # Save initial zstar for basic model with no plasticity or evolution
    if(t == 1){
      zstar_set <- zstar
    }
    
    # Calculate fractional inundation time (FIT)
    fit <-(mhw_vec[t] - z_vec[t]) / (mhw_vec[t] - mlw_vec[t])
    # Restrict FIT to be within or equal to [0,1]
    if(fit > 1) fit <- 1
    if(fit < 0) fit <- 0
    
    # The product of FIT and q should not be greater than 1 or less than 0
    # i.e. capture of sediment in any single tide must be <= 1
    q_new <- q * fit
    if (q_new > 1) q_new <- 1
    if (q_new < 0) q_new <- 0
    
    # Calculate water depth
    d = mhw_vec[t] - z_vec[t]
    
    # Calculate mineral sediment accretion (0.5*d: Flood depth in meters is the
    # same as water volume when in m and assuming a 1m2 area of interest); Note:
    # this is not mediated by plant traits
    sediment_mass <- (0.5 * d * ssc * n_tides * q_new)
    dsdt_out[t+1] <- sediment_mass / rho_m
    
    # Predict aboveground biomass
    zpeakstar <- (zmaxstar + zminstar) / 2
    a <- -((-zminstar * bmax - zmaxstar * bmax) / ((zminstar - zpeakstar) * (-zmaxstar + zpeakstar)))
    b <- -(bmax / ((zminstar - zpeakstar) * (-zmaxstar + zpeakstar)))
    c <- (zminstar * zmaxstar * bmax) / ((zminstar - zpeakstar) * (zmaxstar - zpeakstar))
    biomass <- a*zstar + b*zstar^2 + c
    if(biomass < 0) biomass <- 0
    
    if(no_plastic == TRUE){
      lnrs_store[t+1] <- rs_int + rs_slope*zstar_set
    }else{
      # Calculate mean phenotype at current z*
      if(t == 1){
        lnrs_mean <- rs_int + rs_slope * zstar
        
      }else{
        #lnrs_mean <- lnrs_int_store[t] + rs_slope * zstar
        lnrs_mean <- rs_int + rs_slope * zstar
      }
      
      # Calculate optimal phenotype at current z*
      if(t == 1){
        lnrs_optimal <- rs_int_opt + rs_slope_opt * zstar
      }else{
        #lnrs_optimal <- lnrs_int_opt_store[t] + rs_slope_opt * zstar
        lnrs_optimal <- rs_int_opt + rs_slope_opt * zstar
      }
      
      # Calculate additive genetic variance
      sigma2a <- h2 * sigma2p
      
      # Calculate change in breeding value (i.e. change in reaction norm intercept
      # due to evolution)
      lnrs_int_change <- -abs(strength_selec)*(lnrs_mean - lnrs_optimal)*sigma2a
      
      # Store new intercept value for mean reaction norm
      lnrs_int_store[t+1] <- lnrs_int_store[t] + lnrs_int_change
      
      # Store new intercept value for optimal reaction norm
      #lnrs_int_opt_store[t+1] <- lnrs_int_opt_store[t] + lnrs_int_change
      
      # Calculate new root-to-shoot ratio (add change in intercept to )
      #lnrs_store[t+1] <- lnrs_int_store[t+1] + rs_slope * zstar
      lnrs_store[t+1] <- lnrs_int_store[t+1] + rs_slope * zstar
      
      
    }
    
    # Predict organic accretion
    organic_mass <- kr * exp(lnrs_store[t+1]) * bg_tr * biomass
    dodt_out[t+1] <- organic_mass / rho_o
    
    # Save aboveground biomass
    agb_vec[t] <- biomass
    
    # Calculate carbon sequestration rate
    carbon[t+1] <- organic_mass*0.44
    
    # Calculate total accretion
    dzdt_out[t+1] <- dsdt_out[t+1] + dodt_out[t+1]
    
    # Update elevation given accretion
    z_vec[t+1] <- z_vec[t] + dzdt_out[t+1]
    
    zstar[t] <- zstar
  }
  return(list(lnrs_store = lnrs_store,
              lnrs_int_store = lnrs_int_store,
              lnrs_int_opt_store = lnrs_int_opt_store,
              dsdt_out = dsdt_out,
              dodt_out = dodt_out,
              dzdt_out = dzdt_out,
              carbon = carbon,
              z_vec = z_vec,
              zstar_vec = zstar_vec,
              h2 = h2,
              sigma2p = sigma2p,
              strength_selec = strength_selec,
              biomass = agb_vec))
}