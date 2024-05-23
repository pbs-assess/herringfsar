# Load and wrangle simulated data
sim_dat <- sim_fsar_data("wide") %>%
  tibble() %>%
  rename(
    Year = year, Catch = `Catch-MT`, TAC = `TAC-MT`, SSB_med = `SSB-MT`,
    SSB_min = `SSBlow-MT`, SSB_max = `SSBhigh-MT`, LRP = `SSBlrp-MT`,
    USR = `SSBusr-MT`, F_med = `F-1/yr`, F_min = `Flow-1/yr`,
    F_max = `Fhigh-1/yr`, F_lim = `Flim-1/yr`, M_med = `M-1/yr`,
    R_med = `R-E06`, R_min = `Rlow-E06`, R_max = `Rhigh-E06`
  )
