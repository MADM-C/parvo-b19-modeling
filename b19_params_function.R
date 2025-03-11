load_params <- function() {
  
  l_params <- list(
    ### Transition probabilities
    ## Natural history parameters
    # proportion of population immune
    p_imm     = 0.5,
    # infection rate (0.5-2% during nonepidemic years, 5-10% during epidemic years)
    p_inf    = 0.075,
    #  early infection rate calculation
    #p_einf = ((1/14)*(27-sqrt(729-560*(0.07)))),
    # p_einf = .00770601, 14-20 weeks inf = .0023118
    ## Clinical parameters
    # proportion of maternal infections detected
    p_det_surv = 0.5,
    p_det_sq   = 0.05,
    p_det_inf  = 0.5,  # probability of detection if assumed infected at first visit
    p_det_sus  = 0.5,  # probability of detection if assumed susceptible at first visit
    p_det_imm  = 0.05, # probability of detection if assumed immune at first visit
    # probability of fetus developing severe fetal anemia
    p_sfa      = 0.075,
    # probability of transfusion after development of severe fetal anemia
    p_det_it    = (0.95),  # maternal infection detected
    p_und_it    = (0.80),  # maternal infection undetected
    
    # probability of b19 fetal death being stillbirth
    p_sb = 0.13, 
    
    # probability of live birth
    p_lb        = 0.9,
    # probability of live birth if mother undetected and severe fetal anemia
    p_und_sfa_lb = (0.9)*(0.807),
    # probability of live birth if mother detected and severe fetal anemia
    p_det_sfa_lb = (0.9)*(0.891),
    # probability of live birth if severe fetal anemia is untreated in detected/undetected
    p_det_sfa_lb_nt  = (0.01)*(0.9),
    p_und_sfa_lb_nt = (0.01)*(0.9),
    # probability of miscarriage (fetal death when maternal infection occurs pre 14 weeks)
    p_mc = 0.99,
    
    ## test characteristics
    # IgG - Past Immunity
    p_igg_sens_immune = 0.9576,
    p_igg_spec_immune = 0.9785,
    
    # IgG - Acute Infection
    p_igg_sens_inf = 0.9331,
    p_igg_spec_inf = 0.9785,
    
    # IgM - Past Immunity
    p_igm_sens_immune = 0.9489,
    p_igm_spec_immune = 0.9213,
    
    # IgM - Acute Infection
    p_igm_sens_inf = 0.8723,
    p_igm_spec_inf = 0.9213,
    
    # PCR
    p_pcr_sens = 0.927,
    p_pcr_spec = 1,
    
    # Vaccination
    p_vax     = 0.5,  # maybe 33%?
    p_vax_eff = 0.9,
    
    # QALY
    #qaly_lb = 78.8,
    #qaly_fd = 0,
    
    # Cost
    # average cost of live birth
    c_lb = 11000,
    # average cost of fetal death
    c_fd = 6000,
    # average cost of severe fetal anemia
    c_sfa = 2000,
    # average cost of detecting infection in mother
    c_mid = 130,
    # average cost of monitoring for severe fetal anemia (400 per ultrasound, with average of 6 weeks monitoring)
    c_mon = 2400,
    # average cost of status quo surveillance
    c_sq = 0,
    # average cost of surveillance
    c_surv = 124000,
    
    # Deaths averted
    # number of births
    d_lb = 0,
    # number of deaths
    d_fd = 1,
    # non pvb19 associated birth or death
    d_p_lb = 0,
    # pvb19 associated death
    d_p_fd = 1,
    
    # Expected number of intrauterine transfusions
    # number of transfusions
    n_it = 1,
    # number of no transfusions
    n_nt = 0,
    
    # population parameters
    pop_size = 3600000
  )
  
  return(l_params)
}

get_params_sa <- function() {
  l_params_sa <- list(
    ### Transition probabilities
    ## Natural history parameters
    # proportion of population immune
    p_imm     = c(0.25, 0.75),
    # infection rate (1-4% during nonepidemic years, 10-15% during epidemic years)
    p_inf    = c(0, 0.1),
    ## Clinical parameters
    # proportion of maternal infections detected
    p_det_sq   = 0.05,
    p_det_surv = c(0.05, 0.5),
    p_det_inf  = 0.5,  # probability of detection if assumed infected at first visit
    p_det_sus  = 0.5,  # probability of detection if assumed susceptible at first visit
    p_det_imm  = 0.05, # probability of detection if assumed immune at first visit
    # probability of fetus developing severe fetal anemia
    p_sfa      = c(0.02, 0.1),
    # probability of transfusion given
    p_det_it    = 0.95,  # maternal infection detected
    p_und_it    = c(0.6, 1),  # maternal infection undetected
    
    # probability of live birth
    p_lb         = 0.9,
    # probability of live birth if mother undetected and severe fetal anemia
    p_und_sfa_lb = (0.9)*(0.807),
    # probability of live birth if mother detected and severe fetal anemia
    p_det_sfa_lb = (0.9)*(0.891),
    # probability of live birth if severe fetal anemia is untreated in detected/undetected
    p_det_sfa_lb_nt  = (0.01)*(0.9),
    p_und_sfa_lb_nt = (0.01)*(0.9),
    
    ## test characteristics
    # IgG - Past Immunity
    p_igg_sens_immune = c(qbeta(0.025, 1016, 45), qbeta(0.975, 1016, 45)),
    p_igg_spec_immune = c(qbeta(0.025, 727, 16), qbeta(0.975, 727, 16)),
    
    # IgG - Acute Infection
    p_igg_sens_inf = c(qbeta(0.025, 251, 18), qbeta(0.975, 251, 18)),
    p_igg_spec_inf = c(qbeta(0.025, 727, 16), qbeta(0.975, 727, 16)),
    
    # IgM - Past Immunity
    p_igm_sens_immune = c(qbeta(0.025, 464, 25), qbeta(0.975, 464, 25)),
    p_igm_spec_immune = c(qbeta(0.025, 199, 17), qbeta(0.975, 199, 17)),
    
    # IgM - Acute Infection
    p_igm_sens_inf = c(qbeta(0.025, 123, 18), qbeta(0.975, 123, 18)),
    p_igm_spec_inf = c(qbeta(0.025, 199, 17), qbeta(0.975, 199, 17)),
    
    # PCR
    p_pcr_sens = c(qbeta(0.025, 38, 3), qbeta(0.975, 38, 3)),
    p_pcr_spec = c(qbeta(0.025, 63, 0), qbeta(0.975, 63, 0)),
    
    # Vaccination
    p_vax     = 0.2,  # maybe 33%?
    p_vax_eff = 0.9,
    
    # QALY
    #qaly_lb = 78.8,
    #qaly_fd = 0,
    
    # Cost
    # average cost of live birth
    c_lb = 11000,
    # average cost of fetal death
    c_fd = 6000,
    # average cost of severe fetal anemia
    c_sfa = 2000,
    # average cost of detecting infection in mother
    c_mid = 130,
    # average cost of monitoring for severe fetal anemia (400 per ultrasound, with average of 6 weeks monitoring)
    c_mon = 2400,
    # average cost of status quo surveillance
    c_sq = 0,
    # average cost of surveillance
    c_surv = 124000
  )
}
