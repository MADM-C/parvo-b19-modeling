# Model Function
parvo_model <- function(l_params_all) {
  with(
    as.list(l_params_all),
    {
      # p_inf_10 <- p_inf_20 <- p_inf / 2
      p_inf_10 <- p_inf_20 <- (1 - sqrt(1 - p_inf))
      # vector of weights for each policy
      # status quo vector
      v_w_sq <- c((1 - p_imm) * p_inf       * p_det_sq       * p_sfa       * p_det_it       * p_det_sfa_lb,
                  (1 - p_imm) * p_inf       * p_det_sq       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
                  (1 - p_imm) * p_inf       * p_det_sq       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
                  (1 - p_imm) * p_inf       * p_det_sq       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
                  (1 - p_imm) * p_inf       * p_det_sq       * (1 - p_sfa) * p_lb,
                  (1 - p_imm) * p_inf       * p_det_sq       * (1 - p_sfa) * (1 - p_lb),
                  (1 - p_imm) * p_inf       * (1 - p_det_sq) * p_sfa       * p_und_it       * p_und_sfa_lb,
                  (1 - p_imm) * p_inf       * (1 - p_det_sq) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
                  (1 - p_imm) * p_inf       * (1 - p_det_sq) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
                  (1 - p_imm) * p_inf       * (1 - p_det_sq) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
                  (1 - p_imm) * p_inf       * (1 - p_det_sq) * (1 - p_sfa) * p_lb,
                  (1 - p_imm) * p_inf       * (1 - p_det_sq) * (1 - p_sfa) * (1 - p_lb),
                  (1 - p_imm) * (1 - p_inf) * p_lb,
                  (1 - p_imm) * (1 - p_inf) * (1 - p_lb),
                  p_imm       * p_lb,
                  p_imm       * (1 - p_lb))

      # testing at first prenatal visit vector
      v_w_test <- c(
        ### Susceptible, current or recent infection, TP IgG
        ## TP IgM
        # TP PCR
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * p_det_it       * p_det_sfa_lb,
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * (1 - p_sfa)   * p_lb,
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * (1 - p_sfa)   * (1 - p_lb),
        # FN PCR
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * p_und_it       * p_und_sfa_lb,
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * (1 - p_sfa)   * p_lb,
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * (1 - p_sfa)   * (1 - p_lb),
        ## FN IgM
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * p_und_sfa_lb,
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * (1 - p_sfa) * p_lb,
        (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * (1 - p_sfa) * (1 - p_lb),
        ### Susceptible, current or recent infection, FN IgG
        ## TP IgM
        # TP PCR
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * p_det_it       * p_det_sfa_lb,
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * (1 - p_sfa) * p_lb,
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * (1 - p_sfa) * (1 - p_lb),
        # FN PCR
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * p_und_it       * p_det_sfa_lb,
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * p_und_it       * (1 - p_det_sfa_lb),
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * (1 - p_sfa) * p_lb,
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * (1 - p_sfa) * (1 - p_lb),
        ## FN IgM
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * p_und_sfa_lb,
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * (1 - p_sfa) * p_lb,
        (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * (1 - p_sfa) * (1 - p_lb),
        ### Susceptible, no infection at first visit, FP IgG
        ## FP IgM
        # FP PCR
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * p_det_it       * p_det_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa)   * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa)   * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * p_und_it       * p_und_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa)   * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa)   * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * (1 - p_lb),
        # TN PCR
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * p_det_it       * p_det_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * (1 - p_sfa)   * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * (1 - p_sfa)   * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * p_und_it       * p_und_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * (1 - p_sfa)   * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * (1 - p_sfa)   * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * (1 - p_lb),
        ## TN IgM
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * p_det_it       * p_det_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * (1 - p_sfa)    * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * (1 - p_sfa)    * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * p_und_it       * p_und_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * (1 - p_sfa)    * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * (1 - p_sfa)    * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * (1 - p_inf_20)   * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * (1 - p_inf_20)   * (1 - p_lb),
        ### Susceptible, no infection at first visit, TN IgG
        ## FP IgM
        # FP PCR
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * p_det_it       * p_det_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa) * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa) * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * p_und_it       * p_und_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa) * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa) * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * (1 - p_lb),
        # TN PCR
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * p_det_it       * p_det_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * (1 - p_sfa) * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * (1 - p_sfa) * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * p_und_it       * p_und_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * (1 - p_sfa) * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * (1 - p_sfa) * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * (1 - p_lb),
        ## TN IgM
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * p_det_it       * p_det_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * (1 - p_sfa)    * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * (1 - p_sfa)    * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * p_und_it       * p_und_sfa_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * (1 - p_sfa)    * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * (1 - p_sfa)    * (1 - p_lb),
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * (1 - p_inf_20)   * p_lb,
        (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * (1 - p_inf_20)   * (1 - p_lb),
        ### Immune, TP IgG
        ## FP IgM
        # FP PCR
        p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * p_lb,
        p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * (1- p_lb),
        # TN PCR
        p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * p_pcr_spec * p_lb,
        p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * p_pcr_spec * (1- p_lb),
        ## TN IgM
        p_imm * p_igg_sens_immune * p_igm_spec_immune * p_lb,
        p_imm * p_igg_sens_immune * p_igm_spec_immune * (1- p_lb),
        ### Immune, FN IgG
        ## FP IgM
        # FP PCR
        p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * p_lb,
        p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * (1- p_lb),
        # TN PCR
        p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * p_pcr_spec * p_lb,
        p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * p_pcr_spec * (1- p_lb),
        ## TN IgM
        p_imm * (1 - p_igg_sens_immune) * p_igm_spec_immune * p_lb,
        p_imm * (1 - p_igg_sens_immune) * p_igm_spec_immune * (1- p_lb)
      )

      v_w_vax <- c(p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * p_det_sq       * p_sfa       * p_det_it       * p_det_sfa_lb,
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * p_det_sq       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * p_det_sq       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * p_det_sq       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * p_det_sq       * (1 - p_sfa) * p_lb,
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * p_det_sq       * (1 - p_sfa) * (1 - p_lb),
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * (1 - p_det_sq) * p_sfa       * p_und_it       * p_und_sfa_lb,
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * (1 - p_det_sq) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * (1 - p_det_sq) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * (1 - p_det_sq) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * (1 - p_det_sq) * (1 - p_sfa) * p_lb,
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf       * (1 - p_det_sq) * (1 - p_sfa) * (1 - p_lb),
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf) * p_lb,
                   p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf) * (1 - p_lb),
                   p_vax       * p_vax_eff       * (1 - p_imm) * p_lb,
                   p_vax       * p_vax_eff       * (1 - p_imm) * (1 - p_lb),
                   p_vax       * p_imm           * p_lb,
                   p_vax       * p_imm           * (1 - p_lb),
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * p_det_sq       * p_sfa       * p_det_it       * p_det_sfa_lb,
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * p_det_sq       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * p_det_sq       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * p_det_sq       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * p_det_sq       * (1 - p_sfa) * p_lb,
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * p_det_sq       * (1 - p_sfa) * (1 - p_lb),
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * (1 - p_det_sq) * p_sfa       * p_und_it       * p_und_sfa_lb,
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * (1 - p_det_sq) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * (1 - p_det_sq) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * (1 - p_det_sq) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * (1 - p_det_sq) * (1 - p_sfa) * p_lb,
                   (1 - p_vax) * (1 - p_imm)     * p_inf       * (1 - p_det_sq) * (1 - p_sfa) * (1 - p_lb),
                   (1 - p_vax) * (1 - p_imm)     * (1 - p_inf) * p_lb,
                   (1 - p_vax) * (1 - p_imm)     * (1 - p_inf) * (1 - p_lb),
                   (1 - p_vax) * p_imm           * p_lb,
                   (1 - p_vax) * p_imm           * (1 - p_lb))

      # Vax + testing at first prenatal visit vector
      v_w_vax_test <- c(
        ### Susceptible, current or recent infection, TP IgG
        ## TP IgM
        # TP PCR
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * p_det_it       * p_det_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * p_det_it       * (1 - p_det_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * (1 - p_det_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * (1 - p_sfa)   * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * (1 - p_sfa)   * (1 - p_lb),
        # FN PCR
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * p_und_it       * p_und_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * p_und_it       * (1 - p_und_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * (1 - p_sfa)   * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * (1 - p_sfa)   * (1 - p_lb),
        ## FN IgM
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * p_und_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * (1 - p_sfa) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * (1 - p_sfa) * (1 - p_lb),
        ### Susceptible, current or recent infection, FN IgG
        ## TP IgM
        # TP PCR
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * p_det_it       * p_det_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * (1 - p_sfa) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * (1 - p_sfa) * (1 - p_lb),
        # FN PCR
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * p_und_it       * p_det_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * p_und_it       * (1 - p_det_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * (1 - p_sfa) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * (1 - p_sfa) * (1 - p_lb),
        ## FN IgM
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * p_und_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * (1 - p_sfa) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * (1 - p_sfa) * (1 - p_lb),
        ### Susceptible, no infection at first visit, FP IgG
        ## FP IgM
        # FP PCR
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * p_det_it       * p_det_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * p_det_it       * (1 - p_det_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * (1 - p_det_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa)   * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa)   * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * p_und_it       * p_und_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * p_und_it       * (1 - p_und_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa)   * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa)   * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * (1 - p_lb),
        # TN PCR
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * p_det_it       * p_det_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * p_det_it       * (1 - p_det_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * (1 - p_det_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * (1 - p_sfa)   * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * (1 - p_sfa)   * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * p_und_it       * p_und_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * p_und_it       * (1 - p_und_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * (1 - p_sfa)   * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * (1 - p_sfa)   * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * (1 - p_lb),
        ## TN IgM
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * p_det_it       * p_det_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * p_det_it       * (1 - p_det_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * (1 - p_det_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * (1 - p_sfa)    * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * (1 - p_sfa)    * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * p_und_it       * p_und_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * p_und_it       * (1 - p_und_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * (1 - p_sfa)    * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * (1 - p_sfa)    * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * (1 - p_inf_20)   * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * (1 - p_inf_20)   * (1 - p_lb),
        ### Susceptible, no infection at first visit, TN IgG
        ## FP IgM
        # FP PCR
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * p_det_it       * p_det_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa) * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * p_und_it       * p_und_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa) * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * (1 - p_lb),
        # TN PCR
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * p_det_it       * p_det_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * (1 - p_sfa) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * (1 - p_sfa) * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * p_und_it       * p_und_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * (1 - p_sfa) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * (1 - p_sfa) * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * (1 - p_lb),
        ## TN IgM
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * p_det_it       * p_det_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * p_det_it       * (1 - p_det_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * (1 - p_det_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * (1 - p_sfa)    * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * (1 - p_sfa)    * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * p_und_it       * p_und_sfa_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * p_und_it       * (1 - p_und_sfa_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * (1 - p_und_it) * p_sfa_lb_nt,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * (1 - p_sfa)    * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * (1 - p_sfa)    * (1 - p_lb),
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * (1 - p_inf_20)   * p_lb,
        p_vax       * (1 - p_vax_eff) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * (1 - p_inf_20)   * (1 - p_lb),
        ### Immune, TP IgG
        ## FP IgM
        # FP PCR
        p_vax * p_vax_eff * (1 - p_imm) * p_igg_sens_immune * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * (1- p_lb),
        p_vax * p_vax_eff * (1 - p_imm) * p_igg_sens_immune * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * p_lb,
        # TN PCR
        p_vax * p_vax_eff * (1 - p_imm) * p_igg_sens_immune * (1 - p_igm_spec_immune) * p_pcr_spec * p_lb,
        p_vax * p_vax_eff * (1 - p_imm) * p_igg_sens_immune * (1 - p_igm_spec_immune) * p_pcr_spec * (1- p_lb),
        ## TN IgM
        p_vax * p_vax_eff * (1 - p_imm) * p_igg_sens_immune * p_igm_spec_immune * p_lb,
        p_vax * p_vax_eff * (1 - p_imm) * p_igg_sens_immune * p_igm_spec_immune * (1- p_lb),
        ### Immune, FN IgG
        ## FP IgM
        # FP PCR
        p_vax * p_vax_eff * (1 - p_imm) * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * p_lb,
        p_vax * p_vax_eff * (1 - p_imm) * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * (1- p_lb),
        # TN PCR
        p_vax * p_vax_eff * (1 - p_imm) * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * p_pcr_spec * p_lb,
        p_vax * p_vax_eff * (1 - p_imm) * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * p_pcr_spec * (1- p_lb),
        ## TN IgM
        p_vax * p_vax_eff * (1 - p_imm) * (1 - p_igg_sens_immune) * p_igm_spec_immune * p_lb,
        p_vax * p_vax_eff * (1 - p_imm) * (1 - p_igg_sens_immune) * p_igm_spec_immune * (1- p_lb),
        ### Immune, TP IgG
        ## FP IgM
        # FP PCR
        p_vax * p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * p_lb,
        p_vax * p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * (1- p_lb),
        # TN PCR
        p_vax * p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * p_pcr_spec * p_lb,
        p_vax * p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * p_pcr_spec * (1- p_lb),
        ## TN IgM
        p_vax * p_imm * p_igg_sens_immune * p_igm_spec_immune * p_lb,
        p_vax * p_imm * p_igg_sens_immune * p_igm_spec_immune * (1- p_lb),
        ### Immune, FN IgG
        ## FP IgM
        # FP PCR
        p_vax * p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * p_lb,
        p_vax * p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * (1- p_lb),
        # TN PCR
        p_vax * p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * p_pcr_spec * p_lb,
        p_vax * p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * p_pcr_spec * (1- p_lb),
        ## TN IgM
        p_vax * p_imm * (1 - p_igg_sens_immune) * p_igm_spec_immune * p_lb,
        p_vax * p_imm * (1 - p_igg_sens_immune) * p_igm_spec_immune * (1- p_lb),

        ### Susceptible, current or recent infection, TP IgG
        ## TP IgM
        # TP PCR
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * p_det_it       * p_det_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * p_sfa         * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * (1 - p_sfa)   * p_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * p_pcr_sens  * (1 - p_sfa)   * (1 - p_lb),
        # FN PCR
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * p_und_it       * p_und_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * p_sfa         * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * (1 - p_sfa)   * p_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * p_igm_sens_inf       * (1 - p_pcr_sens)  * (1 - p_sfa)   * (1 - p_lb),
        ## FN IgM
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * p_und_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * (1 - p_sfa) * p_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * p_igg_sens_inf * (1 - p_igm_sens_inf) * (1 - p_sfa) * (1 - p_lb),
        ### Susceptible, current or recent infection, FN IgG
        ## TP IgM
        # TP PCR
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * p_det_it       * p_det_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * (1 - p_sfa) * p_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * p_pcr_sens       * (1 - p_sfa) * (1 - p_lb),
        # FN PCR
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * p_und_it       * p_det_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * p_und_it       * (1 - p_det_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * (1 - p_sfa) * p_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * p_igm_sens_inf * (1 - p_pcr_sens) * (1 - p_sfa) * (1 - p_lb),
        ## FN IgM
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * p_und_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * (1 - p_sfa) * p_lb,
        (1 - p_vax) * (1 - p_imm) * p_inf_10 * (1 - p_igg_sens_inf) * (1 - p_igm_sens_inf) * (1 - p_sfa) * (1 - p_lb),
        ### Susceptible, no infection at first visit, FP IgG
        ## FP IgM
        # FP PCR
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * p_det_it       * p_det_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa         * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa)   * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa)   * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * p_und_it       * p_und_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa         * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa)   * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa)   * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * (1 - p_lb),
        # TN PCR
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * p_det_it       * p_det_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * p_sfa         * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * (1 - p_sfa)   * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_imm       * (1 - p_sfa)   * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * p_und_it       * p_und_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * p_sfa         * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * (1 - p_sfa)   * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_imm) * (1 - p_sfa)   * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * (1 - p_lb),
        ## TN IgM
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * p_det_it       * p_det_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * p_sfa          * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * (1 - p_sfa)    * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * p_det_imm       * (1 - p_sfa)    * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * p_und_it       * p_und_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * p_sfa          * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * (1 - p_sfa)    * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * p_inf_20         * (1 - p_det_imm) * (1 - p_sfa)    * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * (1 - p_inf_20)   * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * (1 - p_igg_spec_immune) * p_igm_spec_inf       * (1 - p_inf_20)   * (1 - p_lb),
        ### Susceptible, no infection at first visit, TN IgG
        ## FP IgM
        # FP PCR
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * p_det_it       * p_det_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa) * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * p_det_inf       * (1 - p_sfa) * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * p_und_it       * p_und_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa) * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * p_inf_20       * (1 - p_det_inf) * (1 - p_sfa) * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * (1 - p_pcr_spec) * (1 - p_inf_20) * (1 - p_lb),
        # TN PCR
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * p_det_it       * p_det_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * p_sfa       * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * (1 - p_sfa) * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * p_det_sus       * (1 - p_sfa) * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * p_und_it       * p_und_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * p_sfa       * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * (1 - p_sfa) * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * p_inf_20       * (1 - p_det_sus) * (1 - p_sfa) * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * (1 - p_igm_spec_inf) * p_pcr_spec       * (1 - p_inf_20) * (1 - p_lb),
        ## TN IgM
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * p_det_it       * p_det_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * p_det_it       * (1 - p_det_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * (1 - p_det_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * p_sfa          * (1 - p_det_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * (1 - p_sfa)    * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * p_det_sus       * (1 - p_sfa)    * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * p_und_it       * p_und_sfa_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * p_und_it       * (1 - p_und_sfa_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * (1 - p_und_it) * p_sfa_lb_nt,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * p_sfa          * (1 - p_und_it) * (1 - p_sfa_lb_nt),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * (1 - p_sfa)    * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * p_inf_20         * (1 - p_det_sus) * (1 - p_sfa)    * (1 - p_lb),
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * (1 - p_inf_20)   * p_lb,
        (1 - p_vax) * (1 - p_imm) * (1 - p_inf_10) * p_igg_spec_immune * p_igm_spec_inf       * (1 - p_inf_20)   * (1 - p_lb),
        ### Immune, TP IgG
        ## FP IgM
        # FP PCR
        (1 - p_vax) * p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * p_lb,
        (1 - p_vax) * p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * (1- p_lb),
        # TN PCR
        (1 - p_vax) * p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * p_pcr_spec * p_lb,
        (1 - p_vax) * p_imm * p_igg_sens_immune * (1 - p_igm_spec_immune) * p_pcr_spec * (1- p_lb),
        ## TN IgM
        (1 - p_vax) * p_imm * p_igg_sens_immune * p_igm_spec_immune * p_lb,
        (1 - p_vax) * p_imm * p_igg_sens_immune * p_igm_spec_immune * (1- p_lb),
        ### Immune, FN IgG
        ## FP IgM
        # FP PCR
        (1 - p_vax) * p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * p_lb,
        (1 - p_vax) * p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * (1 - p_pcr_spec) * (1- p_lb),
        # TN PCR
        (1 - p_vax) * p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * p_pcr_spec * p_lb,
        (1 - p_vax) * p_imm * (1 - p_igg_sens_immune) * (1 - p_igm_spec_immune) * p_pcr_spec * (1- p_lb),
        ## TN IgM
        (1 - p_vax) * p_imm * (1 - p_igg_sens_immune) * p_igm_spec_immune * p_lb,
        (1 - p_vax) * p_imm * (1 - p_igg_sens_immune) * p_igm_spec_immune * (1- p_lb)
      )

      ## Reward vectors
      # vector of deaths - status quo
      v_death_sq <- rep(c(d_lb, d_fd), 8)

      # vector of transfusions - status quo
      v_transfusion_sq <- c(rep(c(n_it, n_it, n_nt, n_nt, n_nt, n_nt), 2),
                            n_nt, n_nt, n_nt, n_nt)

      # vector of deaths - testing
      v_death_test <- rep(c(d_lb, d_fd), 66)

      # vector of transfusions - testing
      v_transfusion_test <- c(rep(c(n_it, n_it, n_nt, n_nt, n_nt, n_nt), 6),
                              rep(c(rep(c(n_it, n_it, n_nt, n_nt, n_nt, n_nt), 2), n_nt, n_nt), 6),
                              rep(n_nt, 12))

      # vector of deaths - vaccination
      v_death_vax <- rep(c(d_lb, d_fd), 17)

      # vector of transfusions - vaccination
      v_transfusion_vax <- c(rep(c(n_it, n_it, n_nt, n_nt, n_nt, n_nt), 2),
                             rep(n_nt, 6),
                             rep(c(n_it, n_it, n_nt, n_nt, n_nt, n_nt), 2),
                             rep(n_nt, 4))

      # vector of deaths - vaccination + testing
      v_death_vax_test <- rep(c(d_lb, d_fd), 138)

      # vector of transfusions - vaccination + testing
      v_transfusion_vax_test <- c(rep(c(n_it, n_it, n_nt, n_nt, n_nt, n_nt), 6),
                                  rep(c(rep(c(n_it, n_it, n_nt, n_nt, n_nt, n_nt), 2), n_nt, n_nt), 6),
                                  rep(n_nt, 24),
                                  rep(c(n_it, n_it, n_nt, n_nt, n_nt, n_nt), 6),
                                  rep(c(rep(c(n_it, n_it, n_nt, n_nt, n_nt, n_nt), 2), n_nt, n_nt), 6),
                                  rep(n_nt, 12))

      # total deaths per policy
      total_death_sq       <- v_w_sq %*% v_death_sq %*% pop_size
      total_death_test     <- v_w_test %*% v_death_test %*% pop_size
      total_death_vax      <- v_w_vax %*% v_death_vax %*% pop_size
      total_death_vax_test <- v_w_vax_test %*% v_death_vax_test %*% pop_size

      # total transfusions per policy
      total_it_sq       <- v_w_sq %*% v_transfusion_sq %*% pop_size
      total_it_test     <- v_w_test %*% v_transfusion_test %*% pop_size
      total_it_vax      <- v_w_vax %*% v_transfusion_vax %*% pop_size
      total_it_vax_test <- v_w_vax_test %*% v_transfusion_vax_test %*% pop_size

      # Incremental Results
      # inc_deaths <- total_death_sq - total_death_test
      # inc_it     <- total_it_sq - total_it_test

      # print(sum(v_w_sq))
      # print(sum(v_w_test))
      # print(sum(v_w_vax))
      # print(sum(v_w_vax_test))

      df_output <- data.frame(Strategy =  v_names_str,
                              Deaths  =  c(total_death_sq, total_death_test,
                                           total_death_vax,
                                           total_death_vax_test),
                              Transfusions = c(total_it_sq, total_it_test,
                                               total_it_vax,
                                               total_it_vax_test))
      return(df_output)})
}
