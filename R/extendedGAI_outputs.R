# Plot outputs from the extended GAI applied to two butterfly species

library(data.table)
library(ggplot2)
library(ggpubr)

# Species
splist <- c("Chalk Hill Blue", "Gatekeeper")

years <- 1976:2022

gaiplots <- list()
for(spp in splist){

  # Extended model output
  mda <- readRDS(paste0("Output/", spp, "_myGAI_annual",
                        "_1976to2022.rds"))
  # Extended model with linear beta output
  mdaal <- readRDS(paste0("Output/", spp, "_myGAI_aal",
                          "_1976to2022.rds"))
  # Estimate posthoc linear trend from exp_beta and use that for plotting linear extended GAI output
  syears <- years - min(years) + 1
  lm_beta <- lm(mda$betas.out$beta ~ syears)
  temp <- data.table(YEAR = years, p = predict(lm_beta),
                    z = predict(lm_beta)[1] + mdaal$betas.out$beta - mdaal$betas.out$beta[1])
  temp[, z_lower := z - 1.96*mdaal$betas.out$beta_SE*sqrt(mdaal$D)]
  temp[, z_upper := z + 1.96*mdaal$betas.out$beta_SE*sqrt(mdaal$D)]
  temp[, z_mid := predict(lm_beta)[median(syears)] + mdaal$betas.out$beta -
                                                    mdaal$betas.out$beta[median(syears)]]
  temp[, z_mid_lower := z_mid - 1.96*mdaal$betas.out$beta_SE*sqrt(mdaal$D)]
  temp[, z_mid_upper := z_mid + 1.96*mdaal$betas.out$beta_SE*sqrt(mdaal$D)]


  beta_plot <- ggplot(mda$betas.out, aes(YEAR, beta))+
    theme_bw()+
    geom_point()+
    geom_line()+
    geom_ribbon(aes(ymin = beta - 1.96*beta_SE*sqrt(mda$D),
                    ymax = beta + 1.96*beta_SE*sqrt(mda$D)), alpha = .3)+
    #geom_line(data = temp, aes(YEAR, z), col = "#48a464")+
    geom_line(data = temp, aes(YEAR, z_mid), col = "#48a464")+
    geom_ribbon(aes(y= z_mid, ymin = z_mid_lower,
                    ymax = z_mid_upper), data = temp, alpha = .3, fill = "#48a464")+
   # geom_ribbon(aes(y= z, ymin = z_lower,
    #                ymax = z_upper), data = temp, alpha = .3, fill = "#48a464")+
    geom_line(data = temp, aes(YEAR, p), col = "#A44888", linetype="dashed", linewidth = 1)+
    theme(text = element_text(size = 22),
          plot.margin=margin(1.5,.5,.25,.5,"cm"))+
    xlab("Year (r)")+
    ylab(expression(beta[r]))+
    ggtitle("(i)")

  # Get extended GAI index on log 10 scale
  ci_mda <- mda$betas.out
  ci_mda[, beta := log(exp_beta)]
  ci_mda[, TRMOBS := 2 + (beta - mean(beta))/log(10)]

  # Calculate SE on log20 scale using delta method
  betai <- which(names(mda$par)=="beta")
  # For first year
  mda_trmobs_se <- msm::deltamethod(as.formula(paste("~((",
                                                     paste("x",1:length(betai), sep="", collapse="+"),
                                                     ")/", length(betai), ")/log(10)", sep="")),
                                    mda$modelfit$par[betai],
                                    solve(mda$modelfit$hessian)[betai, betai])
  # For all other years
  mda_trmobs_se <- c(mda_trmobs_se,
                     sapply(betai,
                            function(i)
                              msm::deltamethod(as.formula(paste("~(x1 - (",
                                                                paste("x",1:length(betai), sep="", collapse="+"),
                                                                ")/", length(betai), ")/log(10)", sep="")),
                                               mda$modelfit$par[c(i, betai[-i])],
                                               solve(mda$modelfit$hessian)[c(i, betai[-i]), c(i, betai[-i])])))

  ci_mda[, TRMOBS_SE := mda_trmobs_se]

  # Abundance index outputs using two-stage GAI and bootstrap
  twostage_output <- readRDS("Data/UKBMS_twostageGAI_output_2sp.rds")[SPECIES == spp]



  index_plot <- ggplot(ci_mda, aes(YEAR, TRMOBS))+
    theme_bw()+
    geom_line()+geom_point()+
    geom_line(data = twostage_output, col="blue")+
    geom_point(data = twostage_output, col="blue")+
    geom_ribbon(aes(ymin = LOW_TRMOBS,
                    ymax = UPP_TRMOBS),  data = twostage_output,
                alpha = .2, fill = "blue")+
    geom_ribbon(aes(ymin = TRMOBS - 1.96*TRMOBS_SE*sqrt(mda$D),
                    ymax = TRMOBS + 1.96*TRMOBS_SE*sqrt(mda$D)), alpha = .4)+
    theme(text = element_text(size = 22),
          plot.margin=margin(1.5,.5,.25,.5,"cm"))+
    xlab("Year (r)")+
    ylab("Abundance index (log10)")+
    ggtitle("(ii)")

  gaiplots[[spp]]$beta_plot <- beta_plot
  gaiplots[[spp]]$index_plot <- index_plot

}




rssplot <- ggarrange(ggarrange(gaiplots[[splist[1]]]$beta_plot, gaiplots[[splist[1]]]$index_plot, nrow = 1),
                     ggarrange(gaiplots[[splist[2]]]$beta_plot, gaiplots[[splist[2]]]$index_plot, nrow = 1),
                     nrow = 2,
                     labels = c(paste("A -", splist[1]),paste("B -", splist[2])),
                     font.label = list(size = 24, color = "black", face = "bold", family = NULL))

ggsave(plot = rssplot, paste0("Figures/eGAI_RSS_plot.png"),
       width = 15, height = 13)
