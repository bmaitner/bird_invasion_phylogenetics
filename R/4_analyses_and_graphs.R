#4 analyses and graphs for ms

library(lme4)

intros <- readRDS("data/mean_phy_metrics_output.rds")
intros <- as.data.frame(intros)
intros[3:33] <- apply(X = intros[3:33],MARGIN = 2,FUN = function(x){as.numeric(as.character(x))})





##############################

#Is it worth including source range size and richness, or are they too correlated?

cor.test(intros$s_range_size,intros$s_richness)

#############################################################################
#Comparison of native range phylogenetic metrics
  #We compared the relative predictive ability of three native range phylogenetic metrics (PD, MPD, NND) by comparing the fits of generalized linear models. 
#The models contained the native region phylogenetic metric,native region range size as predictor variables, 
  #region of introduction (Hawaii, Florida, New Zealand) as a random effect and 
  #introduction success as a binary response (successful vs failed) variable with a logit link. 
#The predictor variable native range size was significantly correlated with PD (Pearson correlation correlation = 0.38, p < 0.05), but not MPD or NND (p > 0.05). 
#In our dataset, establishment success showed a weak, non-significant trend of positive association with PD of the native (source) range (Fig. 2, Supplementary Table S2). In contrast, establishment success showed significant, negative relationships with both native (source) range MPD and NND (Fig. 2, Table S2), with MPD showing a slightly better fit (1 < ΔAIC < 2; Table S2). 
#We focused all further analyses on these latter two metrics.
  # phy correction


source_pd <- glmer(failure ~ s_pd + (1|site_i))
intro

  glm(formula = )
  
  mpd_region_rs<-glmer(focal_success~focal_nat_mpd+(1|focal_region),data = rescaled_combined_nat_inv_data_output,family = "binomial")
  nnd_region_rs<-glmer(focal_success~focal_nat_nnd+(1|focal_region),data = rescaled_combined_nat_inv_data_output,family = "binomial")
  pd_region_rs<-glmer(focal_success~focal_nat_pd+(1|focal_region),data = rescaled_combined_nat_inv_data_output,family = "binomial")
  


#############################################################################
# Graph: (Fig 2)Comparison of phylogenetic metrics in tests of the Evolutionary Imbalance Hypothesis. Shown are estimated model coefficients and standard errors. Full models included native region phylogenetic metric (MPD, NND, PD) and the area of the focal species’ native range as predictors, region of introduction as a random effect and introduction success as a binary response variable. ** indicates p < 0.01.

# Graph(s) (fig 3, possible with other relative metrics (vpd,spd,kpd)) The effect of relative MPD on introduced success. Relative MPD was calculated as [Recipient Region MPD - Native Region MPD] / Native Region MPD. The recipient region here includes only native species within that region. Positive relative MPD values indicate species moving from relatively closely-related native communities into relatively distantly-related recipient communities, while negative values indicate the converse. Bin width is proportional to the number of observations in that bin.

# Figure 4/ model: Models of introduction success based on Mean Phylogenetic Distance (MPD) of the invader to community members. Models included the area of the species’ native range as a covariate, region of introduction as a random effect, and either 1) MPD to a species’ native source community, 2) MPD to the recipient community, or 3) MPD for both native and introduced regions as predictors. The native or recipient (invaded) community was circumscribed 3 ways: 1) only native species, 2) only successfully established introduced species, and 3) native species and successfully established introduced species. 

# Figure/ model: As above but NND

# Graph: Comparison of R2 values and partial R2 values for models incorporating propagule pressure. All models included establishment success as a binary response variable. Predictor variables included: (A) propagule number and size; (B) propagule number and size and recipient community MPD; and (C) propagule number and size, habitat generalism, brood value and brain size residuals. 