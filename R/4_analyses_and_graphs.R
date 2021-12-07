#4 analyses and graphs for ms

#Load required or useful libraries
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
library(lme4)
library(MuMIn)#model competition
library(bbmle)
library(stargazer)#table creation
library(arm)#like stargazer, use "display(glm.inv)"
library(glmulti)
library(leaps)
library(MASS)
library(lattice)
library(Hmisc)
library(usdm)
library(car)
library(rr2)
library(lmSupport)
library(ape)
library(phyr)

source("R/functions/replicate_models.R")

#Load intro data and format columns
intros <- readRDS("data/mean_phy_metrics_output.rds")
intros <- as.data.frame(intros)
intros[3:33] <- apply(X = intros[3:33],
                      MARGIN = 2,
                      FUN = function(x){as.numeric(as.character(x))})

#drop 2 species with source ranges so small they weren't rasterized
intros <- na.omit(intros)


#rescale predictors
intros[3:10] <- scale(intros[3:10])
intros[12:33] <- scale(intros[12:33])

#historically, for some reason, the column corresponding to whether a species was established (1) or not (0), was called "failures".  Let's fix that to avoid confusion
colnames(intros)[which(colnames(intros) == "r_failure")]<-"r_success"
intros$r_success <- as.integer(intros$r_success)

#save a csv version to share
write.csv(x = intros,file = "data/mean_phy_metrics_output.csv",row.names = F)
##############################

#Is it worth including source range size and richness, or are they too correlated?

cor.test(intros$s_range_size,intros$s_richness)  # cor = 0.66, so may be worth including both
cor.test(intros$s_pd,intros$s_richness) #cor = .99


#############################################################################

#Native range size vs. source phy metric
cor.test(x = intros$s_range_size,y = intros$s_pd)#0.66 , p <0.05
cor.test(x = intros$s_range_size,y = intros$s_nnd)#-0.09 , p =0.054
cor.test(x = intros$s_range_size,y = intros$s_mpd)#-0.01, p =0.91
cor.test(x = intros$s_range_size,y = intros$s_vpd)#0.04, p =0.39
cor.test(x = intros$s_range_size,y = intros$s_spd)#0.09, p =0.07
cor.test(x = intros$s_range_size,y = intros$s_kpd)#-0.02, p =0.74
colnames(intros)

#Native range richness vs. source phy metric
cor.test(x = intros$s_richness,y = intros$s_pd)#0.99 , p < 0.05 *
cor.test(x = intros$s_richness,y = intros$s_nnd)#-0.10 , p < 0.05 *
cor.test(x = intros$s_richness,y = intros$s_mpd)#-0.02, p > 0.05
cor.test(x = intros$s_richness,y = intros$s_vpd)#-0.12, p  < 0.05 *
cor.test(x = intros$s_richness,y = intros$s_spd)#0.09, p > 0.05
cor.test(x = intros$s_richness,y = intros$s_kpd)#0.00, p > 0.05


#calculate success vs metric
metric.pd <- glmer(formula = r_success ~ s_pd + s_range_size + s_richness + (1|site_i),
                   data = intros,
                   family = "binomial")
metric.pd_no_richness <- glmer(formula = r_success ~ s_pd + s_range_size + (1|site_i),
                               data = intros,
                               family = "binomial")
metric.nnd <- glmer(formula = r_success ~ s_nnd + s_range_size + s_richness + (1|site_i),
                    data = intros,
                    family = "binomial")
metric.mpd <- glmer(formula = r_success ~ s_mpd + s_range_size + s_richness + (1|site_i),
                    data = intros,
                    family = "binomial")
metric.vpd <- glmer(formula = r_success ~ s_vpd + s_range_size + s_richness + (1|site_i),
                    data = intros,
                    family = "binomial")
metric.spd <- glmer(formula = r_success ~ s_spd + s_range_size + s_richness + (1|site_i),
                    data = intros,
                    family = "binomial")
metric.kpd <- glmer(formula = r_success ~ s_kpd + s_range_size + s_richness + (1|site_i),
                    data = intros,
                    family = "binomial")

#look at models individually
summary(metric.pd_no_richness) #ns, +
summary(metric.nnd) # s, -
summary(metric.mpd) # s, -
summary(metric.vpd) # s, +
summary(metric.spd) # s, +
summary(metric.kpd) # s, -

selection.metric <- model.sel(metric.pd,metric.nnd,metric.mpd,metric.vpd,metric.spd,metric.kpd) #interestingly, variance seems to win

stargazer(metric.pd,metric.nnd,metric.mpd,metric.vpd,metric.spd,metric.kpd,type="html",out="figures_and_tables/source_range_metrics.htm")

# Code for calculating model with phy corrections

phylo.metric.pd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                       formula = r_success ~ s_pd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                       data = intros,
                       family = "binomial")

phylo.metric.pd_no_richness <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                        formula = r_success ~ s_pd + s_range_size + (1|site_i) + (1|species_i__),
                                        data = intros,
                                        family = "binomial")


phylo.metric.nnd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                         formula = r_success ~ s_nnd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                         data = intros,
                                         family = "binomial")


phylo.metric.mpd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                        formula = r_success ~ s_mpd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                        data = intros,
                                        family = "binomial")

phylo.metric.vpd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                         formula = r_success ~ s_vpd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                         data = intros,
                                         family = "binomial")

phylo.metric.spd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                         formula = r_success ~ s_spd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                         data = intros,
                                         family = "binomial")

phylo.metric.kpd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                         formula = r_success ~ s_kpd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                         data = intros,
                                         family = "binomial")

#pull out mean phy corrections parameters

phylo.metric.pd.mean <- summarize_replicates(phylo_reps_output = phylo.metric.pd)
phylo.metric.pd_no_richness.mean <- summarize_replicates(phylo_reps_output = phylo.metric.pd_no_richness)
phylo.metric.nnd.mean <- summarize_replicates(phylo_reps_output = phylo.metric.nnd)
phylo.metric.mpd.mean <- summarize_replicates(phylo_reps_output = phylo.metric.mpd)
phylo.metric.vpd.mean <- summarize_replicates(phylo_reps_output = phylo.metric.vpd)
phylo.metric.spd.mean <- summarize_replicates(phylo_reps_output = phylo.metric.spd)
phylo.metric.kpd.mean <- summarize_replicates(phylo_reps_output = phylo.metric.kpd)

# Cut the replicated stuff to save memory if need be
  # rm(phylo.metric.pd,
  #    phylo.metric.pd_no_richness,
  #    phylo.metric.nnd,
  #    phylo.metric.mpd,
  #    phylo.metric.vpd,
  #    phylo.metric.spd,
  #    phylo.metric.kpd)

estimates_df<-as.data.frame(matrix(nrow = 14,ncol = 6))
colnames(estimates_df)<-c("metric","estimate","se","p_value","phy_corrected","AIC")

estimates_df[1,1]<- "PD_w_richness"
estimates_df[1,2:4]<-phylo.metric.pd.mean$means[2,]
estimates_df[1,5]<- "Y"
estimates_df$AIC[1] <- phylo.metric.pd.mean$aic_mean

estimates_df[2,1]<- "PD"
estimates_df[2,2:4]<-phylo.metric.pd_no_richness.mean$means[2,]
estimates_df[2,5]<- "Y"
estimates_df$AIC[2] <- phylo.metric.pd_no_richness.mean$aic_mean

estimates_df[3,1]<- "NND"
estimates_df[3,2:4]<-phylo.metric.nnd.mean$means[2,]
estimates_df[3,5]<- "Y"
estimates_df$AIC[3] <- phylo.metric.nnd.mean$aic_mean

estimates_df[4,1]<- "MPD"
estimates_df[4,2:4]<-phylo.metric.mpd.mean$means[2,]
estimates_df[4,5]<- "Y"
estimates_df$AIC[4] <- phylo.metric.mpd.mean$aic_mean


estimates_df[5,1]<- "VPD"
estimates_df[5,2:4]<-phylo.metric.vpd.mean$means[2,]
estimates_df[5,5]<- "Y"
estimates_df$AIC[5] <- phylo.metric.vpd.mean$aic_mean

estimates_df[6,1]<- "SPD"
estimates_df[6,2:4]<-phylo.metric.spd.mean$means[2,]
estimates_df[6,5]<- "Y"
estimates_df$AIC[6] <- phylo.metric.spd.mean$aic_mean

estimates_df[7,1]<- "KPD"
estimates_df[7,2:4]<-phylo.metric.kpd.mean$means[2,]
estimates_df[7,5]<- "Y"
estimates_df$AIC[7] <- phylo.metric.kpd.mean$aic_mean



estimates_df[8,1]<- "PD_w_richness"
estimates_df[8,2:4]<-summary(metric.pd)$coefficients[2,c(1,2,4)]
estimates_df[8,5]<- "N"
estimates_df$AIC[8] <- summary(metric.pd)$AIC[1]

estimates_df[9,1]<- "PD"
estimates_df[9,2:4]<-summary(metric.pd_no_richness)$coefficients[2,c(1,2,4)]
estimates_df[9,5]<- "N"
estimates_df$AIC[9] <- summary(metric.pd_no_richness)$AIC[1]

estimates_df[10,1]<- "NND"
estimates_df[10,2:4]<-summary(metric.nnd)$coefficients[2,c(1,2,4)]
estimates_df[10,5]<- "N"
estimates_df$AIC[10] <- summary(metric.nnd)$AIC[1]

estimates_df[11,1]<- "MPD"
estimates_df[11,2:4]<-summary(metric.mpd)$coefficients[2,c(1,2,4)]
estimates_df[11,5]<- "N"
estimates_df$AIC[11] <- summary(metric.mpd)$AIC[1]

estimates_df[12,1]<- "VPD"
estimates_df[12,2:4]<-summary(metric.vpd)$coefficients[2,c(1,2,4)]
estimates_df[12,5]<- "N"
estimates_df$AIC[12] <- summary(metric.vpd)$AIC[1]

estimates_df[13,1]<- "SPD"
estimates_df[13,2:4]<-summary(metric.spd)$coefficients[2,c(1,2,4)]
estimates_df[13,5]<- "N"
estimates_df$AIC[13] <- summary(metric.spd)$AIC[1]

estimates_df[14,1]<- "KPD"
estimates_df[14,2:4]<-summary(metric.kpd)$coefficients[2,c(1,2,4)]
estimates_df[14,5]<- "N"
estimates_df$AIC[14] <- summary(metric.kpd)$AIC[1]


estimates_df$significant<-"No"
estimates_df$significant[which(estimates_df$p_value<0.05)]<-"Yes"


estimates_df$phy_corrected <- gsub(pattern = "Y",
                                   replacement = "Yes",
                                   x = estimates_df$phy_corrected)
estimates_df$phy_corrected <- gsub(pattern = "N",
                                   replacement = "No",
                                   x = estimates_df$phy_corrected)

estimates_df <- estimates_df[which(estimates_df$metric != "PD_w_richness"),]
estimates_df$metric <- factor(estimates_df$metric,
                              levels = c("PD","NND", "MPD", "VPD", "SPD", "KPD" ))


library(tidyverse)
estimates_df %>%
  mutate(
    hypothesis_supported = case_when(
      metric == "PD" & estimate > 0 ~ "EIH",
      metric == "NND" & estimate < 0 ~ "EIH",
      metric == "MPD" & estimate < 0 ~ "EIH",
      metric == "VPD" & estimate > 0 ~ "EIH",
      metric == "SPD" & estimate  > 0 ~ "EIH",
      metric == "KPD" & estimate < 0 ~ "EIH",
      
    )) -> estimates_df

#Make dummy stuff so that CCH is shows up on figure legend
estimates_df <- rbind(estimates_df,
                      replicate(n = ncol(estimates_df),
                                expr = NA))
estimates_df$hypothesis_supported[13] <- "CCH"
estimates_df$phy_corrected[13]<-"No"
estimates_df$significant[13] <- "No"
estimates_df$metric[13] <- "PD"


hyp.colors <- c(EIH = "blue", CCH = "red")
corr.patterns <- c(Yes = NA, No= "stripe")

source_metric_figure <-
ggplot(estimates_df, aes(x = metric,
                         y = estimate,
                         fill = hypothesis_supported,
                         alpha = significant,
                         color = hypothesis_supported,
                         pattern_ = phy_corrected))+
  geom_bar_pattern(mapping = aes(pattern = phy_corrected),
                   stat="identity",
                   position = position_dodge(),
                   pattern_density=0.5,
                   pattern_color = "black",
                   pattern_alpha=1,
                   pattern_fill = "white")+
  scale_fill_manual(values = hyp.colors)+
  scale_color_manual(values=hyp.colors)+
  scale_pattern_manual(values = corr.patterns)+
  geom_errorbar(aes(ymin = estimate-se,
                    ymax = estimate+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9),
                color="black",
                alpha=1,
                size=1
  )+
  theme(text = element_text(size=20))+
  theme_bw(base_size = 30)+
  labs(fill="Hypothesis \n supported",
       alpha="p < 0.05?",
       pattern="Phylogenetic \n correction?")+
  ylab("model coefficient")+
  guides(color=FALSE)


# ggsave(filename = "figures_and_tables/source_metric_figure.pdf",
#        plot = source_metric_figure,
#        width = 14,
#        height = 10)
# 
# ggsave(filename = "figures_and_tables/source_metric_figure.jpg",
#        plot = source_metric_figure,
#        width = 14,
#        height = 10)

#Get AIC estimates  
aictable_source_metrics<-estimates_df[order(estimates_df$AIC,decreasing = F),]

write.csv(x = aictable_source_metrics,file = "figures_and_tables/source_metric_aic_table.csv",row.names = F)

estimates_df<-aictable_source_metrics

#############################################################################




# Graph(s) (fig 3, possible with other relative metrics (vpd,spd,kpd)) 
#The effect of relative MPD on introduced success. 
#Relative MPD was calculated as [Recipient Region MPD - Native Region MPD] / Native Region MPD. 
#The recipient region here includes only native species within that region. 
#Positive relative MPD values indicate species moving from relatively closely-related native communities into relatively distantly-related recipient communities, 
#while negative values indicate the converse. Bin width is proportional to the number of observations in that bin.


intros$delta_nnd_s_r_n <- (intros$r_n_nnd  - intros$s_nnd)
intros$delta_mpd_s_r_n <- (intros$r_n_mpd  - intros$s_mpd)
intros$delta_vpd_s_r_n <- (intros$r_n_vpd  - intros$s_vpd)
intros$delta_spd_s_r_n <- (intros$r_n_spd  - intros$s_spd)
intros$delta_kpd_s_r_n <- (intros$r_n_kpd  - intros$s_kpd)

delta.nnd <- glmer(formula = r_success ~ r_n_nnd + s_range_size + s_richness + (1|site_i),
                   data = intros,
                   family = "binomial")

delta.mpd <- glmer(formula = r_success ~ r_n_mpd + s_range_size + s_richness + (1|site_i),
                   data = intros,
                   family = "binomial")

delta.vpd <- glmer(formula = r_success ~ r_n_vpd + s_range_size + s_richness + (1|site_i),
                   data = intros,
                   family = "binomial")

delta.spd <- glmer(formula = r_success ~ r_n_spd + s_range_size + s_richness + (1|site_i),
                   data = intros,
                   family = "binomial")

delta.kpd <- glmer(formula = r_success ~ r_n_kpd + s_range_size + s_richness + (1|site_i),
                   data =intros,
                   family = "binomial")

delta.nnd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:10],
                                        formula = r_success ~ delta_nnd_s_r_n + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                        data = intros,
                                        family = "binomial")

delta.mpd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:10],
                                  formula = r_success ~ delta_mpd_s_r_n + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                  data = intros,
                                  family = "binomial")

delta.vpd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:10],
                                  formula = r_success ~ delta_vpd_s_r_n + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                  data = intros,
                                  family = "binomial")

delta.spd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:10],
                                  formula = r_success ~ delta_spd_s_r_n + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                  data = intros,
                                  family = "binomial")

delta.kpd <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:10],
                                  formula = r_success ~ delta_kpd_s_r_n + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                  data = intros,
                                  family = "binomial")


selection.delta <- model.sel(delta.nnd,delta.mpd,delta.vpd,delta.spd,delta.kpd) #interestingly, variance seems to win
#models all equally useless

#none of the deltas are significant
stargazer(delta.nnd,delta.mpd,delta.vpd,delta.spd,delta.kpd,type="html",out="figures_and_tables/delta_range_metrics.htm")

intros$word_success <- "Failure"
intros$word_success[which(intros$r_success==1)]<-"Established"

plot(as.factor(intros$word_success)~intros$delta_nnd_s_r_n,xlab = "NND(recipient) - NND(source)",ylab = "Establishment Success")
plot(as.factor(intros$word_success)~intros$delta_mpd_s_r_n,xlab = "MPD(recipient) - MPD(source)",ylab = "Establishment Success")
plot(as.factor(intros$word_success)~intros$delta_vpd_s_r_n,xlab = "VPD(recipient) - VPD(source)",ylab = "Establishment Success")
plot(as.factor(intros$word_success)~intros$delta_spd_s_r_n,xlab = "SPD(recipient) - SPD(source)",ylab = "Establishment Success")
plot(as.factor(intros$word_success)~intros$delta_kpd_s_r_n,xlab = "KPD(recipient) - KPD(source)",ylab = "Establishment Success")

###########################################################################

#NND and MPD (phylo/not phylo)

#NND
#natives only in recipient
full_nnd_s <-     glmer(formula = r_success ~ s_nnd + s_range_size + s_richness + (1|site_i),
                        data = intros,
                        family = "binomial")

full_nnd_s_r_n <- glmer(formula = r_success ~ s_nnd + s_range_size + s_richness + r_n_nnd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_nnd_s_r_n_int <- glmer(formula = r_success ~ s_nnd * r_n_nnd + s_nnd + s_range_size + s_richness + r_n_nnd  +(1|site_i),
                            data = intros,
                            family = "binomial")

full_nnd_r_n <- glmer(formula = r_success ~  s_range_size + s_richness + r_n_nnd  + (1|site_i),
                      data = intros,
                      family = "binomial")

#established only in recipient
full_nnd_s #already done

full_nnd_s_r_e <- glmer(formula = r_success ~ s_nnd + s_range_size + s_richness + r_e_nnd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_nnd_s_r_e_int <- glmer(formula = r_success ~ s_nnd * r_e_nnd + s_nnd + s_range_size + s_richness + r_e_nnd  +(1|site_i),
                            data = intros,
                            family = "binomial")

full_nnd_r_e <- glmer(formula = r_success ~  s_range_size + s_richness + r_e_nnd  + (1|site_i),
                      data = intros,
                      family = "binomial")

#established + native  in recipient
full_nnd_s #already done

full_nnd_s_r_ne <- glmer(formula = r_success ~ s_nnd + s_range_size + s_richness + r_ne_nnd  +(1|site_i),
                         data = intros,
                         family = "binomial")

full_nnd_s_r_ne_int <- glmer(formula = r_success ~ s_nnd * r_ne_nnd + s_nnd + s_range_size + s_richness + r_ne_nnd  +(1|site_i),
                             data = intros,
                             family = "binomial")

full_nnd_r_ne <- glmer(formula = r_success ~  s_range_size + s_richness + r_ne_nnd  + (1|site_i),
                       data = intros,
                       family = "binomial")

summary(full_nnd_s)

summary(full_nnd_r_n)
summary(full_nnd_r_e)
summary(full_nnd_r_ne)

summary(full_nnd_s_r_n)
summary(full_nnd_s_r_e)
summary(full_nnd_s_r_ne)
summary(full_nnd_s_r_ne_int)


#phy corrected

#natives only in recipeint
phylo.full_nnd_s <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                  formula = r_success ~ s_nnd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                  data = intros,
                                  family = "binomial")

phylo.full_nnd_s_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                         formula = r_success ~ s_nnd + s_range_size + s_richness + r_n_nnd  +(1|site_i) + (1|species_i__),
                                         data = intros,
                                         family = "binomial")

phylo.full_nnd_s_r_n_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_nnd*r_n_nnd + s_nnd + s_range_size + s_richness + r_n_nnd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")


phylo.full_nnd_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~  s_range_size + s_richness + r_n_nnd  + (1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")




#established only in recipient
phylo.full_nnd_s #done

phylo.full_nnd_s_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_nnd + s_range_size + s_richness + r_e_nnd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_nnd_s_r_e_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_nnd*r_e_nnd + s_nnd + s_range_size + s_richness + r_e_nnd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_nnd_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_e_nnd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")


#native and established in recipient
phylo.full_nnd_s #done
phylo.full_nnd_s_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_nnd + s_range_size + s_richness + r_ne_nnd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_nnd_s_r_ne_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                              formula = r_success ~ s_nnd * r_ne_nnd +s_nnd + s_range_size + s_richness + r_ne_nnd  +(1|site_i) + (1|species_i__),
                                              data = intros,
                                              family = "binomial")



phylo.full_nnd_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_ne_nnd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")


#natives only in recipient

phylo.full_nnd_s.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_s)#507.2
phylo.full_nnd_s_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_s_r_n)#508.6
phylo.full_nnd_s_r_n_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_s_r_n_int)#
phylo.full_nnd_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_r_n)#510.7

#established only in recipient

phylo.full_nnd_s.mean # done
phylo.full_nnd_s_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_s_r_e) #508.3
phylo.full_nnd_s_r_e_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_s_r_e_int) #
phylo.full_nnd_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_r_e)#510.3

#native and established in recipient

phylo.full_nnd_s.mean # done
phylo.full_nnd_s_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_s_r_ne)#508.7
phylo.full_nnd_s_r_ne_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_s_r_ne_int)#
phylo.full_nnd_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_nnd_r_ne)#510.4

######################################
#MPD

#natives only in recipeint
full_mpd_s <-     glmer(formula = r_success ~ s_mpd + s_range_size + s_richness + (1|site_i),
                        data = intros,
                        family = "binomial")

full_mpd_s_r_n <- glmer(formula = r_success ~ s_mpd + s_range_size + s_richness + r_n_mpd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_mpd_s_r_n_int <- glmer(formula = r_success ~ s_mpd*r_n_mpd + s_mpd + s_range_size + s_richness + r_n_mpd  +(1|site_i),
                            data = intros,
                            family = "binomial")

full_mpd_r_n <- glmer(formula = r_success ~  s_range_size + s_richness + r_n_mpd  + (1|site_i),
                      data = intros,
                      family = "binomial")

#established only in recipeint
full_mpd_s #already done

full_mpd_s_r_e <- glmer(formula = r_success ~ s_mpd + s_range_size + s_richness + r_e_mpd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_mpd_s_r_e_int <- glmer(formula = r_success ~ s_mpd*r_e_mpd + s_mpd + s_range_size + s_richness + r_e_mpd  +(1|site_i),
                            data = intros,
                            family = "binomial")

full_mpd_r_e <- glmer(formula = r_success ~  s_range_size + s_richness + r_e_mpd  + (1|site_i),
                      data = intros,
                      family = "binomial")

#established + native  in recipeint
full_mpd_s #already done

full_mpd_s_r_ne <- glmer(formula = r_success ~ s_mpd + s_range_size + s_richness + r_ne_mpd  + (1|site_i),
                         data = intros,
                         family = "binomial")

full_mpd_s_r_ne_int <- glmer(formula = r_success ~ s_mpd*r_ne_mpd + s_mpd + s_range_size + s_richness + r_ne_mpd  +(1|site_i),
                             data = intros,
                             family = "binomial")

full_mpd_r_ne <- glmer(formula = r_success ~  s_range_size + s_richness + r_ne_mpd  + (1|site_i),
                       data = intros,
                       family = "binomial")

summary(full_mpd_s)

summary(full_mpd_r_n)
summary(full_mpd_r_e)
summary(full_mpd_r_ne)

summary(full_mpd_s_r_n)
summary(full_mpd_s_r_e)
summary(full_mpd_s_r_ne)



#phy corrected

#natives only in recipeint
phylo.full_mpd_s <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                         formula = r_success ~ s_mpd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                         data = intros,
                                         family = "binomial")

phylo.full_mpd_s_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_mpd + s_range_size + s_richness + r_n_mpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_mpd_s_r_n_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_mpd*r_n_mpd +s_mpd + s_range_size + s_richness + r_n_mpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")


phylo.full_mpd_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_n_mpd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")




#established only in recipient
phylo.full_mpd_s #done

phylo.full_mpd_s_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_mpd + s_range_size + s_richness + r_e_mpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_mpd_s_r_e_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_mpd*r_e_mpd + s_mpd + s_range_size + s_richness + r_e_mpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_mpd_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_e_mpd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")


#native and established in recipient
phylo.full_mpd_s #done
phylo.full_mpd_s_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                              formula = r_success ~ s_mpd + s_range_size + s_richness + r_ne_mpd  +(1|site_i) + (1|species_i__),
                                              data = intros,
                                              family = "binomial")

phylo.full_mpd_s_r_ne_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                              formula = r_success ~ s_mpd*r_ne_mpd + s_mpd + s_range_size + s_richness + r_ne_mpd  +(1|site_i) + (1|species_i__),
                                              data = intros,
                                              family = "binomial")


phylo.full_mpd_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                            formula = r_success ~  s_range_size + s_richness + r_ne_mpd  + (1|site_i) + (1|species_i__),
                                            data = intros,
                                            family = "binomial")


#natives only in recipient

phylo.full_mpd_s.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_s)#506.4
phylo.full_mpd_s_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_s_r_n) #490.5
phylo.full_mpd_s_r_n_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_s_r_n_int) #
phylo.full_mpd_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_r_n)#495.3

#established only in recipient

phylo.full_mpd_s.mean # done
phylo.full_mpd_s_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_s_r_e)#507.9
phylo.full_mpd_s_r_e_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_s_r_e_int)#
phylo.full_mpd_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_r_e)#510.1

#native and established in recipient

phylo.full_mpd_s.mean # done
phylo.full_mpd_s_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_s_r_ne)#504.5
phylo.full_mpd_s_r_ne_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_s_r_ne_int)#
phylo.full_mpd_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_mpd_r_ne)#510.2


#Cut full models to save space
# rm(phylo.full_mpd_r_e,
#    phylo.full_mpd_r_n,
#    phylo.full_mpd_r_ne,
#    phylo.full_mpd_s,
#    phylo.full_mpd_s_r_n,
#    phylo.full_mpd_s_r_e,
#    phylo.full_mpd_s_r_ne,
#    phylo.full_mpd_s_r_n_int,
#    phylo.full_mpd_s_r_e_int,
#    phylo.full_mpd_s_r_ne_int)



# Figure 4/ model: Models of introduction success based on Mean Phylogenetic Distance (MPD) of the invader to community members. 
# Models included the area of the species’ native range as a covariate, region of introduction as a random effect, and either 
  #1) MPD to a species’ native source community, 
  #2) MPD to the recipient community, or 
  #3) MPD for both native and introduced regions as predictors.
#The native or recipient (invaded) community was circumscribed 3 ways: 
  #1) only native species, 
  #2) only successfully established introduced species, and 
  #3) native species and successfully established introduced species. 

# Figure/ model: As above but NND
############################################################################

#Overall (NND,MPD,VPD,SPD,KPD) 
  #each has 7 versions:
    # source v recipient v source + recipient
    # native v establshed v native +established

#natives only in recipeint
full_vpd_s <-     glmer(formula = r_success ~ s_vpd + s_range_size + s_richness + (1|site_i),
                        data = intros,
                        family = "binomial")

full_vpd_s_r_n <- glmer(formula = r_success ~ s_vpd + s_range_size + s_richness + r_n_vpd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_vpd_s_r_n_int <- glmer(formula = r_success ~ s_vpd*r_n_vpd + s_vpd + s_range_size + s_richness + r_n_vpd  +(1|site_i),
                            data = intros,
                            family = "binomial")

full_vpd_r_n <- glmer(formula = r_success ~  s_range_size + s_richness + r_n_vpd  + (1|site_i),data=intros,family = "binomial")

#established only in recipeint
full_vpd_s #already done

full_vpd_s_r_e <- glmer(formula = r_success ~ s_vpd + s_range_size + s_richness + r_e_vpd  +(1|site_i),
                        data=intros,
                        family = "binomial")

full_vpd_s_r_e_int <- glmer(formula = r_success ~ s_vpd*r_e_vpd + s_vpd + s_range_size + s_richness + r_e_vpd  +(1|site_i),
                        data=intros,
                        family = "binomial")

full_vpd_r_e <- glmer(formula = r_success ~  s_range_size + s_richness + r_e_vpd  + (1|site_i),
                      data=intros,
                      family = "binomial")

#established + native  in recipeint
full_vpd_s #already done
full_vpd_s_r_ne <- glmer(formula = r_success ~ s_vpd + s_range_size + s_richness + r_ne_vpd  +(1|site_i),
                         data = intros,
                         family = "binomial")

full_vpd_s_r_ne_int <- glmer(formula = r_success ~ s_vpd*r_ne_vpd + s_vpd + s_range_size + s_richness + r_ne_vpd  +(1|site_i),
                         data = intros,
                         family = "binomial")

full_vpd_r_ne <- glmer(formula = r_success ~  s_range_size + s_richness + r_ne_vpd  + (1|site_i),
                       data = intros,
                       family = "binomial")


#VPD

#phy corrected

#natives only in recipeint
phylo.full_vpd_s <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                         formula = r_success ~ s_vpd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                         data = intros,
                                         family = "binomial")

phylo.full_vpd_s_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_vpd + s_range_size + s_richness + r_n_vpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_vpd_s_r_n_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_vpd*r_n_vpd + s_vpd + s_range_size + s_richness + r_n_vpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")


phylo.full_vpd_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_n_vpd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")




#established only in recipient
phylo.full_vpd_s #done

phylo.full_vpd_s_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_vpd + s_range_size + s_richness + r_e_vpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_vpd_s_r_e_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_vpd*r_e_vpd + s_vpd + s_range_size + s_richness + r_e_vpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_vpd_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_e_vpd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")


#native and established in recipient
phylo.full_vpd_s #done
phylo.full_vpd_s_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                              formula = r_success ~ s_vpd + s_range_size + s_richness + r_ne_vpd  +(1|site_i) + (1|species_i__),
                                              data = intros,
                                              family = "binomial")

phylo.full_vpd_s_r_ne_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                              formula = r_success ~ s_vpd*r_ne_vpd + s_vpd + s_range_size + s_richness + r_ne_vpd  +(1|site_i) + (1|species_i__),
                                              data = intros,
                                              family = "binomial")


phylo.full_vpd_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                            formula = r_success ~  s_range_size + s_richness + r_ne_vpd  + (1|site_i) + (1|species_i__),
                                            data = intros,
                                            family = "binomial")


#natives only in recipient

phylo.full_vpd_s.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_s)#
phylo.full_vpd_s_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_s_r_n) #
phylo.full_vpd_s_r_n_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_s_r_n_int) #
phylo.full_vpd_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_r_n)#
#established only in recipient

phylo.full_vpd_s.mean # done
phylo.full_vpd_s_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_s_r_e)#
phylo.full_vpd_s_r_e_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_s_r_e_int)#
phylo.full_vpd_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_r_e)#

#native and established in recipient

phylo.full_vpd_s.mean # done
phylo.full_vpd_s_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_s_r_ne)#
phylo.full_vpd_s_r_ne_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_s_r_ne_int)#
phylo.full_vpd_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_vpd_r_ne)#

#Cut full models to save space
# rm(phylo.full_vpd_r_e,
#    phylo.full_vpd_r_n,
#    phylo.full_vpd_r_ne,
#    phylo.full_vpd_s,
#    phylo.full_vpd_s_r_n,
#    phylo.full_vpd_s_r_e,
#    phylo.full_vpd_s_r_ne,
#    phylo.full_vpd_s_r_n_int,
#    phylo.full_vpd_s_r_e_int,
#    phylo.full_vpd_s_r_ne_int)


#SPD

#natives only in recipeint
full_spd_s <-     glmer(formula = r_success ~ s_spd + s_range_size + s_richness + (1|site_i),
                        data = intros,
                        family = "binomial")

full_spd_s_r_n <- glmer(formula = r_success ~ s_spd + s_range_size + s_richness + r_n_spd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_spd_s_r_n_int <- glmer(formula = r_success ~ s_spd*r_n_spd + s_spd + s_range_size + s_richness + r_n_spd  +(1|site_i),
                            data = intros,
                            family = "binomial")

full_spd_r_n <- glmer(formula = r_success ~  s_range_size + s_richness + r_n_spd  + (1|site_i),
                      data = intros,
                      family = "binomial")

#established only in recipeint
full_spd_s #already done

full_spd_s_r_e <- glmer(formula = r_success ~ s_spd + s_range_size + s_richness + r_e_spd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_spd_s_r_e_int <- glmer(formula = r_success ~ s_spd*r_e_spd + s_spd + s_range_size + s_richness + r_e_spd  +(1|site_i),
                            data = intros,
                            family = "binomial")

full_spd_r_e <- glmer(formula = r_success ~  s_range_size + s_richness + r_e_spd  + (1|site_i),
                      data = intros,
                      family = "binomial")

#established + native  in recipeint
full_spd_s #already done

full_spd_s_r_ne <- glmer(formula = r_success ~ s_spd + s_range_size + s_richness + r_ne_spd  +(1|site_i),
                         data = intros,
                         family = "binomial")

full_spd_s_r_ne_int <- glmer(formula = r_success ~ s_spd*r_ne_spd + s_spd + s_range_size + s_richness + r_ne_spd  +(1|site_i),
                             data = intros,
                             family = "binomial")

full_spd_r_ne <- glmer(formula = r_success ~  s_range_size + s_richness + r_ne_spd  + (1|site_i),
                       data = intros,
                       family = "binomial")

#phy corrected

#natives only in recipeint
phylo.full_spd_s <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                         formula = r_success ~ s_spd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                         data = intros,
                                         family = "binomial")

phylo.full_spd_s_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_spd + s_range_size + s_richness + r_n_spd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_spd_s_r_n_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_spd*r_n_spd + s_spd + s_range_size + s_richness + r_n_spd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_spd_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_n_spd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")


#established only in recipient
phylo.full_spd_s #done

phylo.full_spd_s_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_spd + s_range_size + s_richness + r_e_spd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_spd_s_r_e_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_spd*r_e_spd + s_spd + s_range_size + s_richness + r_e_spd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_spd_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_e_spd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")


#native and established in recipient
phylo.full_spd_s #done
phylo.full_spd_s_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                              formula = r_success ~ s_spd + s_range_size + s_richness + r_ne_spd  +(1|site_i) + (1|species_i__),
                                              data = intros,
                                              family = "binomial")

phylo.full_spd_s_r_ne_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                              formula = r_success ~ s_spd*r_ne_spd + s_spd + s_range_size + s_richness + r_ne_spd  +(1|site_i) + (1|species_i__),
                                              data = intros,
                                              family = "binomial")

phylo.full_spd_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                            formula = r_success ~  s_range_size + s_richness + r_ne_spd  + (1|site_i) + (1|species_i__),
                                            data = intros,
                                            family = "binomial")


#natives only in recipient

phylo.full_spd_s.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_s)#
phylo.full_spd_s_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_s_r_n) #
phylo.full_spd_s_r_n_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_s_r_n_int) #
phylo.full_spd_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_r_n)#
#established only in recipient

phylo.full_spd_s.mean # done
phylo.full_spd_s_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_s_r_e)#
phylo.full_spd_s_r_e_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_s_r_e_int)
phylo.full_spd_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_r_e)#

#native and established in recipient

phylo.full_spd_s.mean # done
phylo.full_spd_s_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_s_r_ne)#
phylo.full_spd_s_r_ne_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_s_r_ne_int)#
phylo.full_spd_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_spd_r_ne)#


#Cut full models to save space
# rm(phylo.full_spd_r_e,
#    phylo.full_spd_r_n,
#    phylo.full_spd_r_ne,
#    phylo.full_spd_s,
#    phylo.full_spd_s_r_n,
#    phylo.full_spd_s_r_e,
#    phylo.full_spd_s_r_ne,
#    phylo.full_spd_s_r_n_int,
#    phylo.full_spd_s_r_e_int,
#    phylo.full_spd_s_r_ne_int)



#KPD

#natives only in recipient
full_kpd_s <-     glmer(formula = r_success ~ s_kpd + s_range_size + s_richness + (1|site_i),
                        data = intros,
                        family = "binomial")

full_kpd_s_r_n <- glmer(formula = r_success ~ s_kpd + s_range_size + s_richness + r_n_kpd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_kpd_s_r_n_int <- glmer(formula = r_success ~ s_kpd*r_n_kpd + s_kpd + s_range_size + s_richness + r_n_kpd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_kpd_r_n <- glmer(formula = r_success ~  s_range_size + s_richness + r_n_kpd  + (1|site_i),
                      data = intros,
                      family = "binomial")

#established only in recipeint
full_kpd_s #already done

full_kpd_s_r_e <- glmer(formula = r_success ~ s_kpd + s_range_size + s_richness + r_e_kpd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_kpd_s_r_e_int <- glmer(formula = r_success ~ s_kpd*r_e_kpd + s_kpd + s_range_size + s_richness + r_e_kpd  +(1|site_i),
                        data = intros,
                        family = "binomial")

full_kpd_r_e <- glmer(formula = r_success ~  s_range_size + s_richness + r_e_kpd  + (1|site_i),data=intros,family = "binomial")

#established + native  in recipeint
full_kpd_s #already done

full_kpd_s_r_ne <- glmer(formula = r_success ~ s_kpd + s_range_size + s_richness + r_ne_kpd  +(1|site_i),
                         data = intros,
                         family = "binomial")

full_kpd_s_r_ne_int <- glmer(formula = r_success ~ s_kpd*r_ne_kpd + s_kpd + s_range_size + s_richness + r_ne_kpd  +(1|site_i),
                         data = intros,
                         family = "binomial")

full_kpd_r_ne <- glmer(formula = r_success ~  s_range_size + s_richness + r_ne_kpd  + (1|site_i),
                       data = intros,
                       family = "binomial")


#phy corrected

#natives only in recipeint
phylo.full_kpd_s <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                         formula = r_success ~ s_kpd + s_range_size + s_richness + (1|site_i) + (1|species_i__),
                                         data = intros,
                                         family = "binomial")

phylo.full_kpd_s_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_kpd + s_range_size + s_richness + r_n_kpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_kpd_s_r_n_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_kpd*r_n_kpd + s_kpd + s_range_size + s_richness + r_n_kpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_kpd_r_n <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_n_kpd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")


#established only in recipient
phylo.full_kpd_s #done

phylo.full_kpd_s_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_kpd + s_range_size + s_richness + r_e_kpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_kpd_s_r_e_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                             formula = r_success ~ s_kpd*r_e_kpd + s_kpd + s_range_size + s_richness + r_e_kpd  +(1|site_i) + (1|species_i__),
                                             data = intros,
                                             family = "binomial")

phylo.full_kpd_r_e <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~  s_range_size + s_richness + r_e_kpd  + (1|site_i) + (1|species_i__),
                                           data = intros,
                                           family = "binomial")


#native and established in recipient
phylo.full_kpd_s #done

phylo.full_kpd_s_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                              formula = r_success ~ s_kpd + s_range_size + s_richness + r_ne_kpd  +(1|site_i) + (1|species_i__),
                                              data = intros,
                                              family = "binomial")

phylo.full_kpd_s_r_ne_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                              formula = r_success ~ s_kpd*r_ne_kpd + s_kpd + s_range_size + s_richness + r_ne_kpd  +(1|site_i) + (1|species_i__),
                                              data = intros,
                                              family = "binomial")

phylo.full_kpd_r_ne <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                            formula = r_success ~  s_range_size + s_richness + r_ne_kpd  + (1|site_i) + (1|species_i__),
                                            data = intros,
                                            family = "binomial")


#natives only in recipient

phylo.full_kpd_s.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_s)#
phylo.full_kpd_s_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_s_r_n) #
phylo.full_kpd_s_r_n_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_s_r_n_int) #
phylo.full_kpd_r_n.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_r_n)#
#established only in recipient

phylo.full_kpd_s.mean # done
phylo.full_kpd_s_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_s_r_e)#
phylo.full_kpd_s_r_e_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_s_r_e_int)#
phylo.full_kpd_r_e.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_r_e)#

#native and established in recipient

phylo.full_kpd_s.mean # done
phylo.full_kpd_s_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_s_r_ne)#
phylo.full_kpd_s_r_ne_int.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_s_r_ne_int)#
phylo.full_kpd_r_ne.mean <- summarize_replicates(phylo_reps_output = phylo.full_kpd_r_ne)#

#Cut full models to save space
# rm(phylo.full_kpd_r_e,
#    phylo.full_kpd_r_n,
#    phylo.full_kpd_r_ne,
#    phylo.full_kpd_s,
#    phylo.full_kpd_s_r_n,
#    phylo.full_kpd_s_r_e,
#    phylo.full_kpd_s_r_ne,
#    phylo.full_kpd_s_r_n_int,
#    phylo.full_kpd_s_r_e_int,
#    phylo.full_kpd_s_r_ne_int)


#############################################################################
#############################################################################

#summarize phylogenetic models
source("R/functions/summarize_phylo_models.R")
phylo_model_summary <- summarize_phyl0_models(digits_to_round = 3)
nonphylo_model_summary <- summarize_nonphylo_models(digits_to_round = 3)
setdiff(colnames(phylo_model_summary),colnames(nonphylo_model_summary))
setdiff(colnames(nonphylo_model_summary),colnames(phylo_model_summary))
write.csv(x = phylo_model_summary,file = "figures_and_tables/phylo_model_summary.csv")
write.csv(x = nonphylo_model_summary,file = "figures_and_tables/nonphylo_model_summary.csv")

#############################################################################
#############################################################################

#Model selection, non phylogenetically corrected


#Load Sol data

sol<-read.csv("Data/1221523Database1.csv",stringsAsFactors = F)
sol_nz<-subset(sol,subset = sol$Location_of_introduction == "New_Zealand")
us_nz<-subset(intros,intros$site_i=="N")

sol_hi<-subset(sol,subset = sol$Location_of_introduction == "Hawaiian_Islands")
us_hi<-subset(intros,intros$site_i=="H")


#The data in NZ is broken down by different introductions, need to convert to #events, # individuals

summarized_nz_output<-NULL
for(i in 1:length(unique(sol_nz$Species))){
  species_i<-as.character(unique(sol_nz$Species)[i] ) 
  data_i<-subset(sol_nz,sol_nz$Species==species_i)  
  
  n_events<-nrow(data_i)
  n_individuals<-sum(data_i$Propagule_size,na.rm = T)
  
  if(sum(data_i$Outcome_binary)>0){success<-1}else{success<-0}  
  
  species_data_i<-unique(data_i[,c(9:32,34:35)])  
  climate_match_mean<-mean(data_i$Climate_matching,na.rm = T)  
  
  output_i<-cbind(species_i,success,n_individuals,n_events,species_data_i,climate_match_mean)
  summarized_nz_output<-rbind(summarized_nz_output,output_i)
  
  
}

rm(species_i,data_i,n_events,n_individuals,success,species_data_i,climate_match_mean,output_i,i)



#The data in HI is broken down by different introductions, need to convert to #events, # individuals

summarized_hi_output<-NULL
for(i in 1:length(unique(sol_hi$Species))){
  species_i<-as.character(unique(sol_hi$Species)[i] ) 
  data_i<-subset(sol_hi,sol_hi$Species==species_i)  
  
  n_events<-nrow(data_i)
  n_individuals<-sum(data_i$Propagule_size,na.rm = T)
  
  if(sum(data_i$Outcome_binary)>0){success<-1}else{success<-0}  
  
  species_data_i<-unique(data_i[,c(9:32,34:35)])  
  climate_match_mean<-mean(data_i$Climate_matching,na.rm = T)  
  
  output_i<-cbind(species_i,success,n_individuals,n_events,species_data_i,climate_match_mean)
  summarized_hi_output<-rbind(summarized_hi_output,output_i)
  
  
}

rm(species_i,data_i,n_events,n_individuals,success,species_data_i,climate_match_mean,output_i,i)



#Combine with our data:
full_data_set_nz <- merge(x = us_nz,y = summarized_nz_output,by = "species_i")
full_data_set_nz <- full_data_set_nz[which(!is.na(full_data_set_nz$site_i)),]
full_data_set_nz <- full_data_set_nz[which(!is.na(full_data_set_nz$success)),]

full_data_set_hi <- merge(x = us_hi,y = summarized_hi_output,by = "species_i")
full_data_set_hi <- full_data_set_hi[which(!is.na(full_data_set_hi$site_i)),]
full_data_set_hi <- full_data_set_hi[which(!is.na(full_data_set_hi$success)),]

full_data_set<- rbind(full_data_set_hi,full_data_set_nz)
rm(full_data_set_hi,full_data_set_nz,us_nz,us_hi,sol_hi,sol_nz,sol)






# Traits included in model selection were those identified by Sol et al. (2012) in their best-fit model and included 
# body mass, the residuals of brain mass against body mass, the value of a brood relative to expected lifetime reproductive output, and habitat generalism (the number of habitat types used by a species). 
# We used stepwise model selection (R Core Team 2016) to test hypotheses based on 
  # 1) Propagule pressure only; 
  # 2) Propagule pressure and MPD metrics; and 
  # 3) Propagule pressure and species’ traits. 
# Stepwise model selection was based on AIC and included both forward and backward selection.  
# Propagule pressure metrics included both the number of individuals and number of introduction events. 
# MPD metrics included both source and recipient community MPD.


# 1) Propagule pressure only;


# 2) Propagule pressure and phylo metrics


# 3) Propagule pressure and species’ traits


# 4) Propagule pressure and species’ traits and phylo metrics

##############################################################################
#Model selection in a non-phylo context
full_data_set$site_i<- as.character(full_data_set$site_i)
full_data_set$site_i<- as.factor(full_data_set$site_i)
str(full_data_set)

full_data_set<-full_data_set[c("r_success","s_range_size","s_richness","s_pd","s_nnd","s_mpd","s_vpd","s_spd","s_kpd","r_n_nnd","r_n_mpd","r_n_vpd","r_n_spd","r_n_kpd","r_e_nnd","r_e_mpd","r_e_vpd",
  "r_e_spd","r_e_kpd","r_ne_nnd","r_ne_mpd","r_ne_vpd","r_ne_spd","r_ne_kpd","n_individuals","n_events","Body_mass","Brain_residual","Brood_value","Habitat_generalism","site_i")]

full_data_set <- na.omit(full_data_set)
full_data_set$site_numeric<- as.numeric(full_data_set$site_i)




full_model = r_success ~ s_range_size + s_richness + s_pd + s_nnd + s_mpd + s_vpd + s_spd + s_kpd + r_n_nnd + r_n_mpd + r_n_vpd + r_n_spd + r_n_kpd + 
  r_e_nnd + r_e_mpd + r_e_vpd + r_e_spd + r_e_kpd + r_ne_nnd + r_ne_mpd + r_ne_vpd + r_ne_spd + r_ne_kpd + 
  n_individuals + n_events + Body_mass + Brain_residual + Brood_value + Habitat_generalism + (1|site_i) + (1|species_i__)

?glmer
prop_pressure_full <-lme4::glmer(formula = r_success ~ s_range_size + s_richness + s_pd + s_nnd + s_mpd + s_vpd + s_spd + s_kpd + r_n_nnd + r_n_mpd + r_n_vpd + r_n_spd + r_n_kpd + 
                      r_e_nnd + r_e_mpd + r_e_vpd + r_e_spd + r_e_kpd + r_ne_nnd + r_ne_mpd + r_ne_vpd + r_ne_spd + r_ne_kpd + 
                      n_individuals + n_events + Body_mass + Brain_residual + Brood_value + Habitat_generalism + (1|site_i),
                    family = "binomial",
                    data = full_data_set)

LMERConvenienceFunctions::fitLMER.fnc(model = prop_pressure_full,method = "AIC")



prop_pressure_full <-glm(formula = r_success ~ s_range_size + s_richness + s_pd + s_nnd + s_mpd + s_vpd + s_spd + s_kpd + r_n_nnd + r_n_mpd + r_n_vpd + r_n_spd + r_n_kpd + 
                                   r_e_nnd + r_e_mpd + r_e_vpd + r_e_spd + r_e_kpd + r_ne_nnd + r_ne_mpd + r_ne_vpd + r_ne_spd + r_ne_kpd + 
                                   n_individuals + n_events + Body_mass + Brain_residual + Brood_value + Habitat_generalism +(1|site_numeric),
                                 family = "binomial",
                                 data = full_data_set)

prop_pressure_full <-glm(formula = r_success ~ 1,
                         family = "binomial",
                         data = full_data_set)


step_non_phylo <- MASS::stepAIC(object = prop_pressure_full,direction = "both",
                  scope = list(upper = ~s_range_size + s_richness + s_pd + s_nnd + s_mpd + s_vpd + s_spd + s_kpd + r_n_nnd + r_n_mpd + r_n_vpd + r_n_spd + r_n_kpd + 
                             r_e_nnd + r_e_mpd + r_e_vpd + r_e_spd + r_e_kpd + r_ne_nnd + r_ne_mpd + r_ne_vpd + r_ne_spd + r_ne_kpd + 
                             n_individuals + n_events + Body_mass + Brain_residual + Brood_value + Habitat_generalism +site_i,
                           lower = ~ 1   ) )

summary(step_non_phylo)
r2glmm::glmPQL(glm.mod = step_non_phylo)
r2glmm::r2beta(step_non_phylo)



?MASS::stepAIC
lmerTest::step(object = prop_pressure_full,reduce.random=F)
?lmerTest::step

?stepAIC

?step
prop_only_w_rs<-step(prop_only,direction = "both")
summary(prop_only_w_rs)
r.squaredGLMM(prop_only_w_rs)



##############################################################################
##############################################################################
#Model selection in a pglmm context


  #use a priori models
full_model = r_success ~ s_range_size + s_richness + s_pd + s_nnd + s_mpd + s_vpd + s_spd + s_kpd + r_n_nnd + r_n_mpd + r_n_vpd + r_n_spd + r_n_kpd + 
  r_e_nnd + r_e_mpd + r_e_vpd + r_e_spd + r_e_kpd + r_ne_nnd + r_ne_mpd + r_ne_vpd + r_ne_spd + r_ne_kpd + 
  n_individuals + n_events + Body_mass + Brain_residual + Brood_value + Habitat_generalism + (1|site_i) + (1|species_i__)


# 1) Propagule pressure only;

prop_only <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                           formula = r_success ~ n_individuals + n_events + (1|site_i) + (1|species_i__),
                                           data = full_data_set,
                                           family = "binomial")

prop_only_mean <- summarize_replicates(phylo_reps_output = prop_only)



# 2) Propagule pressure and phylo metrics

prop_and_phylo <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                  formula = r_success ~ s_mpd + r_n_mpd + n_individuals + n_events + (1|site_i) + (1|species_i__),
                                  data = full_data_set,
                                  family = "binomial")

prop_and_phylo_mean <- summarize_replicates(phylo_reps_output = prop_and_phylo)


#2.1) Propagule pressure and phylo metrics with interaction

prop_and_phylo_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                       formula = r_success ~ s_mpd*r_n_mpd+ s_mpd + r_n_mpd + n_individuals + n_events + (1|site_i) + (1|species_i__),
                                       data = full_data_set,
                                       family = "binomial")

prop_and_phylo_int_mean <- summarize_replicates(phylo_reps_output = prop_and_phylo_int)




# 3) Propagule pressure and species’ traits

prop_and_traits <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                       formula = r_success ~ Body_mass + Brain_residual + Brood_value + Habitat_generalism + n_individuals + n_events + (1|site_i) + (1|species_i__),
                                       data = full_data_set,
                                       family = "binomial")

prop_and_traits_mean <- summarize_replicates(phylo_reps_output = prop_and_traits)

# 4) Propagule pressure and species’ traits and phylo metrics

prop_and_phylo_and_traits <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                       formula = r_success ~ s_mpd + r_n_mpd + Body_mass + Brain_residual + Brood_value + Habitat_generalism + 
                                         n_individuals + n_events + (1|site_i) + (1|species_i__),
                                       data = full_data_set,
                                       family = "binomial")

prop_and_phylo__and_traits_mean <- summarize_replicates(phylo_reps_output = prop_and_phylo_and_traits)

# 4.1) Propagule pressure and species’ traits and phylo metrics with int

prop_and_phylo_and_traits_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                                  formula = r_success ~ s_mpd*r_n_mpd + s_mpd + r_n_mpd + Body_mass + Brain_residual + Brood_value + Habitat_generalism + 
                                                    n_individuals + n_events + (1|site_i) + (1|species_i__),
                                                  data = full_data_set,
                                                  family = "binomial")

prop_and_phylo__and_traits_int_mean <- summarize_replicates(phylo_reps_output = prop_and_phylo_and_traits_int)




# 2) Phylo only

phylo_only <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                       formula = r_success ~ s_mpd + r_n_mpd + (1|site_i) + (1|species_i__),
                                       data = full_data_set,
                                       family = "binomial")

phylo_only_mean <- summarize_replicates(phylo_reps_output = phylo_only)


# 2.1) Phylo only int

phylo_only_int <- replicate_phylo_glmm(tree_list = list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)[1:100],
                                   formula = r_success ~ s_mpd*r_n_mpd + s_mpd + r_n_mpd + (1|site_i) + (1|species_i__),
                                   data = full_data_set,
                                   family = "binomial")

phylo_only_int_mean <- summarize_replicates(phylo_reps_output = phylo_only_int)





prop_only_mean

prop_and_phylo_mean

#prop_and_phylo_int_mean #lower AIC than without int

prop_only_mean$aic_mean - prop_and_phylo_mean$aic_mean #AIC diff when adding phylogeny to propagule model

prop_and_traits_mean

prop_and_phylo__and_traits_mean

#prop_and_phylo__and_traits_mean$aic_mean-prop_and_phylo__and_traits_int_mean$aic_mean
#prop_and_phylo__and_traits_int_mean #lower AIC than without int
#prop_and_phylo_mean$aic_mean - prop_and_phylo_int_mean$aic_mean


prop_and_phylo_mean$aic_mean - prop_and_traits_mean$aic_mean

phylo_only_mean

#phylo_only_int_mean #lower AIC than without int
#phylo_only_mean$aic_mean - phylo_only_int_mean$aic_mean



############################################################################
# Graph: Comparison of R2 values and partial R2 values for models incorporating propagule pressure. All models included establishment success as a binary response variable. Predictor variables included: (A) propagule number and size; (B) propagule number and size and recipient community MPD; and (C) propagule number and size, habitat generalism, brood value and brain size residuals. 

############################################################################

cor(intros$s_mpd,intros$s_richness)
cor(intros$s_richness,intros$s_nnd)
cor(intros$s_richness,intros$s_mpd)
cor(intros$s_richness,intros$s_vpd)
cor(intros$s_richness,intros$s_spd)
cor(intros$s_richness,intros$s_kpd)

cor(intros$s_mpd,intros$s_v)

