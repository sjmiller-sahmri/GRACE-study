#Load libraries
library(vegan)

#Modified cox
ModifiedCox <- function(x){
  n <- length(x)
  y <- log(x)
  y.m <- mean(y)
  y.var <- var(y)
  
  my.t <- qt(0.975, df = n-1) # 95% Confidence interval
  
  my.mean <- mean(x)
  upper <- y.m + y.var/2 + my.t*sqrt(y.var/n + y.var^2/(2*(n - 1)))
  lower <- y.m + y.var/2 - my.t*sqrt(y.var/n + y.var^2/(2*(n - 1)))
  
  return(list(upper = exp(upper), mean = my.mean, lower = exp(lower)))
 }

Singleton <- c(141, 129, 108, 41, 39, 19, 14, 13, 11, 10, 9, 8, 7, 6, 4, 4, 1, 1, 1, 1)
ModifiedCox(Singleton)

#PERMANOVA
names <- data[,c(1,588:603)]

taxa <- data[,1:587]
rownames(taxa) <- taxa[,1]
taxa <- taxa[,2:587]
sqrt_taxa <- sqrt(taxa)

perm_exposureevents_adj <- adonis2(sqrt_taxa ~ num_exposure_events + depression + gastro_ref_disease + hypertension + pain + age_at_enrolment_years + sex + days_since_entry + site, data=data, permutations=999, method="bray")
summary(perm_exposuredays)

perm_dayssinceexposure_adj <- adonis2(sqrt_taxa ~ days_since_most_recent_abx + depression + gastro_ref_disease + hypertension + pain + age_at_enrolment_years + sex + days_since_entry + site, data=data, permutations=999, method="bray")
summary(perm_dayssinceexposure)

perm_exposuredays_adj <- adonis2(sqrt_taxa ~ total_exposed_days + depression + gastro_ref_disease + hypertension + pain + age_at_enrolment_years + sex + days_since_entry + site, data=data, permutations=999, method="bray")
summary(perm_exposuredays)

perm_uniqueclasses_adj <- adonis2(sqrt_taxa ~ num_unique_classes + depression + gastro_ref_disease + hypertension + pain + age_at_enrolment_years + sex + days_since_entry + site, data=data, permutations=999, method="bray")
summary(perm_exposuredays)

perm_Cephalosporin_adj <- adonis2(sqrt_taxa ~ any_cephalosporin_used + depression + gastro_ref_disease + hypertension + pain + age_at_enrolment_years + sex + days_since_entry + site, data=data, permutations=999, method="bray")
summary(perm_Cephalosporin_adj)

perm_Diaminopyrimidine_adj <- adonis2(sqrt_taxa ~ any_diaminopyrimidine_used + depression + gastro_ref_disease + hypertension + pain + age_at_enrolment_years + sex + days_since_entry + site, data=data, permutations=999, method="bray")
summary(perm_Diaminopyrimidine_adj)

perm_Penicillin_adj <- adonis2(sqrt_taxa ~ any_penam_used + depression + gastro_ref_disease + hypertension + pain + age_at_enrolment_years + sex + days_since_entry + site, data=data, permutations=999, method="bray")
summary(perm_Penicillin_adj)

perm_Tetracycline_adj <- adonis2(sqrt_taxa ~ any_tetracycline_used + depression + gastro_ref_disease + hypertension + pain + age_at_enrolment_years + sex + days_since_entry + site, data=data, permutations=999, method="bray")
summary(perm_Tetracycline_adj)
