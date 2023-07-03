#################
# EHRs as complex surveys: Holloway CCHS
# Citation: Goldstein ND, Jones J, Kahal D, Burstyn I. Inferring population HIV viral load from a single HIV clinic's electronic health record. Manuscript in preparation.
# 2/27/23 -- Neal Goldstein
#################


### FUNCTIONS ###

#library(survey)
library(psych) #describeby
library(WeightIt) #trim excessive weights
#library(Hmisc) #weighted variance


### READ DATA  ###

plwh = read.csv("holloway.csv", na.strings="", stringsAsFactors=F)

#recodes
plwh$Age = ifelse(plwh$Age<35, 0, ifelse(plwh$Age>=35 & plwh$Age<45, 1, ifelse(plwh$Age>=45 & plwh$Age<55, 2, 3)))
plwh$Gender = ifelse(plwh$Gender=="Female" | plwh$Gender=="Transgender MtF", "F", "M")
plwh$Race = ifelse(plwh$RaceEthnicity=="White (non-Hispanic)", "W", ifelse(plwh$RaceEthnicity=="Black or African-American", "B", ifelse(plwh$RaceEthnicity=="Hispanic", "H", NA)))
plwh$VL = plwh$LastVL

#remove extraneous variables
plwh$RaceEthnicity = NULL
plwh$LastVLDate = NULL
plwh$LastVL = NULL

#complete records only
plwh = plwh[complete.cases(plwh), ]
n_clinic_a = nrow(plwh)

plwh$Clinic_A = T

#census populations for denominators (based on Delaware catchment)
#https://data.census.gov/table?q=Delaware&tid=ACSDP1Y2021.DP05 (2021 ACS 1-yr)
n_pop = 1003384
n_white = 690900
n_black = 243431
n_hispanic = 101213
n_male = 485908
n_female = 517476
n_35less = 53039+58410+59908+63288+57635+128564
n_35to44 = 122088
n_45to54 = 115300
n_55more = 67381+76220+125407+57183+18961

#plwh populations for denominators (based on state-wide AIDS report, assuming exchangeability)
#https://www.dhss.delaware.gov/dph/dpc/files/comphivplan.pdf, pg38 HIV care continuum
n_plwh = 3841
n_incare = 2984
n_notincare = n_plwh-n_incare
n_incare_white = 958
n_incare_black = 1725
n_incare_hispanic = 222
n_incare_male = 2125
n_incare_female = 859
n_incare_35less = 7+57+330
n_incare_35to44 = 432
n_incare_45to54 = 703
n_incare_55more = 1006+449


### WEIGHTED ANALYSIS ###

#weights are specified as:
#K = source population (n_pop)
#P = PLWH (n_plwh)
#N = PLWH and in care (n_incare)
#S = clinic sample (n_clinic_a...h)
#weight ~ 1 / Beta((S + 1), (N + 1 - S))
#where N ~ Binom(pi, K) and pi = N/K (prevalence of HIV+ and in care in source population)

#overall for each clinic (if selection forces unrelated to characteristics)
#wt_pi_overall = n_incare / n_pop
#wt_hiv_pop_overall = rbinom(1, n_pop, wt_pi_overall)
#wt_a_overall = 1/rbeta(1, (n_clinic_a+1), (wt_hiv_pop_overall + 1 - n_clinic_a))
#wt_b_overall = 1/rbeta(1, (n_clinic_b+1), (wt_hiv_pop_overall + 1 - n_clinic_b))

#simulate 1000 of each weight
clinic_a = NA
clinic_a_biased = NA
clinic_a_bayes_n1 = NA
clinic_a_bayes_n2 = NA
clinic_a_bayes_n3 = NA

i=1
for (j in 1:1000) {
  
  cat("\n\n************** ","Simulation: ",i," Weight: ",j," **************\n",sep="")
  
  #weights by race
  wt_pi_white = n_incare_white / n_pop
  wt_hiv_pop_white = rbinom(1, n_pop, wt_pi_white)
  while (wt_hiv_pop_white < sum(plwh$Race=="W")) { wt_hiv_pop_white = rbinom(1, n_pop, wt_pi_white) } #resample if N<S
  wt_a_white = 1/rbeta(1, (sum(plwh$Race=="W")+1), (wt_hiv_pop_white + 1 - sum(plwh$Race=="W")))
  
  wt_pi_black = n_incare_black / n_pop
  wt_hiv_pop_black = rbinom(1, n_pop, wt_pi_black)
  while (wt_hiv_pop_black < sum(plwh$Race=="B")) { wt_hiv_pop_black = rbinom(1, n_pop, wt_pi_black) } #resample if N<S
  wt_a_black = 1/rbeta(1, (sum(plwh$Race=="B")+1), (wt_hiv_pop_black + 1 - sum(plwh$Race=="B")))
  
  wt_pi_hispanic = n_incare_hispanic / n_pop
  wt_hiv_pop_hispanic = rbinom(1, n_pop, wt_pi_hispanic)
  while (wt_hiv_pop_hispanic < sum(plwh$Race=="H")) { wt_hiv_pop_hispanic = rbinom(1, n_pop, wt_pi_hispanic) } #resample if N<S
  wt_a_hispanic = 1/rbeta(1, (sum(plwh$Race=="H")+1), (wt_hiv_pop_hispanic + 1 - sum(plwh$Race=="H")))
  
  #weights by gender
  wt_pi_male = n_incare_male / n_pop
  wt_hiv_pop_male = rbinom(1, n_pop, wt_pi_male)
  while (wt_hiv_pop_male < sum(plwh$Gender=="M")) { wt_hiv_pop_male = rbinom(1, n_pop, wt_pi_male) } #resample if N<S
  wt_a_male = 1/rbeta(1, (sum(plwh$Gender=="M")+1), (wt_hiv_pop_male + 1 - sum(plwh$Gender=="M")))
  
  wt_pi_female = n_incare_female / n_pop
  wt_hiv_pop_female = rbinom(1, n_pop, wt_pi_female)
  while (wt_hiv_pop_female < sum(plwh$Gender=="F")) { wt_hiv_pop_female = rbinom(1, n_pop, wt_pi_female) } #resample if N<S
  wt_a_female = 1/rbeta(1, (sum(plwh$Gender=="F")+1), (wt_hiv_pop_female + 1 - sum(plwh$Gender=="F")))
  
  #weights by age category
  wt_pi_35less = n_incare_35less / n_pop
  wt_hiv_pop_35less = rbinom(1, n_pop, wt_pi_35less)
  while (wt_hiv_pop_35less < sum(plwh$Age==0)) { wt_hiv_pop_35less = rbinom(1, n_pop, wt_pi_35less) } #resample if N<S
  wt_a_35less = 1/rbeta(1, (sum(plwh$Age==0)+1), (wt_hiv_pop_35less + 1 - sum(plwh$Age==0)))
  
  wt_pi_35to44 = n_incare_35to44 / n_pop
  wt_hiv_pop_35to44 = rbinom(1, n_pop, wt_pi_35to44)
  while (wt_hiv_pop_35to44 < sum(plwh$Age==1)) { wt_hiv_pop_35to44 = rbinom(1, n_pop, wt_pi_35to44) } #resample if N<S
  wt_a_35to44 = 1/rbeta(1, (sum(plwh$Age==1)+1), (wt_hiv_pop_35to44 + 1 - sum(plwh$Age==1)))
  
  wt_pi_45to54 = n_incare_45to54 / n_pop
  wt_hiv_pop_45to54 = rbinom(1, n_pop, wt_pi_45to54)
  while (wt_hiv_pop_45to54 < sum(plwh$Age==2)) { wt_hiv_pop_45to54 = rbinom(1, n_pop, wt_pi_45to54) } #resample if N<S
  wt_a_45to54 = 1/rbeta(1, (sum(plwh$Age==2)+1), (wt_hiv_pop_45to54 + 1 - sum(plwh$Age==2)))
  
  wt_pi_55more = n_incare_55more / n_pop
  wt_hiv_pop_55more = rbinom(1, n_pop, wt_pi_55more)
  while (wt_hiv_pop_55more < sum(plwh$Age==3)) { wt_hiv_pop_55more = rbinom(1, n_pop, wt_pi_55more) } #resample if N<S
  wt_a_55more = 1/rbeta(1, (sum(plwh$Age==3)+1), (wt_hiv_pop_55more + 1 - sum(plwh$Age==3)))
  
  #assign strata specific weights to population
  plwh$Weight_a = ifelse(plwh$Clinic_A==T & plwh$Age==0, wt_a_35less, ifelse(plwh$Clinic_A==T & plwh$Age==1, wt_a_35to44, ifelse(plwh$Clinic_A==T & plwh$Age==2, wt_a_45to54, ifelse(plwh$Clinic_A==T & plwh$Age==3, wt_a_55more, NA)))) *
    ifelse(plwh$Clinic_A==T & plwh$Gender=="M", wt_a_male, ifelse(plwh$Clinic_A==T & plwh$Gender=="F", wt_a_female, NA)) *
    ifelse(plwh$Clinic_A==T & plwh$Race=="W", wt_a_white, ifelse(plwh$Clinic_A==T & plwh$Race=="B", wt_a_black, ifelse(plwh$Clinic_A==T & plwh$Race=="H", wt_a_hispanic, NA)))
  
  #weight trimming parameter: see trim function in WeightIt
  wt_trim = 0.98
  
  #weighted clinics using geometric mean: exp(mean(log(X)))
  clinic_a = c(clinic_a, exp(weighted.mean(log(plwh$VL[plwh$Clinic_A==T]), suppressMessages(trim(plwh$Weight_a[plwh$Clinic_A==T], at=wt_trim, lower=T)))))
  
  #weighted clinics, weights misspecified using geometric mean 
  #misspecification form: logit_pO=log(pT/(1-pT))+b*logVL, where pT is the observed probability (inverse of weight) and b is the bias factor (0=no bias)
  bias_factor = 0.1
  sample_probs_observed = 1/plwh$Weight_a[plwh$Clinic_A==T]
  biased_probs = log(sample_probs_observed/(1-sample_probs_observed))+bias_factor*log(plwh$VL[plwh$Clinic_A==T])
  sample_probs_biased = exp(biased_probs)/(1+exp(biased_probs))
  clinic_a_biased = c(clinic_a_biased, exp(weighted.mean(log(plwh$VL[plwh$Clinic_A==T]), suppressMessages(trim((1/sample_probs_biased), at=wt_trim, lower=T)))))
  
  #bayesian weighted, prior: informed by VL measures in population (see population initialization)
  m0 = log10(10^(log10(23348) - 0.5*(1.2^2))) #geometric mean
  #s20 = 1.2^2 #variance, not used
  #n0 = clinic_n*adjustment #prior sample size
  
  #bayesian weighted, posterior, https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf (eq 86, prior mean and variance unknown)
  clinic_xbar = log10(exp(weighted.mean(log(plwh$VL[plwh$Clinic_A==T]), suppressMessages(trim(plwh$Weight_a[plwh$Clinic_A==T], at=wt_trim, lower=T)))))
  clinic_n = sum(suppressMessages(trim(plwh$Weight_a[plwh$Clinic_A==T], at=wt_trim, lower=T)))
  n1 = clinic_n*0.25 #sample size relative to clinic size
  n2 = clinic_n*0.5 #sample size relative to clinic size
  n3 = clinic_n*2 #sample size relative to clinic size
  clinic_a_bayes_n1 = c(clinic_a_bayes_n1, 10^(((n1*m0) + (clinic_n*clinic_xbar))/(n1 + clinic_n)))
  clinic_a_bayes_n2 = c(clinic_a_bayes_n2, 10^(((n2*m0) + (clinic_n*clinic_xbar))/(n2 + clinic_n)))
  clinic_a_bayes_n3 = c(clinic_a_bayes_n3, 10^(((n3*m0) + (clinic_n*clinic_xbar))/(n3 + clinic_n)))
  
  rm(m0,n1,n2,n3,clinic_xbar,clinic_n,wt_trim,sample_probs_observed,biased_probs,sample_probs_biased,bias_factor)
}
rm(j)
rm(i)


### VISUALIZE RESULTS ###

#naive to weighted
boxplot(log10(clinic_a), at=2, xlim = c(0.7, 6.3), ylim = range(1.5, 3.0), xaxt = "n", ylab="Mean Viral Load (log10 copies/mL)")
boxplot(log10(clinic_a_biased), at=3, xaxt="n", add=T)
boxplot(log10(clinic_a_bayes_n1), at=4, xaxt="n", add=T)
boxplot(log10(clinic_a_bayes_n2), at=5, xaxt="n", add=T)
boxplot(log10(clinic_a_bayes_n3), at=6, xaxt="n", add=T)
points(1,log10(exp(mean(log(plwh$VL)))), pch=18, cex=2)
lines(x=c(0,20), y=rep(log10(exp(mean(log(plwh$VL)))),2), lty=2)
#lines(x=c(0,20), y=rep(log10(10^(log10(23348) - 0.5*(1.2^2))),2), lty=3)
axis(1, at=1:6, labels=c("Observed","Weighted","Biased",expression(Bayesian^1),expression(Bayesian^2),expression(Bayesian^3)), tick=T, las=3, cex.axis=1)


### PAPER DESCRIPTIVES ###

nrow(plwh)
log10(exp(mean(log(plwh$VL)))) #geometric mean
log10(sd(plwh$VL))
sum(plwh$VL<=20)
sum(plwh$VL<=20)/nrow(plwh)

median(log10(clinic_a_biased), na.rm=T)
median(log10(clinic_a_bayes_n1), na.rm=T)
median(log10(clinic_a_bayes_n2), na.rm=T)
median(log10(clinic_a_bayes_n3), na.rm=T)

