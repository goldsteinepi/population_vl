#################
# EHRs as complex surveys: simulation
# Citation: Goldstein ND, Jones J, Kahal D, Burstyn I. Inferring population HIV viral load from a single HIV clinic's electronic health record. Manuscript in preparation.
# 1/11/23 -- Neal Goldstein
#################


### FUNCTIONS ###

#library(survey)
library(psych) #describeby
library(WeightIt) #trim excessive weights
#library(Hmisc) #weighted variance
library(gmodels) #CrossTable


### CREATE DATA and SIMULATE ###

#simulate 1000 of each clinic
set.seed(777)
#sim_results = data.frame("Simulation"=1:1000,"VL_population"=NA,"VL_population_nodx"=NA,"VL_population_dx"=NA,"VL_population_notincare"=NA,"VL_population_incare"=NA,"VL_clinic_A"=NA,"VL_clinic_B"=NA,"VL_clinic_C"=NA,"VL_clinic_D"=NA,"VL_clinic_E"=NA,"VL_clinic_F"=NA,"VL_clinic_G"=NA,"VL_clinic_A_weighted"=NA,"VL_clinic_B_weighted"=NA,"VL_clinic_C_weighted"=NA,"VL_clinic_D_weighted"=NA,"VL_clinic_E_weighted"=NA,"VL_clinic_F_weighted"=NA,"VL_clinic_G_weighted"=NA,"VL_clinic_A_rmse"=NA,"VL_clinic_B_rmse"=NA,"VL_clinic_C_rmse"=NA,"VL_clinic_D_rmse"=NA,"VL_clinic_E_rmse"=NA,"VL_clinic_F_rmse"=NA,"VL_clinic_G_rmse"=NA,"VL_clinic_A_weighted_under"=NA,"VL_clinic_B_weighted_under"=NA,"VL_clinic_C_weighted_under"=NA,"VL_clinic_D_weighted_under"=NA,"VL_clinic_E_weighted_under"=NA,"VL_clinic_F_weighted_under"=NA,"VL_clinic_G_weighted_under"=NA,"VL_clinic_A_rmse_under"=NA,"VL_clinic_B_rmse_under"=NA,"VL_clinic_C_rmse_under"=NA,"VL_clinic_D_rmse_under"=NA,"VL_clinic_E_rmse_under"=NA,"VL_clinic_F_rmse_under"=NA,"VL_clinic_G_rmse_under"=NA,"VL_clinic_A_weighted_over"=NA,"VL_clinic_B_weighted_over"=NA,"VL_clinic_C_weighted_over"=NA,"VL_clinic_D_weighted_over"=NA,"VL_clinic_E_weighted_over"=NA,"VL_clinic_F_weighted_over"=NA,"VL_clinic_G_weighted_over"=NA,"VL_clinic_A_rmse_over"=NA,"VL_clinic_B_rmse_over"=NA,"VL_clinic_C_rmse_over"=NA,"VL_clinic_D_rmse_over"=NA,"VL_clinic_E_rmse_over"=NA,"VL_clinic_F_rmse_over"=NA,"VL_clinic_G_rmse_over"=NA,"VL_clinic_A_weighted_mixed"=NA,"VL_clinic_B_weighted_mixed"=NA,"VL_clinic_C_weighted_mixed"=NA,"VL_clinic_D_weighted_mixed"=NA,"VL_clinic_E_weighted_mixed"=NA,"VL_clinic_F_weighted_mixed"=NA,"VL_clinic_G_weighted_mixed"=NA,"VL_clinic_A_rmse_mixed"=NA,"VL_clinic_B_rmse_mixed"=NA,"VL_clinic_C_rmse_mixed"=NA,"VL_clinic_D_rmse_mixed"=NA,"VL_clinic_E_rmse_mixed"=NA,"VL_clinic_F_rmse_mixed"=NA,"VL_clinic_G_rmse_mixed"=NA,"VL_clinic_A_weighted_bayes"=NA,"VL_clinic_B_weighted_bayes"=NA,"VL_clinic_C_weighted_bayes"=NA,"VL_clinic_D_weighted_bayes"=NA,"VL_clinic_E_weighted_bayes"=NA,"VL_clinic_F_weighted_bayes"=NA,"VL_clinic_G_weighted_bayes"=NA,"VL_clinic_A_rmse_bayes"=NA,"VL_clinic_B_rmse_bayes"=NA,"VL_clinic_C_rmse_bayes"=NA,"VL_clinic_D_rmse_bayes"=NA,"VL_clinic_E_rmse_bayes"=NA,"VL_clinic_F_rmse_bayes"=NA,"VL_clinic_G_rmse_bayes"=NA,stringsAsFactors=F)
#sim_results = data.frame("Simulation"=1:1000, "VL_population"=NA, "VL_population_nodx"=NA, "VL_population_dx"=NA, "VL_population_notincare"=NA, "VL_population_incare"=NA, "VL_clinic_A"=NA, "VL_clinic_B"=NA, "VL_clinic_C"=NA, "VL_clinic_A_weighted"=NA, "VL_clinic_B_weighted"=NA, "VL_clinic_C_weighted"=NA, "VL_clinic_A_rmse"=NA, "VL_clinic_B_rmse"=NA, "VL_clinic_C_rmse"=NA, "VL_clinic_A_weighted_biased"=NA, "VL_clinic_B_weighted_biased"=NA, "VL_clinic_C_weighted_biased"=NA, "VL_clinic_A_rmse_biased"=NA, "VL_clinic_B_rmse_biased"=NA, "VL_clinic_C_rmse_biased"=NA, "VL_clinic_A_weighted_over"=NA, "VL_clinic_B_weighted_over"=NA, "VL_clinic_C_weighted_over"=NA, "VL_clinic_A_rmse_over"=NA, "VL_clinic_B_rmse_over"=NA, "VL_clinic_C_rmse_over"=NA, "VL_clinic_A_weighted_mixed"=NA, "VL_clinic_B_weighted_mixed"=NA, "VL_clinic_C_weighted_mixed"=NA, "VL_clinic_A_rmse_mixed"=NA, "VL_clinic_B_rmse_mixed"=NA, "VL_clinic_C_rmse_mixed"=NA, "VL_clinic_A_weighted_bayes_n1"=NA, "VL_clinic_B_weighted_bayes_n1"=NA, "VL_clinic_C_weighted_bayes_n1"=NA, "VL_clinic_A_weighted_bayes_n2"=NA, "VL_clinic_B_weighted_bayes_n2"=NA, "VL_clinic_C_weighted_bayes_n2"=NA,"VL_clinic_A_weighted_bayes_n3"=NA, "VL_clinic_B_weighted_bayes_n3"=NA, "VL_clinic_C_weighted_bayes_n3"=NA, "VL_clinic_C_weighted_bayes_m1"=NA, "VL_clinic_A_rmse_bayes_n1"=NA, "VL_clinic_B_rmse_bayes_n1"=NA, "VL_clinic_C_rmse_bayes_n1"=NA, "VL_clinic_A_rmse_bayes_n2"=NA, "VL_clinic_B_rmse_bayes_n2"=NA, "VL_clinic_C_rmse_bayes_n2"=NA, "VL_clinic_A_rmse_bayes_n3"=NA, "VL_clinic_B_rmse_bayes_n3"=NA, "VL_clinic_C_rmse_bayes_n3"=NA, "VL_clinic_C_rmse_bayes_m1"=NA, stringsAsFactors=F)
sim_results = data.frame("Simulation"=1:1000, "VL_population"=NA, "VL_population_nodx"=NA, "VL_population_dx"=NA, "VL_population_notincare"=NA, "VL_population_incare"=NA, "VL_clinic_A"=NA, "VL_clinic_B"=NA, "VL_clinic_C"=NA, "VL_clinic_A_weighted"=NA, "VL_clinic_B_weighted"=NA, "VL_clinic_C_weighted"=NA, "VL_clinic_A_rmse"=NA, "VL_clinic_B_rmse"=NA, "VL_clinic_C_rmse"=NA, "VL_clinic_A_weighted_biased"=NA, "VL_clinic_B_weighted_biased"=NA, "VL_clinic_C_weighted_biased"=NA, "VL_clinic_A_rmse_biased"=NA, "VL_clinic_B_rmse_biased"=NA, "VL_clinic_C_rmse_biased"=NA, "VL_clinic_A_weighted_bayes_n1"=NA, "VL_clinic_B_weighted_bayes_n1"=NA, "VL_clinic_C_weighted_bayes_n1"=NA, "VL_clinic_A_weighted_bayes_n2"=NA, "VL_clinic_B_weighted_bayes_n2"=NA, "VL_clinic_C_weighted_bayes_n2"=NA,"VL_clinic_A_weighted_bayes_n3"=NA, "VL_clinic_B_weighted_bayes_n3"=NA, "VL_clinic_C_weighted_bayes_n3"=NA, "VL_clinic_C_weighted_bayes_m1"=NA, "VL_clinic_A_rmse_bayes_n1"=NA, "VL_clinic_B_rmse_bayes_n1"=NA, "VL_clinic_C_rmse_bayes_n1"=NA, "VL_clinic_A_rmse_bayes_n2"=NA, "VL_clinic_B_rmse_bayes_n2"=NA, "VL_clinic_C_rmse_bayes_n2"=NA, "VL_clinic_A_rmse_bayes_n3"=NA, "VL_clinic_B_rmse_bayes_n3"=NA, "VL_clinic_C_rmse_bayes_n3"=NA, "VL_clinic_C_rmse_bayes_m1"=NA, stringsAsFactors=F)
sim_pop = data.frame("Simulation"=1:1000, "Population_age0"=NA, "Population_age1"=NA, "Population_age2"=NA, "Population_age3"=NA, "Population_genderM"=NA, "Population_genderF"=NA, "Population_raceW"=NA, "Population_raceB"=NA, "Population_raceH"=NA, "Clinic_A_age0"=NA, "Clinic_A_age1"=NA, "Clinic_A_age2"=NA, "Clinic_A_age3"=NA, "Clinic_A_genderM"=NA, "Clinic_A_genderF"=NA, "Clinic_A_raceW"=NA, "Clinic_A_raceB"=NA, "Clinic_A_raceH"=NA, "Clinic_B_age0"=NA, "Clinic_B_age1"=NA, "Clinic_B_age2"=NA, "Clinic_B_age3"=NA, "Clinic_B_genderM"=NA, "Clinic_B_genderF"=NA, "Clinic_B_raceW"=NA, "Clinic_B_raceB"=NA, "Clinic_B_raceH"=NA, "Clinic_C_age0"=NA, "Clinic_C_age1"=NA, "Clinic_C_age2"=NA, "Clinic_C_age3"=NA, "Clinic_C_genderM"=NA, "Clinic_C_genderF"=NA, "Clinic_C_raceW"=NA, "Clinic_C_raceB"=NA, "Clinic_C_raceH"=NA, stringsAsFactors=F)

for (i in 1:1000) {
  
  cat("\n\n************** ","Simulation: ",i," **************\n",sep="")
  
  
  ### CREATE HIV+ POPULATION ###
  
  n_pop = 1000000
  n_plwh = 10000 #1% of population is HIV+
  plwh = data.frame("ID"=1:n_plwh,"Age"=NA,"Gender"=NA,"Race"=NA,"VL"=NA,stringsAsFactors=F)
  
  #demographics: matched to "prevalence" https://www.cdc.gov/hiv/library/reports/hiv-surveillance/vol-33/content/national-profile.html
  plwh$Age = sample(c(0,0,1,1,2,2,2,3,3,3,3), n_plwh, replace=T) #0=<35, 1=35-44, 2=45-54, 3=>54
  plwh$Gender = sample(c("M","M","M","F"), n_plwh, replace=T)
  plwh$Race = sample(c("W","W","W","B","B","B","B","H","H"), n_plwh, replace=T)
  
  prop.table(table(plwh$Age))
  prop.table(table(plwh$Gender))
  prop.table(table(plwh$Race))
  
  #VL tuned to San Francisco study
  #arithmetic mean = 23348: https://pubmed.ncbi.nlm.nih.gov/20548786/
  #standard deviation (log10) = 1.2: https://stacks.cdc.gov/view/cdc/28147
  #calculate geometric mean = 10^(log10(23348) - 0.5*(1.2^2)): https://towardsdatascience.com/log-normal-distribution-a-simple-explanation-7605864fb67c
  #rlnorm expects input on natural log scale, not log 10 scale
  plwh$VL = rlnorm(n_plwh, meanlog=(log(10^(log10(23348) - 0.5*(1.2^2)))), sdlog=log(10^1.2))
  
  #race adjustments: W (21087/23348=0.90), B (26404/23348=1.13), H (26774/23348=1.15)
  #gender adjustments: M (23062/23348=0.99), F (27614/23348=1.18)
  #age adjustments: <35 ((33077+23348)/2)/23348=1.21), 35-44 ((23348+18717)/2)/23348=0.90), 45-54 ((18717+15725)/2)/23348=0.74), 55+ ((15725+18824)/2)/23348=0.74)
  #sample adjustments from a random normal (mean=u, sd=0.1)
  plwh$VL = round(plwh$VL * 
                    ifelse(plwh$Age==0, rnorm(n_plwh, mean=1.21, sd=0.1), ifelse(plwh$Age==1, rnorm(n_plwh, mean=0.90, sd=0.1), rnorm(n_plwh, mean=0.74, sd=0.1))) *
                    ifelse(plwh$Race=="W", rnorm(n_plwh, mean=0.90, sd=0.1), ifelse(plwh$Race=="B", rnorm(n_plwh, mean=1.13, sd=0.1), rnorm(n_plwh, mean=1.15, sd=0.1))) *
                    ifelse(plwh$Gender=="M", rnorm(n_plwh, mean=0.99, sd=0.1), rnorm(n_plwh, mean=1.18, sd=0.1)))
  
  #prevent 0 VL
  plwh$VL = ifelse(plwh$VL==0, 1, plwh$VL)
  
  #checks
  # hist(plwh$VL, breaks="fd")
  # describe(plwh$VL)
  # exp(mean(log(plwh$VL))) #geometric mean
  # describeBy(plwh$VL,plwh$Age)
  # describeBy(plwh$VL,plwh$Race)
  # describeBy(plwh$VL,plwh$Gender)
  
  
  ### SAMPLE THOSE NOT DX ###
  
  #10% of population unaware
  plwh$Dx = T
  plwh$Dx_prob = 1
  n_nodx = n_plwh * 0.10
  
  while (sum(plwh$Dx==F)<n_nodx) {
    
    #obtain potential from the population
    potential = sample(plwh$ID[plwh$Dx==T], 1)
    
    rand_in = runif(1)
    rand_notin = runif(1)
    
    #oversample by young <=40, male, black (95% chance of being unaware): https://pubmed.ncbi.nlm.nih.gov/17031318/
    if (plwh$Age[plwh$ID==potential]==0 && plwh$Gender[plwh$ID==potential]=="M" && plwh$Race[plwh$ID==potential]=="B" && rand_in<=0.95) {
      
      plwh$Dx[plwh$ID==potential]=F
      plwh$Dx_prob[plwh$ID==potential]=rand_in
      
    } else if (rand_notin<=0.05) {
      
      #otherwise 5% chance a given person is unaware
      plwh$Dx[plwh$ID==potential]=F
      plwh$Dx_prob[plwh$ID==potential]=rand_notin
      
    } else {
      #do not sample
    }
    
    rm(potential,rand_in,rand_notin)
  }
  
  #checks
  #prop.table(table(plwh$Age[plwh$Dx==F]))
  #prop.table(table(plwh$Gender[plwh$Dx==F]))
  #prop.table(table(plwh$Race[plwh$Dx==F]))
  
  
  ### SAMPLE THOSE IN CARE ###
  
  rand_in = runif(n_plwh)
  rand_notin = runif(n_plwh)
  
  #% suppressed in clinic (target=72%): https://sfdph.org/dph/files/reports/RptsHIVAIDS/AnnualReport2021-Red.pdf
  plwh$Clinic = ifelse(plwh$VL<200 & plwh$Dx==T & rand_in<0.9, T, ifelse(plwh$VL>200 & plwh$Dx==T & rand_notin<0.06, T, F))
  plwh$Clinic_prob = ifelse(plwh$VL<200 & plwh$Dx==T & rand_in<0.9, rand_in, ifelse(plwh$VL>200 & plwh$Dx==T & rand_notin<0.06, rand_notin, 0))
  n_incare = sum(plwh$Clinic)
  
  #prop.table(table(plwh$VL[plwh$Clinic==T]<200))
  #describe(plwh$VL[plwh$Clinic==T])
  #describe(plwh$VL[plwh$Clinic==F])
  
  
  ### SAMPLE CLINICS ###
  
  #sample an EHR from population of active patients engaged in care
  #Clinic A: 250 patients, oversampled by male, white, >=45 from the source
  #Clinic B: 250 patients, oversampled by male, white, >=45 from those dx 
  #Clinic C: 250 patients, oversampled by male, white, >=45 from those in care
  
  plwh$Clinic_A = F
  plwh$Clinic_A_prob = 0
  plwh$Clinic_B = F
  plwh$Clinic_B_prob = 0
  plwh$Clinic_C = F
  plwh$Clinic_C_prob = 0
  n_clinic_a = 250
  n_clinic_b = 250
  n_clinic_c = 250
  
  #sample clinic a
  while (sum(plwh$Clinic_A)<n_clinic_a) {
    
    #obtain potential from the source population
    potential = sample(plwh$ID[plwh$Clinic_A==F], 1)
    
    rand_in = runif(1)
    rand_notin = runif(1)
    
    #oversample by male, white, >=45 (no initial target)
    if (plwh$Age[plwh$ID==potential]>=2 && plwh$Gender[plwh$ID==potential]=="M" && plwh$Race[plwh$ID==potential]=="W" && rand_in<=0.95) {
      
      plwh$Clinic_A[plwh$ID==potential]=T
      plwh$Clinic_A_prob[plwh$ID==potential]=rand_in
      
    } else if (rand_notin<=0.05) {
      
      plwh$Clinic_A[plwh$ID==potential]=T
      plwh$Clinic_A_prob[plwh$ID==potential]=rand_notin
      
    } else {
      #do not sample
    }
    
    rm(potential,rand_in,rand_notin)
  }
  
  #checks
  #prop.table(table(plwh$Age[plwh$Clinic_A==T]))
  ## 0     1     2     3 
  ## 0.056 0.068 0.360 0.516 
  #prop.table(table(plwh$Gender[plwh$Clinic_A==T]))
  ## F     M 
  ## 0.052 0.948 
  #prop.table(table(plwh$Race[plwh$Clinic_A==T]))
  ## B     H     W 
  ## 0.144 0.068 0.788 
  
  #sample clinic b
  while (sum(plwh$Clinic_B)<n_clinic_b) {
    
    #obtain potential from the population of those aware of their status
    potential = sample(plwh$ID[plwh$Dx==T & plwh$Clinic_B==F], 1)
    
    rand_in = runif(1)
    rand_notin = runif(1)
    
    #oversample by male, white, >=45 (target balance of clinic A)
    if (plwh$Age[plwh$ID==potential]>=2 && plwh$Gender[plwh$ID==potential]=="M" && plwh$Race[plwh$ID==potential]=="W" && rand_in<=0.95) {
      
      plwh$Clinic_B[plwh$ID==potential]=T
      plwh$Clinic_B_prob[plwh$ID==potential]=rand_in
      
    } else if (rand_notin<=0.05) {
      
      plwh$Clinic_B[plwh$ID==potential]=T
      plwh$Clinic_B_prob[plwh$ID==potential]=rand_notin
      
      
    } else {
      #do not sample
    }
    
    rm(potential,rand_in,rand_notin)
  }
  
  #checks
  #prop.table(table(plwh$Age[plwh$Clinic_B==T]))
  #prop.table(table(plwh$Gender[plwh$Clinic_B==T]))
  #prop.table(table(plwh$Race[plwh$Clinic_B==T]))
  
  #sample clinic c
  while (sum(plwh$Clinic_C)<n_clinic_c) {
    
    #obtain potential from the population of those in care
    potential = sample(plwh$ID[plwh$Clinic==T & plwh$Clinic_C==F], 1)
    
    rand_in = runif(1)
    rand_notin = runif(1)
    
    #oversample by male, white, >=45 (target balance of clinic A)
    if (plwh$Age[plwh$ID==potential]>=2 && plwh$Gender[plwh$ID==potential]=="M" && plwh$Race[plwh$ID==potential]=="W" && rand_in<=0.95) {
      
      plwh$Clinic_C[plwh$ID==potential]=T
      plwh$Clinic_C_prob[plwh$ID==potential]=rand_in
      
    } else if (rand_notin<=0.05) {
      
      plwh$Clinic_C[plwh$ID==potential]=T
      plwh$Clinic_C_prob[plwh$ID==potential]=rand_notin
      
    } else {
      #do not sample
    }
    
    rm(potential,rand_in,rand_notin)
  }
  
  #checks
  #prop.table(table(plwh$Age[plwh$Clinic_C==T]))
  #prop.table(table(plwh$Gender[plwh$Clinic_C==T]))
  #prop.table(table(plwh$Race[plwh$Clinic_C==T]))
  
  # #checks
  # sum(plwh$Clinic_A)
  # sum(plwh$Clinic_B)
  # sum(plwh$Clinic_C)
  # sum(plwh$Dx)
  # sum(plwh$Clinic)
  
  
  ### WEIGHTS ###
  
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
  clinic_b = NA
  clinic_c = NA
  clinic_a_biased = NA
  clinic_b_biased = NA
  clinic_c_biased = NA
  clinic_a_bayes_n1 = NA
  clinic_a_bayes_n2 = NA
  clinic_a_bayes_n3 = NA
  clinic_b_bayes_n1 = NA
  clinic_b_bayes_n2 = NA
  clinic_b_bayes_n3 = NA
  clinic_c_bayes_n1 = NA
  clinic_c_bayes_n2 = NA
  clinic_c_bayes_n3 = NA
  clinic_c_bayes_m1 = NA
  
  for (j in 1:1000) {
    
    cat("\n\n************** ","Simulation: ",i," Weight: ",j," **************\n",sep="")
    
    #weights by race
    wt_pi_white = sum(plwh$Race=="W" & plwh$Clinic==T) / n_pop
    wt_hiv_pop_white = rbinom(1, n_pop, wt_pi_white)
    while (sum(wt_hiv_pop_white < c(sum(plwh$Race=="W" & plwh$Clinic_A==T),sum(plwh$Race=="W" & plwh$Clinic_B==T),sum(plwh$Race=="W" & plwh$Clinic_C==T),sum(plwh$Race=="W" & plwh$Clinic_D==T),sum(plwh$Race=="W" & plwh$Clinic_E==T),sum(plwh$Race=="W" & plwh$Clinic_F==T),sum(plwh$Race=="W" & plwh$Clinic_G==T),sum(plwh$Race=="W" & plwh$Clinic_H==T))) > 0) { wt_hiv_pop_white = rbinom(1, n_pop, wt_pi_white) } #resample if N<S
    wt_a_white = 1/rbeta(1, (sum(plwh$Race=="W" & plwh$Clinic_A==T)+1), (wt_hiv_pop_white + 1 - sum(plwh$Race=="W" & plwh$Clinic_A==T)))
    wt_b_white = 1/rbeta(1, (sum(plwh$Race=="W" & plwh$Clinic_B==T)+1), (wt_hiv_pop_white + 1 - sum(plwh$Race=="W" & plwh$Clinic_B==T)))
    wt_c_white = 1/rbeta(1, (sum(plwh$Race=="W" & plwh$Clinic_C==T)+1), (wt_hiv_pop_white + 1 - sum(plwh$Race=="W" & plwh$Clinic_C==T)))
    
    wt_pi_black = sum(plwh$Race=="B" & plwh$Clinic==T) / n_pop
    wt_hiv_pop_black = rbinom(1, n_pop, wt_pi_black)
    while (sum(wt_hiv_pop_black < c(sum(plwh$Race=="B" & plwh$Clinic_A==T),sum(plwh$Race=="B" & plwh$Clinic_B==T),sum(plwh$Race=="B" & plwh$Clinic_C==T),sum(plwh$Race=="B" & plwh$Clinic_D==T),sum(plwh$Race=="B" & plwh$Clinic_E==T),sum(plwh$Race=="B" & plwh$Clinic_F==T),sum(plwh$Race=="B" & plwh$Clinic_G==T),sum(plwh$Race=="B" & plwh$Clinic_H==T))) > 0) { wt_hiv_pop_black = rbinom(1, n_pop, wt_pi_black) } #resample if N<S
    wt_a_black = 1/rbeta(1, (sum(plwh$Race=="B" & plwh$Clinic_A==T)+1), (wt_hiv_pop_black + 1 - sum(plwh$Race=="B" & plwh$Clinic_A==T)))
    wt_b_black = 1/rbeta(1, (sum(plwh$Race=="B" & plwh$Clinic_B==T)+1), (wt_hiv_pop_black + 1 - sum(plwh$Race=="B" & plwh$Clinic_B==T)))
    wt_c_black = 1/rbeta(1, (sum(plwh$Race=="B" & plwh$Clinic_C==T)+1), (wt_hiv_pop_black + 1 - sum(plwh$Race=="B" & plwh$Clinic_C==T)))
    
    wt_pi_hispanic = sum(plwh$Race=="H" & plwh$Clinic==T) / n_pop
    wt_hiv_pop_hispanic = rbinom(1, n_pop, wt_pi_hispanic)
    while (sum(wt_hiv_pop_hispanic < c(sum(plwh$Race=="H" & plwh$Clinic_A==T),sum(plwh$Race=="H" & plwh$Clinic_B==T),sum(plwh$Race=="H" & plwh$Clinic_C==T),sum(plwh$Race=="H" & plwh$Clinic_D==T),sum(plwh$Race=="H" & plwh$Clinic_E==T),sum(plwh$Race=="H" & plwh$Clinic_F==T),sum(plwh$Race=="H" & plwh$Clinic_G==T),sum(plwh$Race=="H" & plwh$Clinic_H==T))) > 0) { wt_hiv_pop_hispanic = rbinom(1, n_pop, wt_pi_hispanic) } #resample if N<S
    wt_a_hispanic = 1/rbeta(1, (sum(plwh$Race=="H" & plwh$Clinic_A==T)+1), (wt_hiv_pop_hispanic + 1 - sum(plwh$Race=="H" & plwh$Clinic_A==T)))
    wt_b_hispanic = 1/rbeta(1, (sum(plwh$Race=="H" & plwh$Clinic_B==T)+1), (wt_hiv_pop_hispanic + 1 - sum(plwh$Race=="H" & plwh$Clinic_B==T)))
    wt_c_hispanic = 1/rbeta(1, (sum(plwh$Race=="H" & plwh$Clinic_C==T)+1), (wt_hiv_pop_hispanic + 1 - sum(plwh$Race=="H" & plwh$Clinic_C==T)))
    
    #weights by gender
    wt_pi_male = sum(plwh$Gender=="M" & plwh$Clinic==T) / n_pop
    wt_hiv_pop_male = rbinom(1, n_pop, wt_pi_male)
    while (sum(wt_hiv_pop_male < c(sum(plwh$Gender=="M" & plwh$Clinic_A==T),sum(plwh$Gender=="M" & plwh$Clinic_B==T),sum(plwh$Gender=="M" & plwh$Clinic_C==T),sum(plwh$Gender=="M" & plwh$Clinic_D==T),sum(plwh$Gender=="M" & plwh$Clinic_E==T),sum(plwh$Gender=="M" & plwh$Clinic_F==T),sum(plwh$Gender=="M" & plwh$Clinic_G==T),sum(plwh$Gender=="M" & plwh$Clinic_H==T))) > 0) { wt_hiv_pop_male = rbinom(1, n_pop, wt_pi_male) } #resample if N<S
    wt_a_male = 1/rbeta(1, (sum(plwh$Gender=="M" & plwh$Clinic_A==T)+1), (wt_hiv_pop_male + 1 - sum(plwh$Gender=="M" & plwh$Clinic_A==T)))
    wt_b_male = 1/rbeta(1, (sum(plwh$Gender=="M" & plwh$Clinic_B==T)+1), (wt_hiv_pop_male + 1 - sum(plwh$Gender=="M" & plwh$Clinic_B==T)))
    wt_c_male = 1/rbeta(1, (sum(plwh$Gender=="M" & plwh$Clinic_C==T)+1), (wt_hiv_pop_male + 1 - sum(plwh$Gender=="M" & plwh$Clinic_C==T)))
    
    wt_pi_female = sum(plwh$Gender=="F" & plwh$Clinic==T) / n_pop
    wt_hiv_pop_female = rbinom(1, n_pop, wt_pi_female)
    while (sum(wt_hiv_pop_female < c(sum(plwh$Gender=="F" & plwh$Clinic_A==T),sum(plwh$Gender=="F" & plwh$Clinic_B==T),sum(plwh$Gender=="F" & plwh$Clinic_C==T),sum(plwh$Gender=="F" & plwh$Clinic_D==T),sum(plwh$Gender=="F" & plwh$Clinic_E==T),sum(plwh$Gender=="F" & plwh$Clinic_F==T),sum(plwh$Gender=="F" & plwh$Clinic_G==T),sum(plwh$Gender=="F" & plwh$Clinic_H==T))) > 0) { wt_hiv_pop_female = rbinom(1, n_pop, wt_pi_female) } #resample if N<S
    wt_a_female = 1/rbeta(1, (sum(plwh$Gender=="F" & plwh$Clinic_A==T)+1), (wt_hiv_pop_female + 1 - sum(plwh$Gender=="F" & plwh$Clinic_A==T)))
    wt_b_female = 1/rbeta(1, (sum(plwh$Gender=="F" & plwh$Clinic_B==T)+1), (wt_hiv_pop_female + 1 - sum(plwh$Gender=="F" & plwh$Clinic_B==T)))
    wt_c_female = 1/rbeta(1, (sum(plwh$Gender=="F" & plwh$Clinic_C==T)+1), (wt_hiv_pop_female + 1 - sum(plwh$Gender=="F" & plwh$Clinic_C==T)))
    
    #weights by age category
    wt_pi_35less = sum(plwh$Age==0 & plwh$Clinic==T) / n_pop
    wt_hiv_pop_35less = rbinom(1, n_pop, wt_pi_35less)
    while (sum(wt_hiv_pop_35less < c(sum(plwh$Age==0 & plwh$Clinic_A==T),sum(plwh$Age==0 & plwh$Clinic_B==T),sum(plwh$Age==0 & plwh$Clinic_C==T),sum(plwh$Age==0 & plwh$Clinic_D==T),sum(plwh$Age==0 & plwh$Clinic_E==T),sum(plwh$Age==0 & plwh$Clinic_F==T),sum(plwh$Age==0 & plwh$Clinic_G==T),sum(plwh$Age==0 & plwh$Clinic_H==T))) > 0) { wt_hiv_pop_35less = rbinom(1, n_pop, wt_pi_35less) } #resample if N<S
    wt_a_35less = 1/rbeta(1, (sum(plwh$Age==0 & plwh$Clinic_A==T)+1), (wt_hiv_pop_35less + 1 - sum(plwh$Age==0 & plwh$Clinic_A==T)))
    wt_b_35less = 1/rbeta(1, (sum(plwh$Age==0 & plwh$Clinic_B==T)+1), (wt_hiv_pop_35less + 1 - sum(plwh$Age==0 & plwh$Clinic_B==T)))
    wt_c_35less = 1/rbeta(1, (sum(plwh$Age==0 & plwh$Clinic_C==T)+1), (wt_hiv_pop_35less + 1 - sum(plwh$Age==0 & plwh$Clinic_C==T)))
    
    wt_pi_35to44 = sum(plwh$Age==1 & plwh$Clinic==T) / n_pop
    wt_hiv_pop_35to44 = rbinom(1, n_pop, wt_pi_35to44)
    while (sum(wt_hiv_pop_35to44 < c(sum(plwh$Age==1 & plwh$Clinic_A==T),sum(plwh$Age==1 & plwh$Clinic_B==T),sum(plwh$Age==1 & plwh$Clinic_C==T),sum(plwh$Age==1 & plwh$Clinic_D==T),sum(plwh$Age==1 & plwh$Clinic_E==T),sum(plwh$Age==1 & plwh$Clinic_F==T),sum(plwh$Age==1 & plwh$Clinic_G==T),sum(plwh$Age==1 & plwh$Clinic_H==T))) > 0) { wt_hiv_pop_35to44 = rbinom(1, n_pop, wt_pi_35to44) } #resample if N<S
    wt_a_35to44 = 1/rbeta(1, (sum(plwh$Age==1 & plwh$Clinic_A==T)+1), (wt_hiv_pop_35to44 + 1 - sum(plwh$Age==1 & plwh$Clinic_A==T)))
    wt_b_35to44 = 1/rbeta(1, (sum(plwh$Age==1 & plwh$Clinic_B==T)+1), (wt_hiv_pop_35to44 + 1 - sum(plwh$Age==1 & plwh$Clinic_B==T)))
    wt_c_35to44 = 1/rbeta(1, (sum(plwh$Age==1 & plwh$Clinic_C==T)+1), (wt_hiv_pop_35to44 + 1 - sum(plwh$Age==1 & plwh$Clinic_C==T)))
    
    wt_pi_45to54 = sum(plwh$Age==2 & plwh$Clinic==T) / n_pop
    wt_hiv_pop_45to54 = rbinom(1, n_pop, wt_pi_45to54)
    while (sum(wt_hiv_pop_45to54 < c(sum(plwh$Age==2 & plwh$Clinic_A==T),sum(plwh$Age==2 & plwh$Clinic_B==T),sum(plwh$Age==2 & plwh$Clinic_C==T),sum(plwh$Age==2 & plwh$Clinic_D==T),sum(plwh$Age==2 & plwh$Clinic_E==T),sum(plwh$Age==2 & plwh$Clinic_F==T),sum(plwh$Age==2 & plwh$Clinic_G==T),sum(plwh$Age==2 & plwh$Clinic_H==T))) > 0) { wt_hiv_pop_45to54 = rbinom(1, n_pop, wt_pi_45to54) } #resample if N<S
    wt_a_45to54 = 1/rbeta(1, (sum(plwh$Age==2 & plwh$Clinic_A==T)+1), (wt_hiv_pop_45to54 + 1 - sum(plwh$Age==2 & plwh$Clinic_A==T)))
    wt_b_45to54 = 1/rbeta(1, (sum(plwh$Age==2 & plwh$Clinic_B==T)+1), (wt_hiv_pop_45to54 + 1 - sum(plwh$Age==2 & plwh$Clinic_B==T)))
    wt_c_45to54 = 1/rbeta(1, (sum(plwh$Age==2 & plwh$Clinic_C==T)+1), (wt_hiv_pop_45to54 + 1 - sum(plwh$Age==2 & plwh$Clinic_C==T)))
    
    wt_pi_55more = sum(plwh$Age==3 & plwh$Clinic==T) / n_pop
    wt_hiv_pop_55more = rbinom(1, n_pop, wt_pi_55more)
    while (sum(wt_hiv_pop_55more < c(sum(plwh$Age==3 & plwh$Clinic_A==T),sum(plwh$Age==3 & plwh$Clinic_B==T),sum(plwh$Age==3 & plwh$Clinic_C==T),sum(plwh$Age==3 & plwh$Clinic_D==T),sum(plwh$Age==3 & plwh$Clinic_E==T),sum(plwh$Age==3 & plwh$Clinic_F==T),sum(plwh$Age==3 & plwh$Clinic_G==T),sum(plwh$Age==3 & plwh$Clinic_H==T))) > 0) { wt_hiv_pop_55more = rbinom(1, n_pop, wt_pi_55more) } #resample if N<S
    wt_a_55more = 1/rbeta(1, (sum(plwh$Age==3 & plwh$Clinic_A==T)+1), (wt_hiv_pop_55more + 1 - sum(plwh$Age==3 & plwh$Clinic_A==T)))
    wt_b_55more = 1/rbeta(1, (sum(plwh$Age==3 & plwh$Clinic_B==T)+1), (wt_hiv_pop_55more + 1 - sum(plwh$Age==3 & plwh$Clinic_B==T)))
    wt_c_55more = 1/rbeta(1, (sum(plwh$Age==3 & plwh$Clinic_C==T)+1), (wt_hiv_pop_55more + 1 - sum(plwh$Age==3 & plwh$Clinic_C==T)))
    
    #assign strata specific weights to population
    plwh$Weight_a = ifelse(plwh$Clinic_A==T & plwh$Age==0, wt_a_35less, ifelse(plwh$Clinic_A==T & plwh$Age==1, wt_a_35to44, ifelse(plwh$Clinic_A==T & plwh$Age==2, wt_a_45to54, ifelse(plwh$Clinic_A==T & plwh$Age==3, wt_a_55more, NA)))) *
      ifelse(plwh$Clinic_A==T & plwh$Gender=="M", wt_a_male, ifelse(plwh$Clinic_A==T & plwh$Gender=="F", wt_a_female, NA)) *
      ifelse(plwh$Clinic_A==T & plwh$Race=="W", wt_a_white, ifelse(plwh$Clinic_A==T & plwh$Race=="B", wt_a_black, ifelse(plwh$Clinic_A==T & plwh$Race=="H", wt_a_hispanic, NA)))
    
    plwh$Weight_b = ifelse(plwh$Clinic_B==T & plwh$Age==0, wt_b_35less, ifelse(plwh$Clinic_B==T & plwh$Age==1, wt_b_35to44, ifelse(plwh$Clinic_B==T & plwh$Age==2, wt_b_45to54, ifelse(plwh$Clinic_B==T & plwh$Age==3, wt_b_55more, NA)))) *
      ifelse(plwh$Clinic_B==T & plwh$Gender=="M", wt_b_male, ifelse(plwh$Clinic_B==T & plwh$Gender=="F", wt_b_female, NA)) *
      ifelse(plwh$Clinic_B==T & plwh$Race=="W", wt_b_white, ifelse(plwh$Clinic_B==T & plwh$Race=="B", wt_b_black, ifelse(plwh$Clinic_B==T & plwh$Race=="H", wt_b_hispanic, NA)))
    
    plwh$Weight_c = ifelse(plwh$Clinic_C==T & plwh$Age==0, wt_c_35less, ifelse(plwh$Clinic_C==T & plwh$Age==1, wt_c_35to44, ifelse(plwh$Clinic_C==T & plwh$Age==2, wt_c_45to54, ifelse(plwh$Clinic_C==T & plwh$Age==3, wt_c_55more, NA)))) *
      ifelse(plwh$Clinic_C==T & plwh$Gender=="M", wt_c_male, ifelse(plwh$Clinic_C==T & plwh$Gender=="F", wt_c_female, NA)) *
      ifelse(plwh$Clinic_C==T & plwh$Race=="W", wt_c_white, ifelse(plwh$Clinic_C==T & plwh$Race=="B", wt_c_black, ifelse(plwh$Clinic_C==T & plwh$Race=="H", wt_c_hispanic, NA)))
    
    #weight trimming parameter: see trim function in WeightIt
    wt_trim = 0.98
    
    #weighted clinics using geometric mean: exp(mean(log(X)))
    clinic_a = c(clinic_a, exp(weighted.mean(log(plwh$VL[plwh$Clinic_A==T]), suppressMessages(trim(plwh$Weight_a[plwh$Clinic_A==T], at=wt_trim, lower=T)))))
    clinic_b = c(clinic_b, exp(weighted.mean(log(plwh$VL[plwh$Clinic_B==T]), suppressMessages(trim(plwh$Weight_b[plwh$Clinic_B==T], at=wt_trim, lower=T)))))
    clinic_c = c(clinic_c, exp(weighted.mean(log(plwh$VL[plwh$Clinic_C==T]), suppressMessages(trim(plwh$Weight_c[plwh$Clinic_C==T], at=wt_trim, lower=T)))))
    
    #weighted clinics, weights misspecified using geometric mean 
    #misspecification form: logit_pO=log(pT/(1-pT))+b*logVL, where pT is the observed probability (inverse of weight) and b is the bias factor (0=no bias)
    bias_factor = 0.1
    sample_probs_observed = 1/plwh$Weight_a[plwh$Clinic_A==T]
    biased_probs = log(sample_probs_observed/(1-sample_probs_observed))+bias_factor*log(plwh$VL[plwh$Clinic_A==T])
    sample_probs_biased = exp(biased_probs)/(1+exp(biased_probs))
    clinic_a_biased = c(clinic_a_biased, exp(weighted.mean(log(plwh$VL[plwh$Clinic_A==T]), suppressMessages(trim((1/sample_probs_biased), at=wt_trim, lower=T)))))
    
    sample_probs_observed = 1/plwh$Weight_b[plwh$Clinic_B==T]
    biased_probs = log(sample_probs_observed/(1-sample_probs_observed))+bias_factor*log(plwh$VL[plwh$Clinic_B==T])
    sample_probs_biased =exp(biased_probs)/(1+exp(biased_probs))
    clinic_b_biased = c(clinic_b_biased, exp(weighted.mean(log(plwh$VL[plwh$Clinic_B==T]), suppressMessages(trim((1/sample_probs_biased), at=wt_trim, lower=T)))))
    
    sample_probs_observed = 1/plwh$Weight_c[plwh$Clinic_C==T]
    biased_probs = log(sample_probs_observed/(1-sample_probs_observed))+bias_factor*log(plwh$VL[plwh$Clinic_C==T])
    sample_probs_biased =exp(biased_probs)/(1+exp(biased_probs))
    clinic_c_biased = c(clinic_c_biased, exp(weighted.mean(log(plwh$VL[plwh$Clinic_C==T]), suppressMessages(trim((1/sample_probs_biased), at=wt_trim, lower=T)))))
    
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
    
    clinic_xbar = log10(exp(weighted.mean(log(plwh$VL[plwh$Clinic_B==T]), suppressMessages(trim(plwh$Weight_b[plwh$Clinic_B==T], at=wt_trim, lower=T)))))
    clinic_n = sum(suppressMessages(trim(plwh$Weight_b[plwh$Clinic_B==T], at=wt_trim, lower=T)))
    n1 = clinic_n*0.25 #sample size relative to clinic size
    n2 = clinic_n*0.5 #sample size relative to clinic size
    n3 = clinic_n*2 #sample size relative to clinic size
    clinic_b_bayes_n1 = c(clinic_b_bayes_n1, 10^(((n1*m0) + (clinic_n*clinic_xbar))/(n1 + clinic_n)))
    clinic_b_bayes_n2 = c(clinic_b_bayes_n2, 10^(((n2*m0) + (clinic_n*clinic_xbar))/(n2 + clinic_n)))
    clinic_b_bayes_n3 = c(clinic_b_bayes_n3, 10^(((n3*m0) + (clinic_n*clinic_xbar))/(n3 + clinic_n)))
    
    clinic_xbar = log10(exp(weighted.mean(log(plwh$VL[plwh$Clinic_C==T]), suppressMessages(trim(plwh$Weight_c[plwh$Clinic_C==T], at=wt_trim, lower=T)))))
    clinic_n = sum(suppressMessages(trim(plwh$Weight_c[plwh$Clinic_C==T], at=wt_trim, lower=T)))
    n1 = clinic_n*0.25 #sample size relative to clinic size
    n2 = clinic_n*0.5 #sample size relative to clinic size
    n3 = clinic_n*2 #sample size relative to clinic size
    clinic_c_bayes_n1 = c(clinic_c_bayes_n1, 10^(((n1*m0) + (clinic_n*clinic_xbar))/(n1 + clinic_n)))
    clinic_c_bayes_n2 = c(clinic_c_bayes_n2, 10^(((n2*m0) + (clinic_n*clinic_xbar))/(n2 + clinic_n)))
    clinic_c_bayes_n3 = c(clinic_c_bayes_n3, 10^(((n3*m0) + (clinic_n*clinic_xbar))/(n3 + clinic_n)))
    
    #clinic C special case, adjust prior mean m0 to use weighted clinic B
    m0 = log10(exp(weighted.mean(log(plwh$VL[plwh$Clinic_B==T]), suppressMessages(trim(plwh$Weight_b[plwh$Clinic_B==T], at=wt_trim, lower=T)))))
    n1 = clinic_n*0.25
    clinic_c_bayes_m1 = c(clinic_c_bayes_m1, 10^(((n1*m0) + (clinic_n*clinic_xbar))/(n1 + clinic_n)))
    
    rm(m0,n1,n2,n3,clinic_xbar,clinic_n,wt_trim,sample_probs_observed,biased_probs,sample_probs_biased,bias_factor)
  }
  rm(j)
  
  
  ### SUMMARY MEASURES ### 
  
  #population characteristics
  sim_pop$Population_age0[i] = sum(plwh$Age==0)
  sim_pop$Population_age1[i] = sum(plwh$Age==1)
  sim_pop$Population_age2[i] = sum(plwh$Age==2)
  sim_pop$Population_age3[i] = sum(plwh$Age==3)
  sim_pop$Population_genderM[i] = sum(plwh$Gender=="M")
  sim_pop$Population_genderF[i] = sum(plwh$Gender=="F")
  sim_pop$Population_raceW[i] = sum(plwh$Race=="W")
  sim_pop$Population_raceB[i] = sum(plwh$Race=="B")
  sim_pop$Population_raceH[i] = sum(plwh$Race=="H")
  
  sim_pop$Clinic_A_age0[i] = sum(plwh$Age[plwh$Clinic_A==T]==0)
  sim_pop$Clinic_A_age1[i] = sum(plwh$Age[plwh$Clinic_A==T]==1)
  sim_pop$Clinic_A_age2[i] = sum(plwh$Age[plwh$Clinic_A==T]==2)
  sim_pop$Clinic_A_age3[i] = sum(plwh$Age[plwh$Clinic_A==T]==3)
  sim_pop$Clinic_A_genderM[i] = sum(plwh$Gender[plwh$Clinic_A==T]=="M")
  sim_pop$Clinic_A_genderF[i] = sum(plwh$Gender[plwh$Clinic_A==T]=="F")
  sim_pop$Clinic_A_raceW[i] = sum(plwh$Race[plwh$Clinic_A==T]=="W")
  sim_pop$Clinic_A_raceB[i] = sum(plwh$Race[plwh$Clinic_A==T]=="B")
  sim_pop$Clinic_A_raceH[i] = sum(plwh$Race[plwh$Clinic_A==T]=="H")
  
  sim_pop$Clinic_B_age0[i] = sum(plwh$Age[plwh$Clinic_B==T]==0)
  sim_pop$Clinic_B_age1[i] = sum(plwh$Age[plwh$Clinic_B==T]==1)
  sim_pop$Clinic_B_age2[i] = sum(plwh$Age[plwh$Clinic_B==T]==2)
  sim_pop$Clinic_B_age3[i] = sum(plwh$Age[plwh$Clinic_B==T]==3)
  sim_pop$Clinic_B_genderM[i] = sum(plwh$Gender[plwh$Clinic_B==T]=="M")
  sim_pop$Clinic_B_genderF[i] = sum(plwh$Gender[plwh$Clinic_B==T]=="F")
  sim_pop$Clinic_B_raceW[i] = sum(plwh$Race[plwh$Clinic_B==T]=="W")
  sim_pop$Clinic_B_raceB[i] = sum(plwh$Race[plwh$Clinic_B==T]=="B")
  sim_pop$Clinic_B_raceH[i] = sum(plwh$Race[plwh$Clinic_B==T]=="H")
  
  sim_pop$Clinic_C_age0[i] = sum(plwh$Age[plwh$Clinic_C==T]==0)
  sim_pop$Clinic_C_age1[i] = sum(plwh$Age[plwh$Clinic_C==T]==1)
  sim_pop$Clinic_C_age2[i] = sum(plwh$Age[plwh$Clinic_C==T]==2)
  sim_pop$Clinic_C_age3[i] = sum(plwh$Age[plwh$Clinic_C==T]==3)
  sim_pop$Clinic_C_genderM[i] = sum(plwh$Gender[plwh$Clinic_C==T]=="M")
  sim_pop$Clinic_C_genderF[i] = sum(plwh$Gender[plwh$Clinic_C==T]=="F")
  sim_pop$Clinic_C_raceW[i] = sum(plwh$Race[plwh$Clinic_C==T]=="W")
  sim_pop$Clinic_C_raceB[i] = sum(plwh$Race[plwh$Clinic_C==T]=="B")
  sim_pop$Clinic_C_raceH[i] = sum(plwh$Race[plwh$Clinic_C==T]=="H")
  
  #true VL using geometric mean for all VLs: exp(mean(log(X)))
  sim_results$VL_population[i] = exp(mean(log(plwh$VL)))
  
  #undiagnosed
  sim_results$VL_population_nodx[i] = exp(mean(log(plwh$VL[plwh$Dx==F])))
  
  #diagnosed
  sim_results$VL_population_dx[i] = exp(mean(log(plwh$VL[plwh$Dx==T])))
  
  #not in care
  sim_results$VL_population_notincare[i] = exp(mean(log(plwh$VL[plwh$Clinic==F])))
  
  #in care
  sim_results$VL_population_incare[i] = exp(mean(log(plwh$VL[plwh$Clinic==T])))
  
  #clinic VLs
  sim_results$VL_clinic_A[i] = exp(mean(log(plwh$VL[plwh$Clinic_A==T])))
  sim_results$VL_clinic_B[i] = exp(mean(log(plwh$VL[plwh$Clinic_B==T])))
  sim_results$VL_clinic_C[i] = exp(mean(log(plwh$VL[plwh$Clinic_C==T])))
  
  #weighted clinics
  sim_results$VL_clinic_A_weighted[i] = mean(clinic_a, na.rm=T)
  sim_results$VL_clinic_B_weighted[i] = mean(clinic_b, na.rm=T)
  sim_results$VL_clinic_C_weighted[i] = mean(clinic_c, na.rm=T)
  
  #weighted clinics, misspecified weights
  sim_results$VL_clinic_A_weighted_biased[i] = mean(clinic_a_biased, na.rm=T)
  sim_results$VL_clinic_B_weighted_biased[i] = mean(clinic_b_biased, na.rm=T)
  sim_results$VL_clinic_C_weighted_biased[i] = mean(clinic_c_biased, na.rm=T)
  
  #rmse, weighted clinics
  sim_results$VL_clinic_A_rmse[i] = sqrt(mean((clinic_a - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_B_rmse[i] = sqrt(mean((clinic_b - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_C_rmse[i] = sqrt(mean((clinic_c - exp(mean(log(plwh$VL))))^2, na.rm=T))
  
  #rmse, misspecified weights
  sim_results$VL_clinic_A_rmse_biased[i] = sqrt(mean((clinic_a_biased - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_B_rmse_biased[i] = sqrt(mean((clinic_b_biased - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_C_rmse_biased[i] = sqrt(mean((clinic_c_biased - exp(mean(log(plwh$VL))))^2, na.rm=T))
  
  #bayes, weighted clinics
  sim_results$VL_clinic_A_weighted_bayes_n1[i] = mean(clinic_a_bayes_n1, na.rm=T)
  sim_results$VL_clinic_A_weighted_bayes_n2[i] = mean(clinic_a_bayes_n2, na.rm=T)
  sim_results$VL_clinic_A_weighted_bayes_n3[i] = mean(clinic_a_bayes_n3, na.rm=T)
  sim_results$VL_clinic_B_weighted_bayes_n1[i] = mean(clinic_b_bayes_n1, na.rm=T)
  sim_results$VL_clinic_B_weighted_bayes_n2[i] = mean(clinic_b_bayes_n2, na.rm=T)
  sim_results$VL_clinic_B_weighted_bayes_n3[i] = mean(clinic_b_bayes_n3, na.rm=T)
  sim_results$VL_clinic_C_weighted_bayes_n1[i] = mean(clinic_c_bayes_n1, na.rm=T)
  sim_results$VL_clinic_C_weighted_bayes_n2[i] = mean(clinic_c_bayes_n2, na.rm=T)
  sim_results$VL_clinic_C_weighted_bayes_n3[i] = mean(clinic_c_bayes_n3, na.rm=T)
  sim_results$VL_clinic_C_weighted_bayes_m1[i] = mean(clinic_c_bayes_m1, na.rm=T)
  
  #rmse, bayes weighted clinics
  sim_results$VL_clinic_A_rmse_bayes_n1[i] = sqrt(mean((clinic_a_bayes_n1 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_A_rmse_bayes_n2[i] = sqrt(mean((clinic_a_bayes_n2 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_A_rmse_bayes_n3[i] = sqrt(mean((clinic_a_bayes_n3 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_B_rmse_bayes_n1[i] = sqrt(mean((clinic_b_bayes_n1 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_B_rmse_bayes_n2[i] = sqrt(mean((clinic_b_bayes_n2 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_B_rmse_bayes_n3[i] = sqrt(mean((clinic_b_bayes_n3 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_C_rmse_bayes_n1[i] = sqrt(mean((clinic_c_bayes_n1 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_C_rmse_bayes_n2[i] = sqrt(mean((clinic_c_bayes_n2 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_C_rmse_bayes_n3[i] = sqrt(mean((clinic_c_bayes_n3 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  sim_results$VL_clinic_C_rmse_bayes_m1[i] = sqrt(mean((clinic_c_bayes_m1 - exp(mean(log(plwh$VL))))^2, na.rm=T))
  
}
rm(i)

#save
save.image("simulation.RData")


### VISUALIZE RESULTS ###

#load
load("simulation.RData")

#naive to weighted, linear scale
boxplot((sim_results$VL_clinic_A), at=2, xlim=c(0.7, 17.3), ylim=range(0, 20000), xaxt = "n", ylab="Mean Viral Load (copies/mL)")
boxplot((sim_results$VL_clinic_B), at=7, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C), at=12, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_A_weighted), at=3, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_B_weighted), at=8, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_weighted), at=13, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_A_weighted_bayes_n1), at=4, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_B_weighted_bayes_n1), at=9, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_weighted_bayes_n1), at=14, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_A_weighted_bayes_n2), at=5, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_B_weighted_bayes_n2), at=10, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_weighted_bayes_n2), at=15, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_A_weighted_bayes_n3), at=6, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_B_weighted_bayes_n3), at=11, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_weighted_bayes_n3), at=16, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_weighted_bayes_m1), at=17, xaxt="n", add=T)
points(1,(mean(sim_results$VL_population)), pch=18, cex=2)
lines(x=c(0,20), y=rep((mean(sim_results$VL_population)),2), lty=2)
polygon(x=c(1.5,1.5,6.5,6.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
polygon(x=c(11.5,11.5,17.5,17.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
axis(1, at=1:17, labels=c("Population",rep(c("Observed","Weighted",expression(Bayesian^1),expression(Bayesian^2),expression(Bayesian^3)),3),expression(Bayesian^4)), tick=T, las=3, cex.axis=1)
axis(3, at=c(4,9,14.5), labels=c("Clinic A\n(all plwh)","Clinic B\n(diagnosed)","Clinic C\n(in care)"), tick=F, cex.axis=0.8)

#naive to weighted, log10 scale
boxplot(log10(sim_results$VL_clinic_A), at=2, xlim=c(0.7, 17.3), ylim=range(1.8, 4.4), xaxt = "n", ylab="Mean Viral Load (log10 copies/mL)")
boxplot(log10(sim_results$VL_clinic_B), at=7, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C), at=12, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_A_weighted), at=3, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_B_weighted), at=8, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_weighted), at=13, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_A_weighted_bayes_n1), at=4, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_B_weighted_bayes_n1), at=9, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_weighted_bayes_n1), at=14, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_A_weighted_bayes_n2), at=5, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_B_weighted_bayes_n2), at=10, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_weighted_bayes_n2), at=15, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_A_weighted_bayes_n3), at=6, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_B_weighted_bayes_n3), at=11, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_weighted_bayes_n3), at=16, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_weighted_bayes_m1), at=17, xaxt="n", add=T)
points(1,log10(mean(sim_results$VL_population)), pch=18, cex=2)
lines(x=c(0,20), y=rep(log10(mean(sim_results$VL_population)),2), lty=2)
polygon(x=c(1.5,1.5,6.5,6.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
polygon(x=c(11.5,11.5,17.5,17.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
axis(1, at=1:17, labels=c("Population",rep(c("Observed","Weighted",expression(Bayesian^1),expression(Bayesian^2),expression(Bayesian^3)),3),expression(Bayesian^4)), tick=T, las=3, cex.axis=1)
axis(3, at=c(4,9,14.5), labels=c("Clinic A\n(all plwh)","Clinic B\n(diagnosed)","Clinic C\n(in care)"), tick=F, cex.axis=0.8)

#RMSE, linear scale
boxplot((sim_results$VL_clinic_A_rmse), at=1, xlim=c(0.7, 13.3), ylim=range(0, 20000), xaxt = "n", ylab="RMSE Viral Load (copies/mL)")
boxplot((sim_results$VL_clinic_B_rmse), at=5, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_rmse), at=9, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_A_rmse_bayes_n1), at=2, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_B_rmse_bayes_n1), at=6, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_rmse_bayes_n1), at=10, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_A_rmse_bayes_n2), at=3, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_B_rmse_bayes_n2), at=7, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_rmse_bayes_n2), at=11, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_A_rmse_bayes_n3), at=4, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_B_rmse_bayes_n3), at=8, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_rmse_bayes_n3), at=12, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_rmse_bayes_m1), at=13, xaxt="n", add=T)
axis(1, at=1:13, labels=c(rep(c("Weighted",expression(Bayesian^1),expression(Bayesian^2),expression(Bayesian^3)),3),expression(Bayesian^4)), tick=T, las=3, cex.axis=1)
polygon(x=c(0.5,0.5,4.5,4.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
polygon(x=c(8.5,8.5,13.5,13.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
axis(3, at=c(2.5,6.5,11), labels=c("Clinic A\n(all plwh)","Clinic B\n(diagnosed)","Clinic C\n(in care)"), tick=F, cex.axis=0.8)

#RMSE, log10 scale
boxplot(log10(sim_results$VL_clinic_A_rmse), at=1, xlim=c(0.7, 13.3), ylim=range(1.8, 4.4), xaxt = "n", ylab="RMSE Viral Load (log10 copies/mL)")
boxplot(log10(sim_results$VL_clinic_B_rmse), at=5, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_rmse), at=9, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_A_rmse_bayes_n1), at=2, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_B_rmse_bayes_n1), at=6, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_rmse_bayes_n1), at=10, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_A_rmse_bayes_n2), at=3, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_B_rmse_bayes_n2), at=7, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_rmse_bayes_n2), at=11, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_A_rmse_bayes_n3), at=4, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_B_rmse_bayes_n3), at=8, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_rmse_bayes_n3), at=12, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_rmse_bayes_m1), at=13, xaxt="n", add=T)
axis(1, at=1:13, labels=c(rep(c("Weighted",expression(Bayesian^1),expression(Bayesian^2),expression(Bayesian^3)),3),expression(Bayesian^4)), tick=T, las=3, cex.axis=1)
polygon(x=c(0.5,0.5,4.5,4.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
polygon(x=c(8.5,8.5,13.5,13.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
axis(3, at=c(2.5,6.5,11), labels=c("Clinic A\n(all plwh)","Clinic B\n(diagnosed)","Clinic C\n(in care)"), tick=F, cex.axis=0.8)

#weight misspecification, linear scale
boxplot((sim_results$VL_clinic_A_weighted), at=2, xlim=c(0.7, 7.3), ylim=range(0, 20000), xaxt = "n", ylab="Mean Viral Load (copies/mL)")
boxplot((sim_results$VL_clinic_B_weighted), at=4, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_weighted), at=6, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_A_weighted_biased), at=3, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_B_weighted_biased), at=5, xaxt="n", add=T)
boxplot((sim_results$VL_clinic_C_weighted_biased), at=7, xaxt="n", add=T)
points(1,(mean(sim_results$VL_population)), pch=18, cex=2)
lines(x=c(0,20), y=rep((mean(sim_results$VL_population)),2), lty=2)
polygon(x=c(1.5,1.5,3.5,3.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
polygon(x=c(5.5,5.5,7.5,7.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
axis(1, at=1:7, labels=c("Population",rep(c("Weighted","Biased"),3)), tick=T, las=3, cex.axis=1)
axis(3, at=c(2.5,4.5,6.5), labels=c("Clinic A\n(all plwh)","Clinic B\n(diagnosed)","Clinic C\n(in care)"), tick=F, cex.axis=0.8)

#weight misspecification, log10 scale
boxplot(log10(sim_results$VL_clinic_A_weighted), at=2, xlim=c(0.7, 7.3), ylim=range(1.8, 4.4), xaxt = "n", ylab="Mean Viral Load (log10 copies/mL)")
boxplot(log10(sim_results$VL_clinic_B_weighted), at=4, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_weighted), at=6, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_A_weighted_biased), at=3, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_B_weighted_biased), at=5, xaxt="n", add=T)
boxplot(log10(sim_results$VL_clinic_C_weighted_biased), at=7, xaxt="n", add=T)
points(1,log10(mean(sim_results$VL_population)), pch=18, cex=2)
lines(x=c(0,20), y=rep(log10(mean(sim_results$VL_population)),2), lty=2)
polygon(x=c(1.5,1.5,3.5,3.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
polygon(x=c(5.5,5.5,7.5,7.5),y=c(0,60000,60000,0), col=grey(0.5,0.3), border=NA)
axis(1, at=1:7, labels=c("Population",rep(c("Weighted","Biased"),3)), tick=T, las=3, cex.axis=1)
axis(3, at=c(2.5,4.5,6.5), labels=c("Clinic A\n(all plwh)","Clinic B\n(diagnosed)","Clinic C\n(in care)"), tick=F, cex.axis=0.8)


# ### COMPLEX SURVEY ###
# 
# pop.stratum = as.data.frame(table(plwh$Strata))
# plwh$Stratum = interaction(plwh$Age, plwh$Gender, plwh$Race)
# plwh2 = merge(plwh,pop.stratum,by.x="Stratum", by.y="Var1")
# clinica=plwh2[plwh2$Clinic_A==T,]
# 
# #define survey
# ehr_survey = svydesign(id=~1, strata=~Stratum, data=clinica, fpc=~Freq)
# options(survey.adjust.domain.lonely=TRUE)
# options(survey.lonely.psu="adjust")
# 
# #weighted mean
# svymean(~log(VL), ehr_survey)


### PAPER DESCRIPTIVES ###

#note: these are all based on the individual simulation's GM (actual VL values not saved)
median(log10(sim_results$VL_population))
median((sim_results$VL_population))
quantile((sim_results$VL_population), probs=c(0.025,0.975))

median(log10(sim_results$VL_clinic_A))
median((sim_results$VL_clinic_A))
quantile((sim_results$VL_clinic_A), probs=c(0.025,0.975))
median(log10(sim_results$VL_clinic_B))
median((sim_results$VL_clinic_B))
quantile((sim_results$VL_clinic_B), probs=c(0.025,0.975))
median(log10(sim_results$VL_clinic_C))
median((sim_results$VL_clinic_C))
quantile((sim_results$VL_clinic_C), probs=c(0.025,0.975))

median((sim_results$VL_clinic_A_weighted))
quantile((sim_results$VL_clinic_A_weighted), probs=c(0.025,0.975))
median((sim_results$VL_clinic_B_weighted))
quantile((sim_results$VL_clinic_B_weighted), probs=c(0.025,0.975))
median((sim_results$VL_clinic_C_weighted))
quantile((sim_results$VL_clinic_C_weighted), probs=c(0.025,0.975))

median((sim_results$VL_clinic_C_weighted_bayes_n1))
quantile((sim_results$VL_clinic_C_weighted_bayes_n1), probs=c(0.025,0.975))
median((sim_results$VL_clinic_C_weighted_bayes_n2))
quantile((sim_results$VL_clinic_C_weighted_bayes_n2), probs=c(0.025,0.975))
median((sim_results$VL_clinic_C_weighted_bayes_n3))
quantile((sim_results$VL_clinic_C_weighted_bayes_n3), probs=c(0.025,0.975))
median((sim_results$VL_clinic_C_weighted_bayes_m1))
quantile((sim_results$VL_clinic_C_weighted_bayes_m1), probs=c(0.025,0.975))

median((sim_results$VL_clinic_A_rmse))
quantile((sim_results$VL_clinic_A_rmse), probs=c(0.025,0.975))
median((sim_results$VL_clinic_B_rmse))
quantile((sim_results$VL_clinic_B_rmse), probs=c(0.025,0.975))
median((sim_results$VL_clinic_C_rmse))
quantile((sim_results$VL_clinic_C_rmse), probs=c(0.025,0.975))

median((sim_results$VL_clinic_A_rmse_bayes_n3))
quantile((sim_results$VL_clinic_A_rmse_bayes_n3), probs=c(0.025,0.975))
median((sim_results$VL_clinic_B_rmse_bayes_n3))
quantile((sim_results$VL_clinic_B_rmse_bayes_n3), probs=c(0.025,0.975))
median((sim_results$VL_clinic_C_rmse_bayes_n3))
quantile((sim_results$VL_clinic_C_rmse_bayes_n3), probs=c(0.025,0.975))

#table
median(sim_pop$Population_age0); median(sim_pop$Population_age0)/n_plwh*100
median(sim_pop$Population_age1); median(sim_pop$Population_age1)/n_plwh*100
median(sim_pop$Population_age2); median(sim_pop$Population_age2)/n_plwh*100
median(sim_pop$Population_age3); median(sim_pop$Population_age3)/n_plwh*100
median(sim_pop$Population_genderF); median(sim_pop$Population_genderF)/n_plwh*100
median(sim_pop$Population_genderM); median(sim_pop$Population_genderM)/n_plwh*100
median(sim_pop$Population_raceB); median(sim_pop$Population_raceB)/n_plwh*100
median(sim_pop$Population_raceW); median(sim_pop$Population_raceW)/n_plwh*100
median(sim_pop$Population_raceH); median(sim_pop$Population_raceH)/n_plwh*100

median(sim_pop$Clinic_A_age0); median(sim_pop$Clinic_A_age0)/n_clinic_a*100
median(sim_pop$Clinic_A_age1); median(sim_pop$Clinic_A_age1)/n_clinic_a*100
median(sim_pop$Clinic_A_age2); median(sim_pop$Clinic_A_age2)/n_clinic_a*100
median(sim_pop$Clinic_A_age3); median(sim_pop$Clinic_A_age3)/n_clinic_a*100
median(sim_pop$Clinic_A_genderF); median(sim_pop$Clinic_A_genderF)/n_clinic_a*100
median(sim_pop$Clinic_A_genderM); median(sim_pop$Clinic_A_genderM)/n_clinic_a*100
median(sim_pop$Clinic_A_raceB); median(sim_pop$Clinic_A_raceB)/n_clinic_a*100
median(sim_pop$Clinic_A_raceW); median(sim_pop$Clinic_A_raceW)/n_clinic_a*100
median(sim_pop$Clinic_A_raceH); median(sim_pop$Clinic_A_raceH)/n_clinic_a*100

median(sim_pop$Clinic_B_age0); median(sim_pop$Clinic_B_age0)/n_clinic_b*100
median(sim_pop$Clinic_B_age1); median(sim_pop$Clinic_B_age1)/n_clinic_b*100
median(sim_pop$Clinic_B_age2); median(sim_pop$Clinic_B_age2)/n_clinic_b*100
median(sim_pop$Clinic_B_age3); median(sim_pop$Clinic_B_age3)/n_clinic_b*100
median(sim_pop$Clinic_B_genderF); median(sim_pop$Clinic_B_genderF)/n_clinic_b*100
median(sim_pop$Clinic_B_genderM); median(sim_pop$Clinic_B_genderM)/n_clinic_b*100
median(sim_pop$Clinic_B_raceB); median(sim_pop$Clinic_B_raceB)/n_clinic_b*100
median(sim_pop$Clinic_B_raceW); median(sim_pop$Clinic_B_raceW)/n_clinic_b*100
median(sim_pop$Clinic_B_raceH); median(sim_pop$Clinic_B_raceH)/n_clinic_b*100

median(sim_pop$Clinic_C_age0); median(sim_pop$Clinic_C_age0)/n_clinic_c*100
median(sim_pop$Clinic_C_age1); median(sim_pop$Clinic_C_age1)/n_clinic_c*100
median(sim_pop$Clinic_C_age2); median(sim_pop$Clinic_C_age2)/n_clinic_c*100
median(sim_pop$Clinic_C_age3); median(sim_pop$Clinic_C_age3)/n_clinic_c*100
median(sim_pop$Clinic_C_genderF); median(sim_pop$Clinic_C_genderF)/n_clinic_c*100
median(sim_pop$Clinic_C_genderM); median(sim_pop$Clinic_C_genderM)/n_clinic_c*100
median(sim_pop$Clinic_C_raceB); median(sim_pop$Clinic_C_raceB)/n_clinic_c*100
median(sim_pop$Clinic_C_raceW); median(sim_pop$Clinic_C_raceW)/n_clinic_c*100
median(sim_pop$Clinic_C_raceH); median(sim_pop$Clinic_C_raceH)/n_clinic_c*100
