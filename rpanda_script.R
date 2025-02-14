#RPANDA
library(RPANDA)
library(ggplot2)
library(phytools)
library(pspline)
library(gridExtra)
library(cowplot)

setwd('/Users/marcosqueiroz1/Library/CloudStorage/GoogleDrive-marvin.danque@gmail.com/My Drive/PD/Biogeo/rpanda')
ctree <- read.tree('/Users/marcosqueiroz1/Library/CloudStorage/GoogleDrive-marvin.danque@gmail.com/My Drive/PD/Biogeo/0_trees/2024/parental_v_r_h_ML_calibrated_2024_wo_outgroup_edgelengths_01.nwk')


#functions of constant and variable speciation (lambda) and extinction (mu)
#The variable x represents time, running from the present to the past, 
#while the variable y is a vector containing the different parameters involved 
#in the definition of the temporal dependency. The parameters in y are therefore 
#the parameters that will be estimated by maximum likelihood.

lambda.cst <- function(x,y){y}
lambda.var <- function(x,y){y[1]*exp(y[2]*x)}
mu.cst <- function(x,y){y[1]}
mu.var <- function(x,y){y[1]*exp(y[2]*x)}


lamb_par_init<-c(0.05,0.01)
mu_par_init<-c(0.005)

#function of constant and variable speciation, extinction and environment
#We fit a simple model with an exponential dependence of the speciation rate on
#the environmental variable, no time dependence, and no extinction. We thus 
#define the following:

lambda.cst.env <- function(t, x, y){y[1]*exp(y[2]*x)}
mu.cst.env <- function(t,x,y){0}

lamb_par_init<-c(0.10,0.01)
mu_par_init<-c()

#temperature dataset
#temp_means <- read.csv('./sd_means.csv', header = T)
#temp_mins <- read.csv('./sd_mins.csv', header = T)
#temp_maxs <- read.csv('./sd_maxs.csv', header = T)

temp_means_asia <- read.csv('./sd_means_asia.csv', header = T)
temp_means_asia = temp_means_asia[order(temp_means_asia$time), ]
temp_mins_asia <- read.csv('./sd_mins_asia.csv', header = T)
temp_mins_asia = temp_mins_asia[order(temp_mins_asia$time), ]
temp_maxs_asia <- read.csv('./sd_maxs_asia.csv', header = T)
temp_maxs_asia = temp_maxs_asia[order(temp_maxs_asia$time), ]

#Extract the total time
tot_time<-max(node.age(ctree)$ages)

#Fitting the model by maximum likelihood
#res<-fit_bd(ctree,tot_time,lambda.cst,mu.cst,lamb_par_init,mu_par_init,f=67/130,expo.lamb=TRUE,cst.mu=TRUE)

res<-fit_env(ctree,temp_means_asia,tot_time,lambda.cst.env,mu.cst.env,lamb_par_init,mu_par_init,f=62/130,fix.mu=TRUE,dt=1e-3)

#parameter estimate of y[1],which is the speciation rate that would correspond
#to a temperature o f0°C
res$lamb_par[1]

#returns the maximum parameter estimate of y[2], which is the rate of change in
#speciation rate with temperature
res$lamb_par[2]

#A positive value thus suggests a positive effect of the environmental variable
#(here temperature) on speciation rates.

plot_fit_env(res,temp_means_asia,tot_time)

res_min = fit_env(ctree,temp_mins_asia,tot_time,lambda.cst.env,mu.cst.env,lamb_par_init,mu_par_init,f=62/130,fix.mu=TRUE,dt=1e-3)

#parameter estimate of y[1],which is the speciation rate that would correspond
#to a temperature o f0°C
res_min$lamb_par[1]

#returns the maximum parameter estimate of y[2], which is the rate of change in
#speciation rate with temperature
res_min$lamb_par[2]

plot_fit_env(res_min,temp_mins_asia,tot_time)

re_max = fit_env(ctree,temp_maxs_asia,tot_time,lambda.cst.env,mu.cst.env,lamb_par_init,mu_par_init,f=62/130,fix.mu=TRUE,dt=1e-3)

plot_fit_env(re_max,temp_maxs_asia,tot_time)


########################
#### 2024 ##############
########################

setwd('/Users/marcosqueiroz1/Library/CloudStorage/GoogleDrive-marvin.danque@gmail.com/My Drive/PD/Biogeo/rpanda')
ctree <- read.tree('/Users/marcosqueiroz1/Library/CloudStorage/GoogleDrive-marvin.danque@gmail.com/My Drive/PD/Biogeo/0_trees/2024/parental_v_r_h_ML_calibrated_2024_wo_outgroup_edgelengths_01.nwk')

#Extract the total time
tot_time<-max(node.age(ctree)$ages)

#Temperature data from Oreocharis (Gesneriaceae) paper
#Moonson data from Oreocharis (Gesneriaceae) paper

temp_oreo = read.table('CenozoicTemperature.txt', header = T)
temp_oreo = temp_oreo[temp_oreo$Age <= tot_time + 1,]
moon_oreo = read.table('AsianMonsoonsfrom_15_Ma.txt', header = T)
#moon_oreo = moon_oreo[moon_oreo$Age <= tot_time + 1,]




#sampling proportion
f=62/120

#conditioning on a speciation event at the crown age and survival of the 2 daugther lineages (use when the stem age is not known, in this case tot_time should be the crown age).
cond='crown'

################################################
##### Birth-Death models (no env dependency) ###
################################################


#Scenario 1 - Constant Birth-Death
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(0.1)
mu_par<-c(0.01)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

BcstDcst = fit_bd(ctree, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=62/120,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)

#Scenario 2 - Yule (constant Birth)
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(0.1)
mu_par<-c(0)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

BcstYule<-fit_bd(ctree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=62/120,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)

################################################
## Env Dependence (exponential variation) ######
################################################


#Scenario 3 - Speciation (lambda) rate exponentially depends to the environmental data, extinction (mu) is NULL

#t - time of the tree
#x - environmental data
#y - a numeric vector of the parameters controlling the time and env variation
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.1,0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

BEnvVar_EXPO_temp <-fit_env(ctree,temp_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)
BEnvVar_EXPO_moon <-fit_env(ctree,moon_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)



#Scenario 4 - Speciation (lambda) rate exponentially depends to the env data and extinction (mu) is constant
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}
# lamb_par<-c(abs(treei_BEnvVar_EXPO$lamb_par[1]),treei_BEnvVar_EXPO$lamb_par[2]) # difference 
lamb_par<-c(0.1,0)
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

BEnvVarDCST_EXPO_temp <-fit_env(ctree,temp_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)
BEnvVarDCST_EXPO_moon <-fit_env(ctree,moon_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)

#Scenario 5 - Speciation is constant and extinction exponentianly depends to the env data 
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
# lamb_par<-c(treei_BCSTDCST$lamb_par[1])
lamb_par<-c(0.1)
mu_par<-c(0.01,0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

BCSTDEnvVar_EXPO_temp <- fit_env(ctree,temp_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)
BCSTDEnvVar_EXPO_moon <- fit_env(ctree,moon_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)

#Scenario 6 - Speciation and extinction both exponentially depend of the env data
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(0.1,0)
mu_par<-c(0.01,0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

BEnvVarDEnvVar_EXPO_temp <- fit_env(ctree,temp_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)
BEnvVarDEnvVar_EXPO_moon <- fit_env(ctree,moon_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)

################################################
###### Env Dependence (linear variation) ######
################################################

#Scenario 7 - Speciation linearly depends of env data and extinction is NULL
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){0}
# lamb_par<-c(abs(treei_BEnvVar_EXPO$lamb_par[1]),0)
lamb_par<-c(0.1,0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

BEnvVar_LIN_temp <- fit_env(ctree,temp_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)
BEnvVar_LIN_moon <- fit_env(ctree,moon_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)

#Scenario 8 - Speciation linearly depends of env data and extinction is constant

f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]}
# lamb_par<-c(abs(treei_BEnvVar_LIN$lamb_par[1]),treei_BEnvVar_LIN$lamb_par[2])
lamb_par<-c(0.1,0)
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

BEnvVarDCST_LIN_temp <- fit_env(ctree,temp_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)
BEnvVarDCST_LIN_moon <- fit_env(ctree,moon_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)

#Scenario 9 - Speciation is constant and extinction linearly depends of env data

f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]+y[2]*x}
# lamb_par<-c(treei_BCSTDCST$lamb_par[1])
lamb_par<-c(0.1)
mu_par<-c(0.02,0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

BCSTDEnvVar_LIN_temp <- fit_env(ctree,temp_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)
BCSTDEnvVar_LIN_moon <- fit_env(ctree,moon_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)

#Scenario 10 - Speciation and extinction both linearly depends of env data

f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]+y[2]*x}
# lamb_par<-c(abs(treei_BEnvVarDCST_LIN$lamb_par[1]),treei_BEnvVarDCST_LIN$lamb_par[2])
lamb_par<-c(0.1,0)
mu_par<-c(0.02,0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

BEnvVarDEnvVar_LIN_temp <- fit_env(ctree,temp_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)
BEnvVarDEnvVar_LIN_moon <- fit_env(ctree,moon_oreo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond, dt = 1e-3)


############# RESULTS ###########################################

results<-matrix(NA,18,9)
colnames(results)<-c("Models","Parameters","logL","AICc","Lambda","Alpha","Mu","Beta", "AICcWt")

#Models
results[,1]<-c("BCST_DCST",
               "BCST_Yule",
               "BEnvVar_EXPO_temp",
               "BEnvVar_EXPO_moon",
               "BEnvVarDCST_EXPO_temp",
               "BEnvVarDCST_EXPO_moon",
               "BCSTDEnvVar_EXPO_temp",
               "BCSTDEnvVar_EXPO_moon",
               "BEnvVarDEnvVar_EXPO_temp",
               "BEnvVarDEnvVar_EXPO_moon",
               "BEnvVar_LIN_temp",
               "BEnvVar_LIN_moon",
               "BEnvVarDCST_LIN_temp",
               "BEnvVarDCST_LIN_moon",
               "BCSTDEnvVar_LIN_temp",
               "BCSTDEnvVar_LIN_moon",
               "BEnvVarDEnvVar_LIN_temp",
               "BEnvVarDEnvVar_LIN_moon")

#Parameters
results[1,2] = 2
results[2,2] = 1
results[3,2] = 2
results[4,2] = 2
results[5,2] = 3
results[6,2] = 3
results[7,2] = 3
results[8,2] = 3
results[9,2] = 4
results[10,2] = 4
results[11,2] = 2
results[12,2] = 2
results[13,2] = 3
results[14,2] = 3
results[15,2] = 3
results[16,2] = 3
results[17,2] = 4
results[18,2] = 4

#logL
results[1,3] = Bcst_Dcst$LH
results[2,3] = BcstYule$LH
results[3,3] = BEnvVar_EXPO_temp$LH
results[4,3] = BEnvVar_EXPO_moon$LH
results[5,3] = BEnvVarDCST_EXPO_temp$LH
results[6,3] = BEnvVarDCST_EXPO_moon$LH
results[7,3] = BCSTDEnvVar_EXPO_temp$LH
results[8,3] = BCSTDEnvVar_EXPO_moon$LH
results[9,3] = BEnvVarDEnvVar_EXPO_temp$LH
results[10,3] = BEnvVarDEnvVar_EXPO_moon$LH
results[11,3] = BEnvVar_LIN_temp$LH
results[12,3] = BEnvVar_LIN_moon$LH
results[13,3] = BEnvVarDCST_LIN_temp$LH
results[14,3] = BEnvVarDCST_LIN_moon$LH
results[15,3] = BCSTDEnvVar_LIN_temp$LH
results[16,3] = BCSTDEnvVar_LIN_moon$LH
results[17,3] = BEnvVarDEnvVar_LIN_temp$LH
results[18,3] = BEnvVarDEnvVar_LIN_moon$LH

#AICc
results[1,4] = BcstDcst$aic
results[2,4] = BcstYule$aic
results[3,4] = BEnvVar_EXPO_temp$aic
results[4,4] = BEnvVar_EXPO_moon$aic
results[5,4] = BEnvVarDCST_EXPO_temp$aic
results[6,4] = BEnvVarDCST_EXPO_moon$aic
results[7,4] = BCSTDEnvVar_EXPO_temp$aic
results[8,4] = BCSTDEnvVar_EXPO_moon$aic
results[9,4] = BEnvVarDEnvVar_EXPO_temp$aic
results[10,4] = BEnvVarDEnvVar_EXPO_moon$aic
results[11,4] = BEnvVar_LIN_temp$aic
results[12,4] = BEnvVar_LIN_moon$aic
results[13,4] = BEnvVarDCST_LIN_temp$aic
results[14,4] = BEnvVarDCST_LIN_moon$aic
results[15,4] = BCSTDEnvVar_LIN_temp$aic
results[16,4] = BCSTDEnvVar_LIN_moon$aic
results[17,4] = BEnvVarDEnvVar_LIN_temp$aic
results[18,4] = BEnvVarDEnvVar_LIN_moon$aic

#Lambda0
results[1,5] = abs(BcstDcst$lamb_par[1])
results[2,5] = abs(BcstYule$lamb_par[1])
results[3,5] = abs(BEnvVar_EXPO_temp$lamb_par[1])
results[4,5] = abs(BEnvVar_EXPO_moon$lamb_par[1])
results[5,5] = abs(BEnvVarDCST_EXPO_temp$lamb_par[1])
results[6,5] = abs(BEnvVarDCST_EXPO_moon$lamb_par[1])
results[7,5] = abs(BCSTDEnvVar_EXPO_temp$lamb_par[1])
results[8,5] = abs(BCSTDEnvVar_EXPO_moon$lamb_par[1])
results[9,5] = abs(BEnvVarDEnvVar_EXPO_temp$lamb_par[1])
results[10,5] = abs(BEnvVarDEnvVar_EXPO_moon$lamb_par[1])
results[11,5] = abs(BEnvVar_LIN_temp$lamb_par[1])
results[12,5] = abs(BEnvVar_LIN_moon$lamb_par[1])
results[13,5] = abs(BEnvVarDCST_LIN_temp$lamb_par[1])
results[14,5] = abs(BEnvVarDCST_LIN_moon$lamb_par[1])
results[15,5] = abs(BCSTDEnvVar_LIN_temp$lamb_par[1])
results[16,5] = abs(BCSTDEnvVar_LIN_moon$lamb_par[1])
results[17,5] = abs(BEnvVarDEnvVar_LIN_temp$lamb_par[1])
results[18,5] = abs(BEnvVarDEnvVar_LIN_moon$lamb_par[1])

#Alpha
# - no alpha for commented lines
# results[1,6] = BcstDcst$lamb_par[2]
# results[2,6] = BcstYule$lamb_par[2]
results[3,6] = BEnvVar_EXPO_temp$lamb_par[2]
results[4,6] = BEnvVar_EXPO_moon$lamb_par[2]
results[5,6] = BEnvVarDCST_EXPO_temp$lamb_par[2]
results[6,6] = BEnvVarDCST_EXPO_moon$lamb_par[2]
# results[7,6] = BCSTDEnvVar_EXPO_temp$lamb_par[2]
# results[8,6] = BCSTDEnvVar_EXPO_moon$lamb_par[2]
results[9,6] = BEnvVarDEnvVar_EXPO_temp$lamb_par[2]
results[10,6] = BEnvVarDEnvVar_EXPO_moon$lamb_par[2]
results[11,6] = BEnvVar_LIN_temp$lamb_par[2]
results[12,6] = BEnvVar_LIN_moon$lamb_par[2]
results[13,6] = BEnvVarDCST_LIN_temp$lamb_par[2]
results[14,6] = BEnvVarDCST_LIN_moon$lamb_par[2]
# results[15,6] = BCSTDEnvVar_LIN_temp$lamb_par[2]
# results[16,6] = BCSTDEnvVar_LIN_moon$lamb_par[2]
results[17,6] = BEnvVarDEnvVar_LIN_temp$lamb_par[2]
results[18,6] = BEnvVarDEnvVar_LIN_moon$lamb_par[2]

#Mu0
# - no mu for commented lines
results[1,7] = abs(BcstDcst$mu_par[1])
#results[2,7] = abs(BcstYule$mu_par[1])
#results[3,7] = abs(BEnvVar_EXPO_temp$mu_par[1])
#results[4,7] = abs(BEnvVar_EXPO_moon$mu_par[1])
results[5,7] = abs(BEnvVarDCST_EXPO_temp$mu_par[1])
results[6,7] = abs(BEnvVarDCST_EXPO_moon$mu_par[1])
results[7,7] = abs(BCSTDEnvVar_EXPO_temp$mu_par[1])
results[8,7] = abs(BCSTDEnvVar_EXPO_moon$mu_par[1])
results[9,7] = abs(BEnvVarDEnvVar_EXPO_temp$mu_par[1])
results[10,7] = abs(BEnvVarDEnvVar_EXPO_moon$mu_par[1])
#results[11,7] = abs(BEnvVar_LIN_temp$mu_par[1])
#results[12,7] = abs(BEnvVar_LIN_moon$mu_par[1])
results[13,7] = abs(BEnvVarDCST_LIN_temp$mu_par[1])
results[14,7] = abs(BEnvVarDCST_LIN_moon$mu_par[1])
results[15,7] = abs(BCSTDEnvVar_LIN_temp$mu_par[1])
results[16,7] = abs(BCSTDEnvVar_LIN_moon$mu_par[1])
results[17,7] = abs(BEnvVarDEnvVar_LIN_temp$mu_par[1])
results[18,7] = abs(BEnvVarDEnvVar_LIN_moon$mu_par[1])

#Beta
# - no beta for commented lines
#results[1,8] = BcstDcst$mu_par[2]
#results[2,8] = BcstYule$mu_par[2]
#results[3,8] = BEnvVar_EXPO_temp$mu_par[2]
#results[4,8] = BEnvVar_EXPO_moon$mu_par[2]
#results[5,8] = BEnvVarDCST_EXPO_temp$mu_par[2]
#results[6,8] = BEnvVarDCST_EXPO_moon$mu_par[2]
results[7,8] = BCSTDEnvVar_EXPO_temp$mu_par[2]
results[8,8] = BCSTDEnvVar_EXPO_moon$mu_par[2]
results[9,8] = BEnvVarDEnvVar_EXPO_temp$mu_par[2]
results[10,8] = BEnvVarDEnvVar_EXPO_moon$mu_par[2]
#results[11,8] = BEnvVar_LIN_temp$mu_par[2]
#results[12,8] = BEnvVar_LIN_moon$mu_par[2]
#results[13,8] = BEnvVarDCST_LIN_temp$mu_par[2]
#results[14,8] = BEnvVarDCST_LIN_moon$mu_par[2]
results[15,8] = BCSTDEnvVar_LIN_temp$mu_par[2]
results[16,8] = BCSTDEnvVar_LIN_moon$mu_par[2]
results[17,8] = BEnvVarDEnvVar_LIN_temp$mu_par[2]
results[18,8] = BEnvVarDEnvVar_LIN_moon$mu_par[2]

#AICc weights
results[,9] = aic.w(as.numeric(results[,4]))

#Save the results
write.table(results, file = 'RPANDA_Curcuma_results.txt', row.names = F)
save(BcstDcst, file = 'BcstDcst.Rdata')
save(BcstYule, file = 'BcstYule.Rdata')
save(BEnvVar_EXPO_temp, file = 'BEnvVar_EXPO_temp.Rdata')
save(BEnvVar_EXPO_moon, file = 'BEnvVar_EXPO_moon.Rdata')
save(BEnvVarDCST_EXPO_temp, file = 'BEnvVarDCST_EXPO_temp.Rdata')
save(BEnvVarDCST_EXPO_moon, file = 'BEnvVarDCST_EXPO_moon.Rdata')
save(BCSTDEnvVar_EXPO_temp, file = 'BCSTDEnvVar_EXPO_temp.Rdata')
save(BCSTDEnvVar_EXPO_moon, file = 'BCSTDEnvVar_EXPO_moon.Rdata')
save(BEnvVarDEnvVar_EXPO_temp, file = 'BEnvVarDEnvVar_EXPO_temp.Rdata')
save(BEnvVarDEnvVar_EXPO_moon, file = 'BEnvVarDEnvVar_EXPO_moon.Rdata')
save(BEnvVar_LIN_temp, file = 'BEnvVar_LIN_temp.Rdata')
save(BEnvVar_LIN_moon, file = 'BEnvVar_LIN_moon.Rdata')
save(BEnvVarDCST_LIN_temp, file = 'BEnvVarDCST_LIN_temp.Rdata')
save(BEnvVarDCST_LIN_moon, file = 'BEnvVarDCST_LIN_moon.Rdata')
save(BCSTDEnvVar_LIN_temp, file = 'BCSTDEnvVar_LIN_temp.Rdata')
save(BCSTDEnvVar_LIN_moon, file = 'BCSTDEnvVar_LIN_moon.Rdata')
save(BEnvVarDEnvVar_LIN_temp, file = 'BEnvVarDEnvVar_LIN_temp.Rdata')
save(BEnvVarDEnvVar_LIN_moon, file = 'BEnvVarDEnvVar_LIN_moon.Rdata')

resi<-list("Clade_age" = tot_time,
           "Taxon_sampling" = Ntip(ctree),
           "Sampling_fraction", f,
           "BcstDcst" = BcstDcst,
           "BcstYule" = BcstYule,
           "BEnvVar_EXPO_temp" = BEnvVar_EXPO_temp,
           "BEnvVar_EXPO_moon" = BEnvVar_EXPO_moon,
           "BEnvVarDCST_EXPO_temp" = BEnvVarDCST_EXPO_temp,
           "BEnvVarDCST_EXPO_moon" = BEnvVarDCST_EXPO_moon,
           "BCSTDEnvVar_EXPO_temp" = BCSTDEnvVar_EXPO_temp,
           "BCSTDEnvVar_EXPO_moon" = BCSTDEnvVar_EXPO_moon,
           "BEnvVarDEnvVar_EXPO_temp" = BEnvVarDEnvVar_EXPO_temp,
           "BEnvVarDEnvVar_EXPO_moon" = BEnvVarDEnvVar_EXPO_moon,
           "BEnvVar_LIN_temp" = BEnvVar_LIN_temp,
           "BEnvVar_LIN_moon" = BEnvVar_LIN_moon,
           "BEnvVarDCST_LIN_temp" = BEnvVarDCST_LIN_temp,
           "BEnvVarDCST_LIN_moon" = BEnvVarDCST_LIN_moon,
           "BCSTDEnvVar_LIN_temp" = BCSTDEnvVar_LIN_temp,
           "BCSTDEnvVar_LIN_moon" = BCSTDEnvVar_LIN_moon,
           "BEnvVarDEnvVar_LIN_temp" = BEnvVarDEnvVar_LIN_temp,
           "BEnvVarDEnvVar_LIN_moon" = BEnvVarDEnvVar_LIN_moon)


###### Plots #########


#Temperature vs Monsoon regimes

tm_mn = ggplot() +
  geom_line(data = temp_oreo, aes(x = Age, y = Temperature), color = 'blue') +
  geom_line(data = moon_oreo, aes(x = Age, y = Monsoons), color = 'magenta') +
  scale_y_continuous(
    name = 'Temperature',
    sec.axis = sec_axis(~./4, name = 'Monsoons')
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(color = 'blue'),
    axis.title.y.right = element_text(color = 'magenta')
  )

tm_mn

tmp_plot = ggplot(data = temp_oreo, aes(x = -Age, y = Temperature)) +
  geom_line(alpha = 0.5, color = 'blue') +
  theme_classic() +
  stat_smooth(method = loess, color = 'black') +
  xlab('Time (Mya)') +
  ylab('Mean temperature (°C)')

tmp_plot 

# Lineage through time
lineage_time <- as.data.frame(ltt.plot.coords(ctree))

lineage_time$time <- lineage_time$time * -1

lineage_plot = ggplot(data = lineage_time, aes(x = time, y = N)) +
  geom_line() +
  theme_classic()

lineage_plot

#Temperature and lineage plot
tmp_lin_plot = ggplot() +
  geom_line(data = temp_oreo, aes(x = -Age, y = Temperature), color = 'blue', alpha = 0.5) +
  geom_smooth(data = temp_oreo, aes(x = -Age, y = Temperature), color = 'black') +
   geom_line(data = lineage_time, aes(x = -time, y = N/4), color = 'orange', linetype = 2, size = 1.5) +
   scale_y_continuous(
     name = 'Mean Temperature (°C)',
     sec.axis = sec_axis(~.*4, name = 'Lineages')
   ) +
  labs(x = 'Time (Mya)') +
  xlim(-15, 0) +
   theme_classic()

tmp_lin_plot


#Speciation, extinction and diversification plots

best_model = BEnvVarDEnvVar_LIN_temp
t = seq(0, tot_time, length.out = 100)

#Speciation vs time
sps_t = sapply(t, best_model$f.lamb)
t_vs_sps = as.data.frame(cbind(t, sps_t))
colnames(t_vs_sps) = c('time', 'sps_t')

t_vs_sps_plot = ggplot(data = t_vs_sps, aes(x = -time, y = sps_t)) +
  geom_line(color = 'blue') +
  ylim(0,10)+
  theme_minimal() +
  labs(x = 'Time (Mya)',
       y = expression('Speciation Rate '~(lambda)))

t_vs_sps_plot

#Speciation vs environment
df = smooth.spline(temp_oreo[, 1], temp_oreo[, 2])$df
spline_result = sm.spline(temp_oreo[, 1], temp_oreo[, 2], df = df)
spline_predict = predict(spline_result, t)

sps_e = sapply(t, best_model$f.lamb)

e_vs_sps_e = as.data.frame(cbind(spline_predict, sps_e))
colnames(e_vs_sps_e) = c('spline_pred', 'sps_e')

e_vs_sps_e_plot = ggplot(data = e_vs_sps_e, aes(x = spline_pred, y = sps_e)) +
  geom_point(color = 'blue') + 
  ylim(0,10)+
  theme_minimal() +
  labs(x = 'Temperature (°C)',
       y = expression('Speciation Rate '~(lambda)))


e_vs_sps_e_plot

#Extinction vs time
ext_t = sapply(t, best_model$f.mu)
t_vs_ext_t = as.data.frame(cbind(t, ext_t))
colnames(t_vs_ext_t) = c('time', 'ext_t')

t_vs_ext_plot = ggplot(data = t_vs_ext_t, aes(x = -time, y = ext_t)) +
  geom_line(color = 'red') +
  theme_minimal() +
  labs(x = 'Time (Mya)',
       y = expression('Extinction Rate '~(mu)))

t_vs_ext_plot

#Extinction vs environment
ext_e = sapply(t, best_model$f.mu)

e_vs_ext_e = as.data.frame(cbind(spline_predict, ext_e))
colnames(e_vs_ext_e) = c('spline_pred', 'ext_e')

e_vs_ext_e_plot = ggplot(data = e_vs_ext_e, aes(x = spline_pred, y = ext_e)) +
  geom_point(color = 'red') + 
  theme_minimal() +
  labs(x = 'Temperature (°C)',
       y = expression('Extinction Rate '~(mu)))


e_vs_ext_e_plot

#Net Diversification (time)
r = best_model$f.lamb(t) - best_model$f.mu(t)

net_div_time = as.data.frame(cbind(t,r))

net_div_plot_time = ggplot(data = net_div_time, aes(x = -t, y = r)) +
  geom_line(color = 'black') +
  theme_minimal() +
  labs(x = 'Time (Mya)',
       y = 'Diversification rate') +
  xlim(-15, 0)

net_div_plot_time

#Net Diversifiation (env)

net_div_env = as.data.frame(cbind(spline_predict, r))
colnames(net_div_env) = c('spline_pred', 'r')

net_div_plot_env = ggplot(data = net_div_env, aes(x = spline_pred, y = r)) +
  geom_point(color = 'black') +
  theme_minimal() +
  labs(x = 'Temperature (°C)',
       y = 'Diversification rate')

net_div_plot_env


#Grid
grid.arrange(t_vs_sps_plot, t_vs_ext_plot, ncol=2) #speciation and extinction vs time
grid.arrange(e_vs_sps_e_plot, e_vs_ext_e_plot, ncol = 2) #speciation and extinction vs env
grid.arrange(net_div_plot_time, net_div_plot_env, ncol = 2) #net diversification

grid.arrange(t_vs_sps_plot, e_vs_sps_e_plot, ncol=2) #speciation vs time and env
grid.arrange(t_vs_ext_plot, e_vs_ext_e_plot, ncol = 2) #extinction vs time and env

svglite('speciation_extinction_diversification_plots.svg', height = 11.7, width = 8.3)
plot_grid(t_vs_sps_plot, e_vs_sps_e_plot,
             t_vs_ext_plot, e_vs_ext_e_plot,
             net_div_plot_time, net_div_plot_env,
             nrow=3, ncol=2, align = 'v',
          labels = c('(a)', '(b)',
                     '(c)', '(d)',
                     '(e)', '(f)'))

dev.off()

svglite('temp_vs_time_and_netdiversification.svg', height = 11.7, width = 4.2)
plot_grid(tmp_lin_plot, net_div_plot_time, net_div_plot_env, align = 'v', ncol = 1) #temperature and lineage vs time and netdiversification plots
dev.off()




