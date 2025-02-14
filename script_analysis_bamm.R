library(BAMMtools)
library(coda)
library(svglite)

setwd('/Users/marcosqueiroz1/Library/CloudStorage/GoogleDrive-marvin.danque@gmail.com/My Drive/PD/Biogeo/test_bamm')

ctree <- read.tree('./Curcuma_v_r_h_without_v_ultrametric.nwk')
cdata <- getEventData(ctree, eventdata = './default_params/event_data.txt', burnin = 0.1)

#Checking MCMC convergence
mcmcout <- read.csv("./default_params/Curcuma_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]


effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
#Should be at least 200 (my was 8075 and 8488)

####
#How many rate shifts?
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
post_probs[1] / post_probs[2]

#model with name ‘0’ is model M_0, or a model with no rate shifts
#The probability of model ‘0’ is the posterior probability of a model with just 
#a single evolutionary rate dynamic (no rate shifts)
#model with name '1' is the one with a single diversification shift, and so one...


shift_probs <- summary(cdata)
#Our case shows that there's 84% chance of no diversification shifts in our tree.

#### MCMC
#We can compare the models in a Bayesian approach
bfmat <- computeBayesFactors(mcmcout, expectedNumberOfShifts=1, burnin=0.1)

#The output here is a pairwise table comparing all models with each other.
#In the first column '0', check the Bayes factor in favor for the other models
#In our case, the evidence in favor of a model with a rate shift, relative to the null model, is low (0.31)
#Bayes factors greater than 20 generally imply strong evidence for one model over another; 
#values greater than 50 are very strong evidence in favor of the numerator model

plotPrior(mcmcout, expectedNumberOfShifts = 1)

#Plotting
plot.bammdata(cdata, lwd=2)

css <- credibleShiftSet(cdata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)
plot.credibleshiftset(css)

#Only one shift....

best <- getBestShiftConfiguration(cdata, expectedNumberOfShifts=1)
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=2.5)
#If there is only one rate regime, then you have no rate shifts: 
#the single rate regime starts at the root and describes the entire tree.

###Clade specific

allrates <- getCladeRates(cdata)
mean_allrates <- mean(allrates$lambda)
mean_allrates
quantile_allrates <- quantile(allrates$lambda, c(0.05, 0.95))
quantile_allrates

plotTree(ctree, ftype = 'i', fsize = 0.6, offset = 0.5)
nodelabels()

#Clade Curcuma node 104

curcuma_rates <- getCladeRates(cdata, node = 104)
mean_curcuma_rates <- mean(curcuma_rates$lambda)
mean_curcuma_rates
quantile_curcuma_rates <- quantile(curcuma_rates$lambda, c(0.05, 0.95))
quantile_curcuma_rates

#non hybrids rates
non_curcuma_rates <- getCladeRates(cdata, node = 104, nodetype = "exclude")
mean_non_curcuma_rates <- mean(non_curcuma_rates$lambda)
mean_non_curcuma_rates
quantile_non_curcuma_rates <- quantile(non_curcuma_rates$lambda, c(0.05, 0.95))
quantile_non_curcuma_rates

#plot the specific clade rates

svglite('./with_incomplete_taxon_sampling/rates_clade_specific.svg', height = 4, width = 6)
par(mfrow = c(1,3))
st <- max(branching.times(ctree))
plotRateThroughTime(cdata, intervalCol = 'red', avgCol="red", start.time=st, ylim=c(0,1), cex.axis=1)
text(x=4, y= 0.8, label="All Curcuma", font=.5, cex=1.0, pos=4)

plotRateThroughTime(cdata, intervalCol = 'blue', avgCol="blue", start.time=st, node = 104, ylim=c(0,1), cex.axis=1)
text(x=4, y= 0.8, label="Hybrid Curcuma", font=.5, cex=1.0, pos=4)

plotRateThroughTime(cdata, intervalCol="darkgreen", avgCol="darkgreen", start.time=st, node=104, nodetype = "exclude", ylim=c(0,1), cex.axis=1)
text(x=5, y= 0.8, label="Non-hybrid Curcuma", font=.5, cex=1.0, pos=4)

dev.off()

#plot speciation, extinction and netdiversification for all tree
#speciation
svglite('./rates_plots.svg', width = 4.2, height = 11.6)
par(mfrow = c(3,1))
plotRateThroughTime(cdata, intervalCol = 'blue', avgCol="blue", start.time=st, ylim=c(0,1), cex.axis=1)

#extinction
extinction_rate_df = plotRateThroughTime(cdata, intervalCol = 'red', avgCol="red", start.time=st, ylim=c(0,1), cex.axis=1, ratetype = 'extinction')

#netdiv
netdiv_rate_df = plotRateThroughTime(cdata, intervalCol = 'magenta', avgCol="magenta", start.time=st, ylim=c(0,1), cex.axis=1, ratetype = 'netdiv')

dev.off()

#plot speciation, extinction and netdiversification for all tree (black and white)
svglite('./rates_plots_bw.svg', width = 4.2, height = 11.6)
par(mfrow = c(3,1))
plotRateThroughTime(cdata, intervalCol = 'gray80', avgCol="black", start.time=st, ylim=c(0,1), cex.axis=1)

#extinction
plotRateThroughTime(cdata, intervalCol = 'gray80', avgCol="black", start.time=st, ylim=c(0,1), cex.axis=1, ratetype = 'extinction')

#netdiv
plotRateThroughTime(cdata, intervalCol = 'gray80', avgCol="black", start.time=st, ylim=c(0,1), cex.axis=1, ratetype = 'netdiv')

dev.off()

svglite('./bamm_tree.svg')
plot.bammdata(cdata, lwd=3, legend = T)
axisPhylo()
dev.off()
