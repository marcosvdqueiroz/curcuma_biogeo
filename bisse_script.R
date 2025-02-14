#BISSE

library(geiger)
library(phytools)
library(diversitree)
library(ggplot2)
library(tidyr)
library(dplyr)

setwd('/Users/marcosqueiroz1/Library/CloudStorage/GoogleDrive-marvin.danque@gmail.com/My Drive/PD/Biogeo/test_bisse')

#read the tree and the traits files
ctree <- read.tree('/Users/marcosqueiroz1/Library/CloudStorage/GoogleDrive-marvin.danque@gmail.com/My Drive/PD/Biogeo/0_trees/2024/parental_v_r_h_ML_calibrated_2024_wo_outgroup_edgelengths_01.nwk')
ctraits <- read.csv('./Curcuma_v_r_h_withough_v_traits_2024.txt', row.names = 1, sep = '\t', stringsAsFactors = T)

#extract trait data
traits <- ctraits[,1]

#set names
names(traits) <- rownames(ctraits)

#plot the tree
plotTree(ctree, ftype = 'i', fsize = 0.6, offset = 0.5)

#add tip labels
tiplabels(pie = to.matrix(traits, 0:1)[ctree$tip.label,],
          piecol = c('white', 'black'), cex = 0.4)
#create legend
legend('topleft', c('non-hybrid', 'hybrid'),
       pch = 21, pt.cex = 1.6,
       cex = 0.8, bty = 'n',
       pt.bg = c('white', 'black'))


#### make BiSSE likelihood function
bisse.model <- make.bisse(ctree, traits)

#find reasonable parameter values for optimization
p <- starting.point.bisse(ctree)
p

#optimize BiSSE model
bisse.mle <- find.mle(bisse.model, p)
bisse.mle

#lambda 0 is the speciation rate of the non-hybrids 2.447086e-01
#lambda 1 is the speciation rate of the hybrids 3.522379e-01

#Let's compare if this model significantly better explains our data compared to a null model in which
#the speciation and exticntion rates are constant through the tree.

#create a constrained null model
bissenull.model <- constrain(bisse.model, lambda1 ~ lambda0, mu1 ~ mu0)

#optimize null model
bissenulll.mle <- find.mle(bissenull.model, p[c(-2,-4)])
coef(bissenulll.mle)
logLik(bissenulll.mle)

#run likelihood-ratio test

bisseAnova <- anova(bisse.mle, null = bissenulll.mle)
bisseAnova

aicw(setNames(bisseAnova$AIC, rownames(bisseAnova)))


#Bayesian MCMC
prior <- make.prior.exponential(1/2*0.4)
prior

bisse.mcmc <- mcmc(bisse.model, bisse.mle$par, nsteps = 10000, prior = prior, w = 0.1, print.every = 100)

par(mfrow = c(1,2), mar = c(5.1, 4.1, 3.1, 2.1))
col <- setNames(c('#d8b365', '#5ab4ac'),
                c('non-hybrid', 'hybrid'))
profiles.plot(bisse.mcmc[,c('lambda0', 'lambda1')],
              col.line = col, las = 1, bty = 'n',
              xlab = expression(lambda), cex.axis = 0.7)
#legend('topleft', names(col), pch = 15, col = col,
#       pt.cex = 1.5, bty = 'n', cex = 0.7)
profiles.plot(bisse.mcmc[,c('mu0', 'mu1')],
              col.line = col, las = 1, bty = 'n',
              xlab = expression(mu), cex.axis = 0.7)
legend('topright', names(col), pch = 15, col = col,
       pt.cex = 1.5, bty = 'n', cex = 0.7)

#sums
#lambda 1
sum(bisse.mcmc$lambda1 > bisse.mcmc$lambda0) / length(bisse.mcmc$lambda1)

#lambda 0
sum(bisse.mcmc$lambda0 > bisse.mcmc$lambda1) / length(bisse.mcmc$lambda1)

#0.7719 lambda 1. 
#0.226 lambda 0

#High probability that the diversification rate is higher in hybrids than non-hybrids


#mu
sum(bisse.mcmc$mu1 > bisse.mcmc$mu0) / length(bisse.mcmc$mu1)
#0.536


###ggplot2 histogram

lambda_df = as.data.frame(bisse.mcmc[,c('lambda0', 'lambda1')])
lambda_df_long = pivot_longer(lambda_df, cols = c(lambda0, lambda1), names_to = 'variable', values_to = 'value')
          

# Create the  density plots
ggplot(lambda_df_long, aes(x = value, fill = variable, color = variable)) +
  geom_density(alpha = 0.7) +
  labs(title = "Density Plots of lambda0 and lambda1",
       x = "Value",
       y = "Density",
       fill = "Variable",
       color = "Variable") +
  theme_minimal()

lambda_df %>% summarise(across(everything(), list(mean = mean, sd = sd))) 

