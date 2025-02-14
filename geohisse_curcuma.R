#### GeoHiSSE code ####

# For more information, see:
# Caetano, D. S., O'Meara, B. C., & Beaulieu, J. M. (2018). 
# Hidden state models improve state‐dependent diversification approaches, 
# including biogeographical models. Evolution, 72(11), 2308-2324.

library(ape)
library(hisse, lib.loc= '../RPackages')
library(parallel)



### Import data
tree <- read.tree('parental_v_r_h_ML_calibrated_2024_wo_outgroup_edgelengths_01.nwk') # load tree file 
dist <- read.table('curcuma_areas_high_vs_low_seasonality.txt', header = T) # load distribution file 

# Preparing data - areas have to be as 0 (11 - widespread), 
# 1 (10, endemic of first area - low seasonality) 
# and 2 (01, endemic of second area - high seasonality)

areas <- as.data.frame(rep(1, nrow(dist)))
dist <- cbind(dist, areas)
colnames(dist)[4] <- "area"

for (i in 1:length(dist$area)){
  if (dist[i, "high_seasonality"] == 1 && dist[i, "low_seasonality"]  == 1){
    dist[i, "area"] = 0 
  }
  if (dist[i, "high_seasonality"] == 0 && dist[i, "low_seasonality"]  == 1){
    dist[i, "area"] = 1
  }
  if (dist[i, "high_seasonality"] == 1 && dist[i, "low_seasonality"]  == 0){
    dist[i, "area"] = 2
  }
}
states<-dist[,c("species", "area")]

table(states$area) # check if species-richness in each range make sense

# 2 - high seasonality"
# 1 - low seasonality"
# 0 - widespread

# Load sampling fraction for the group
#sf<-c(1,1,1) # e.g. if it's fully sampled 
sf <- c(1, 0.7, 1) #widespread 100%, low seasonality 70%, high seasonality 100%
  
#
phy=tree
dat=states
outname='Curcuma'

# We used the same 18 models of Caetano et al. (2018) (plus a second set of models including jump dispersal) - see their original publication for more information

###############################################################################
## Block of GeoSSE-like models.
## Here extirpation is linked to range reduction.
###############################################################################

## Model 1 - Dispersal parameters vary only, no range-dependent diversification. 
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE, separate.extirpation=FALSE)
mod1 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                      hidden.states=FALSE, trans.rate=trans.rate,
                      assume.cladogenetic=TRUE)
saveRDS(mod1, file=paste0(outname, "_mod1.rds"))
print('Model 1 is done!')


## Model 2. Canonical GeoSSE model, range effect on diversification 
speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE, separate.extirpation=FALSE)
mod2 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                      hidden.states=FALSE, trans.rate=trans.rate,
                      assume.cladogenetic=TRUE)
saveRDS(mod2, file=paste0(outname, "_mod2.rds"))
print('Model 2 is done!')

## Model 3. Heterogeneous diversification, not tied to range evolution.
## Assumes three distinct diversification rates.
## Dispersion parameters across hidden areas are the same.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,2,2,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
mod3 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                      hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)
saveRDS(mod3, file=paste0(outname, "_mod3.rds"))
print('Model 3 is done!')


## Model 4. Heterogeneous diversification, tied to range evolution. 
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE, separate.extirpation=FALSE)
mod4 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                      hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)
saveRDS(mod4, file=paste0(outname, "_mod4.rds"))
print('Model 4 is done!')


## Model 5. Heterogeneous diversification, not tied to range evolution. 
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2))
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
mod5 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                      hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)
saveRDS(mod5, file=paste0(outname, "_mod5.rds"))
print('Model 5 is done!')

## Model 6. Heterogeneous diversification, not tied to range evolution. 
## Assumes two distinct diversification rates.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
mod6 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                 hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE)
saveRDS(mod6, file=paste0(outname, "_mod6.rds"))
print('Model 6 is done!')


###############################################################################
## Block of cladogenetic models not GeoSSE-like.
## Here extirpation is NOT linked to range reduction.
## So range reduction is different from the extinction of an endemic lineage.
## Jumps between endemic areas are not allowed.
###############################################################################

## Model 7 - Dispersal parameters vary only, no range-dependent diversification. 
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE, separate.extirpation=TRUE)
mod7 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                 hidden.states=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE)
saveRDS(mod7, file=paste0(outname, "_mod7.rds"))
print('Model 7 is done!')
     
## Model 8. GeoSSE model, with range effect on diversification turnover <- c(1,2,3)
speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE, separate.extirpation=TRUE)
mod8 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                 hidden.states=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE)
saveRDS(mod8, file=paste0(outname, "_mod8.rds"))
print('Model 8 is done!')
        

## Model 9. Heterogeneous diversification, not tied to range evolution.
## Assumes three distinct diversification rates.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,2,2,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, make.null=TRUE, include.jumps=FALSE, separate.extirpation=TRUE)
mod9 <- try( GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                      hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
saveRDS(mod9, file=paste0(outname, "_mod9.rds"))
print('Model 9 is done!')

## Model 10. Heterogeneous diversification, tied to range evolution.
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE, separate.extirpation=TRUE)
mod10 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                       hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod10, file=paste0(outname, "_mod10.rds"))
print('Model 10 is done!')

## Model 11. Heterogeneous diversification, not tied to range evolution. 
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2))
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, make.null=TRUE, include.jumps=FALSE, separate.extirpation=TRUE)
mod11 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                       hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod11, file=paste0(outname, "_mod11.rds"))
print('Model 11 is done!')
        

## Model 12. Heterogeneous diversification, not tied to range evolution. 
## Assumes two distinct diversification rates.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE, include.jumps=FALSE, separate.extirpation=TRUE)
mod12 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod12, file=paste0(outname, "_mod12.rds"))
print('Model 12 is done!')
  

###############################################################################
## Block of anagenetic geographic models (MuSSE).
## Here models emulate GeoSSE (or GeoHiSSE) but changes only happen along branches.
###############################################################################
## Model 13. Transitions only. No character effect on diversification
speciation <- c(1,1,1)
extirpation <- c(1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE, separate.extirpation=TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod13 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                  hidden.states=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod13, file=paste0(outname, "_mod13.rds"))
print('Model 13 is done!')

## Model 14. Character effect on diversification.
speciation <- c(1,2,3)
extirpation <- c(1,2,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE, separate.extirpation=TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod14 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                  hidden.states=FALSE, trans.rate=trans.rate.mod, assume.cladogenetic=FALSE))
saveRDS(mod14, file=paste0(outname, "_mod14.rds"))
print('Model 14 is done!')
        
## Model 15. No character effect on diversification.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,1,2,2,2,3,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
mod15 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation, 
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod15, file=paste0(outname, "_mod15.rds"))
print('Model 15 is done!')

## Model 16. Character effect on diversification, with a hidden state
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4,5,6)
trans.rate <- trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE, separate.extirpation=TRUE)
mod16 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation, hidden.states=TRUE,
                  trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod16, file=paste0(outname, "_mod16.rds"))
print('Model 16 is done!')


## Model 17. No character effect on diversification, multiple shifts
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
trans.rate <- trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
mod17 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation, hidden.states=TRUE,
                  trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod17, file=paste0(outname, "_mod17.rds"))
print('Model 17 is done!')
   

## Model 18. No character effect on diversification, multiple shifts.
speciation <- c(rep(1,3), rep(2,3))
extirpation <- c(rep(1,3), rep(2,3))
trans.rate <- trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
mod18 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation, hidden.states=TRUE,
                  trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod18, file=paste0(outname, "_mod18.rds"))
print('Model 18 is done!')
     

#################### JUMP MODELS ##########################
#################### JUMP MODELS ##########################
#################### JUMP MODELS ##########################
#################### JUMP MODELS ##########################
#################### JUMP MODELS ##########################

## argument "include.jumps" set to TRUE

###############################################################################
## Block of GeoSSE-like models.
## Here extirpation is linked to range reduction.
###############################################################################
## Model 19 - Dispersal parameters vary only, no range-dependent diversification. 
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE 
                                    , separate.extirpation=FALSE) 
mod19 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE)) 
saveRDS(mod19, file=paste0(outname, "_mod19.rds"))
print('Model 19 is done!')

## Model 20. Canonical GeoSSE model, range effect on diversification 
speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                    , separate.extirpation=FALSE) 
mod20 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                  hidden.states=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod20, file=paste0(outname, "_mod20.rds"))
print('Model 20 is done!')

## Model 21. Heterogeneous diversification, not tied to range evolution.
## Assumes three distinct diversification rates.
## Dispersion parameters across hidden areas are the same.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,2,2,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, make.null=TRUE
                                    , include.jumps=TRUE, separate.extirpation=FALSE) 
mod21 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod21, file=paste0(outname, "_mod21.rds"))
print('Model 21 is done!')


## Model 22. Heterogeneous diversification, tied to range evolution. 
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE , 
                                    separate.extirpation=FALSE)
mod22 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation,
                  eps=extirpation, hidden.states=TRUE, 
                  trans.rate=trans.rate, assume.cladogenetic=TRUE)
saveRDS(mod22, file=paste0(outname, "_mod22.rds"))
print('Model 22 is done!')


## Model 23. Heterogeneous diversification, not tied to range evolution. 
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)) 
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)) 
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, make.null=TRUE, 
                                    include.jumps=TRUE, separate.extirpation=FALSE) 
mod23 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod23, file=paste0(outname, "_mod23.rds"))
print('Model 23 is done!')


## Model 24. Heterogeneous diversification, not tied to range evolution. 
## Assumes two distinct diversification rates.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE , 
                                    separate.extirpation=FALSE)
mod24 <- GeoHiSSE(phy, dat, f=sf, turnover=speciation
                  , eps=extirpation, hidden.states=TRUE
                  , trans.rate=trans.rate, assume.cladogenetic=TRUE)
saveRDS(mod24, file=paste0(outname, "_mod24.rds"))
print('Model 24 is done!')


###############################################################################
## Block of GeoSSE+extinction models.
## Here extirpation is NOT linked to range reduction.
## Range reduction is different from the extinction of an endemic lineage.
###############################################################################
## Model 25 - Dispersal parameters vary only, no range-dependent diversification. 
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                    , separate.extirpation=TRUE) 
mod25 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=FALSE , trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod25, file=paste0(outname, "_mod25.rds"))
print('Model 25 is done!')

## Model 26. GeoSSE model, with range effect on diversification speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                    , separate.extirpation=TRUE) 
mod26 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod26, file=paste0(outname, "_mod26.rds"))
print('Model 26 is done!')

## Model 27. Heterogeneous diversification, not tied to range evolution.
## Assumes three distinct diversification rates.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,2,2,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, make.null=TRUE,include.jumps=TRUE,
                                    separate.extirpation=TRUE)
mod27 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod27, file=paste0(outname, "_mod27.rds"))
print('Model 27 is done!')

## Model 28. Heterogeneous diversification, tied to range evolution.
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE, 
                                    separate.extirpation=TRUE) 
mod28 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod28, file=paste0(outname, "_mod28.rds"))
print('Model 28 is done!')

## Model 29. Heterogeneous diversification, not tied to range evolution. 
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)) 
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)) 
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, make.null=TRUE, include.jumps=TRUE,
                                    separate.extirpation=TRUE) 
mod29 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod29, file=paste0(outname, "_mod29.rds"))
print('Model 29 is done!')

## Model 30. Heterogeneous diversification, not tied to range evolution. 
## Assumes two distinct diversification rates.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE , 
                                    separate.extirpation=TRUE)
mod30 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE))
saveRDS(mod30, file=paste0(outname, "_mod30.rds"))
print('Model 30 is done!')

###############################################################################
## Block of anagenetic geographic models (MuSSE).
## Here models emulate GeoSSE (or GeoHiSSE) but changes only happen along branches.
###############################################################################
## Model 31. Transitions only. No character effect on diversification
speciation <- c(1,1,1)
extirpation <- c(1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                    , separate.extirpation=TRUE) 
mod31 <- try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, 
                  eps=extirpation, hidden.states=FALSE, trans.rate=trans.rate, 
                  assume.cladogenetic=FALSE))
saveRDS(mod31, file=paste0(outname, "_mod31.rds"))
print('Model 31 is done!')

## Model 32. Character effect on diversification.
speciation <- c(1,2,3)
extirpation <- c(1,2,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                    , separate.extirpation=TRUE) 
mod32 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation , eps=extirpation,
                  hidden.states=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod32, file=paste0(outname, "_mod32.rds"))
print('Model 32 is done!')

## Model 33. No character effect on diversification.
speciation <- c(1,1,1,2,2,2,3,3,3)
extirpation <- c(1,1,1,2,2,2,3,3,3)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, include.jumps=TRUE
                                    , separate.extirpation=TRUE, make.null=TRUE) 
mod33 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod33, file=paste0(outname, "_mod33.rds"))
print('Model 33 is done!')

## Model 34. Character effect on diversification, with a hidden state
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4,5,6)
trans.rate <- trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE
                                                  , separate.extirpation=TRUE)
mod34 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod34, file=paste0(outname, "_mod34.rds"))
print('Model 34 is done!')

## Model 35. No character effect on diversification, multiple shifts
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, include.jumps=TRUE,
                                    separate.extirpation=TRUE, make.null=TRUE) 
mod35 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation,
                  hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod35, file=paste0(outname, "_mod35.rds"))
print('Model 35 is done!')

## Model 36. No character effect on diversification, multiple shifts.
speciation <- c(rep(1,3), rep(2,3))
extirpation <- c(rep(1,3), rep(2,3))
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE, 
                                    separate.extirpation=TRUE, make.null=TRUE)
mod36 <-  try(GeoHiSSE(phy, dat, f=sf, turnover=speciation, eps=extirpation
                  , hidden.states=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE))
saveRDS(mod36, file=paste0(outname, "_mod36.rds"))
print('Model 36 is done!')