## ------------------------------------------------------------------------
library(ggplot2); library(rstan); library(doParallel)
#library(cowplot)
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

#Load the data
# load("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/na.sp/all.targets.jac.luh.final.het.con_not.all.census_impute.buffer_NRI.tyler_NPI.tyler_na.sp_HOM.1.3.RData")
load("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/na.sp/all.targets.jac.luh.final.het.con_not.all.census_impute.buffer_NRI.tyler_NPI.tyler_Na.sp_HOM1.3_nri.filter_new.nri.npi_0.5_impute.RData")
# load("/Volumes/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/na.sp/all.targets.jac.luh.final.het.con_not.all.census_impute.buffer_NRI.tyler_NPI.tyler_na.sp_HOM.1.3.RData")

#Rename the data to fit the Stan model
df <- all.targets.final
rm(all.targets.final)

#Subset the columns and rename to fit the Stan model
df <- subset(df, select = c(tag, DBH, Latin, growth, TrioML, nri.ses, Quadrat, Census, nc.con, nc.het, npi.ses))
colnames(df) <- c("tree", "dbh", "species", "growth", "inbr", "nri.ses", "quadrat", "census", "nc.con", "nc.het", "npi.ses")

#Subset to get just some species
# df <- df[df$species %in% c("Jacaranda copaia", "Luehea seemannii", "Cecropia insignis") ,]

#NC is a list for some reason. Transform it to numeric
df$nc.con <- unlist(df$nc.con)
df$nc.het <- unlist(df$nc.het)

# #Scale inbreeding first
# df$inbr2 <- df$inbr
# df$inbr2.scale <- scale(df$inbr2, center=-0.0001, scale=FALSE)
# qplot(df$inbr2.scale)
# #Then transform it
# df$inbr2.scale.logit <- log(na.pass(df$inbr2.scale)/(1-(na.pass(df$inbr2.scale))))
# qplot(df$inbr2.scale.logit)


# #Add a small amount from a gamma distribution
# df$inbr3 <- df$inbr
# df$inbr3 <- na.pass(df$inbr) + rgamma(nrow(df), 1, 10000)
# df$inbr3.logit <-  log(na.pass(df$inbr3)/(1-(na.pass(df$inbr3))))
# qplot(df$inbr3.logit)


#Transform inbreeding
df$inbr4 <- df$inbr
#Grab anything that is 0 and not NA and add a random value drawn from a uniform distribution bounded between 0 and the smallest non-zero inbreeding value
#df$inbr4 <- df$inbr4 + runif(1, 0, sort(unique(na.omit(df$inbr4)))[2])
df$inbr4 <- df$inbr4 + 0.0001
# qplot(df$inbr4)
#Transform it with logit function
df$inbr4.logit <- log(na.pass(df$inbr4)/(1-(na.pass(df$inbr4))))
#Look at the distribution
# qplot(df$inbr4.logit, binwidth=0.1)

#Transform NRI
# df$nri2 <- df$nri
# df$nri2 <- df$nri2 + 0.0001
# df$nri.logit <- log(na.pass(df$nri2)/(1-(na.pass(df$nri2))))

#Exclude targets with 0 and negative growth
df <- df[df$growth > 0 ,]

# #Exclude targets with NaN for NRI
# df <- df[!is.nan(df$nri) ,]

#Exclude targets with NC of infinity or 0
df <- df[df$nc.con != Inf & df$nc.con != -Inf ,]
df <- df[df$nc.het != Inf & df$nc.het != -Inf ,]

#Filter out when row is NA
df <- df[!is.na(df$species) ,]
#Filter out when npi.ses is NA
df <- df[!is.na(df$npi.ses) ,]
#Filter out with nri.ses is NA
df <- df[!is.na(df$nri.ses) ,]

#Add small amount to NC con and NRI
df$nc.con <- df$nc.con + 0.001


#############Log transform everything and check distribution for normality
df$log.dbh <- log(df$dbh)
df$log.nc.con <- log((df$nc.con))
df$log.nc.het <- log(df$nc.het)
df$log.growth <- log(df$growth)
# df$log.nri <- log(df$nri)
# df$log.npi <- log(df$npi)

#Make census a factor
df$census <- as.factor(df$census)

#Make numeric target species ID column
df$species.num[df$species == "Jacaranda copaia"] <- 1
df$species.num[df$species == "Luehea seemannii"] <- 2
df$species.num[df$species == "Cecropia insignis"] <- 3
df$species.num[df$species == "Tetragastris panamensis"] <- 4
df$species.num[df$species == "Triplaris cumingiana"] <- 5
df$species.num[df$species == "Virola sebifera"] <- 6

#Make numeric census column
df$census.num <- as.numeric(df$census)
#Make a census x species column to allow for the census intercept to differ by species
df$census.sp <- as.numeric(as.factor(paste0(df$species.num, df$census)))

#Define internal tree index
df$tree.ind <- 1:nrow(df)

#Quadrat.ind
df$quadrat.ind <- as.factor(df$quadrat)
levels(df$quadrat.ind) <- 1:length(unique(levels(df$quadrat.ind)))
df$quadrat.ind <- as.numeric(df$quadrat.ind)
length(unique(df$quadrat.ind))
#Quadrat.sp
df$quadrat.sp <- as.numeric(as.factor(paste0(df$species.num, df$quadrat.ind)))

#Find means and SDs of data
#Note that we are log transforming the growth, then standardizing it!
#Note that we also have to exlude the NAs in inbreeding b/c we can't take the mean and sd with them
###Standardize the data separately by species
#Jac
log.growthM.jac <- mean(df$log.growth[df$species == "Jacaranda copaia"]); log.growthSD.jac <- sd(df$log.growth[df$species == "Jacaranda copaia"])
dbhM.jac <- mean(df$dbh[df$species == "Jacaranda copaia"]); dbhSD.jac <- sd(df$dbh[df$species == "Jacaranda copaia"])
log.dbhM.jac <- mean(df$log.dbh[df$species == "Jacaranda copaia"]); log.dbhSD.jac <- sd(df$log.dbh[df$species == "Jacaranda copaia"])
inbr.transM.jac <- mean(na.exclude(df$inbr4.logit[df$species == "Jacaranda copaia"])); inbr.transSD.jac <- sd(na.exclude(df$inbr4.logit[df$species == "Jacaranda copaia"]))
inbrM.jac <- mean(na.exclude(df$inbr[df$species == "Jacaranda copaia"])); inbrSD.jac <- sd(na.exclude(df$inbr[df$species == "Jacaranda copaia"]))
# log.nri.ses.sesM.jac <- mean(na.exclude(df$log.nri.ses.ses[df$species == "Jacaranda copaia"])); log.nri.ses.sesSD.jac <- sd(na.exclude(df$log.nri.ses.ses[df$species == "Jacaranda copaia"]))
nri.sesM.jac <- mean(na.exclude(df$nri.ses[df$species == "Jacaranda copaia"])); nri.sesSD.jac <- sd(na.exclude(df$nri.ses[df$species == "Jacaranda copaia"]))
log.nc.conM.jac <- mean(df$log.nc.con[df$species == "Jacaranda copaia"]); log.nc.conSD.jac <- sd(df$log.nc.con[df$species == "Jacaranda copaia"])
log.nc.hetM.jac <- mean(df$log.nc.het[df$species == "Jacaranda copaia"]); log.nc.hetSD.jac <- sd(df$log.nc.het[df$species == "Jacaranda copaia"])
npi.sesM.jac <- mean(df$npi.ses[df$species == "Jacaranda copaia"]); npi.sesSD.jac <- sd(df$npi[df$species == "Jacaranda copaia"])
# log.npi.sesM.jac <- mean(df$log.npi.ses[df$species == "Jacaranda copaia"]); log.npi.sesSD.jac <- sd(df$log.npi.ses[df$species == "Jacaranda copaia"])
#Luh
log.growthM.luh <- mean(df$log.growth[df$species == "Luehea seemannii"]); log.growthSD.luh <- sd(df$log.growth[df$species == "Luehea seemannii"])
dbhM.luh <- mean(df$dbh[df$species == "Luehea seemannii"]); dbhSD.luh <- sd(df$dbh[df$species == "Luehea seemannii"])
log.dbhM.luh <- mean(df$log.dbh[df$species == "Luehea seemannii"]); log.dbhSD.luh <- sd(df$log.dbh[df$species == "Luehea seemannii"])
inbr.transM.luh <- mean(na.exclude(df$inbr4.logit[df$species == "Luehea seemannii"])); inbr.transSD.luh <- sd(na.exclude(df$inbr4.logit[df$species == "Luehea seemannii"]))
inbrM.luh <- mean(na.exclude(df$inbr[df$species == "Luehea seemannii"])); inbrSD.luh <- sd(na.exclude(df$inbr[df$species == "Luehea seemannii"]))
# log.nri.sesM.luh <- mean(na.exclude(df$log.nri.ses[df$species == "Luehea seemannii"])); log.nri.sesSD.luh <- sd(na.exclude(df$log.nri.ses[df$species == "Luehea seemannii"]))
nri.sesM.luh <- mean(na.exclude(df$nri.ses[df$species == "Luehea seemannii"])); nri.sesSD.luh <- sd(na.exclude(df$nri.ses[df$species == "Luehea seemannii"]))
log.nc.conM.luh <- mean(df$log.nc.con[df$species == "Luehea seemannii"]); log.nc.conSD.luh <- sd(df$log.nc.con[df$species == "Luehea seemannii"])
log.nc.hetM.luh <- mean(df$log.nc.het[df$species == "Luehea seemannii"]); log.nc.hetSD.luh <- sd(df$log.nc.het[df$species == "Luehea seemannii"])
npi.sesM.luh <- mean(df$npi.ses[df$species == "Luehea seemannii"]); npi.sesSD.luh <- sd(df$npi.ses[df$species == "Luehea seemannii"])
# log.npi.sesM.luh <- mean(df$log.npi.ses[df$species == "Luehea seemannii"]); log.npi.sesSD.luh <- sd(df$log.npi.ses[df$species == "Luehea seemannii"])
#Cec
log.growthM.cec <- mean(df$log.growth[df$species == "Cecropia insignis"]); log.growthSD.cec <- sd(df$log.growth[df$species == "Cecropia insignis"])
dbhM.cec <- mean(df$dbh[df$species == "Cecropia insignis"]); dbhSD.cec <- sd(df$dbh[df$species == "Cecropia insignis"])
log.dbhM.cec <- mean(df$log.dbh[df$species == "Cecropia insignis"]); log.dbhSD.cec <- sd(df$log.dbh[df$species == "Cecropia insignis"])
inbr.transM.cec <- mean(na.exclude(df$inbr4.logit[df$species == "Cecropia insignis"])); inbr.transSD.cec <- sd(na.exclude(df$inbr4.logit[df$species == "Cecropia insignis"]))
inbrM.cec <- mean(na.exclude(df$inbr[df$species == "Cecropia insignis"])); inbrSD.cec <- sd(na.exclude(df$inbr[df$species == "Cecropia insignis"]))
# log.nri.sesM.cec <- mean(na.exclude(df$log.nri.ses[df$species == "Cecropia insignis"])); log.nri.sesSD.cec <- sd(na.exclude(df$log.nri.ses[df$species == "Cecropia insignis"]))
nri.sesM.cec <- mean(na.exclude(df$nri.ses[df$species == "Cecropia insignis"])); nri.sesSD.cec <- sd(na.exclude(df$nri.ses[df$species == "Cecropia insignis"]))
log.nc.conM.cec <- mean(df$log.nc.con[df$species == "Cecropia insignis"]); log.nc.conSD.cec <- sd(df$log.nc.con[df$species == "Cecropia insignis"])
log.nc.hetM.cec <- mean(df$log.nc.het[df$species == "Cecropia insignis"]); log.nc.hetSD.cec <- sd(df$log.nc.het[df$species == "Cecropia insignis"])
npi.sesM.cec <- mean(df$npi.ses[df$species == "Cecropia insignis"]); npi.sesSD.cec <- sd(df$npi.ses[df$species == "Cecropia insignis"])
# log.npi.sesM.cec <- mean(df$log.npi.ses[df$species == "Cecropia insignis"]); log.npi.sesSD.cec <- sd(df$log.npi.ses[df$species == "Cecropia insignis"])
#Tet
log.growthM.tet <- mean(df$log.growth[df$species == "Tetragastris panamensis"]); log.growthSD.tet <- sd(df$log.growth[df$species == "Tetragastris panamensis"])
dbhM.tet <- mean(df$dbh[df$species == "Tetragastris panamensis"]); dbhSD.tet <- sd(df$dbh[df$species == "Tetragastris panamensis"])
log.dbhM.tet <- mean(df$log.dbh[df$species == "Tetragastris panamensis"]); log.dbhSD.tet <- sd(df$log.dbh[df$species == "Tetragastris panamensis"])
inbr.transM.tet <- mean(na.exclude(df$inbr4.logit[df$species == "Tetragastris panamensis"])); inbr.transSD.tet <- sd(na.exclude(df$inbr4.logit[df$species == "Tetragastris panamensis"]))
inbrM.tet <- mean(na.exclude(df$inbr[df$species == "Tetragastris panamensis"])); inbrSD.tet <- sd(na.exclude(df$inbr[df$species == "Tetragastris panamensis"]))
# log.nri.sesM.tet <- mean(na.exclude(df$log.nri.ses[df$species == "Tetragastris panamensis"])); log.nri.sesSD.tet <- sd(na.exclude(df$log.nri.ses[df$species == "Tetragastris panamensis"]))
nri.sesM.tet <- mean(na.exclude(df$nri.ses[df$species == "Tetragastris panamensis"])); nri.sesSD.tet <- sd(na.exclude(df$nri.ses[df$species == "Tetragastris panamensis"]))
log.nc.conM.tet <- mean(df$log.nc.con[df$species == "Tetragastris panamensis"]); log.nc.conSD.tet <- sd(df$log.nc.con[df$species == "Tetragastris panamensis"])
log.nc.hetM.tet <- mean(df$log.nc.het[df$species == "Tetragastris panamensis"]); log.nc.hetSD.tet <- sd(df$log.nc.het[df$species == "Tetragastris panamensis"])
npi.sesM.tet <- mean(df$npi.ses[df$species == "Tetragastris panamensis"]); npi.sesSD.tet <- sd(df$npi.ses[df$species == "Tetragastris panamensis"])
# log.npi.sesM.tet <- mean(df$log.npi.ses[df$species == "Tetragastris panamensis"]); log.npi.sesSD.tet <- sd(df$log.npi.ses[df$species == "Tetragastris panamensis"])
#Tri
log.growthM.tri <- mean(df$log.growth[df$species == "Triplaris cumingiana"]); log.growthSD.tri <- sd(df$log.growth[df$species == "Triplaris cumingiana"])
dbhM.tri <- mean(df$dbh[df$species == "Triplaris cumingiana"]); dbhSD.tri <- sd(df$dbh[df$species == "Triplaris cumingiana"])
log.dbhM.tri <- mean(df$log.dbh[df$species == "Triplaris cumingiana"]); log.dbhSD.tri <- sd(df$log.dbh[df$species == "Triplaris cumingiana"])
inbr.transM.tri <- mean(na.exclude(df$inbr4.logit[df$species == "Triplaris cumingiana"])); inbr.transSD.tri <- sd(na.exclude(df$inbr4.logit[df$species == "Triplaris cumingiana"]))
inbrM.tri <- mean(na.exclude(df$inbr[df$species == "Triplaris cumingiana"])); inbrSD.tri <- sd(na.exclude(df$inbr[df$species == "Triplaris cumingiana"]))
# log.nri.sesM.tri <- mean(na.exclude(df$log.nri.ses[df$species == "Triplaris cumingiana"])); log.nri.sesSD.tri <- sd(na.exclude(df$log.nri.ses[df$species == "Triplaris cumingiana"]))
nri.sesM.tri <- mean(na.exclude(df$nri.ses[df$species == "Triplaris cumingiana"])); nri.sesSD.tri <- sd(na.exclude(df$nri.ses[df$species == "Triplaris cumingiana"]))
log.nc.conM.tri <- mean(df$log.nc.con[df$species == "Triplaris cumingiana"]); log.nc.conSD.tri <- sd(df$log.nc.con[df$species == "Triplaris cumingiana"])
log.nc.hetM.tri <- mean(df$log.nc.het[df$species == "Triplaris cumingiana"]); log.nc.hetSD.tri <- sd(df$log.nc.het[df$species == "Triplaris cumingiana"])
npi.sesM.tri <- mean(df$npi.ses[df$species == "Triplaris cumingiana"]); npi.sesSD.tri <- sd(df$npi.ses[df$species == "Triplaris cumingiana"])
# log.npi.sesM.tri <- mean(df$log.npi.ses[df$species == "Triplaris cumingiana"]); log.npi.sesSD.tri <- sd(df$log.npi.ses[df$species == "Triplaris cumingiana"])
#Vir
log.growthM.vir <- mean(df$log.growth[df$species == "Virola sebifera"]); log.growthSD.vir <- sd(df$log.growth[df$species == "Virola sebifera"])
dbhM.vir <- mean(df$dbh[df$species == "Virola sebifera"]); dbhSD.vir <- sd(df$dbh[df$species == "Virola sebifera"])
log.dbhM.vir <- mean(df$log.dbh[df$species == "Virola sebifera"]); log.dbhSD.vir <- sd(df$log.dbh[df$species == "Virola sebifera"])
inbr.transM.vir <- mean(na.exclude(df$inbr4.logit[df$species == "Virola sebifera"])); inbr.transSD.vir <- sd(na.exclude(df$inbr4.logit[df$species == "Virola sebifera"]))
inbrM.vir <- mean(na.exclude(df$inbr[df$species == "Virola sebifera"])); inbrSD.vir <- sd(na.exclude(df$inbr[df$species == "Virola sebifera"]))
# log.nri.sesM.vir <- mean(na.exclude(df$log.nri.ses[df$species == "Virola sebifera"])); log.nri.sesSD.vir <- sd(na.exclude(df$log.nri.ses[df$species == "Virola sebifera"]))
nri.sesM.vir <- mean(na.exclude(df$nri.ses[df$species == "Virola sebifera"])); nri.sesSD.vir <- sd(na.exclude(df$nri.ses[df$species == "Virola sebifera"]))
log.nc.conM.vir <- mean(df$log.nc.con[df$species == "Virola sebifera"]); log.nc.conSD.vir <- sd(df$log.nc.con[df$species == "Virola sebifera"])
log.nc.hetM.vir <- mean(df$log.nc.het[df$species == "Virola sebifera"]); log.nc.hetSD.vir <- sd(df$log.nc.het[df$species == "Virola sebifera"])
npi.sesM.vir <- mean(df$npi.ses[df$species == "Virola sebifera"]); npi.sesSD.vir <- sd(df$npi.ses[df$species == "Virola sebifera"])
# log.npi.sesM.vir <- mean(df$log.npi.ses[df$species == "Virola sebifera"]); log.npi.sesSD.vir <- sd(df$log.npi.ses[df$species == "Virola sebifera"])


#Create standardized data columns
#Jac
df$log.growth.std[df$species == "Jacaranda copaia"] <- (df$log.growth[df$species == "Jacaranda copaia"] - log.growthM.jac) / log.growthSD.jac
df$dbh.std[df$species == "Jacaranda copaia"] <- (df$dbh[df$species == "Jacaranda copaia"] - dbhM.jac) / dbhSD.jac
df$log.dbh.std[df$species == "Jacaranda copaia"] <- (df$log.dbh[df$species == "Jacaranda copaia"] - log.dbhM.jac) / log.dbhSD.jac
df$inbr.trans.std[df$species == "Jacaranda copaia"] <- (df$inbr4.logit[df$species == "Jacaranda copaia"] - inbr.transM.jac) / inbr.transSD.jac
df$inbr.std[df$species == "Jacaranda copaia"] <- (df$inbr[df$species == "Jacaranda copaia"] - inbrM.jac) / inbrSD.jac
# df$log.nri.ses.std[df$species == "Jacaranda copaia"] <- (df$log.nri.ses[df$species == "Jacaranda copaia"] - log.nri.sesM.jac) / log.nri.sesSD.jac
df$log.nc.con.std[df$species == "Jacaranda copaia"] <- (df$log.nc.con[df$species == "Jacaranda copaia"] - log.nc.conM.jac) / log.nc.conSD.jac
df$log.nc.het.std[df$species == "Jacaranda copaia"] <- (df$log.nc.het[df$species == "Jacaranda copaia"] - log.nc.hetM.jac) / log.nc.hetSD.jac
df$npi.ses.std[df$species == "Jacaranda copaia"] <- (df$npi.ses[df$species == "Jacaranda copaia"] - npi.sesM.jac) / npi.sesSD.jac
# df$log.npi.ses.std[df$species == "Jacaranda copaia"] <- (df$log.npi.ses[df$species == "Jacaranda copaia"] - log.npi.sesM.jac) / log.npi.sesSD.jac
df$nri.ses.std[df$species == "Jacaranda copaia"] <- (df$nri.ses[df$species == "Jacaranda copaia"] - nri.sesM.jac) / nri.sesSD.jac

#Luh
df$log.growth.std[df$species == "Luehea seemannii"] <- (df$log.growth[df$species == "Luehea seemannii"] - log.growthM.luh) / log.growthSD.luh
df$dbh.std[df$species == "Luehea seemannii"] <- (df$dbh[df$species == "Luehea seemannii"] - dbhM.luh) / dbhSD.luh
df$log.dbh.std[df$species == "Luehea seemannii"] <- (df$log.dbh[df$species == "Luehea seemannii"] - log.dbhM.luh) / log.dbhSD.luh
df$inbr.trans.std[df$species == "Luehea seemannii"] <- (df$inbr4.logit[df$species == "Luehea seemannii"] - inbr.transM.luh) / inbr.transSD.luh
df$inbr.std[df$species == "Luehea seemannii"] <- (df$inbr[df$species == "Luehea seemannii"] - inbrM.luh) / inbrSD.luh
# df$log.nri.ses.std[df$species == "Luehea seemannii"] <- (df$log.nri.ses[df$species == "Luehea seemannii"] - log.nri.sesM.luh) / log.nri.sesSD.luh
df$log.nc.con.std[df$species == "Luehea seemannii"] <- (df$log.nc.con[df$species == "Luehea seemannii"] - log.nc.conM.luh) / log.nc.conSD.luh
df$log.nc.het.std[df$species == "Luehea seemannii"] <- (df$log.nc.het[df$species == "Luehea seemannii"] - log.nc.hetM.luh) / log.nc.hetSD.luh
df$npi.ses.std[df$species == "Luehea seemannii"] <- (df$npi.ses[df$species == "Luehea seemannii"] - npi.sesM.luh) / npi.sesSD.luh
# df$log.npi.ses.std[df$species == "Luehea seemannii"] <- (df$log.npi.ses[df$species == "Luehea seemannii"] - log.npi.sesM.luh) / log.npi.sesSD.luh
df$nri.ses.std[df$species == "Luehea seemannii"] <- (df$nri.ses[df$species == "Luehea seemannii"] - nri.sesM.luh) / nri.sesSD.luh

#Cec
df$log.growth.std[df$species == "Cecropia insignis"] <- (df$log.growth[df$species == "Cecropia insignis"] - log.growthM.cec) / log.growthSD.cec
df$dbh.std[df$species == "Cecropia insignis"] <- (df$dbh[df$species == "Cecropia insignis"] - dbhM.cec) / dbhSD.cec
df$log.dbh.std[df$species == "Cecropia insignis"] <- (df$log.dbh[df$species == "Cecropia insignis"] - log.dbhM.cec) / log.dbhSD.cec
df$inbr.trans.std[df$species == "Cecropia insignis"] <- (df$inbr4.logit[df$species == "Cecropia insignis"] - inbr.transM.cec) / inbr.transSD.cec
df$inbr.std[df$species == "Cecropia insignis"] <- (df$inbr[df$species == "Cecropia insignis"] - inbrM.cec) / inbrSD.cec
# df$log.nri.ses.std[df$species == "Cecropia insignis"] <- (df$log.nri.ses[df$species == "Cecropia insignis"] - log.nri.sesM.cec) / log.nri.sesSD.cec
df$log.nc.con.std[df$species == "Cecropia insignis"] <- (df$log.nc.con[df$species == "Cecropia insignis"] - log.nc.conM.cec) / log.nc.conSD.cec
df$log.nc.het.std[df$species == "Cecropia insignis"] <- (df$log.nc.het[df$species == "Cecropia insignis"] - log.nc.hetM.cec) / log.nc.hetSD.cec
df$npi.ses.std[df$species == "Cecropia insignis"] <- (df$npi.ses[df$species == "Cecropia insignis"] - npi.sesM.cec) / npi.sesSD.cec
# df$log.npi.ses.std[df$species == "Cecropia insignis"] <- (df$log.npi.ses[df$species == "Cecropia insignis"] - log.npi.sesM.cec) / log.npi.sesSD.cec
df$nri.ses.std[df$species == "Cecropia insignis"] <- (df$nri.ses[df$species == "Cecropia insignis"] - nri.sesM.cec) / nri.sesSD.cec

#Tet
df$log.growth.std[df$species == "Tetragastris panamensis"] <- (df$log.growth[df$species == "Tetragastris panamensis"] - log.growthM.tet) / log.growthSD.tet
df$dbh.std[df$species == "Tetragastris panamensis"] <- (df$dbh[df$species == "Tetragastris panamensis"] - dbhM.tet) / dbhSD.tet
df$log.dbh.std[df$species == "Tetragastris panamensis"] <- (df$log.dbh[df$species == "Tetragastris panamensis"] - log.dbhM.tet) / log.dbhSD.tet
df$inbr.trans.std[df$species == "Tetragastris panamensis"] <- (df$inbr4.logit[df$species == "Tetragastris panamensis"] - inbr.transM.tet) / inbr.transSD.tet
df$inbr.std[df$species == "Tetragastris panamensis"] <- (df$inbr[df$species == "Tetragastris panamensis"] - inbrM.tet) / inbrSD.tet
# df$log.nri.ses.std[df$species == "Tetragastris panamensis"] <- (df$log.nri.ses[df$species == "Tetragastris panamensis"] - log.nri.sesM.tet) / log.nri.sesSD.tet
df$log.nc.con.std[df$species == "Tetragastris panamensis"] <- (df$log.nc.con[df$species == "Tetragastris panamensis"] - log.nc.conM.tet) / log.nc.conSD.tet
df$log.nc.het.std[df$species == "Tetragastris panamensis"] <- (df$log.nc.het[df$species == "Tetragastris panamensis"] - log.nc.hetM.tet) / log.nc.hetSD.tet
df$npi.ses.std[df$species == "Tetragastris panamensis"] <- (df$npi.ses[df$species == "Tetragastris panamensis"] - npi.sesM.tet) / npi.sesSD.tet
# df$log.npi.ses.std[df$species == "Tetragastris panamensis"] <- (df$log.npi.ses[df$species == "Tetragastris panamensis"] - log.npi.sesM.tet) / log.npi.sesSD.tet
df$nri.ses.std[df$species == "Tetragastris panamensis"] <- (df$nri.ses[df$species == "Tetragastris panamensis"] - nri.sesM.tet) / nri.sesSD.tet

#Tri
df$log.growth.std[df$species == "Triplaris cumingiana"] <- (df$log.growth[df$species == "Triplaris cumingiana"] - log.growthM.tri) / log.growthSD.tri
df$dbh.std[df$species == "Triplaris cumingiana"] <- (df$dbh[df$species == "Triplaris cumingiana"] - dbhM.tri) / dbhSD.tri
df$log.dbh.std[df$species == "Triplaris cumingiana"] <- (df$log.dbh[df$species == "Triplaris cumingiana"] - log.dbhM.tri) / log.dbhSD.tri
df$inbr.trans.std[df$species == "Triplaris cumingiana"] <- (df$inbr4.logit[df$species == "Triplaris cumingiana"] - inbr.transM.tri) / inbr.transSD.tri
df$inbr.std[df$species == "Triplaris cumingiana"] <- (df$inbr[df$species == "Triplaris cumingiana"] - inbrM.tri) / inbrSD.tri
# df$log.nri.ses.std[df$species == "Triplaris cumingiana"] <- (df$log.nri.ses[df$species == "Triplaris cumingiana"] - log.nri.sesM.tri) / log.nri.sesSD.tri
df$log.nc.con.std[df$species == "Triplaris cumingiana"] <- (df$log.nc.con[df$species == "Triplaris cumingiana"] - log.nc.conM.tri) / log.nc.conSD.tri
df$log.nc.het.std[df$species == "Triplaris cumingiana"] <- (df$log.nc.het[df$species == "Triplaris cumingiana"] - log.nc.hetM.tri) / log.nc.hetSD.tri
df$npi.ses.std[df$species == "Triplaris cumingiana"] <- (df$npi.ses[df$species == "Triplaris cumingiana"] - npi.sesM.tri) / npi.sesSD.tri
# df$log.npi.ses.std[df$species == "Triplaris cumingiana"] <- (df$log.npi.ses[df$species == "Triplaris cumingiana"] - log.npi.sesM.tri) / log.npi.sesSD.tri
df$nri.ses.std[df$species == "Triplaris cumingiana"] <- (df$nri.ses[df$species == "Triplaris cumingiana"] - nri.sesM.tri) / nri.sesSD.tri

#Vir
df$log.growth.std[df$species == "Virola sebifera"] <- (df$log.growth[df$species == "Virola sebifera"] - log.growthM.vir) / log.growthSD.vir
df$dbh.std[df$species == "Virola sebifera"] <- (df$dbh[df$species == "Virola sebifera"] - dbhM.vir) / dbhSD.vir
df$log.dbh.std[df$species == "Virola sebifera"] <- (df$log.dbh[df$species == "Virola sebifera"] - log.dbhM.vir) / log.dbhSD.vir
df$inbr.trans.std[df$species == "Virola sebifera"] <- (df$inbr4.logit[df$species == "Virola sebifera"] - inbr.transM.vir) / inbr.transSD.vir
df$inbr.std[df$species == "Virola sebifera"] <- (df$inbr[df$species == "Virola sebifera"] - inbrM.vir) / inbrSD.vir
# df$log.nri.ses.std[df$species == "Virola sebifera"] <- (df$log.nri.ses[df$species == "Virola sebifera"] - log.nri.sesM.vir) / log.nri.sesSD.vir
df$log.nc.con.std[df$species == "Virola sebifera"] <- (df$log.nc.con[df$species == "Virola sebifera"] - log.nc.conM.vir) / log.nc.conSD.vir
df$log.nc.het.std[df$species == "Virola sebifera"] <- (df$log.nc.het[df$species == "Virola sebifera"] - log.nc.hetM.vir) / log.nc.hetSD.vir
df$npi.ses.std[df$species == "Virola sebifera"] <- (df$npi.ses[df$species == "Virola sebifera"] - npi.sesM.vir) / npi.sesSD.vir
# df$log.npi.ses.std[df$species == "Virola sebifera"] <- (df$log.npi.ses[df$species == "Virola sebifera"] - log.npi.sesM.vir) / log.npi.sesSD.vir
df$nri.ses.std[df$species == "Virola sebifera"] <- (df$nri.ses[df$species == "Virola sebifera"] - nri.sesM.vir) / nri.sesSD.vir


######################################
#Split off data with missing inbreeding
######################################
# df.mis <- df[is.na(df$inbr.trans) ,]
# #Filter out missing from regular data
df.not.mis <- df[!is.na(df$inbr.trans.std) ,]
# #Make a df.all with first non-missing and then missing trees for the quad and census effects to index. This will be accompanied by adding the number of non-missing trees to the index of the missing trees in the likelihood.
# df.all <- rbind(df.not.mis, df.mis)

#Drop Luehea and Virola b/c they have so few individuals
df.not.mis <- df.not.mis[!(df.not.mis$species %in% c("Luehea seemannii", "Virola sebifera")) ,]

#Make census a factor
df.not.mis$census <- as.factor(df.not.mis$census)

#Make numeric target species ID column
df.not.mis$species.num[df.not.mis$species == "Jacaranda copaia"] <- 1
df.not.mis$species.num[df.not.mis$species == "Cecropia insignis"] <- 2
df.not.mis$species.num[df.not.mis$species == "Triplaris cumingiana"] <- 3

#Make numeric census column
df.not.mis$census.num <- as.numeric(df.not.mis$census)
#Make a census x species column to allow for the census intercept to differ by species
df.not.mis$census.sp <- as.numeric(as.factor(paste0(df.not.mis$species.num, df.not.mis$census)))

#Define internal tree index
df.not.mis$tree.ind <- 1:nrow(df.not.mis)

#Make numeric tree tag column for random tree intercept
df.not.mis$tree.num[df.not.mis$species == "Jacaranda copaia"] <- as.numeric(as.factor(df.not.mis$tree[df.not.mis$species == "Jacaranda copaia"]))
df.not.mis$tree.num[df.not.mis$species == "Cecropia insignis"] <- (as.numeric(as.factor(df.not.mis$tree[df.not.mis$species == "Cecropia insignis"])))+length(unique(df.not.mis$tree[df.not.mis$species == "Jacaranda copaia"]))
df.not.mis$tree.num[df.not.mis$species == "Triplaris cumingiana"] <- (as.numeric(as.factor(df.not.mis$tree[df.not.mis$species == "Triplaris cumingiana"])))+length(unique(df.not.mis$tree[df.not.mis$species == "Jacaranda copaia"]))+length(unique(df.not.mis$tree[df.not.mis$species == "Cecropia insignis"]))

#Make a tree observation-level column for individual error terms
df.not.mis$obs.num <- 1:nrow(df.not.mis)

#Quadrat.ind
df.not.mis$quadrat.ind <- as.factor(df.not.mis$quadrat)
levels(df.not.mis$quadrat.ind) <- 1:length(unique(levels(df.not.mis$quadrat.ind)))
df.not.mis$quadrat.ind <- as.numeric(df.not.mis$quadrat.ind)
length(unique(df.not.mis$quadrat.ind))
#Quadrat.sp
df.not.mis$quadrat.sp <- as.numeric(as.factor(paste0(df.not.mis$species.num, df.not.mis$quadrat.ind)))

jac.sub <- df.not.mis[df.not.mis$species == "Jacaranda copaia" ,]

qplot(log.dbh.std, log.growth, color=species, data=df.not.mis)+geom_smooth(method = "lm")
qplot(log.nc.het.std, log.growth, color=species, data=df.not.mis)+geom_smooth(method="lm")
qplot(log.nc.con.std, log.growth, color=species, data=df.not.mis)+geom_smooth(method="lm")
qplot(nri.ses.std, log.growth, color=species, data=df.not.mis)+geom_smooth(method = "lm")
qplot(npi.ses.std, log.growth, color=species, data=df.not.mis)+geom_smooth(method = "lm")


#############################################
#Split up data by levels for multilevel model
#############################################

###Observation regressors
X <- as.matrix(subset(df.not.mis, select = c(log.dbh.std, log.nc.het.std, log.nc.con.std, nri.ses.std, npi.ses.std)))

###Individual tree regressors
#U <- as.matrix(cbind(rep(1, length(unique(df.not.mis$tree))), rep(1, length(unique(df.not.mis$tree)))))
U <- NULL
for (i in 1:length(unique(df.not.mis$tree))) {
  U[i] <- as.vector(unique(df.not.mis$inbr.trans.std[df.not.mis$tree == unique(df.not.mis$tree)[i]]))
}
#Make species vector for individual tree-level data
sp.tree <- NULL
for (i in 1:length(unique(df.not.mis$tree))) {
  sp.tree[i] <- as.vector(unique(df.not.mis$species.num[df.not.mis$tree == unique(df.not.mis$tree)[i]]))
}

###Quadrat regressors
library(gdata)
env <- read.xls("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/Soils\ bci.block20.data_TS.xls", sheet="Data")
# env <- read.xls("/Volumes/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/Soils\ bci.block20.data_TS.xls", sheet="Data")
# env <- read.csv("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/multilevel/Soils\ bci.block20.data_TS.csv")
census.2010 <- read.table("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/Census2010.txt", sep="\t", header=TRUE)
# census.2010 <- read.table("/Volumes/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/Census2010.txt", sep="\t", header=TRUE)
# census.2010 <- read.table("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/Census2010.txt", sep="\t", header=TRUE)

#Get the sequential quadrat numbers in df
df.not.mis.quad.seq <- NULL
for (i in 1:nrow(df.not.mis)) {
  df.not.mis.quad.seq[i] <- env$subquadrat[env$non.seq.quadrat == df.not.mis$quadrat[i]]
}
df.not.mis$quad.seq <- df.not.mis.quad.seq
#Add 1 because these quad.seqs start at 0, but we want them to match the row number in the V matrix below
df.not.mis$quad.seq.row <- df.not.mis$quad.seq + 1

#Log transform the nutrient variables (columns 26:36, 38) -- drop N.min
log.nutrient <- log(subset(env, select = c(Al, B, Ca, Cu, Fe, K, Mg, Mn, P, Zn, N, pH)))
#Do PCA
nutrient.pca <- prcomp(log.nutrient,
                       center = TRUE,
                       scale. = TRUE)
#Grab the PCA results and put it into the target tree dataframe
env <- cbind(env, nutrient.pca$x[,1:3])
#Subset to get only those quadrats in df
env.sub <- env[env$non.seq.quadrat %in% df.not.mis$quadrat ,]
#Number these quadrats from 1 to whatever to be able to index them
env.sub$quad.index <- 1:nrow(env.sub)
#Create matrix of predictors with first column as intercept
V <- as.matrix(subset(env.sub, select = c(PC1, PC2)))
#Make rownames as quadrat names from raw census data (not sequential)
rownames(V) <- env$quadrat
#Note: The PC values for these are going to be different because the PC analysis uses all of the quadrats instead of just a subset.

###Quadrat indicator for each tree
#Need to get the q vector as the quadrat association for each unique tree
#Grab the unsequential quadrat names that each unique tree is located in
q.temp <- NULL
for (i in 1:length(unique(df.not.mis$tree))) {
  q.temp[i] <- unique(df.not.mis$quadrat[df.not.mis$tree == unique(df.not.mis$tree)[i]])
}
#Now grab the quadrat index from the env.sub (V matrix) using these unsequential quadrat names
q.temp2 <- NULL
for (i in 1:length(q.temp)) {
  q.temp2[i] <- env.sub$quad.index[env.sub$non.seq.quadrat == q.temp[i]]
}
q <- q.temp2

###Separate quadrat indicators for each tree by species
#Jac
#Grab the unsequential quadrat names that each unique tree is located in
q_1.temp <- NULL
for (i in 1:length(unique(df.not.mis$tree[df.not.mis$species == "Jacaranda copaia"]))) {
  q_1.temp[i] <- unique(df.not.mis$quadrat[df.not.mis$tree == unique(df.not.mis$tree[df.not.mis$species == "Jacaranda copaia"])[i]])
}
#Now grab the quadrat index from the env.sub (V matrix) using these unsequential quadrat names
q_1.temp2 <- NULL
for (i in 1:length(q_1.temp)) {
  q_1.temp2[i] <- env.sub$quad.index[env.sub$non.seq.quadrat == q_1.temp[i]]
}
q_1 <- q_1.temp2

#Cec
q_2.temp <- NULL
for (i in 1:length(unique(df.not.mis$tree[df.not.mis$species == "Cecropia insignis"]))) {
  q_2.temp[i] <- unique(df.not.mis$quadrat[df.not.mis$tree == unique(df.not.mis$tree[df.not.mis$species == "Cecropia insignis"])[i]])
}
#Now grab the quadrat index from the env.sub (V matrix) using these unsequential quadrat names
q_2.temp2 <- NULL
for (i in 1:length(q_2.temp)) {
  q_2.temp2[i] <- env.sub$quad.index[env.sub$non.seq.quadrat == q_2.temp[i]]
}
q_2 <- q_2.temp2

#Tri
q_3.temp <- NULL
for (i in 1:length(unique(df.not.mis$tree[df.not.mis$species == "Triplaris cumingiana"]))) {
  q_3.temp[i] <- unique(df.not.mis$quadrat[df.not.mis$tree == unique(df.not.mis$tree[df.not.mis$species == "Triplaris cumingiana"])[i]])
}
#Now grab the quadrat index from the env.sub (V matrix) using these unsequential quadrat names
q_3.temp2 <- NULL
for (i in 1:length(q_3.temp)) {
  q_3.temp2[i] <- env.sub$quad.index[env.sub$non.seq.quadrat == q_3.temp[i]]
}
q_3 <- q_3.temp2


# test <- subset(df.not.mis, select = c(species, species.num))
# 
# min(df.not.mis$tree.num[df.not.mis$species == "Jacaranda copaia"])
# max(df.not.mis$tree.num[df.not.mis$species == "Jacaranda copaia"])
# max(df.not.mis$tree.num[df.not.mis$species == "Jacaranda copaia"])-min(df.not.mis$tree.num[df.not.mis$species == "Jacaranda copaia"])
# min(df.not.mis$tree.num[df.not.mis$species == "Luehea seemannii"])
# max(df.not.mis$tree.num[df.not.mis$species == "Luehea seemannii"])
# max(df.not.mis$tree.num[df.not.mis$species == "Luehea seemannii"])-min(df.not.mis$tree.num[df.not.mis$species == "Luehea seemannii"])
# min(df.not.mis$tree.num[df.not.mis$species == "Cecropia insignis"])
# max(df.not.mis$tree.num[df.not.mis$species == "Cecropia insignis"])
# max(df.not.mis$tree.num[df.not.mis$species == "Cecropia insignis"])-min(df.not.mis$tree.num[df.not.mis$species == "Cecropia insignis"])
# min(df.not.mis$tree.num[df.not.mis$species == "Tetragastris panamensis"])
# max(df.not.mis$tree.num[df.not.mis$species == "Tetragastris panamensis"])
# max(df.not.mis$tree.num[df.not.mis$species == "Tetragastris panamensis"])-min(df.not.mis$tree.num[df.not.mis$species == "Tetragastris panamensis"])

### Define the data for stan
data_list <- list(
  #Observation level data
  N <- nrow(X),
  K <- ncol(X),
  y <- df.not.mis$log.growth,
  X <- X,
  j <- df.not.mis$tree.num,
  
  Nsp <- length(unique(df.not.mis$species.num)),
  sp_obs <- df.not.mis$species.num,
  
  #Individual tree level data
  J <- length(unique(df.not.mis$tree.num)),
  J_1 <- length(unique(df.not.mis$tree.num[df.not.mis$species == "Jacaranda copaia"])),
  J_2 <- length(unique(df.not.mis$tree.num[df.not.mis$species == "Cecropia insignis"])),
  J_3 <- length(unique(df.not.mis$tree.num[df.not.mis$species == "Triplaris cumingiana"])),
  J_1_end <- length(unique(df.not.mis$tree.num[df.not.mis$species == "Jacaranda copaia"])),
  J_2_end <- length(unique(df.not.mis$tree.num[df.not.mis$species == "Jacaranda copaia"])) + length(unique(df.not.mis$tree.num[df.not.mis$species == "Cecropia insignis"])),
  J_3_end <- length(unique(df.not.mis$tree.num[df.not.mis$species == "Jacaranda copaia"])) + length(unique(df.not.mis$tree.num[df.not.mis$species == "Cecropia insignis"])) + length(unique(df.not.mis$tree.num[df.not.mis$species == "Triplaris cumingiana"])),
  J_2_start <- J_1_end + 1,
  J_3_start <- J_2_end + 1,
  L <- 1,
  q <- q,
  q_1 <- q_1,
  q_2 <- q_2,
  q_3 <- q_3,
  U <- U,
  
  sp_tree <- sp.tree,
  
  #Quadrat level data
  Q <- length(unique(df.not.mis$quad.seq.row)),
  M <- ncol(V),
  V <- (V)
)


#################################
##Run the model
#################################
# fit <- stan(file = '/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/multilevel/growth/na.sp/multilevel_growth_1-6_NRI_log.lik_NPI_na.sp.stan', data = data_list, iter = 4000, chains = 3, control = list(adapt_delta = 0.9999, stepsize = 0.003))
fit <- stan(file = '/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/na.sp/multilevel_growth_1-6_NRI_log.lik_NPI_na.sp_0.5.stan', data = data_list, iter = 1000, chains = 3, control = list(adapt_delta = 0.9999, stepsize = 0.003, max_treedepth = 15))
# fit <- stan(file = '/Volumes/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/na.sp/multilevel_growth_1-6_NRI_log.lik_NPI_na.sp.stan', data = data_list, iter = 500, chains = 3, control = list(adapt_delta = 0.9999, stepsize = 0.003, max_treedepth = 15))

#Current run:  No intercepts with means parameterized with errors explicitly
#Note that Gelman and Hill recommend using an intercept at one level of the model

#Random run number
run.number <- sample(1:1000, 1)

## ------------------------------------------------------------------------
# save(fit, file=paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/na.sp/growth_1-6_NRI_2000iter_log.lik_NRI.tyler_NPI.tyler_na.sp_HOM1.3_new.nri.npi_0.2", run.number, ".RData"))
# save(fit, file=paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/multilevel/growth/na.sp/growth_na.sp_HOM.1.3_tree.num.fix_4000iter_", run.number, ".RData"))
# load("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/na.sp/growth_1-6_NRI_2000iter_log.lik_NRI.tyler_NPI.tyler_na.sp_HOM1.3_new.nri.npi_0.2935.RData")

# #Save traceplot
# #png(width=1000, height=800, file=paste0("/Volumes/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_missing_growth.no.std_not.centered4_not.chen_figs/trace.plot", run.number, ".png"))
# png(width=1000, height=800, file=paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/trace.plot", run.number, ".png"))
traceplot(fit, pars = c("B"))
traceplot(fit, pars = c("G"))
traceplot(fit, pars = c("D"))

plot(fit, pars = c("B"))

# dev.off()
#
# #Save all plot
# #png(width=1000, height=800, file=paste0("/Volumes/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_missing_growth.no.std_not.centered4_not.chen_figs/all.plot", run.number, ".png"))
# png(width=1000, height=800, file=paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/all.plot", run.number, ".png"))
# plot(fit, pars = c("beta"))
# dev.off()
#
print(fit, pars = c("B"))
print(fit, pars = c("A"))
print(fit, pars = c("G"))
print(fit, pars = c("theta"))
print(fit, pars = c("D"))

#
# png(width=1000, height=800, file=paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/pairs.plot", run.number, ".png"))
pairs(fit, pars = c("lp__", "B"))
pairs(fit, pars = c("lp__", "G"))
pairs(fit, pars = c("lp__", "D"))
# dev.off()


############################
###NEW extract fitted values
############################
HDIofMCMC = function( sampleVec , credMass ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

B.extract <- extract(fit, pars = "B", inc_warmup=FALSE)
A.extract <- extract(fit, pars = "A", inc_warmup = FALSE)
G.extract <- extract(fit, pars = "G", inc_warmup = FALSE)
theta.extract <- extract(fit, pars = "theta", inc_warmup = FALSE)
D.extract <- extract(fit, pars = "D", inc_warmup = FALSE)
y_pred.extract <- extract(fit, pars = "y_pred", inc_warmup = FALSE)

###Betas
#Estimates
B.est <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  c(mean(B.extract$B[,1,i]), mean(B.extract$B[,2,i]), mean(B.extract$B[,3,i]), mean(B.extract$B[,4,i]), mean(B.extract$B[,5,i]))
}
#Lower CI
B.lower <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  c(quantile(B.extract$B[,1,i], 0.05), quantile(B.extract$B[,2,i], 0.05), quantile(B.extract$B[,3,i], 0.05), quantile(B.extract$B[,4,i], 0.05), quantile(B.extract$B[,5,i], 0.05))
}
#Upper CI
B.upper <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  c(quantile(B.extract$B[,1,i], 0.95), quantile(B.extract$B[,2,i], 0.95), quantile(B.extract$B[,3,i], 0.95), quantile(B.extract$B[,4,i], 0.95), quantile(B.extract$B[,5,i], 0.95))
}

###Alphas
#Mean, lower and upper in one DF
A <- foreach (i = 1:J, .combine = 'rbind') %do% {
  c(quantile(A.extract$A[,i], 0.05), mean(A.extract$A[,i]), quantile(A.extract$A[,i], 0.95))
}
rownames(A) <- unique(df.not.mis$tree) #This may not give the right tree names

###Gammmas
#Mean, lower and upper in one DF
G <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  c(quantile(G.extract$G[,i], 0.05), mean(G.extract$G[,i]), quantile(G.extract$G[,i], 0.95))
}
rownames(G) <- c("Jacaranda copaia", "Cecropia insignis", "Triplaris cumingiana")

###Thetas
#Mean, lower and upper in one DF
theta <- foreach (i = 1:Q, .combine = 'rbind') %do% {
  c(quantile(theta.extract$theta[,i], 0.05), mean(theta.extract$theta[,i]), quantile(theta.extract$theta[,i], 0.95))
}
rownames(theta) <- unique(q) #This may not give the right quadrat names

###Deltas
#Mean, lower and upper in one DF
D <- foreach (i = 1:M, .combine = 'rbind') %do% {
  c(quantile(D.extract$D[,i], 0.05), mean(D.extract$D[,i]), quantile(D.extract$D[,i], 0.95))
}
rownames(D) <- c("PC1", "PC2")

#Put betas and gammas together
results <- cbind(B.lower[,1], B.est[,1], B.upper[,1], B.lower[,2], B.est[,2], B.upper[,2], B.lower[,3], B.est[,3], B.upper[,3], B.lower[,4], B.est[,4], B.upper[,4], B.lower[,5], B.est[,5], B.upper[,5], G[,1], G[,2], G[,3])
colnames(results) <- c("DBH_lower", "DBH", "DBH_upper", "NC.het_lower", "NC.het", "NC.het_upper", "NC.con_lower", "NC.con", "NC.con_upper", "NRI_lower", "NRI", "NRI_upper", "NPI_lower", "NPI", "NPI_upper", "Inbr_lower", "Inbr", "Inbr_upper")
rownames(results) <- c("Jacaranda copaia", "Cecropia insignis", "Triplaris cumingiana")
# write.csv(results, file="/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/na.sp/growth_1-6_NRI.tyler_NPI.tyler_na.sp_HOM.1.3.csv")

###################
#Get HDIs instead
###################
###Betas
B.hdi.lower <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  cbind(HDIofMCMC(B.extract$B[,1,i], credMass = 0.95)[1], HDIofMCMC(B.extract$B[,2,i], credMass = 0.95)[1], HDIofMCMC(B.extract$B[,3,i], credMass = 0.95)[1], HDIofMCMC(B.extract$B[,4,i], credMass = 0.95)[1], HDIofMCMC(B.extract$B[,5,i], credMass = 0.95)[1])
}
B.hdi.upper <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  cbind(HDIofMCMC(B.extract$B[,1,i], credMass = 0.95)[2], HDIofMCMC(B.extract$B[,2,i], credMass = 0.95)[2], HDIofMCMC(B.extract$B[,3,i], credMass = 0.95)[2], HDIofMCMC(B.extract$B[,4,i], credMass = 0.95)[2], HDIofMCMC(B.extract$B[,5,i], credMass = 0.95)[2])
}

###Alphas
#Mean, lower and upper in one DF
A.hdi <- foreach (i = 1:J, .combine = 'rbind') %do% {
  c(HDIofMCMC(A.extract$A[,i], credMass = 0.95)[1], mean(A.extract$A[,i]), HDIofMCMC(A.extract$A[,i], credMass = 0.95)[2])
}
rownames(A) <- unique(df.not.mis$tree) #This may not give the right tree names

###Gammas
#Mean, lower and upper in one DF
G.hdi <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  c(HDIofMCMC(G.extract$G[,i], credMass = 0.95)[1], mean(G.extract$G[,i]), HDIofMCMC(G.extract$G[,i], credMass = 0.95)[2])
}
rownames(G) <- c("Jacaranda copaia", "Cecropia insignis", "Triplaris cumingiana")

###Thetas
#Mean, lower and upper in one DF
theta.hdi <- foreach (i = 1:Q, .combine = 'rbind') %do% {
  c(HDIofMCMC(theta.extract$theta[,i], credMass = 0.95)[1], mean(theta.extract$theta[,i]), HDIofMCMC(theta.extract$theta[,i], credMass = 0.95)[2])
}
rownames(theta.hdi) <- unique(q) #This may not give the right quadrat names

###Deltas
#Mean, lower and upper in one DF
D.hdi <- foreach (i = 1:M, .combine = 'rbind') %do% {
  c(HDIofMCMC(D.extract$D[,i], credMass = 0.95)[1], mean(D.extract$D[,i]), HDIofMCMC(D.extract$D[,i], credMass = 0.95)[2])
}
rownames(D.hdi) <- c("PC1", "PC2")

##############
#Get 80% confidence intervals as well
##############

###Betas
B.hdi.lower.80 <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  cbind(HDIofMCMC(B.extract$B[,1,i], credMass = 0.8)[1], HDIofMCMC(B.extract$B[,2,i], credMass = 0.8)[1], HDIofMCMC(B.extract$B[,3,i], credMass = 0.8)[1], HDIofMCMC(B.extract$B[,4,i], credMass = 0.8)[1], HDIofMCMC(B.extract$B[,5,i], credMass = 0.8)[1])
}
B.hdi.upper.80 <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  cbind(HDIofMCMC(B.extract$B[,1,i], credMass = 0.8)[2], HDIofMCMC(B.extract$B[,2,i], credMass = 0.8)[2], HDIofMCMC(B.extract$B[,3,i], credMass = 0.8)[2], HDIofMCMC(B.extract$B[,4,i], credMass = 0.8)[2], HDIofMCMC(B.extract$B[,5,i], credMass = 0.8)[2])
}

###Alphas
#Mean, lower and upper in one DF
A.hdi.80 <- foreach (i = 1:J, .combine = 'rbind') %do% {
  c(HDIofMCMC(A.extract$A[,i], credMass = 0.8)[1], mean(A.extract$A[,i]), HDIofMCMC(A.extract$A[,i], credMass = 0.8)[2])
}
rownames(A.hdi.80) <- unique(df.not.mis$tree) #This may not give the right tree names

###Gammas
#Mean, lower and upper in one DF
G.hdi.80 <- foreach (i = 1:Nsp, .combine = 'rbind') %do% {
  c(HDIofMCMC(G.extract$G[,i], credMass = 0.8)[1], mean(G.extract$G[,i]), HDIofMCMC(G.extract$G[,i], credMass = 0.8)[2])
}
rownames(G.hdi.80) <- c("Jacaranda copaia","Cecropia insignis", "Triplaris cumingiana")

###Thetas
#Mean, lower and upper in one DF
theta.hdi.80 <- foreach (i = 1:Q, .combine = 'rbind') %do% {
  c(HDIofMCMC(theta.extract$theta[,i], credMass = 0.8)[1], mean(theta.extract$theta[,i]), HDIofMCMC(theta.extract$theta[,i], credMass = 0.8)[2])
}
rownames(theta.hdi.80) <- unique(q) #This may not give the right quadrat names

###Deltas
#Mean, lower and upper in one DF
D.hdi.80 <- foreach (i = 1:M, .combine = 'rbind') %do% {
  c(HDIofMCMC(D.extract$D[,i], credMass = 0.8)[1], mean(D.extract$D[,i]), HDIofMCMC(D.extract$D[,i], credMass = 0.8)[2])
}
rownames(D.hdi.80) <- c("PC1", "PC2")

hdi.results <- cbind(B.hdi.lower[,1], B.hdi.lower.80[,1], B.est[,1], B.hdi.upper.80[,1], B.hdi.upper[,1], B.hdi.lower[,2], B.hdi.lower.80[,2], B.est[,2], B.hdi.upper.80[,2], B.hdi.upper[,2], B.hdi.lower[,3], B.hdi.lower.80[,3], B.est[,3], B.hdi.upper.80[,3], B.hdi.upper[,3], B.hdi.lower[,4], B.hdi.lower.80[,4], B.est[,4], B.hdi.upper.80[,4], B.hdi.upper[,4], B.hdi.lower[,5], B.hdi.lower.80[,5], B.est[,5], B.hdi.upper.80[,5], B.hdi.upper[,5], G.hdi[,1], G.hdi.80[,1], G.hdi[,2], G.hdi.80[,3], G.hdi[,3])
colnames(hdi.results) <- c("DBH_lower", "DBH_lower.80", "DBH", "DBH_upper.80", "DBH_upper", "NC.het_lower", "NC.het_lower.80", "NC.het", "NC.het_upper.80", "NC.het_upper", "NC.con_lower", "NC.con_lower.80", "NC.con", "NC.con_upper.80", "NC.con_upper", "NRI_lower", "NRI_lower.80", "NRI", "NRI_upper.80", "NRI_upper", "NPI_lower", "NPI_upper.80", "NPI", "NPI_upper.80", "NPI_upper", "Inbr_lower", "Inbr_lower.80", "Inbr", "Inbr_upper.80", "Inbr_upper")
rownames(hdi.results) <- c("Jacaranda copaia", "Cecropia insignis", "Triplaris cumingiana")
# write.csv(hdi.results, file="/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/growth_1-6_NRI_4000iter_NRI.tyler_NPI.tyler_hdi.csv")

#Save all estimates and HDIs into an Rdata file
# save(B.est, B.hdi.lower, B.hdi.upper, A, G, theta, D, file="/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/growth_1-6_NRI_4000iter_NRI.tyler_NPI.tyler_estimates.RData")


##############################################
#Transform parameters back to original scale
##############################################
#Use this site for help:  https://stats.stackexchange.com/questions/74622/converting-standardized-betas-back-to-original-variables

###Betas
#Estimates
B.est.orig <- matrix(nrow=2,ncol=5)
B.est.orig[,1] <- B.est[,1]*1/log.dbhSD
B.est.orig[,2] <- B.est[,2]*1/log.nc.hetSD
B.est.orig[,3] <- B.est[,3]*1/log.nc.conSD
B.est.orig[,4] <- B.est[,4]*1/log.nriSD
B.est.orig[,5] <- B.est[,5]*1/log.npiSD
#Lower CIs
B.hdi.lower.orig <- matrix(nrow=2,ncol=5)
B.hdi.lower.orig[,1] <- B.hdi.lower[,1]*1/log.dbhSD
B.hdi.lower.orig[,2] <- B.hdi.lower[,2]*1/log.nc.hetSD
B.hdi.lower.orig[,3] <- B.hdi.lower[,3]*1/log.nc.conSD
B.hdi.lower.orig[,4] <- B.hdi.lower[,4]*1/log.nriSD
B.hdi.lower.orig[,5] <- B.hdi.lower[,5]*1/log.npiSD
#Upper CIs
B.hdi.upper.orig <- matrix(nrow=2,ncol=5)
B.hdi.upper.orig[,1] <- B.hdi.upper[,1]*1/log.dbhSD
B.hdi.upper.orig[,2] <- B.hdi.upper[,2]*1/log.nc.hetSD
B.hdi.upper.orig[,3] <- B.hdi.upper[,3]*1/log.nc.conSD
B.hdi.upper.orig[,4] <- B.hdi.upper[,4]*1/log.nriSD
B.hdi.upper.orig[,5] <- B.hdi.upper[,5]*1/log.npiSD

###Gammas
G.est.orig <- G[,2] * 1 / inbr.transSD
G.hdi.lower.orig <- G.hdi[,1] / inbr.transSD
G.hdi.upper.orig <- G.hdi[,3] / inbr.transSD

hdi.results.orig <- cbind(B.hdi.lower.orig[,1], B.est.orig[,1], B.hdi.upper.orig[,1], B.hdi.lower.orig[,2], B.est.orig[,2], B.hdi.upper.orig[,2], B.hdi.lower.orig[,3], B.est.orig[,3], B.hdi.upper.orig[,3], B.hdi.lower.orig[,4], B.est.orig[,4], B.hdi.upper.orig[,4], B.hdi.lower.orig[,5], B.est.orig[,5], B.hdi.upper.orig[,5], G.hdi.lower.orig, G.est.orig, G.hdi.upper.orig)
colnames(hdi.results.orig) <- c("DBH_lower", "DBH", "DBH_upper", "NC.het_lower", "NC.het", "NC.het_upper", "NC.con_lower", "NC.con", "NC.con_upper", "NRI_lower", "NRI", "NRI_upper", "NPI_lower", "NPI", "NPI_upper", "Inbr_lower", "Inbr", "Inbr_upper")
rownames(hdi.results.orig) <- c("Jacaranda copaia", "Luehea seemannii", "Cecropia insignis", "Tetragastris panamensis", "Triplaris cumingiana", "Virola sebifera")
# write.csv(hdi.results.orig, file="/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/growth_1-6_NRI_4000iter_NRI.tyler_NPI.tyler_hdi_orig.csv")

###Make figure of parameter estimates
#Massage hdi.results into correct format
# hdi.results.plot <- data.frame(c(hdi.results[1,1], hdi.results[1,4], hdi.results[1,7], hdi.results[1,10], hdi.results[1,13], hdi.results[1,16], hdi.results[2,1], hdi.results[2,4], hdi.results[2,7], hdi.results[2,10], hdi.results[2,13], hdi.results[2,16], hdi.results[3,1], hdi.results[3,4], hdi.results[3,7], hdi.results[3,10], hdi.results[3,13], hdi.results[3,16], hdi.results[4,1], hdi.results[4,4], hdi.results[4,7], hdi.results[4,10], hdi.results[4,13], hdi.results[4,16], hdi.results[5,1], hdi.results[5,4], hdi.results[5,7], hdi.results[5,10], hdi.results[5,13], hdi.results[5,16], hdi.results[6,1], hdi.results[6,4], hdi.results[6,7], hdi.results[6,10], hdi.results[6,13], hdi.results[6,16]))
# colnames(hdi.results.plot) <- "Lower"
# hdi.results.plot$Estimate <- c(hdi.results[1,2], hdi.results[1,5], hdi.results[1,8], hdi.results[1,11], hdi.results[1,14], hdi.results[1,17], hdi.results[2,2], hdi.results[2,5], hdi.results[2,8], hdi.results[2,11], hdi.results[2,14], hdi.results[2,17], hdi.results[3,2], hdi.results[3,5], hdi.results[3,8], hdi.results[3,11], hdi.results[3,14], hdi.results[3,17], hdi.results[4,2], hdi.results[4,5], hdi.results[4,8], hdi.results[4,11], hdi.results[4,14], hdi.results[4,17], hdi.results[5,2], hdi.results[5,5], hdi.results[5,8], hdi.results[5,11], hdi.results[5,14], hdi.results[5,17], hdi.results[6,2], hdi.results[6,5], hdi.results[6,8], hdi.results[6,11], hdi.results[6,14], hdi.results[6,17])
# hdi.results.plot$Upper <- c(hdi.results[1,3], hdi.results[1,6], hdi.results[1,9], hdi.results[1,12], hdi.results[1,15], hdi.results[1,18], hdi.results[2,3], hdi.results[2,6], hdi.results[2,9], hdi.results[2,12], hdi.results[2,15], hdi.results[2,18], hdi.results[3,3], hdi.results[3,6], hdi.results[3,9], hdi.results[3,12], hdi.results[3,15], hdi.results[3,18], hdi.results[4,3], hdi.results[4,6], hdi.results[4,9], hdi.results[4,12], hdi.results[4,15], hdi.results[4,18], hdi.results[5,3], hdi.results[5,6], hdi.results[5,9], hdi.results[5,12], hdi.results[5,15], hdi.results[5,18], hdi.results[6,3], hdi.results[6,6], hdi.results[6,9], hdi.results[6,12], hdi.results[6,15], hdi.results[6,18])
# hdi.results.plot$species <- c(rep("Jacaranda copaia", 6), rep("Luehea seemannii", 6), rep("Cecropia insignis", 6), rep("Tetragastris panamensis", 6), rep("Triplaris cumingiana", 6), rep("Virola sebifera", 6))
# hdi.results.plot$parameter <- rep(c("DBH", "NC.het", "NC.con", "NRI", "NPI", "Inbreeding"), 6)
# hdi.results.plot$sig <- ifelse((hdi.results.plot$Lower < 0 & hdi.results.plot$Upper < 0) | (hdi.results.plot$Lower > 0 & hdi.results.plot$Upper > 0), "yes", "no")

###Try melting and re-casting the hdi.results dataframe instead of this convoluted crap above. I think something is going wrong with it.
library(reshape2)
hdi.results.plot <- melt(hdi.results)
hdi.results.plot$type <- rep(c(rep("lower", 3), rep("lower.80", 3), rep("estimate", 3), rep("upper.80", 3), rep("upper", 3)), 6)
hdi.results.plot$variable <- c(rep("DBH", 15), rep("NC.het", 15), rep("NC.con", 15), rep("NRI", 15), rep("NPI", 15), rep("Inbreeding", 15))
hdi.results.plot.recast <- dcast(hdi.results.plot, variable+Var1 ~ type)
colnames(hdi.results.plot.recast) <- c("parameter", "species", "estimate", "lower.95", "lower.80", "upper.95", "upper.80")
hdi.results.plot.recast$sig <- ifelse((hdi.results.plot.recast$lower.80 < 0 & hdi.results.plot.recast$upper.80 < 0) | (hdi.results.plot.recast$lower.80 > 0 & hdi.results.plot.recast$upper.80 > 0), "yes", "no")
hdi.results.plot.recast$sig <- as.factor(hdi.results.plot.recast$sig)

point.dodge <- position_dodge(width=0.5)
errorbar.dodge <- position_dodge(width=0.5)
(growth.est.plot <- ggplot(y=estimate, x=parameter, data=hdi.results.plot.recast)+
    geom_point(aes(y=estimate, x=parameter, color=species, shape = sig), size=5, position = point.dodge)+
    geom_errorbar(aes(x=parameter, ymin = lower.80, ymax = upper.80, color=species, linetype = sig), width=0.4, lwd=0.9, position=errorbar.dodge)+
    scale_shape_manual(values = c("yes" = 19, "no" = 1), guide=FALSE)+
    scale_linetype_manual(values = c("yes" = "solid", "no" = "twodash"), guide=FALSE)+
    scale_color_discrete(name = "Species", labels = c("Jacaranda copaia", "Cecropia insignis", "Triplaris cumingiana"))+
    # scale_x_discrete(labels = c("DBH","Inbreeding","Conspecific NC", "Heterospecific NC", "Neighborhood\nphylogenetic index", "Neighborhood\nrelatedness index"))+
    geom_hline(yintercept=0, size=0.9)+
    xlab("Parameter")+
    ylab("Estimate")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 16, angle=90, hjust=1, vjust=0.55, color = "black"),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12)
    ))
# ggsave(growth.est.plot, width=7.5, height=6, file="/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/growth.est.plot_NRI.tyler_NPI.tyler.pdf")

#Try to make a sideways plot
source("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/horizontal_dodge.R")
source("/Volumes/tschappe/Documents/Projects/Wind\ River/Analysis/horizontal_dodge.R")
(growth.est.plot.sideways <- ggplot(y=parameter, x=estimate, data=hdi.results.plot.recast)+
    geom_point(aes(y=parameter, x=estimate, shape=species, color=species, alpha=sig), size=5, position = position_dodgev(height=0.7))+
    geom_errorbarh(aes(x=estimate, y=parameter, xmin = lower.80, xmax = upper.80, color=species, alpha=sig, linetype=sig), lwd=1.5, position = position_dodgev(height=0.7), height=0.7, data=hdi.results.plot.recast)+
    geom_errorbarh(aes(x=estimate, y=parameter, xmin = lower.95, xmax = upper.95, color=species, alpha=sig, linetype=sig), lwd=0.5, position = position_dodgev(height=0.7), height=0.7, data=hdi.results.plot.recast)+
    scale_shape_manual(name = "Species", values = c("Jacaranda copaia" = 16, "Cecropia insignis" = 17, "Triplaris cumingiana" = 18))+
    scale_alpha_manual(name = "Species", values = c("yes" = 1, "no" = 0.2), guide=FALSE)+
    scale_linetype_manual(values = c("yes" = "solid", "no" = "twodash"), guide=FALSE)+
    scale_color_discrete(name = "Species")+
    scale_y_discrete(labels = c("DBH","Inbreeding", "Conspecific NC", "Heterospecific NC", "NPI", "NRI"))+
    scale_x_continuous(breaks=seq(from=-20, to=20, by=0.2))+
    geom_vline(xintercept=0, size=0.9)+
    xlab("Estimate")+
    ylab("Parameter")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 16, hjust=1, vjust=0.55, color = "black"),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size=20),
          axis.title.y = element_blank(),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12)
    ))
# ggsave(growth.est.plot.sideways, width=10, height=8, file="/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/na.sp/growth.est.plot.npi_NRI.tyler_NPI.tyler_na.sp_HOM1.3_new.nri.npi_0.5_impute_sideways.pdf")



###################
###Model validation
###################
###Check out predicted vs true growth values
df.not.mis$growth.pred.log <- foreach(z = 1:nrow(X), .combine = 'c') %do% {
  mean(y_pred.extract$y_pred[,z])
}
df.not.mis$growth.pred <- exp(df.not.mis$growth.pred.log)

###Check residuals on unlogged scale
df.not.mis$resid <- df.not.mis$growth - df.not.mis$growth.pred
#Check normality of residuals
qqnorm(df.not.mis$resid); qqline(df.not.mis$resid)
#Check residuals vs. fitted
(pred.resid <- ggplot(aes(x=growth.pred, y=resid), data=df.not.mis)+
  geom_point(aes(x=growth.pred, y=resid), color="#505050")+
  geom_hline(yintercept=0, size=1.5, color="black")+
  xlab("Predicted growth (mm)")+
  ylab("Residual value")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ))
# ggsave(pred.resid, width=8, height=6, file="/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/pred.resid.nri.pdf")
#Check residuals vs. true values
qplot(df.not.mis$growth, df.not.mis$resid)
##Interesting!!! There is a bit of a funnel showing that the variance increases as the fitted value increases -- may need to allow for different variances among individuals? Or quadrats?

###Check residuals on logged scale
df.not.mis$resid.log <- df.not.mis$log.growth - df.not.mis$growth.pred.log
#Check normality of residuals
qqnorm(df.not.mis$resid.log); qqline(df.not.mis$resid.log)
#Check residuals vs. fitted
qplot(df.not.mis$growth.pred.log, df.not.mis$resid.log)
#Check residuals vs. true values
qplot(df.not.mis$log.growth, df.not.mis$resid.log)

#Plot predicted vs. true growth
qplot(log.growth, growth.pred.log, data=df.not.mis)+geom_smooth(method="lm")
log.growth.lm <- summary(lm(growth.pred.log ~ log.growth, data=df.not.mis))
growth.lm <- summary(lm(growth.pred ~ growth, data=df.not.mis))

growth.r.squared <- growth.lm$adj.r.squared
log.growth.r.squared <- log.growth.lm$adj.r.squared

(pred.true.growth <- ggplot(aes(x=growth, y=growth.pred), data=df.not.mis)+
  geom_point(aes(x=growth, y=growth.pred), color="#707070")+
  geom_smooth(method="lm", fill="#363636", color="black", size=1.2)+
  xlab("Measured growth (mm)")+
  ylab("Estimated growth (mm)")+
  annotate("text", x=25, y=112, label = paste0("R-squared = ", round(growth.r.squared, digits = 2)))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
)
# ggsave(pred.true.growth, width=8, height=6, file="/Users/tschappe/Documents/Projects/Wind\ River/Analysis/new.targets/4.census/real.stan.4census/multilevel/growth/pred.true.growth.plot_NRI.tyler_NPI.tyler.pdf")


###Look at log likelihood to calculate model fit diagnostics
library(loo)
log_lik_1 <- extract_log_lik(fit)
#Look at leave-one-out validation
loo_1 <- loo(log_lik_1)
print(loo_1)
plot(loo_1)
#Look at WAIC
waic_1 <- waic(log_lik_1)
print(waic_1)


#
# ###########Get predicted values
# new_y <- extract(fit,pars="surv_pred")
# df.not.mis$pred <- new_y$surv_pred[1,]
#
# #Grab the quadrat random intercepts
# for (i in 1:length(unique(df.not.mis$tree.ind))) {
#   df.not.mis$quadrat.int[i] <- Quad.est[df.not.mis$quadrat.sp[i]]
# }
# df.not.mis$quadrat.int.exp <- exp(df.not.mis$quadrat.int)



# ###Plot just the non-missing data
# (pred.dbh.plot <- ggplot(aes(x=dbh.std, y=pred), data=df.not.mis)+
#   geom_point(aes(x=dbh.std, y=pred, color="blue")))+
#   geom_point(aes(x=dbh.std, y=status, color="red"))
#
# #ggsave(pred.dbh.plot, width=8, height=6, file = paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_not.missing_growth.no.std_not.centered5_STD_new.nc_all.beta_figs/pred.dbh", run.number, ".pdf"))
# ggsave(pred.dbh.plot, width=8, height=6, file = paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/pred.dbh", run.number, ".pdf"))
#
# pred.growth.lm.sum <- summary(lm(exp(pred) ~ growth, data=df.not.mis))
# r.squared <- pred.growth.lm.sum$adj.r.squared
# (pred.true.growth.plot <- ggplot(data=df.not.mis, aes(x=exp(pred), y=growth))+
#   geom_point(aes(x=exp(pred), y=growth))+
#   geom_abline(intercept = 0, slope=1, color="red")+
#   #geom_smooth(method="lm")+
#   annotate("text", x=100, y=50, label = paste0("r-squared = ", round(r.squared, 3))))
# #ggsave(pred.true.growth.plot, width=8, height=6, file = paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_not.missing_growth.no.std_not.centered5_STD_new.nc_all.beta_figs/pred.true.growth", run.number, ".pdf"))
# ggsave(pred.true.growth.plot, width=8, height=6, file = paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/pred.true.growth", run.number, ".pdf"))
#
# #On log scale
# pred.growth.lm.sum.log <- summary(lm(pred ~ log(growth), data=df.not.mis))
# r.squared.log <- pred.growth.lm.sum.log$adj.r.squared
# (pred.true.growth.plot.log <- ggplot(data=df.not.mis, aes(x=pred, y=log(growth)))+
#     geom_point(aes(x=pred, y=log(growth)))+
#     geom_abline(intercept = 0, slope=1, color="red")+
#     #geom_smooth(method="lm")+
#     annotate("text", x=3, y=5, label = paste0("r-squared = ", round(r.squared.log, 3))))
# #ggsave(pred.true.growth.plot.log, width=8, height=6, file = paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_not.missing_growth.no.std_not.centered5_STD_new.nc_all.beta_figs/pred.true.growth.log", run.number, ".pdf"))
# ggsave(pred.true.growth.plot.log, width=8, height=6, file = paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/pred.true.growth.log", run.number, ".pdf"))



# ############Transform parameters back to original scale
# #Note that since we didn't transform the response variable, we need to back-transform different from the linear regression
# beta1.est.orig <- beta1.est/log.dbhSD
# beta2.est.orig <- beta2.est/pc1SD
# beta3.est.orig <- beta3.est/pc2SD
# beta4.est.orig <- beta4.est/inbr.transSD
# beta5.est.orig <- beta5.est/log.nc.conSD
# beta6.est.orig <- beta6.est/log.nc.hetSD
# beta0.est.orig <- c(0,0)
# beta0.est.orig[1] <- (beta0.est[1]) - sum((beta1.est[1]*log.dbhM/log.dbhSD) + (beta2.est[1]*pc1M/pc1SD) + (beta3.est[1]*pc2M/pc2SD) + (beta4.est[1]*inbr.transM/inbr.transSD) + (beta5.est[1]*log.nc.conM/log.nc.conSD) + (beta6.est[1]*log.nc.hetM/log.nc.hetSD))
# beta0.est.orig[2] <- (beta0.est[2]) - sum((beta1.est[2]*log.dbhM/log.dbhSD) + (beta2.est[2]*pc1M/pc1SD) + (beta3.est[2]*pc2M/pc2SD) + (beta4.est[2]*inbr.transM/inbr.transSD) + (beta5.est[2]*log.nc.conM/log.nc.conSD) + (beta6.est[2]*log.nc.hetM/log.nc.hetSD))
# Census.est.orig <- Census.est
# 
# #Back transform lower CIs
# beta1.lower.orig <- beta1.lower/log.dbhSD
# beta2.lower.orig <- beta2.lower/pc1SD
# beta3.lower.orig <- beta3.lower/pc2SD
# beta4.lower.orig <- beta4.lower/inbr.transSD
# beta5.lower.orig <- beta5.lower/log.nc.conSD
# beta6.lower.orig <- beta6.lower/log.nc.hetSD
# beta0.lower.orig <- c(0,0)
# beta0.lower.orig[1] <- (beta0.lower[1]) - sum((beta1.lower[1]*log.dbhM/log.dbhSD) + (beta2.lower[1]*pc1M/pc1SD) + (beta3.lower[1]*pc2M/pc2SD) + (beta4.lower[1]*inbr.transM/inbr.transSD) + (beta5.lower[1]*log.nc.conM/log.nc.conSD) + (beta6.lower[1]*log.nc.hetM/log.nc.hetSD))
# beta0.lower.orig[2] <- (beta0.lower[2]) - sum((beta1.lower[2]*log.dbhM/log.dbhSD) + (beta2.lower[2]*pc1M/pc1SD) + (beta3.lower[2]*pc2M/pc2SD) + (beta4.lower[2]*inbr.transM/inbr.transSD) + (beta5.lower[2]*log.nc.conM/log.nc.conSD) + (beta6.lower[2]*log.nc.hetM/log.nc.hetSD))
# 
# #Back transform upper CIs
# beta1.upper.orig <- beta1.upper/log.dbhSD
# beta2.upper.orig <- beta2.upper/pc1SD
# beta3.upper.orig <- beta3.upper/pc2SD
# beta4.upper.orig <- beta4.upper/inbr.transSD
# beta5.upper.orig <- beta5.upper/log.nc.conSD
# beta6.upper.orig <- beta6.upper/log.nc.hetSD
# beta0.upper.orig <- c(0,0)
# beta0.upper.orig[1] <- (beta0.upper[1]) - sum((beta1.upper[1]*log.dbhM/log.dbhSD) + (beta2.upper[1]*pc1M/pc1SD) + (beta3.upper[1]*pc2M/pc2SD) + (beta4.upper[1]*inbr.transM/inbr.transSD) + (beta5.upper[1]*log.nc.conM/log.nc.conSD) + (beta6.upper[1]*log.nc.hetM/log.nc.hetSD))
# beta0.upper.orig[2] <- (beta0.upper[2]) - sum((beta1.upper[2]*log.dbhM/log.dbhSD) + (beta2.upper[2]*pc1M/pc1SD) + (beta3.upper[2]*pc2M/pc2SD) + (beta4.upper[2]*inbr.transM/inbr.transSD) + (beta5.upper[2]*log.nc.conM/log.nc.conSD) + (beta6.upper[2]*log.nc.hetM/log.nc.hetSD))
# 
# #Save into a dataframe
# fitted.orig <- data.frame(cbind(c(beta0.lower.orig, beta1.lower.orig, beta2.lower.orig, beta3.lower.orig, beta4.lower.orig, beta5.lower.orig, beta6.lower.orig), c(beta0.est.orig, beta1.est.orig, beta2.est.orig, beta3.est.orig, beta4.est.orig, beta5.est.orig, beta6.est.orig), c(beta0.upper.orig, beta1.upper.orig, beta2.upper.orig, beta3.upper.orig, beta4.upper.orig, beta5.upper.orig, beta6.upper.orig)))
# colnames(fitted.orig) <- c("Lower", "Estimate", "Upper")
# rownames(fitted.orig) <- c("beta0.jac", "beta0.luh", "beta1.jac", "beta1.luh", "beta2.jac", "beta2.luh", "beta3.jac", "beta3.luh", "beta4.jac", "beta4.luh", "beta5.jac", "beta5.luh", "beta6.jac", "beta6.luh")
# 
# ###Make predicted values in original scale for non-missing
# df.not.mis$odds.orig <- foreach(i = 1:nrow(df.not.mis), .combine='c') %do% {
#   beta0.est.orig[df.not.mis$species.num[i]] + (beta1.est.orig[df.not.mis$species.num[i]] * df.not.mis$log.dbh[i]) + (beta2.est.orig[df.not.mis$species.num[i]] * df.not.mis$pc1[i]) + (beta3.est.orig[df.not.mis$species.num[i]] * df.not.mis$pc2[i]) + (beta4.est.orig[df.not.mis$species.num[i]] * df.not.mis$inbr4.logit[i]) + (beta5.est.orig[df.not.mis$species.num[i]] * df.not.mis$log.nc.con[i]) + (beta6.est.orig[df.not.mis$species.num[i]] * df.not.mis$log.nc.het[i]) + Census.est.orig[df.not.mis$census.num[i]]
# }
# df.not.mis$prob.orig <- 1 / (1 + exp(-(df.not.mis$odds.orig)))
# #Make predicted binomial response
# df.not.mis$surv.pred <- rbinom(1, 1, df.not.mis$prob.orig)
# 
# #Plot predicted values over top true growth values
# (pred.dbh.plot.orig <- ggplot(aes(x=dbh, y=prob.orig), data=df.not.mis)+
#     geom_point(aes(x=dbh, y=status, color=species))+
#     geom_line(aes(x=dbh, y=prob.orig, color=species)))
# #ggsave(pred.dbh.plot.orig, width=8, height=6, file = paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_missing_growth.no.std_not.centered4_not.chen_figs/pred.dbh.orig", run.number, ".pdf"))
# #ggsave(pred.dbh.plot.orig, width=8, height=6, file = paste0("/raid1/home/bpp/schappet/wind.river/bayesian.neighborhood.analysis/jac.luh/real/stan_growth/standardized_missing_growth.no.std_not.centered4_not.chen_figs/pred.dbh.orig", run.number, ".pdf"))
# 
# #Regress predicted values in original scale against true growth
# pred.growth.lm.sum.orig <- summary(lm(exp(pred.orig.unlog) ~ growth, data=df))
# r.squared.orig <- pred.growth.lm.sum.orig$adj.r.squared
# (pred.true.growth.plot.orig <- ggplot(data=df, aes(x=exp(pred.orig.unlog), y=growth))+
#     geom_point(aes(x=exp(pred.orig), y=growth))+
#     geom_smooth(method="lm")+
#     annotate("text", x=25, y=100, label = paste0("r-squared = ", round(r.squared.orig, 3))))
# #ggsave(pred.true.growth.plot.orig, width=8, height=6, file = paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_missing_growth.no.std_not.centered4_not.chen_figs/pred.true.growth.orig", run.number, ".pdf"))
# #ggsave(pred.true.growth.plot.orig, width=8, height=6, file = paste0("/raid1/home/bpp/schappet/wind.river/bayesian.neighborhood.analysis/jac.luh/real/stan_growth/standardized_missing_growth.no.std_not.centered4_not.chen_figs/pred.true.growth.orig", run.number, ".pdf"))
# 
# 
# 
# #Look at unlogged residuals
# df.not.mis$resids <- df.not.mis$growth - exp(df.not.mis$pred)
# #Plot them
# #pdf(width = 10, height = 8, file=paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_not.missing_growth.no.std_not.centered5_STD_new.nc_all.beta_figs/resids.", run.number, ".pdf"))
# pdf(width = 10, height = 8, file=paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/resids.", run.number, ".pdf"))
# qplot(df.not.mis$resids, geom="histogram")
# dev.off()
# 
# #pdf(width = 10, height = 8, file=paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_not.missing_growth.no.std_not.centered5_STD_new.nc_all.beta_figs/resids.norm.", run.number, ".pdf"))
# pdf(width = 10, height = 8, file=paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/resids.norm.", run.number, ".pdf"))
# qqnorm(df.not.mis$resids); qqline(df.not.mis$resids)
# dev.off()
# 
# 
# #Look at logged residuals
# df.not.mis$log.resids <- log(df.not.mis$growth) - df.not.mis$pred
# #Plot them
# #pdf(width = 10, height = 8, file=paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_not.missing_growth.no.std_not.centered5_STD_new.nc_all.beta_figs/resids.log.", run.number, ".pdf"))
# pdf(width = 10, height = 8, file=paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/resids.log.", run.number, ".pdf"))
# qplot(df.not.mis$log.resids, geom="histogram")
# dev.off()
# 
# #pdf(width = 10, height = 8, file=paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_not.missing_growth.no.std_not.centered5_STD_new.nc_all.beta_figs/resids.log.norm", run.number, ".pdf"))
# pdf(width = 10, height = 8, file=paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/resids.log.norm", run.number, ".pdf"))
# qqnorm(df.not.mis$log.resids); qqline(df.not.mis$log.resids)
# dev.off()
# 
# #Save targets dataframes
# #save(df, fit, file = paste0("/Users/tschappe/Documents/Projects/Wind\ River/Analysis/bayesian.neighborhood.analysis/jac.luh/real/real_stan/standardized_not.missing_growth.no.std_not.centered5_STD_new.nc_all.beta_figs/output", run.number, ".RData"))
# save(df.not.mis, file = paste0("/raid1/home/bpp/schappet/wind.river/new.targets/4.census/real.stan.4census/4.census_surv_not.missing_het.con_figs/output", run.number, ".RData"))

