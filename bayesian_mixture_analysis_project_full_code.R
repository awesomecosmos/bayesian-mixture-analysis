# Full code for Bayesian Mixture Analysis Project

# importing packages
library(palmerpenguins)
library(ggplot2)
library(tidyverse)
library(MCMCglmm)
library(mixAK)
library(cluster.datasets)
library("mixtools")
library(MASS)

# cleaning data - removing NA values
penguins.clean <- penguins[(!is.na(penguins$bill_length_mm))&
                           (!is.na(penguins$bill_depth_mm))&
                           (!is.na(penguins$species))&
                           (!is.na(penguins$sex)),]

# exploratory plotting of dataset
ggplot(data=penguins.clean,aes(bill_depth_mm,bill_length_mm))+ 
  geom_point(aes(group=species,col=species))

# exploratory data summary
summary(penguins.clean)

# setting variables for future use
my_xlab <- "bill depth (mm)"
my_ylab <- "bill length (mm)"

# plotting
par(bty="n")
layout(matrix(c(0,1,1,1,1,0, 2,2,2,3,3,3), nrow=2, byrow=TRUE))
plot(penguins$bill_depth_mm,penguins$bill_length_mm, col="darkviolet", pch=16, 
     xlab=my_xlab, ylab=my_ylab, main="Bill Depth vs Bill Length")
hist(penguins$bill_depth_mm, prob=TRUE, col="mediumslateblue", 
     xlab=my_xlab, ylab="Density", main="Bill Depth",breaks=100)
hist(penguins$bill_length_mm, prob=TRUE, col="mediumslateblue", 
     xlab=my_ylab, ylab="Density", main="Bill Length",breaks=100)

########################################################
# INITIALISATION
########################################################

bill_depth <- c(penguins.clean$bill_depth_mm)
bill_length <- c(penguins.clean$bill_length_mm)
my_dataset <-  data.frame(bill_depth, bill_length)
my_kmax = 4

# setting initial MCMC values
nMCMC <- c(burn=5000, keep=20000, thin=10, info=10000)
ygrid <- list(bill_depth=seq(12, 22, length=100), 
              bill_length=seq(30, 60, length=100))

########################################################
# MIXTURE MODEL ANALYSIS
########################################################
RJPrior1 <- list(priorK="uniform", Kmax=my_kmax)
RJPrior2 <- list(priorK="uniform", Kmax=my_kmax,
                 delta=1)
parRJMCMC2 <- list(par.u1=c(2, 2),par.u2=c(2, 2),par.u3=c(1, 1))

# Model with three components
RJModel2 <- NMixMCMC(y0=my_dataset, prior=RJPrior2, RJMCMC=parRJMCMC2,
                          nMCMC=nMCMC, scale=list(shift=0, scale=1), PED=TRUE)

# obtaining summary info
summary(RJModel2[[1]]$mu)
summary(RJModel2[[1]]$Sigma)

# setting variables for future use
rnorm_depth_ch1 <- rnorm(10^4, mean=mean(RJModel2[[1]]$pm.y$y1), sd=sd(RJModel2[[1]]$pm.y$y1))
rnorm_depth_ch2 <- rnorm(10^4, mean=mean(RJModel2[[2]]$pm.y$y1), sd=sd(RJModel2[[2]]$pm.y$y1))
rnorm_length_ch1 <- rnorm(10^4, mean=mean(RJModel2[[1]]$pm.y$y2), sd=sd(RJModel2[[1]]$pm.y$y2))
rnorm_length_ch2 <- rnorm(10^4, mean=mean(RJModel2[[2]]$pm.y$y2), sd=sd(RJModel2[[2]]$pm.y$y2))

dat1 = data.frame(x=rnorm_depth_ch1, group="Chain 1")
dat2 = data.frame(x=rnorm_depth_ch2, group="Chain 2")
depth_dat = rbind(dat1, dat2)
dat3 = data.frame(x=rnorm_length_ch1, group="Chain 1")
dat4 = data.frame(x=rnorm_length_ch2, group="Chain 2")
length_dat = rbind(dat3, dat4)

# plotting
ggplot(depth_dat, aes(x, fill=group, colour=group)) + 
   geom_histogram(binwidth=0.1) +
   ggtitle("Bill Depth Normal Posterior Distribution") +
   xlab(my_xlab)

ggplot(depth_dat, aes(x, fill=group, colour=group)) + geom_density(alpha=.2) +
   ggtitle("Bill Depth Normal Posterior Distribution Density") +      
   xlab(my_xlab)

ggplot(length_dat, aes(x, fill=group, colour=group)) + 
   geom_histogram(binwidth=0.1) +
   ggtitle("Bill Length Normal Posterior Distribution") +
   xlab(my_ylab)

ggplot(length_dat, aes(x, fill=group, colour=group)) + geom_density(alpha=.2) +
   ggtitle("Bill Length Normal Posterior Distribution Density") +      
   xlab(my_ylab)

# Model0: Marginal predictive densities
PDensRJ2 <- list()
PDensRJ2[[1]] <- NMixPredDensMarg(RJModel2[[1]], grid = ygrid)
PDensRJ2[[2]] <- NMixPredDensMarg(RJModel2[[2]], grid = ygrid)

# Model0: Plot of the marginal predictive density
plot(PDensRJ2[[1]], K=3, auto.layout=TRUE,
     type="l", col="darkblue", lty=1, lwd=1)

print(RJModel2[[1]]$moves)
print(RJModel2[[2]]$moves)
print(RJModel2)

# plotting
par(mfrow=c(1, 1), bty="n")
plot(PDensRJ2[[1]], K=c(0:my_kmax), xlab=my_xlab,
     lty=c(1, rep(2, 4)), col=c("darkblue", rep("red", 4)),
     lwd=c(2, rep(1, 4)),
     main="Predictive Density",
     sub="[for random number of mixture components (blue) vs K=1:4 (red)]")

# Chain 1 and Chain 2
par(mfrow=c(1, 1))
hist(my_dataset$bill_depth, prob=TRUE, col="darkviolet",
     xlab=my_xlab, ylab="Density", main="Predictive Density for both chains vs data",
     breaks=100)
lines(PDensRJ2[[1]]$x$x1, PDensRJ2[[1]]$dens[[1]],
      col="deeppink", lwd=4)
lines(PDensRJ2[[2]]$x$x1, PDensRJ2[[2]]$dens[[1]],
      col="limegreen", lwd=2, type="l")
legend("topleft", legend = c("Chain 1", "Chain 2"),
       col = c("deeppink", "limegreen"), lty = 1:2, cex = 0.8)

# Model0: Joint predictive density
JPDensModel0 <- list()
JPDensModel0[[1]] <- NMixPredDensJoint2(RJModel2[[1]], grid=ygrid)  
JPDensModel0[[2]] <- NMixPredDensJoint2(RJModel2[[2]], grid=ygrid) 

# Model0: Plot of the joint predictive density
plot(JPDensModel0[[1]])
plot(JPDensModel0[[2]])

# Model0: Plot of the marginal predictive density with the histogram
par(mfrow=c(2, 2), bty="n")
ylimE <- c(0, max(c(PDensRJ2[[1]]$dens[[1]], PDensRJ2[[2]]$dens[[1]])))
ylimW <- c(0, max(c(PDensRJ2[[1]]$dens[[2]], PDensRJ2[[2]]$dens[[2]])))

par(mfrow=c(1, 2), bty="n")
hist(bill_depth, prob=TRUE, col="mediumslateblue", 
     xlab=my_xlab, ylab="Density", main=paste("Predictive Density (Bill Depth)"), 
     ylim=ylimE,breaks=50)
lines(PDensRJ2[[1]]$x$x1, PDensRJ2[[1]]$dens[[1]], col="deeppink", lwd=4)
lines(PDensRJ2[[2]]$x$x1, PDensRJ2[[2]]$dens[[1]], col="limegreen", lwd=2)
legend("topleft", legend = c("Chain 1", "Chain 2"),
       col = c("deeppink", "limegreen"), lty = 1:2, cex = 0.8)
hist(bill_length, prob=TRUE, col="mediumslateblue", 
     xlab=my_ylab, ylab="Density", main=paste("Predictive Density (Bill Length)"), ylim=ylimW,breaks=50)
lines(PDensRJ2[[1]]$x$x2, PDensRJ2[[1]]$dens[[2]], col="deeppink", lwd=4)
lines(PDensRJ2[[2]]$x$x2, PDensRJ2[[2]]$dens[[2]], col="limegreen", lwd=2)

# Model0: Contour plot of the joint predictive density with the scatterplot
par(mfrow=c(2, 1), bty="n")
for (CH in 1:2){
   plot(my_dataset, col="deeppink", xlab=my_xlab, ylab=my_ylab, 
        main=paste("Chain", CH))
   contour(JPDensModel0[[CH]]$x$bill_depth, JPDensModel0[[CH]]$x$bill_length, 
           JPDensModel0[[CH]]$dens[["1-2"]], 
           col="darkviolet", add=TRUE)}

# Model0: Image plot of the joint predictive density with the scatterplot
par(mfrow=c(2, 1), bty="n")
for (CH in 1:2){
   plot(my_dataset, col="darkblue", xlab=my_xlab, ylab=my_ylab,
        main=paste("Chain", CH))
   image(JPDensModel0[[CH]]$x$bill_depth, JPDensModel0[[CH]]$x$bill_length, 
         JPDensModel0[[CH]]$dens[["1-2"]], add=TRUE, 
         col=rev(heat_hcl(33, c=c(80, 30), l=c(30, 90), power=c(1/5, 1.3))))  
   points(my_dataset, col="darkblue")}

# setting variables for future use
adelie <- penguins.clean[(penguins.clean$species=='Adelie'),]
gentoo <- penguins.clean[(penguins.clean$species=='Gentoo'),]
chinstrap <- penguins.clean[(penguins.clean$species=='Chinstrap'),]

# Image plot of the joint predictive density with the scatterplot
par(mfrow=c(1, 1), bty="n")
plot(my_dataset,xlab=my_xlab, ylab=my_ylab,
     main="Joint Predictive Density Countour Map for Chain 1")
image(JPDensModel0[[1]]$x$bill_depth, JPDensModel0[[1]]$x$bill_length, 
         JPDensModel0[[1]]$dens[["1-2"]], add=TRUE, 
         col=rev(heat_hcl(33, c=c(80, 30), l=c(30, 90), power=c(1/5, 1.3)))) 
contour(JPDensModel0[[1]]$x$bill_depth, JPDensModel0[[1]]$x$bill_length, 
           JPDensModel0[[1]]$dens[["1-2"]], 
           col="darkviolet", add=TRUE)
points(adelie$bill_depth_mm, adelie$bill_length_mm, col="deeppink", pch = 20)
points(gentoo$bill_depth_mm, gentoo$bill_length_mm, col="dodgerblue", pch = 20)
points(chinstrap$bill_depth_mm, chinstrap$bill_length_mm, col="limegreen", pch = 20)
legend("topleft", legend = c("Adelie", "Gentoo", "Chinstrap"),
       col = c("deeppink", "dodgerblue", "limegreen"), lty = 1:2, cex = 0.8)

# setting variables for future use
biscoe <- penguins.clean[(penguins.clean$island=='Biscoe'),]
dream <- penguins.clean[(penguins.clean$island=='Dream'),]
torgersen <- penguins.clean[(penguins.clean$island=='Torgersen'),]

# Image plot of the joint predictive density with the scatterplot
par(mfrow=c(1, 1), bty="n")
plot(my_dataset,xlab=my_xlab, ylab=my_ylab,
     main="Joint Predictive Density Countour Map for Chain 1")
image(JPDensModel0[[1]]$x$bill_depth, JPDensModel0[[1]]$x$bill_length, 
      JPDensModel0[[1]]$dens[["1-2"]], add=TRUE, 
      col=rev(heat_hcl(33, c=c(80, 30), l=c(30, 90), power=c(1/5, 1.3)))) 
contour(JPDensModel0[[1]]$x$bill_depth, JPDensModel0[[1]]$x$bill_length, 
        JPDensModel0[[1]]$dens[["1-2"]], 
        col="darkviolet", add=TRUE)
points(biscoe$bill_depth_mm, biscoe$bill_length_mm, col="deeppink", pch = 20)
points(dream$bill_depth_mm, dream$bill_length_mm, col="dodgerblue", pch = 20)
points(torgersen$bill_depth_mm, torgersen$bill_length_mm, col="limegreen", pch = 20)
legend("topleft", legend = c("Biscoe", "Dream", "Torgersen"),
       col = c("deeppink", "dodgerblue", "limegreen"), lty = 1:2, cex = 0.8)

# setting variables for future use
male <- penguins.clean[(penguins.clean$sex=='male'),]
female <- penguins.clean[(penguins.clean$sex=='female'),]

# Image plot of the joint predictive density with the scatterplot
par(mfrow=c(1, 1), bty="n")
plot(my_dataset,xlab=my_xlab, ylab=my_ylab,
     main="Joint Predictive Density Countour Map for Chain 1")
image(JPDensModel0[[1]]$x$bill_depth, JPDensModel0[[1]]$x$bill_length, 
      JPDensModel0[[1]]$dens[["1-2"]], add=TRUE, 
      col=rev(heat_hcl(33, c=c(80, 30), l=c(30, 90), power=c(1/5, 1.3)))) 
contour(JPDensModel0[[1]]$x$bill_depth, JPDensModel0[[1]]$x$bill_length, 
        JPDensModel0[[1]]$dens[["1-2"]], 
        col="darkviolet", add=TRUE)
points(male$bill_depth_mm, male$bill_length_mm, col="royalblue", pch = 20)
points(female$bill_depth_mm, female$bill_length_mm, col="deeppink", pch = 20)
legend("topleft", legend = c("Male", "Female"),
       col = c("royalblue", "deeppink"), lty = 1:2, cex = 0.8)

# Bivariate (for each pair of margins) predictive densities
plot(JPDensModel0[[1]])

par(mfrow=c(1, 1), bty="n")
plot(my_dataset,xlab=my_xlab, ylab=my_ylab,
     main="Joint Predictive Density Countour Map for Chain 1")
image(JPDensModel0[[1]]$x$bill_depth, JPDensModel0[[1]]$x$bill_length, 
      JPDensModel0[[1]]$dens[["1-2"]], add=TRUE, 
      col=rev(heat_hcl(33, c=c(80, 30), l=c(30, 90), power=c(1/5, 1.3)))) 
contour(JPDensModel0[[1]]$x$bill_depth, JPDensModel0[[1]]$x$bill_length, 
        JPDensModel0[[1]]$dens[["1-2"]], 
        col="darkviolet", add=TRUE)

# Plot with data
par(mfrow=c(1, 1), bty="n")
i = 1
j = 2
ICOL <- rev(heat_hcl(20, c=c(80, 30), l=c(30, 90), power=c(1/5, 2)))
VARS <- names(my_dataset) #len=2
SPECIES <- c("Adelie","Gentoo","Chinstrap") #levels(iris[, "Species"])
COLS <- c("red", "darkblue", "darkgreen")
NAME <- paste(i, "-", j, sep="")
MAIN <- paste(VARS[i], "x", VARS[j])
image(JPDensModel0[[1]]$x[[i]], JPDensModel0[[1]]$x[[j]], JPDensModel0[[1]]$dens[[NAME]], col=ICOL,
      xlab=my_xlab, ylab=my_ylab, main=MAIN)
contour(JPDensModel0[[1]]$x$bill_depth, JPDensModel0[[1]]$x$bill_length, 
        JPDensModel0[[1]]$dens[["1-2"]], 
        col="darkviolet", add=TRUE)
counter <- 0
for (spec in SPECIES){
   counter <- counter + 1
   Data <- subset(penguins.clean, species==spec)
   #points(Data[,4], Data[,3], pch=16, col=COLS[spec])}
   x <- Data %>% pull(bill_depth_mm)
   y <- Data %>% pull(bill_length_mm)
   points(x, y, pch=16, col=COLS[counter])}
legend("topleft", legend = c("Adelie", "Gentoo", "Chinstrap"),
       col = COLS, lty = 1:2, cex = 0.8)

# Clustering based on posterior summary statistics of mixture weight, means, 
# variances
Prior <- list(priorK = "fixed", Kmax = 3)
fit <- NMixMCMC(y0=my_dataset, prior=Prior, nMCMC = nMCMC)
p3 <- NMixPlugDA(fit[[1]], my_dataset)
p3 <- p3[, paste("prob", 1:3, sep="")]
apply(p3[1:50,], 2, quantile)

# RJMCMC: Choose chain
CH <- 1

# Model0: Create mcmc objects
start <- RJModel2[[CH]]$nMCMC["burn"] + 1
end <- RJModel2[[CH]]$nMCMC["burn"] + RJModel2[[CH]]$nMCMC["keep"]
chK <- mcmc(RJModel2[[CH]]$K, start=start, end=end)
chgammaInv <- mcmc(RJModel2[[CH]]$gammaInv, start=start, end=end)
chmixture <- mcmc(RJModel2[[CH]]$mixture, start=start, end=end)
chdeviance <- mcmc(RJModel2[[CH]]$deviance, start=start, end=end)

# Model0: Scatterplots and histograms of mixture elements
PCH <- 1;  CEX <- 0.5
set.seed(770328)
SELECT <- sample(end-start+1, size=500, replace=FALSE)

# Fixed K MCMC
FixPrior2 <- list(priorK="fixed",
                  delta=1,
                  priormuQ="independentC")

Keep <- c("iter", "nMCMC", "dim", "prior", "init", "RJMCMC",
               "scale", "state", "freqK", "propK", "DIC", "moves",
               "pm.y", "pm.z", "pm.indDev", "pred.dens", "summ.y.Mean",
               "summ.y.SDCorr", "summ.z.Mean", "summ.z.SDCorr")

set.seed(770328)
FixModel2 <- list()
MPDensModelK <- list()
JPDensModelK <- list()

for (k in 1:5){
   cat(paste("K = ", k, "\n-------------------------------\n", sep=""))
   PriorNow <- FixPrior2
   PriorNow$Kmax <- k
   FixModel2[[k]] <- NMixMCMC(y0=my_dataset, prior=PriorNow, nMCMC=nMCMC,
                                   scale=list(shift=0, scale=1), PED=TRUE)
   cat(paste("\nComputation of pred. densities started on ", date(),
                   "\n", sep=""))
   MPDensModelK[[k]] <- list()
   MPDensModelK[[k]][[1]] <- NMixPredDensMarg(FixModel2[[k]][[1]], grid=ygrid)
   MPDensModelK[[k]][[2]] <- NMixPredDensMarg(FixModel2[[k]][[2]], grid=ygrid)
   cat(paste("Computation of pred. densities finished on ", date(),
                  "\n\n\n", sep=""))
   JPDensModelK[[k]] <- list()
   JPDensModelK[[k]][[1]] <- NMixPredDensJoint2(FixModel2[[k]][[1]], grid=ygrid)
   JPDensModelK[[k]][[2]] <- NMixPredDensJoint2(FixModel2[[k]][[2]], grid=ygrid)
   cat(paste("Computation of joint pred. densities finished on ", date(), 
             "\n\n\n", sep=""))   
   FixModel2[[k]][[1]] <- FixModel2[[k]][[1]][Keep]
   FixModel2[[k]][[2]] <- FixModel2[[k]][[2]][Keep]
   class(FixModel2[[k]][[1]]) <- class(FixModel2[[k]][[2]]) <- "NMixMCMC"}


print(FixModel2[[1]])

NMixCluster(FixModel2)

fix <- NMixEM(my_dataset[, names(my_dataset)[1:2]], K = 3)
y <- seq(12, 22, length=300)
fy <- dMVNmixture(y, weight=fix$weight, mean=fix$mean,
                  Sigma=c(fix$Sigma, matrix(rep(fix$K), nrow = 2, ncol = 2)))
hist(my_dataset$bill_depth, prob=TRUE,
     main="", xlab="Velocity (km/sec)", col="sandybrown")
lines(y, fy, col="darkblue", lwd=2)


# Fixed K: PED and DIC
PED <- RJModel2[[1]]$PED
DIC <- list(Chain1 = RJModel2[[1]]$DIC, Chain2 = RJModel2[[2]]$DIC)

for (k in 1:length(FixModel2)) {
   PED <- rbind(PED, FixModel2[[k]]$PED)
   DIC[[1]] <- rbind(DIC[[1]], FixModel2[[k]][[1]]$DIC)
   DIC[[2]] <- rbind(DIC[[2]], FixModel2[[k]][[2]]$DIC)}
                                                                
# Fixed K: print PED
print(PED)
# Fixed K: print DIC
print(DIC)

# Fixed K MCMC: Eruptions - plot of the marginal predictive density
par(mfrow=c(1, 1), bty="n")
hist(my_dataset$bill_depth, prob=TRUE, col="grey90",
     xlab=my_xlab, ylab="Density", main="")
par(mfrow=c(1, 1), bty="n")

for (k in 1:5){
   lines(MPDensModelK[[k]][[1]]$x$x1, MPDensModelK[[k]][[1]]$dens[[1]], col="red",
         main="Predictive Densities for Fixed K-Values")}

# Fixed K MCMC: Eruptions - plot of the marginal predictive density
for (k in 1:5){
   par(mfrow=c(1, 1), bty="n")
   hist(my_dataset$bill_depth, prob=TRUE, col="darkviolet", breaks=seq(7, 37, by=0.5),
          xlab=my_xlab, ylab="", main=paste("K = ", k, sep=""))
   lines(MPDensModelK[[k]][[1]]$x$x1, MPDensModelK[[k]][[1]]$dens[[1]], col="deeppink", lwd=2)}

# Fixed K MCMC: Waiting - plot of the marginal predictive density
par(mfrow=c(1, 1), bty="n")
hist(my_dataset$bill_length, prob=TRUE, col="grey90", 
     xlab=my_xlab, ylab="Density", main="")
for (k in 1:5){
   lines(MPDensModelK[[k]][[CH]]$x$x2, MPDensModelK[[k]][[CH]]$dens[[2]], 
         col="red")}

# Fixed K MCMC: Waiting - plot of the marginal predictive density
par(mar=c(3, 2, 2, 1)+0.1)
par(mfrow=c(5, 2), bty="n")
for (k in 1:5){
   hist(my_dataset$bill_length, prob=TRUE, col="lightblue", 
        xlab="", ylab="", main=paste("K = ", k, sep=""))
   lines(MPDensModelK[[k]][[CH]]$x$x2, MPDensModelK[[k]][[CH]]$dens[[2]], 
         col="red", lwd=2)} 

# more mixture modelling
plot_mix_comps <- function(x, mu, sigma, lam) {
   lam * dnorm(x, mu, sigma)}

set.seed(1)
wait <- penguins.clean$bill_depth_mm 
mixmdl <- normalmixEM(wait)

# plotting
data.frame(x = mixmdl$x) %>%
   ggplot() +
   geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
                  fill = "darkviolet") +
   stat_function(geom = "line", fun = plot_mix_comps,
                 args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                 colour = "deeppink", lwd = 1.5) +
   stat_function(geom = "line", fun = plot_mix_comps,
                 args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                 colour = "limegreen", lwd = 1.5) +
   ylab("Density")

set.seed(1)
wait <- penguins.clean$bill_length_mm #distance # #faithful$waiting
mixmdl <- normalmixEM(wait)

# plotting
data.frame(x = mixmdl$x) %>%
   ggplot() +
   geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
                  fill = "darkviolet") +
   stat_function(geom = "line", fun = plot_mix_comps,
                 args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                 colour = "deeppink", lwd = 1.5) +
   stat_function(geom = "line", fun = plot_mix_comps,
                 args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                 colour = "limegreen", lwd = 1.5) +
   ylab("Density")


## TRACEPLOTS FROM MIXTURE ANALYSIS
# Choose iters to draw traceplots
CH <-1 
titers <- tstart:tend
chK <- mcmc(RJModel2[[CH]]$K, start=start, end=end)
chgammaInv2 <- mcmc(RJModel2[[CH]]$gammaInv[titers,], start=tstart, end=tend)
chmixture2 <- mcmc(RJModel2[[CH]]$mixture[titers,], start=tstart, end=tend)
chdeviance2 <- mcmc(RJModel2[[CH]]$deviance[titers,], start=tstart, end=tend)

# Model0: Create mcmc objects
start <- RJModel2[[CH]]$nMCMC["burn"] + 1
end <- RJModel2[[CH]]$nMCMC["burn"] + RJModel2[[CH]]$nMCMC["keep"]
chK <- mcmc(RJModel2[[CH]]$K, start=start, end=end)
chgammaInv <- mcmc(RJModel2[[CH]]$gammaInv, start=start, end=end)
chmixture <- mcmc(RJModel2[[CH]]$mixture, start=start, end=end)
chdeviance <- mcmc(RJModel2[[CH]]$deviance, start=start, end=end)

# getting trace plots for both chains
# chain 1
start1 <- RJModel2[[1]]$nMCMC["burn"] + 1
end1 <- RJModel2[[1]]$nMCMC["burn"] + RJModel2[[1]]$nMCMC["keep"]
chK1 <- mcmc(RJModel2[[1]]$K, start=start, end=end)
chgammaInv1 <- mcmc(RJModel2[[1]]$gammaInv, start=start, end=end)
chmixture1 <- mcmc(RJModel2[[1]]$mixture, start=start, end=end)
chdeviance1 <- mcmc(RJModel2[[1]]$deviance, start=start, end=end)

# chain2
start2 <- RJModel2[[2]]$nMCMC["burn"] + 1
end2 <- RJModel2[[2]]$nMCMC["burn"] + RJModel2[[2]]$nMCMC["keep"]
chK2 <- mcmc(RJModel2[[2]]$K, start=start, end=end)
chgammaInv2 <- mcmc(RJModel2[[2]]$gammaInv, start=start, end=end)
chmixture2 <- mcmc(RJModel2[[2]]$mixture, start=start, end=end)
chdeviance2 <- mcmc(RJModel2[[2]]$deviance, start=start, end=end)

lwd <- 0.5   

# plotting trace plots
par(mfrow=c(2, 1), bty="n")
traceplot(chK1, smooth=TRUE, lwd=lwd, main="K Convergence",col="deeppink")
traceplot(chK2, smooth=TRUE, lwd=lwd, col="limegreen")

par(mfrow=c(2, 2), bty="n")
traceplot(chgammaInv1[, "gammaInv1"], smooth=FALSE,
          col="brown", lwd=lwd, main="gamma^{-1}")
traceplot(chgammaInv2[, "gammaInv1"], smooth=FALSE,
          col="brown", lwd=lwd, main="gamma^{-1}")

par(mfrow=c(2, 2), bty="n")
traceplot(chmixture1[, "y.Mean.1"], smooth=FALSE,
          col="deeppink", lwd=lwd, main="Mixture Mean for Bill Depth")
traceplot(chmixture1[, "y.Mean.2"], smooth=FALSE,
          col="deeppink", lwd=lwd, main="Mixture Mean for Bill Length")
traceplot(chmixture2[, "y.Mean.1"], smooth=FALSE,
          col="limegreen", lwd=lwd)
traceplot(chmixture2[, "y.Mean.2"], smooth=FALSE,
          col="limegreen", lwd=lwd)

par(mfrow=c(2, 2), bty="n")
traceplot(chmixture1[, "y.SD.1"], smooth=FALSE,
          col="deeppink", lwd=lwd, main="Standard Deviation for Bill Depth")
traceplot(chmixture1[, "y.SD.2"], smooth=FALSE,
          col="deeppink", lwd=lwd, main="Standard Deviation for Bill Length")
traceplot(chmixture2[, "y.SD.1"], smooth=FALSE,
          col="limegreen", lwd=lwd)
traceplot(chmixture2[, "y.SD.2"], smooth=FALSE,
          col="limegreen", lwd=lwd)

par(mfrow=c(2, 1), bty="n")
traceplot(chdeviance1[, "LogL0"], smooth=FALSE,
          col="red", lwd=lwd, main="Log(L0)")
traceplot(chdeviance2[, "LogL0"], smooth=FALSE,
          col="red", lwd=lwd, main="Log(L0)")

par(mfrow=c(2, 1), bty="n")
traceplot(chdeviance1[, "LogL1"], smooth=FALSE,
          col="red", lwd=lwd, main="Log(L1)")
traceplot(chdeviance2[, "LogL1"], smooth=FALSE,
          col="red", lwd=lwd, main="Log(L1)")

par(mfrow=c(2, 1), bty="n")
traceplot(chdeviance1[, "dev.complete"], smooth=FALSE,
          col="red", lwd=lwd, main="D(complete)")
traceplot(chdeviance2[, "dev.complete"], smooth=FALSE,
          col="red", lwd=lwd, main="D(complete)")

par(mfrow=c(2, 1), bty="n")
traceplot(chdeviance1[, "dev.observed"], smooth=FALSE,
          col="red", lwd=lwd, main="D(observed)")
traceplot(chdeviance2[, "dev.observed"], smooth=FALSE,
          col="red", lwd=lwd, main="D(observed)")

par(mfrow=c(2, 1), bty="n")
densplot(chK1, show.obs=TRUE, col="deeppink", main="Density of Estimated K Components")
densplot(chK2, show.obs=TRUE, col="limegreen")

par(mfrow=c(2, 1), bty="n")
densplot(chgammaInv1[, "gammaInv1"], show.obs=FALSE,
         col="brown", main="gamma^{-1}", xlim=c(0, 30))
densplot(chgammaInv2[, "gammaInv1"], show.obs=FALSE,
         col="brown", main="gamma^{-1}", xlim=c(0, 30))

par(mfrow=c(2, 1), bty="n")
densplot(chmixture1[, "y.Mean.1"], show.obs=FALSE,
         col="darkblue", main="EY", xlim=c(15, 25))
densplot(chmixture2[, "y.Mean.1"], show.obs=FALSE,
         col="darkblue", main="EY", xlim=c(15, 25))

par(mfrow=c(2, 1), bty="n")
densplot(chmixture1[, "y.SD.1"], show.obs=FALSE,
         col="darkblue", main="sd(Y)", xlim=c(0, 12))
densplot(chmixture2[, "y.SD.1"], show.obs=FALSE,
         col="darkblue", main="sd(Y)", xlim=c(0, 12))

par(mfrow=c(2, 1), bty="n")
autocorr.plot(chK1, auto.layout=FALSE, ask=FALSE,
              col="darkgreen", lwd=2, main="K")
autocorr.plot(chK2, auto.layout=FALSE, ask=FALSE,
              col="darkgreen", lwd=2, main="K")

par(mfrow=c(2, 1), bty="n")
autocorr.plot(chgammaInv1[, "gammaInv1"], auto.layout=FALSE, ask=FALSE,
              col="brown", lwd=2, main="gamma^{-1}")
autocorr.plot(chgammaInv2[, "gammaInv1"], auto.layout=FALSE, ask=FALSE,
              col="brown", lwd=2, main="gamma^{-1}")

par(mfrow=c(2, 1), bty="n")
autocorr.plot(chmixture1[, "y.Mean.1"], auto.layout=FALSE, ask=FALSE,
              col="darkblue", lwd=2, main="EY")
autocorr.plot(chmixture2[, "y.Mean.1"], auto.layout=FALSE, ask=FALSE,
              col="darkblue", lwd=2, main="EY")

par(mfrow=c(2, 1), bty="n")
autocorr.plot(chmixture1[, "y.SD.1"], auto.layout=FALSE, ask=FALSE,
              col="darkblue", lwd=2, main="sd(Y)")
autocorr.plot(chmixture2[, "y.SD.1"], auto.layout=FALSE, ask=FALSE,
              col="darkblue", lwd=2, main="sd(Y)")
