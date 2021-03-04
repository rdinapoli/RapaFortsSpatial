#load packages
library(maptools)
library(MuMIn)  
library(raster) 
library(rgdal) 
library(rgeos) 
library(sp) 
library(spatstat)
library(here)

#load data from working directory
pare_points <- readOGR(here('Data', "Pare_Projected.shp"))
shore_line <- readOGR(here('Data', "shore limit.shp"))
DTM <- raster(here('Data', "tingrid_Copy.tif"))
ag_vis <- raster(here('Data',"5m_ag34_vis_Copy.tif"))
tot_vis <- raster(here('Data', "5m_tot_vis2_Copy.tif"))

#Make slope variable
slope <- terrain(DTM, opt='slope', unit='radians', neighbors=8)
#aggregate slope because pare are on flattened areas, want average slope of neighborhood
slope <- aggregate(slope, fact=5)

#convert to spatstat format 
shore <- as.owin(shore_line) #convert to window format
pare <- ppp(pare_points$POINT_X, pare_points$POINT_Y, window=shore) #convert pare to a ppp object/ error? maybe use X_Cor
elev <- as.im(DTM) #convert DTM to a pixel image
ag_vis <- as.im(ag_vis) #convert g visibility to a pixel image
tot_vis <- as.im(tot_vis) #convert total visibility to a pixel image
slope <- as.im(slope)

###########################
### Exploratory Analyses ##
###########################

#Compute nearest neighbor distances for forts 
pare_nn <- nndist(pare)
#mean nn
mean(pare_nn)
#median_nn
median(pare_nn)

#Perform L function test againt 39 realizations of CSR with fixed number of points
set.seed(1234) #set random seed to get reproducible result
pare_L <- envelope(pare, fun=Lest,  fix.n=T,nsim=39) 
#check it
plot(pare_L)

#Examine form of possible covariate effects 
#using relative distribution (nonparametric regression (rhohat))
elev_rh <- rhohat(pare, elev, confidence = 0) #intensity as a function of elevation 
ag_vis_rh <- rhohat(pare, ag_vis, confidence = 0) #intensity as a function of visibility from agricultural plots 
tot_vis_rh <- rhohat(pare, tot_vis, confidence = 0) #intensity as a function of total visibility for the island 
slope_rh <- rhohat(pare, slope, confidence = 0) #intensity as a function of slope


###################
## Basic Maps #####
###################

# par(mfrow=c(2,2))
# plot(elev, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = FALSE), riblab="Elevation")
# plot(pare, pch=16, col='red', add=T)
# plot(shore, add=T)
# plot(slope, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = FALSE), riblab="Slope")
# plot(pare, pch=16, col='red', add=T)
# plot(shore, add=T)
# plot(tot_vis, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = FALSE), riblab="Visibility")
# plot(pare, pch=16, col='red', add=T)
# plot(shore, add=T)
# plot(ag_vis, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = FALSE), riblab="Visibility")
# plot(pare, pch=16, col='red', add=T)
# plot(shore, add=T)
# par(mfrow=c(1,1))

############
## Models ##
############

# Beginning of point process models               
#correction="none" so models do not calculate unobserved points outside the window, 
#since this is an island poinst could not exist in the ocean  
#the default is a "border correction" which uses the fixed threshold value and has 
#an erroded border (not really sure what it means) but either could be valid. 
ppm0 <- ppm(pare, ~1, correction = 'none')
ppm1 <- ppm(pare, ~elev, correction = 'none')
ppm2 <- ppm(pare, ~ag_vis, correction = 'none')
ppm3 <- ppm(pare, ~tot_vis, correction = 'none')
ppm4 <- ppm(pare, ~elev+ag_vis, correction = 'none')
ppm5 <- ppm(pare, ~elev+tot_vis, correction = 'none')
ppmS <- ppm(pare, ~slope, correction= 'none') #added models with slope
ppmSE <- ppm(pare, ~slope+elev, correction= 'none')
ppmSV <- ppm(pare, ~slope+tot_vis, correction= 'none')
ppmSEV <- ppm(pare, ~slope+elev+tot_vis, correction= 'none')

# AICc model comparison
ppm_AICc <- model.sel(ppm0, ppm1, ppm2, ppm3, ppm4, ppm5, ppmS, ppmSE, ppmSV, ppmSEV, rank=AICc)
ppm_AICc# Results
 
#print results of best fitting model
summary(ppm3)

# Compute Residual K-function to evaluate model fit in terms of second-order interaction 
set.seed(1234)
K_sim <- envelope(ppm3, Kres, nsim=99, fix.n=T)
#check fit
plot(K_sim, lwd=3, legend='F')


# Modeling with Gibbs  Hard-core process: forbidden to come close to oneanother within a specified distance   
#Inhomogeneous models include models that involve covariate effects (like elevation)  - 
# Worth considering a hybrid model as well that combines two different processes
# Need to add an edge correction to make sure it does not account for points "outside" of the window. Appropriate when borders are "real" borders like an island          correction="none"
#possibly use a hybrid of Strauss and Area Interaction? 
#Stick with Strauss, hybrid is better for multiple interactions, but this is really just one interaction between forts. 

#Strauss. range of pairwise interaction determined by looking at Kr throguh mean(K_sim$r) = 942  and the mean of nearest neighbor = 971 
ppm6 <- ppm(pare, ~ 1, Strauss(950), correction = 'none')
ppm7 <- ppm(pare, ~ tot_vis, Strauss(950), correction = 'none')
ppm8 <- ppm(pare, ~ elev, Strauss(950), correction = 'none')
ppm9 <- ppm(pare, ~ elev + tot_vis, Strauss(950), correction = 'none')
ppm10 <- ppm(pare, ~ ag_vis, Strauss(950), correction = 'none')
ppm11 <- ppm(pare, ~ elev + ag_vis, Strauss(950), correction = 'none')
ppmS2 <- ppm(pare, ~ slope, Strauss(950), correction= 'none') #added models with slope
ppmSE2 <- ppm(pare, ~ slope+elev, Strauss(950), correction= 'none')
ppmSV2 <- ppm(pare, ~ slope+tot_vis, Strauss(950), correction= 'none')
ppmSEV2 <- ppm(pare, ~ slope+elev+tot_vis, Strauss(950), correction= 'none')

#compare Gibbs models
ppm_AICc2 <- model.sel(ppm6, ppm7, ppm8, ppm9, ppm10, ppm11, ppmS2, ppmSE2, ppmSV2, ppmSEV2,rank=AICc)
ppm_AICc2

#print results of best fitting model
summary(ppm7)

#check residual K
set.seed(1234)
K_sim2 <- envelope(ppm7, Kres, nsim=99, fix.n=T)
#check fit
plot(K_sim2, lwd=3, legend='F')

#plot predicted first-order intensity of best-fitting model
plot(intensity.ppm(ppm7), 
     col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = T),
     main="", riblab="Fitted intensity")
plot(pare, pch=16, col='red', cex=0.25, add=T)

#partial residual plot for elevation
par_res <- parres(ppm7, "tot_vis")
#check fit
plot(par_res)

#######################
#######################
#### Print Figures ####
#######################
#######################
#Nearest neighbor and L function plot
jpeg(file=here('Figures','temporary','NN_and_Lfunction.jpeg'),width = 8, height = 4,units='in',res=300)
par(mfrow=c(1,2))
hist(pare_nn, xlab="meters", breaks=8, col="grey", xlim=c(0,1500), ylim=c(0,5), main="") #creates histogram
mtext(side=3, line=1, at=-1, adj=0, cex=1, "a)")  # creates the textline above the figure
plot(pare_L, lwd=3, xlab="r (meters)", main="", legend=F)
mtext(side=3, line=1, at=-1, adj=0, cex=1, "b)")
par(mfrow=c(1,1))
dev.off()


#plot relative distributions
jpeg(file=here('Figures','temporary','relative_distributions.jpeg'),width = 8, height = 8,units='in',res=300)
par(mfrow=c(2,2))
plot(elev_rh, legend=F, xlab=("Elevation ASL"), main="", xlim=c(0,500), lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "a)")
plot(slope_rh, legend=F, xlab=("Slope (rad)"), main="", xlim=c(0,1.2), lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "b)")
plot(tot_vis_rh, legend=F, xlab=("Visibility"), main="", lwd=3)
mtext(side=3, line=1, at=2550, adj=0, cex=0.9, "c)")
plot(ag_vis_rh, legend=F, xlab=("Ag. Visibility"), main="", lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "d)")
par(mfrow=c(1,1))
dev.off()


#residual diagnostic plots
jpeg(file=here('Figures','temporary','Residual_diagnostics.jpeg'),width = 8, height = 8,units='in',res=300)
par(mfrow=c(2,2))
plot(K_sim, lwd=3, main='', legend=F, xlab="r (meters)") #model 3
mtext(side=3, line=1, at=0, adj=0, cex=1, "a)")
plot(K_sim2, lwd=3, main='', legend=F, xlab="r (meters)") #model 7
mtext(side=3, line=1, at=0, adj=0, cex=1, "b)")
plot(par_res, legend=F, main='', xlab="Visibility")
mtext(side=3, line=1, at=0, adj=0, cex=1, "c)")
par(mar=c(2.25,2.25,2.25,2.25))
plot(intensity.ppm(ppm7), 
     col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = T),
     main="") #riblab="Fitted intensity")
plot(pare, pch=16, col='red', cex=0.5, add=T)
mtext(side=3, line=1, at=759002, adj=0, cex=1, "d)")
par(mfrow=c(1,1))
dev.off()
