# Script for generating age-specific net migration given total net migration
# for the Bayesian FDM and the Deterministic FDM.
# Bayesian FDM: 
#.   It is expected that an estimation has been run (via e.g. estimateBayesFDMwacounty.R)
#.   and MCMC samples of the RC curves have been exported into csv file 
#.   (via e.g. extract_rc_samples.R).
# Deterministic FDM:
#.   This script contains the estimation as well as projection of the deterministic FDM.
#
# For splitting total net migration into total in- and total out-migration,
# we show the use of both, using estimates from a mixed-effects model (MEM), 
# as well as using the heuristic model.
#
# Input: - total net migration of a location (county) to be distributed into ages
#        - trajectories of Rogers-Castro in- and out-curves 
#          extracted from STAN estimation via extract_rc_samples.R
#        - results from a MEM estimation (provided in file MEMest_wacounties.csv)
#        - age-specific population of the same location
#        - age-specific population of a larger geography (state)
#        - For the deterministic FDM age-specific migration over 
#          several time periods is needed.
#
# Hana Sevcikova, 9/24/2025
#
library(data.table)

data.dir <- "."

# county and year to use data for
county <- 53065 # use numerical code, here Stevens
yr <- 2010 # in the pop file 2010 corresponds 2010-2020

# file with the observed data (migration and population)
migpop.file <- file.path(data.dir, "wacounty5x10.csv")
agemigpop <- fread(migpop.file)

# extract total migration and age-specific and total population
G10 <- agemigpop[fips == county & year == yr, sum(migcount)] # total net migration for 10 years
G <- G10 / 10 # average annual total migration
Pxi <- agemigpop[fips == county & year == yr, .(age, pop)] # county age-specific population
Pxw <- agemigpop[name == "State" & year == yr, .(age, pop)] # global age-specific population
P <- Pxi[, sum(pop)] # county total population

# G is our total net migration to distribute into ages using both methods
# P is the corresponding population

###############
# Bayesian FDM
###############
input.dir <- "outputs-2025-09-19/RC" # where are estimation results from the Bayesian FDM

# load input data
samples.file <- file.path(input.dir, "RCcurves_samples5y1000.csv") # file with MCMC samples
mem.file <- file.path(data.dir, "MEMest_wacounties.csv") # MEM results
rcdat <- fread(samples.file)
memdat <- fread(mem.file)

# set random number seed for reproducibility
set.seed(12345)

# split G into A and B using the MEM method (Eq. 18)
A <- memdat[reg_code == county, pmax(b0 * P + b1 * G, G + min_rate * P)]
B <- A - G

# create a dataset with age-specific quantities
agedat <- dcast(rcdat[LocID == county & Parameter != "v"], Age + Trajectory ~ Parameter, value.var = "Value")
setnames(agedat, c("Age", "in", "out"), c("age", "rc_in", "rc_out"))

# add population and the v parameter
agedat[Pxi, p_xi := i.pop, on = "age"] # county pop
agedat[Pxw, p_xw := i.pop, on = "age"] # global pop
agedat[rcdat[LocID == county & Parameter == "v"], v := Value, on = .(Trajectory)]

# compute r^star (population weighted and normalized RC) (Eqs. 27, 28)
agedat[, `:=`(rc_star_in = rc_in * p_xw / sum(rc_in * p_xw), 
              rc_star_out = rc_out * p_xi/ sum(rc_out * p_xi)), by = .(Trajectory)]

# compute age-specific in-migration (iota) and out-migration (o) (Eqs. 27, 28)
agedat[, `:=`(iota_xi = A * rc_star_in, o_xi = B * rc_star_out)]

# compute the standard deviation (Eq. 30)
agedat[, `:=`(sigma = sqrt(pmin((iota_xi + o_xi)/v, (p_xi/2)^2)))]

# draw age-specific net migration from normal distribution (Eq. 29)
agedat[, `:=`(gtilde_xi = rnorm(nrow(agedat), iota_xi - o_xi, sigma))]

# scale so that the total net migration is equal to G (Eq. 31)
agedat[, `:=`(g_xi = gtilde_xi - sigma / sum(sigma) * (sum(gtilde_xi) + B - A)), by = .(Trajectory)]

# check that all trajectories across ages sum to G
all.equal(range(agedat[, sum(g_xi), by = "Trajectory"]$V1), c(G, G))

# compute quantiles
migq <- agedat[, .(med = median(g_xi), low80 = quantile(g_xi, 0.1), 
                   high80 = quantile(g_xi, 0.9)), by = .(age)]


###################
# Deterministic FDM
###################
m <- 0.7 # choice of m on 10-year scale

# Model Rogers-Castro
rogers.castro <- function(x, a1=0.01, alpha1=0.09, a2=0.05, alpha2=0.077, mu2=16.5, 
                          lambda2=0.374, c=0.0003){
    return(a1*exp(-alpha1*x) + a2*exp(-alpha2*(x - mu2) - exp(-lambda2*(x - mu2)
    )) + c)
}

# Estimation
# ----------
# Create dataset with the model RC schedule
ThetaM <- data.table(age = seq(0, 95, by = 5))[, rxM := rogers.castro(age)][, rxM := rxM/sum(rxM)]
#plot(ThetaM$age, ThetaM$rxM, type = "l")

# split G into A and B using the heuristic method
ABdat <- agemigpop[fips == county, .(P = sum(pop), G = sum(migcount)), by = "year"]
ABdat[, `:=`(A = P * m + 0.5 * G, B = P * m - 0.5 * G)]

# extract data for all time periods
agedat.allT <- agemigpop[fips == county, .(age, year, g_xi = migcount)]

# merge with the other datasets
agedat.allT <- merge(merge(agedat.allT, ThetaM, by = "age"), ABdat, by = "year")

# Equation 14
gbar <- agedat.allT[, .(gbar_xi = mean(g_xi)), by = "age"]
agedat.allT[gbar, iota_xit := A * rxM + 0.5 * i.gbar_xi, on = "age"]

# Equation 15
agedat.allT[, iota_star_xit := iota_xit / A]
Rx <- agedat.allT[, .(Rx = mean(iota_star_xit/rxM)), by = "age"]

# Projection for 2020 for one year
# --------------------------------
# Split G decreasing m to annual scale
A <- P * m/10 + 0.5 * G
B <- P * m/10 - 0.5 * G

# Equations 24-26
agedatD <- merge(ThetaM, Rx, by = "age")
agedatD[, `:=`(iota_xi = A * rxM * Rx / sum(rxM * Rx), o_xi = B * rxM * 1/Rx / sum(rxM * 1/Rx))]
agedatD[, `:=`(g_xi = iota_xi - o_xi)]

# check that the sum is equal to the give G
all.equal(agedatD[, sum(g_xi)], G)

################################
# Plot results from both methods
################################
# Bayesian FDM
plot(migq$age, migq$med, type = "l", col = "red", ylim = range(migq$low80, migq$high80),
     xlab = "Age", ylab = "net migration", main = "Age-specific net migration fit")
abline(h = 0, col = "grey")
lines(migq$age, migq$low80, col = "red", lty = 2)
lines(migq$age, migq$high80, col = "red", lty = 2)

# Deterministic FDM
lines(agedatD$age, agedatD$g_xi, col = "green")

# Observed data (per year)
observed <- agemigpop[fips == county & year == yr, .(age, migcount)]
lines(observed$age, observed$migcount/10, col = "blue")

legend("bottomright", legend = c("Bayesian FDM", "Bayesian FDM 80% PI","Deterministic FDM", "Observed"),
       col = c("red", "red", "green", "blue"), lty = c(1,2,1,1), bty = "n")
