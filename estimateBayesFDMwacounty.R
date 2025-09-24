# Script to estimate Bayesian FDM's Rogers-Castro curves via MCMCs using Stan
# Install cmdstanr via (from https://mc-stan.org/cmdstanr/):
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
#
# Hana Sevcikova, 9/19/2025

library(data.table)
library(cmdstanr)

options(cmdstanr_warn_inits = FALSE) # switches off some warnings on missing initial values

# === Settings ====
########################
# where are input data
data.dir <- "."

# name of the input file containing age-specific migration and population 
agemig.file <- file.path(data.dir, "wacounty5x10.csv")

# which file holds the MEM estimates to split totals into in- and out-totals
MEMest.file <- file.path(data.dir, "MEMest_wacounties.csv")

# should the STAN estimation run 
# (can be set to FALSE if results exist and only a plot is desired)
run.stan <- FALSE

# Name of the file containing the STAN model
stan.file <- "bayesian_fdm.stan"

# should plots with results be created
plot.results <- TRUE

# output directory
output.dir <- paste0("outputs-", Sys.Date())

# MCMC setup
nr.chains <- 3 # cannot be larger than what is used in the get.init.values() function
warmup <- 10000
iter <- 100000
thin <- 50

# ==== Functions =====

get.init.values <- function(rc_config, nchains){
    # Function to set initial parameter values 
    # rc_config is a list determining which RC components are included for in- and out-migration,
    # e.g. pre_working_age_in, pre_working_age_out, etc.
    
    # initial values of RC parameters for 3 chains
    # (it's used for both, in- and out-RC)
    pars.all <- list(c(a1 = 0.01, alpha1 = 0.1, a2 = 0.1, alpha2 = 0.1, mu2 = 21, lambda2 = 0.3, a3 = 0.001, alpha3 = 0.05, mu3 = 67, lambda3 = 0.2, a4 = 0.01, lambda4 = 0.02, c = 0.001),
                     c(a1 = 0.05, alpha1 = 0.2, a2 = 0.2, alpha2 = 0.2, mu2 = 30, lambda2 = 0.5, a3 = 0.01, alpha3 = 0.05, mu3 = 60, lambda3 = 0.9, a4 = 0.01, lambda4 = 0.01, c = 0.008),
                     c(a1 = 0.02, alpha1 = 0.45, a2 = 0.25, alpha2 = 0.3, mu2 = 25, lambda2 = 1, a3 = 0.02, alpha3 = 0.25, mu3 = 64, lambda3 = 0.6, a4 = 0.01, lambda4 = 0.01, c = 0.0001)
                     )
    # initial values of the "v" parameter for 3 chains
    pars.com <- list(c(v = 0.5),
                     c(v = 0.2),
                     c(v = 0.9)
                    )
    # check if we have initial values for all chains that we need
    if(length(pars.all) < nchains || length(pars.com) < nchains)
        stop("Initial values defined for less chains than required.")
    
    # exclude RC components based on the settings in rc_config
    # and put everything together
    rc.components <- list( # which parameters correspond to which RC component
        pre_working_age = c("a1", "alpha1"),
        working_age = c("a2", "alpha2", "mu2", "lambda2"),
        retirement = c("a3", "alpha3", "mu3", "lambda3"),
        post_retirement = c("a4", "lambda4")
    )
    allpars.list <- list()
    for(chain in 1:min(length(pars.all), nchains)) { # iterate over chains
        pars.in <- pars.all[[chain]] # parameters for the in-RC
        pars.out <- pars.all[[chain]] # parameters for the out-RC
        # remove parameters based on the config
        remove <- list("in" = c(), "out" = c())
        for(comp in names(rc.components)){
            for(migtype in c("in", "out")){
                comp.name <- paste(comp, migtype, sep = "_")
                if(!is.null(rc_config[[comp.name]]) && !rc_config[[comp.name]])
                    remove[[migtype]] <- c(remove[[migtype]], rc.components[[comp]])
            }
        }
        # collect parameters for this chain
        inpars <- outpars <- list()
        for(par in setdiff(names(pars.in), remove[["in"]]))
            inpars[[paste0(par, "_in")]] <- pars.in[par]
        for(par in setdiff(names(pars.out), remove[["out"]]))
            outpars[[paste0(par, "_out")]] <- pars.out[par]
        allpars.list[[chain]] <- c(inpars, outpars, pars.com[[chain]])
    } # end of chain loop
    return(allpars.list)
}

# ==== Processing =====
#
# load age-specific migration and population
dat.all <- fread(agemig.file)

# compile STAN model
if(run.stan) mod <- cmdstan_model(stan.file)

# Set which region uses which of the four Rogers-Castro's components 
# (1. pre-working age, 2. working age, 3. retirement, 4. post-retirement)
# Each county has 8 values: first 4 are the four components for in-migration, 
#                           the last 4 are for out-migration.
# By default the first three components are used for both, in and out.
# If county is missing, the 'default' item will be used.
settings <- c(list(
                default = c( 1, 1, 1, 0, 1, 1, 1, 0),
                Stevens =  c(1, 1, 1, 0, 1, 1, 0, 0)
             ))

# extract parameters from the estimation of the mixed-effects model (MEM)
memest <- fread(MEMest.file)

# As the MEM was estimated for annual data, we need to convert it 
# to 10 year value
memest[, `:=`(b0 = b0 * 10, min_rate = min_rate * 10)]

# Compute quantities necessary to split total net migration into A (in) and B (out) 
# sum mig and pop over ages 
datAB.all <- dat.all[, .(G = sum(migcount, na.rm = TRUE), P = sum(pop, na.rm = TRUE)), 
                     by = c("name", "year")]
# merge with MEM estimates
datAB.all <- merge(datAB.all, memest, by = "name")
# compute A and B
datAB.all[, `:=`(A = pmax(P * b0 + b1*G, min_rate * P, G + min_rate * P))][
    , `:=`(B = A - G)]

# which counties should this run for
counties <- "Stevens"

# The STAN estimation runs for individual counties.
for(cnty in counties){ 
    # create the output directory if needed
    if(!dir.exists(output.dir))
        dir.create(output.dir, recursive = TRUE)
    
    # extract data for this county
    dat <- dat.all[name == cnty] # age-specific
    datAB <- datAB.all[name == cnty] # totals
    
    # extract which RC components to use
    rxsetting <- if(cnty %in% names(settings)) settings[[cnty]] else settings[["default"]]
    rc.config <- list( # all values below are either 0 or 1
                        pre_working_age_in = rxsetting[1], 
                        pre_working_age_out = rxsetting[5],
                        working_age_in = rxsetting[2],
                        working_age_out = rxsetting[6],
                        retirement_in = rxsetting[3],
                        retirement_out = rxsetting[7],
                        post_retirement_in = rxsetting[4],
                        post_retirement_out = rxsetting[8]
                    )
    
    # extract age and time variables
    ages <- sort(unique(dat$age))
    years <- sort(unique(dat$year))
    N <- length(ages)
    Tm <- length(years)

    # create data matrices and vectors needed for STAN
    mcounts <- popu <- popglob <- matrix(NA, nrow = Tm, ncol = N) # matrices of size times x ages 
    A <- B <- rep(NA, Tm)
    for(iyr in 1:Tm){
        yr <- years[iyr]
        d <- dat[year == yr]
        mcounts[iyr,] <- d$migcount
        popu[iyr,] <- d$pop
        popglob[iyr,] <- dat.all[name == "State" & year == yr, pop] # use State pop as weights for in-migration
        A[iyr] = datAB[year == yr, A]
        B[iyr] = datAB[year == yr, B]
    }
    # fill 0 for missing data
    mcounts[is.na(mcounts)] <- 0
    popu[is.na(popu)] <- 0

    # create the final STAN data object
    stan.data <- c(list(
                    N = N,
                    Tm = Tm,
                    g = round(mcounts),
                    x = ages,
                    pop = round(popu),
                    pop_glob = round(popglob),
                    A = A,
                    B = B
                    ),
                  rc.config
                )
    
    rds.file <- file.path(output.dir, paste0("fit_", cnty, ".rds")) # where the stan results will go
    
    if(run.stan){
        cat("\nEstimating ", cnty, "(", datAB[1, reg_code], ")")
        # run the STAN model
        fit <- mod$sample(data = stan.data, chains = nr.chains, iter_warmup = warmup,
                  save_warmup = FALSE, iter_sampling = iter, thin = thin,
                  output_dir = output.dir, output_basename = paste("samples", cnty, sep = "_"),
                  init = get.init.values(rc.config, nr.chains),
                  refresh = 5000, parallel_chains = nr.chains, seed = 1234,
                  show_exceptions = FALSE
                  )
        # save the resulting object
        fit$save_object(file = rds.file)
    } else { # if not running STAN, but plotting results, load existing results
        if(plot.results){
            if(!file.exists(rds.file))
                stop("File ",  rds.file, " does not exist.")
            fit <- readRDS(rds.file)
        }  
    }
    
    # PLOT results
    ###############
    if (plot.results) {
        # extract samples
        samples <- list()
        for(par in c("g_pred", "rx_in", "rx_out", "rx_star_in", "rx_star_out")){
            samples[[par]] <- fit$draws(variables = par, format = "array")
            # collapse the chains into one dimension
            if(par %in% c("g_pred", "rx_star_in", "rx_star_out")) # 3d -> 2d 
                samples[[par]] <- array(samples[[par]], 
                                        dim = c(dim(samples[[par]])[1] * dim(samples[[par]])[2], 
                                                dim(samples[[par]])[3]/N, N)
                )
            # 2d -> 1d
            else samples[[par]] <- matrix(samples[[par]], 
                                          nrow = dim(samples[[par]])[1] * dim(samples[[par]])[2])
        }
        set.seed(1234)
        nr.traj <- 50 # number of trajectories to plot
        sindex <- sample(1:dim(samples$g_pred)[1], min(nr.traj, dim(samples$g_pred)[1]))
        
        pdf(file.path(output.dir, paste0("results_", cnty, ".pdf")), width = 13, height = 9)
        par(mfrow=c(3,4))
        par(mar=c(3,2,1,0.4)+0.1, mgp=c(1,0.3,0),  oma = c(0, 0, 2, 0))
        
        ## leave one plot frame empty
        plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
        legend("center", legend = paste(cnty, "County"), bty = "n", cex = 2, xjust = 0.5)
        
        # plot r_x
        rx_in <- samples$rx_in
        rx_out <- samples$rx_out
        # rx_in
        plot(ages, rx_in[sindex[1],], col = '#00000020', ylim = range(rx_in[sindex,], -rx_out[sindex,]),
             type = "l", main = "In-/Out RC", ylab = "r_x", xlab = "age")
        for(i in 2:nr.traj)
            lines(ages, rx_in[sindex[i],], col = '#00000020')
        abline(h = 0, col = "grey", lty = 2)
        lines(ages, apply(rx_in, 2, median), col = "blue", lwd = 2)
        # rx_out
        lines(ages, -rx_out[sindex[1],], col = '#00000020', type = "l")
        for(i in 2:nr.traj)
            lines(ages, -rx_out[sindex[i],], col = '#00000020')
        lines(ages, -apply(rx_out, 2, median), col = "green", lwd = 2)
        
        # plot r_x_star (for latest year)
        rx_in <- samples$rx_star_in[, Tm, ]
        rx_out <- samples$rx_star_out[, Tm, ]
        # rx_star_in
        plot(ages, rx_in[sindex[1],], col = '#00000020', ylim = range(rx_in[sindex,], -rx_out[sindex,]),
             type = "l", main = paste0("normalized RC * Pop (", years[Tm], ")"), 
             ylab = "r_x_star", xlab = "age")
        for(i in 2:nr.traj)
            lines(ages, rx_in[sindex[i],], col = '#00000020')
        abline(h = 0, col = "grey", lty = 2)
        lines(ages, apply(rx_in, 2, median), col = "blue", lwd = 2)
        # rx_star_out
        lines(ages, -rx_out[sindex[1],], col = '#00000020', type = "l")
        for(i in 2:nr.traj)
            lines(ages, -rx_out[sindex[i],], col = '#00000020')
        lines(ages, -apply(rx_out, 2, median), col = "green", lwd = 2)
        
        # plot median r_x_star for all years
        rx_star_in <- apply(samples$rx_star_in, c(2,3), median)
        rx_star_out <- apply(samples$rx_star_out, c(2,3), median)
        cols <- topo.colors(Tm)
        plot(ages, rx_star_in[1,], col = cols[1], ylim = range(rx_star_in, -rx_star_out),
             type = "l", main = "normalized RC * Pop by time", ylab = "r_x_star", xlab = "age")
        lines(ages, -rx_star_out[1,], col = cols[1])
        abline(h = 0, col = "grey", lty = 2)
        for(t in 2:Tm){
            lines(ages, rx_star_in[t,], col = cols[t])
            lines(ages, -rx_star_out[t,], col = cols[t])
        }
        legend("topright", legend = years, col = cols, bty = "n", lty = 1)
        
        # plot A & B
        plot(years, rowSums(mcounts), col = 'red', ylim = range(A, -B),
             type = "l", main = "In-/Out-migration totals", ylab = "migration counts", xlab = "", lwd = 2)
        abline(h = 0, col = "grey", lty = 2)
        lines(years, A, col = "blue", lwd = 2)
        lines(years, -B, col = "green", lwd = 2)
        legend("topleft", legend = c(paste0("beta0 = ", round(datAB$b0[1],3)), 
                                     paste0("beta1 = ", round(datAB$b1[1],3))), bty = "n")
        legend("topright", legend = c("G", "A", "-B"), col = c("red", "blue", "green"), bty = "n", lty = 1)
        
        # plot predicted net migration counts
        for(t in 1:Tm) {
            med_pred <- rx_star_in[t,] * A[t] - rx_star_out[t,] * B[t]
            plot(ages, samples$g_pred[sindex[1],t,], col = '#00000020', ylim = range(samples$g_pred[sindex,t,], mcounts[t,], med_pred, 0),
                 type = "l", main = paste(years[t], years[t]+10, sep = " - "), ylab = "net migration counts", xlab = "age")
            for(i in 2:nr.traj){
                lines(ages, samples$g_pred[sindex[i],t,], col = '#00000020')
            }
            abline(h = 0, col = "grey", lty = 2)
            lines(ages, mcounts[t,], col = "blue", lwd = 2)
            lines(ages, med_pred, col = "red", lwd = 2)
            legend("topleft", legend = c("predicted", "observed"), col = c("red", "blue"), lwd = 2, lty = 1, bty = "n")
        }
        dev.off()
    }
}

cat("\nProcess finished.\n")

