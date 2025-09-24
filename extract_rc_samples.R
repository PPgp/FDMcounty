# Script for extracting estimated samples of Rogers-Castro curves
# from the Bayesian FDM while interpolating from 5-year to 1-year age groups.
# It creates two files, one with the median values 
# and one with the RC trajectories.
# 
# Hana Sevcikova, 9/19/2025
#


library(data.table)

data.dir <- "."
input.dir <- "outputs-2025-09-19"
output.dir <- file.path(input.dir, "RC")

# file from which to extract county codes (must have columns "reg_code" and "name")
file.with.county.codes <- file.path(data.dir, "MEMest_wacounties.csv")

nsamples <- 1000
convert.to.1y <- FALSE # should 5-year age groups be converted to 1-year

output.file.samples <- paste0("RCcurves_samples", if(convert.to.1y) "1y" else "5y", nsamples, ".csv")
output.file.medians <- paste0("RCcurves_medians", if(convert.to.1y) "1y" else "5y", ".txt")

# create the output directory if needed
if(!dir.exists(output.dir))
    dir.create(output.dir, recursive = TRUE)

# extract names of counties from file names
files <- list.files(input.dir, pattern = "fit_")
counties <- gsub("fit_|.rds", "", files)

# get county codes
locs <- unique(fread(file.with.county.codes)[, .(reg_code, name)])

# input and output ages
ages <- seq(0, 95, by = 5)
nages <- length(ages)
out.ages <- 0:100

set.seed(1234)

res.samples <- res.medians <- NULL
for(loc in counties){
    cat("\nProcessing", loc)
    code <- locs[name == loc, reg_code]
    rds.file <- file.path(input.dir, paste0("fit_", loc, ".rds"))
    fit <- readRDS(rds.file)
    samples <- fit$draws(variables = c("rx_in", "rx_out", "v"), format = "array")
    samples <- array(samples, dim = c(dim(samples)[1] * dim(samples)[2], # flatten the chain dimension
                                      dim(samples)[3]))
    set.seed(1234)
    sample.index <- sample(1:nrow(samples), nsamples)

    # construct trajectories dataset
    this.res.samples <- NULL
    for(iag in 1:length(ages)){ # iterate over ages
        this.res.samples <- rbind(this.res.samples, 
                                data.table(
                                    LocID = code,
                                    Age = ages[iag], 
                                    Trajectory = 1:nsamples,
                                    Value = samples[sample.index, iag],
                                    Parameter = "rc_in"
                                 ),
                                data.table(
                                    LocID = code,
                                    Age = ages[iag], 
                                    Trajectory = 1:nsamples,
                                    Value = samples[sample.index, iag + nages],
                                    Parameter = "rc_out"
                                )
                            )
    }
    # construct median dataset
    meds <- apply(samples, 2, median)
    this.res.med <- data.table(
                        reg_code = code,
                        age = ages, 
                        rc_in = meds[1:nages],
                        rc_out = meds[(nages + 1):(2*nages)]
                    )
    if(convert.to.1y){
        in.ages <- c(sort(unique(this.res.samples$Age)), 100) # will set RC at age 100 to 0
        interp.dat.samples <- NULL
        interp.dat.med <- data.table(country_code = code, age = out.ages)
        for(attr in c("rc_in", "rc_out")){
            # interpolate trajectories
            for(traj in sort(unique(this.res.samples$Trajectory))){
                this.dat <- this.res.samples[Trajectory == traj & Parameter == attr]
                rc <- approx(in.ages, c(unlist(this.dat[, Value]), 0), xout = out.ages)$y
                if(sum(is.na(rc))) stop("NAs during age interpolation")
                rc <- rc/sum(rc)
                interp.dat.samples <- rbind(interp.dat.samples, 
                                        data.table(LocID = code, Age = out.ages,
                                             Trajectory = traj,
                                             Value = rc, Parameter = attr)
                                    )
            }
            # interpolate medians
            rc <- approx(in.ages, c(unlist(this.res.med[, attr, with = FALSE]), 0), xout = out.ages)$y
            if(sum(is.na(rc))) stop("NAs during age interpolation")
            rc <- rc/sum(rc)
            interp.dat.med <- merge(interp.dat.med, 
                                        data.table(age = out.ages, value = rc))
            setnames(interp.dat.med, "value", attr)
        }
        this.res.samples <- interp.dat.samples
        this.res.med <- interp.dat.med
    }
    # append samples for v
    res.samples <- rbind(res.samples, this.res.samples, 
                        data.table(
                            LocID = code,
                            Age = NA, 
                            Trajectory = 1:nsamples, 
                            Value = samples[sample.index, 2*nages + 1],
                            Parameter = "v"
                         )
                     )
    res.medians <- rbind(res.medians, this.res.med)
}
if(!is.null(res.samples)) {
    res.samples[, Parameter := gsub("rc_", "", Parameter)]
    fname <- file.path(output.dir, output.file.samples)
    fwrite(res.samples, file = fname)
    cat("\nRC trajectories saved into ", fname, "\n")
} else
    warning("No RC trajectories written.")

if(!is.null(res.medians)) {
    setnames(res.medians, c("rc_in", "rc_out"), c("in", "out"))
    fname <- file.path(output.dir, output.file.medians)
    fwrite(res.medians, file = fname, sep = "\t")
    cat("\nRC medians saved into ", fname, "\n")
} else
    warning("No RC medians written.")


