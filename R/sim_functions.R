#' Simulate a single subset correlation of a randomly drawn subset of specified size
#'
#' @param dat data.frame or matrix with two columns which include the data
#' @param fraction Size of the subset as percentage of nrow(dat), e.g. 0.3
#' @param sampling (logical) If TRUE, sample round(nrow(dat) * fraction) rows from
#' the data to calculate the subset correlation. Otherwise take the first
#' round(nrow(dat) * fraction) rows
sim_subcor <- function(dat, fraction, sampling = T){
    stopifnot(ncol(dat) == 2)
    if (sampling) {
        inSubset <- sample(1:nrow(dat), size = round(nrow(dat) * fraction),
                           replace = F)
        subset_data<-dat[inSubset, ]
    } else {
        inSubset <- 1:round(nrow(dat)*fraction)
        subset_data<-dat[inSubset, ]
    }
    cor_full <- cor(dat$x,dat$y)
    cor_subset <- cor(subset_data$x,subset_data$y)
    return(list(cor_full = as.numeric(cor_full),
                cor_subset = as.numeric(cor_subset),
                inSubset = inSubset))
}

#' Simulate independent normally distributed data
#'
#' @param n Sample size
#' @param mu1 Mean of x
#' @param mu2 Mean of y
#' @param sd1 Sd of x
#' @param sd2 Sd of y
getBiCop_indep <- function(N, mu1 = 0, mu2 = 0, sd1 = 1, sd2 = 1) {
    x <- rnorm(N, mu1, sd1)
    y <- rnorm(N, mu2, sd2)
    return(data.frame(x = x, y = y))
}


#' Simulate bivariate normally distributed data with specified correlation, mean and sd
#'
#' The correlation between x and y will be exactly rho.
#'
#' Code from http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
#' @param N Sample size
#' @param rho True correlation coefficient in the sample
#' @param mu1 Expected value of the mean of variable x
#' @param sd1 Expected value of the standard deviation of variable x
getBiCop_exact <- function(N, rho, rho_sub = NA, mu1 = 0, sd1 = 1, fraction = NA) {
    if (sum(is.na(c(rho_sub, fraction))) == 1) stop("Specify both rho_sub and fraction")
    # No differing subset
    if (is.na(rho_sub)) {
        theta <- acos(rho)                            # corresponding angle
        x1    <- rnorm(N, mu1, sd1)                   # fixed given data
        x2    <- rnorm(N, 0, 1)                   # new random data
        X     <- cbind(x1, x2)                        # matrix
        Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

        Id   <- diag(N)                               # identity matrix
        Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
        P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
        x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
        Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
        Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

        x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]
        return(data.frame(x = x1, y = x))
    } else {
        # Different subset
        n_subset <- round(N * fraction)
        n_nonsubset <- N - n_subset
        full_sample <- getBiCop_exact(N = n_nonsubset, rho = rho)
        subset <- getBiCop_exact(N = n_subset, rho = rho_sub)
        df <- rbind(subset, full_sample)
        return(df)
    }
}

#' Simulate bivariate normally distributed data with specified correlation, mean and sd
#'
#' The variables x and y stem from a distribution with correlation rho, but the
#' actual correlation between x and y will vary around that true value.
#'
#' Code from http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
#' @param n Sample size.
#' @param rho True correlation coefficient in the population. If rho_sub is given,
#' this is the correlation in the rest of the sample without the subset.
#' @param rho_sub True correlation coefficient in the subset. Default is no
#' seperate subset.
#' @param fraction The size of the subset as a percentage of the full sample.
#' The subset is then contained in the first round(N * fraction) rows of the
#' returned data frame.
getBiCop <- function(N, rho, rho_sub = NA, fraction = NA, mar.fun=rnorm, x = NULL, ...) {
    if (!is.null(x)) {X1 <- x} else {X1 <- mar.fun(N, ...)}
    if (!is.null(x) & length(x) != N) warning("Variable x does not have the same length as N!")
    if (sum(is.na(c(rho_sub, fraction))) == 1) stop("Specify both rho_sub and fraction")

    # No differing subset
    if (is.na(rho_sub)) {
        C <- matrix(rho, nrow = 2, ncol = 2)
        diag(C) <- 1
        C <- chol(C)
        X2 <- mar.fun(N)
        X <- cbind(X1,X2)
        # induce correlation (does not change X1)
        df <- X %*% C
        df <- data.frame(df)
        colnames(df) <- c("x", "y")
        return(df)
    # Subset with different correlation
    } else {
        n_subset <- round(N * fraction)
        n_nonsubset <- N - n_subset
        full_sample <- getBiCop(N = n_nonsubset, rho = rho)
        subset <- getBiCop(N = n_subset, rho = rho_sub)
        df <- rbind(subset, full_sample)
        return(df)
    }
}

#' Simulate  non-normally distributed data with specified correlation, mean and sd
#'
#' Draws random numbers from distributions with specified skew and kurtosis based on
#' PoisNonNor::RNG_P_NN. Checks the obtained correlation and issues a warning if it
#' deviates more than rhotol from rho.
#'
#' @param N Size of generated data (rows)
#' @param rho Target correlation coefficient
#' @param rho_sub True correlation coefficient in the subset. Default is no
#' seperate subset.
#' @param fraction The size of the subset as a percentage of the full sample.
#' The subset is then contained in the first round(N * fraction) rows of the
#' returned data frame.
#' @param m1 Mean of first variable
#' @param m2 Mean of second variable
#' @param s1 Standard deviation of first variable
#' @param s2 Standard deviation of second variable
#' @param skew1 Skew of first variable
#' @param skew2 Skew of second variable
#' @param kurtosis1 Kurtosis of first variable
#' @param kurtosis2 Kurtosis of second variable
getBiCop_nonnor <- function(N, rho, rho_sub = NA, fraction = NA,
                            mu1 = 0, mu2 = 0, sd1 = 1, sd2 = 1,
                            skew1 = -1, skew2 = 1, kurtosis1 = 1, kurtosis2 = 1,
                            rhotol = 0.1) {

    if (sum(is.na(c(rho_sub, fraction))) == 1) stop("Specify both rho_sub and fraction")
    require(PoisNonNor)
    if (is.na(rho_sub)) {
        cmat <- matrix(c(1, rho, rho, 1), nrow = 2) # Korrelation
        rmat <- matrix(c(skew1, skew2, kurtosis1, kurtosis2), ncol = 2, byrow = F)
        mean_vec <- c(mu1, mu2)
        variance_vec <- c(sd1, sd2)
        dat <- RNG_P_NN(cmat = cmat, rmat = rmat, norow = N, mean.vec = mean_vec,
                        variance.vec = variance_vec)
        colnames(dat) <- c("x", "y")
        dat <- data.frame(dat)
        obtained_cor <- cor(dat)[2]
        if (abs(obtained_cor - rho) > rhotol) warning("Obtained correlation > rhotol")
        return(dat)
    } else {
        n_subset <- round(N * fraction)
        n_nonsubset <- N - n_subset
        full_sample <- getBiCop_nonnor(N = n_nonsubset, rho = rho,
                                       mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2,
                                       skew1 = skew1, skew2 = skew2,
                                       kurtosis1 = kurtosis1, kurtosis2 = kurtosis2)
        subset <- getBiCop_nonnor(N = n_subset, rho = rho_sub,
                                       mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2,
                                       skew1 = skew1, skew2 = skew2,
                                       kurtosis1 = kurtosis1, kurtosis2 = kurtosis2)
        df <- rbind(subset, full_sample)
        return(df)
    }
}

#' Create sequence of integers with increasing step size
#'
#' The sequence starts at initialSamplesize and every next element
#' is (at least, due to rounding) increased by stepPerc percent. The last
#' element is always nmax.
#'
#' @param nmax largest and last value of sequence
#' @param stepPerc every element of the sequence is at least stepPerc percent larger
#' than the preceding one
#' @param initialSamplesize start of the sequence
#' @param manValues manually given values that are added to the sequence even if they
#' do not fit into the step size that is defined by stepPerc
makeSteps <- function(nmax, stepPerc = 0.01, initialSamplesize = 10,
                      manValues = NULL) {
    end <- nmax
    steps <- round(initialSamplesize)
    while(max(steps) < end) {
        steps <- c(steps, ceiling((1 + stepPerc) * steps[length(steps)]))
    }
    steps[length(steps)] <- end
    if (length(manValues) > 0) {
        steps <- sort(unique(c(steps, manValues)))
    }
    return(steps)
}

#' Simulate repetitions of drawing a subset from data of a specified type,
#'  calculate the (subset) correlation and compute coverage probabilities and
#'  type 1 errors for specified inference type
#'
#' @param N (integer) Size of the full sample
#' @param fraction Size of the subset as a percentage of the full set
#' @param distribution One of "normal" or "nonnormal". The type of distribution the
#' variables are simulated from. In case of "nonnormal" the skewness and kurtosis
#' can deviate from a normal distribution
#' @param repetitions The simulation will be repeated this many times
#' @param skew1 Skewness of the distribution of x
#' @param skew2 Skewness of the distribution of y
#' @param kurtosis1 Kurtosis of the distribution of x
#' @param kurtosis2 Kurtosis of the distribution of y
run_sim_fisher <- function(N, fraction, rho = 0, rho_sub = NA, distribution = "normal",
                          kurtosis1 = 1, kurtosis2 = 1, skew1 = 1, skew2 = 1,
                          repetitions = 100, progress_bar = T, summarize = F) {
    stopifnot(distribution %in% c("normal", "nonnormal"))
    stopifnot(repetitions > 0 & length(repetitions) == 1)
    repetitions = round(repetitions)
    stopifnot(fraction <= 1 & fraction > 0)
    stopifnot(rho <= 1 & rho > -1)
    require(foreach)
    # exportfuncs <- c("sim_subcor", "getBiCop", "getBiCop_nonnor", "rtoz", "Sd_z_finite",
    #                  "ztor", "Sd_z", "E_z", "ci_cor_fisher", "pval_cor_fisher")
    if (progress_bar) {
        require(tcltk)
        pb <- txtProgressBar(max=repetitions, style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress=progress)
    } else opts = NULL
    stats <- foreach(i = 1:repetitions,
                     # .export = exportfuncs,
                     .options.redis = list(chunkSize = 5, ftinterval = 100),
                     .options.snow = opts, .packages = c("e1071", "subsetcor")) %dorng% {
        if (is.na(rho_sub)) {
            if (distribution == "normal") {
                dat <- getBiCop(N = N, rho = rho)
            } else if (distribution == "nonnormal") {
                dat <- getBiCop_nonnor(N = N, rho = rho, mu1 = 0, mu2 = 0,
                                       sd1 = 1, sd2 = 1,
                                       kurtosis1 = kurtosis1, kurtosis2 = kurtosis2,
                                       skew1 = skew1, skew2 = skew2)
            }
        } else {
            if (distribution == "normal") {
                dat <- getBiCop(N = N, rho = rho, rho_sub = rho_sub, fraction = fraction)
            } else if (distribution == "nonnormal") {
                dat <- getBiCop_nonnor(N = N, rho = rho, mu1 = 0, mu2 = 0,
                                       sd1 = 1, sd2 = 1,
                                       kurtosis1 = kurtosis1, kurtosis2 = kurtosis2,
                                       skew1 = skew1, skew2 = skew2)
            }
        }
        # Als Test (für den Output) skew und kurtosis einer Variablen
        tempkurtosis <- kurtosis(dat$x)
        tempskew <- skewness(dat$x)
        #
        n_full <- nrow(dat)
        n_subset <- round(n_full * fraction)
        if (is.na(rho_sub)) {
            correlations <- sim_subcor(dat = dat, fraction = fraction)
        } else {
            correlations <- list()
            correlations$cor_full <- cor(dat[(n_subset + 1):n_full, ])[2, 1]
            correlations$cor_subset <- cor(dat[1:n_subset, ])[2, 1]
            correlations$inSubset <- 1:n_subset
        }
        actual_cor <- correlations$cor_full
        ci95 <- ci_cor_fisher(alpha = 0.05, rho = actual_cor, n = n_subset, Npop = n_full)
        ci95_upper <- ci95["cor_ci_upper"]
        ci95_lower <- ci95["cor_ci_lower"]
        ci99 <- ci_cor_fisher(alpha = 0.01, rho = actual_cor, n = n_subset, Npop = n_full)
        ci99_upper <- ci99["cor_ci_upper"]
        ci99_lower <- ci99["cor_ci_lower"]
        pgreater <- pval_cor_fisher(r = correlations$cor_subset, rho = correlations$cor_full,
                                    n = n_subset,
                                    alternative = "greater", Npop = n_full)
        pless <- pval_cor_fisher(r = correlations$cor_subset, rho = correlations$cor_full,
                                 n = n_subset,
                                 alternative = "less", Npop = n_full)
        ptwosided <- pval_cor_fisher(r = correlations$cor_subset, rho = correlations$cor_full,
                                     n = n_subset,
                                     alternative = "two.sided", Npop = n_full)
        res <- c(n_full = n_full,
                 n_subset = n_subset,
                 rho = rho, rho_sub = rho_sub,
                 cor_subset = correlations$cor_subset,
                 cor_full = correlations$cor_full,
                 fraction = fraction,
                 ci95 = ci95_upper, ci95 = ci95_lower,
                 ci99 = ci99_upper, ci99 = ci99_lower,
                 pless = pless, pgreater = pgreater, ptwosided = ptwosided,
                 skew = tempskew, kurtosis = tempkurtosis)
    }
    if (progress_bar) close(pb)
    stats <- do.call(rbind, stats)
    stats <- data.frame(stats)
    if (summarize) {
        return(data.frame(
            inCI95 = stats$cor_subset >= stats$ci95.cor_ci_lower &
                stats$cor_subset <= stats$ci95.cor_ci_upper,
            inCI99 = stats$cor_subset >= stats$ci99.cor_ci_lower &
                stats$cor_subset <= stats$ci99.cor_ci_upper,
            pgreater_frac_significant = mean(stats$pgreater <= 0.05),
            pless_frac_significant = mean(stats$pless <= 0.05),
            ptwosided_frac_significant = mean(stats$ptwosided <= 0.05),
            skew = stats$skew, kurtosis = stats$kurtosis,
            cor_full = stats$cor_full, cor_subset = stats$cor_subset,
            n_full = stats$n_full, n_subset = stats$n_subset,
            rho = stats$rho, rho_sub = rho_sub, fraction = stats$fraction
        ))
    } else {
        return(data.frame(
            CI95_lower = stats$ci95.cor_ci_lower,
            CI95_upper = stats$ci95.cor_ci_upper,
            CI99_lower = stats$ci99.cor_ci_lower,
            CI99_upper = stats$ci99.cor_ci_upper,
            pgreater = stats$pgreater,
            pless = stats$pless,
            ptwosided = stats$ptwosided,
            skew = stats$skew, kurtosis = stats$kurtosis,
            cor_full = stats$cor_full, cor_subset = stats$cor_subset,
            n_full = stats$n_full, n_subset = stats$n_subset,
            rho = stats$rho, rho_sub = stats$rho_sub, fraction = stats$fraction
        ))
    }
}


run_sim_permut <- function(N, fraction, rho = 0, rho_sub = NA, distribution = "normal",
                           permuttype = "simple", n_permut = 1000,
                           kurtosis1 = 1, kurtosis2 = 1, skew1 = 1, skew2 = 1,
                           repetitions = 100, progress_bar = T, summarize = F) {
    stopifnot(distribution %in% c("normal", "nonnormal"))
    stopifnot(repetitions > 0 & length(repetitions) == 1)
    repetitions = round(repetitions)
    stopifnot(fraction <= 1 & fraction > 0)
    stopifnot(rho <= 1 & rho > -1)
    require(foreach)
    # exportfuncs <- c("sim_subcor", "getBiCop", "getBiCop_nonnor", "rtoz", "Sd_z_finite",
    #                  "ztor", "Sd_z", "E_z", "ci_cor_fisher", "pval_cor_fisher",
    #                  "inference_cor_perm", "permute_ecdf", "test_cor_perm_student")
    if (progress_bar) {
        require(tcltk)
        pb <- txtProgressBar(max=repetitions, style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress=progress)
    } else opts = NULL
    stats <- foreach(i = 1:repetitions,
                     # .export = exportfuncs,
                     .options.redis = list(chunkSize = 5, ftinterval = 100),
                     .options.snow = opts, .packages = c("e1071", "subsetcor")) %dorng% {
        if (is.na(rho_sub)) {
            if (distribution == "normal") {
                dat <- getBiCop(N = N, rho = rho)
            } else if (distribution == "nonnormal") {
                dat <- getBiCop_nonnor(N = N, rho = rho, mu1 = 0, mu2 = 0,
                                       sd1 = 1, sd2 = 1,
                                       kurtosis1 = kurtosis1, kurtosis2 = kurtosis2,
                                       skew1 = skew1, skew2 = skew2)
            }
        } else {
            if (distribution == "normal") {
                dat <- getBiCop(N = N, rho = rho, rho_sub = rho_sub, fraction = fraction)
            } else if (distribution == "nonnormal") {
                dat <- getBiCop_nonnor(N = N, rho = rho, mu1 = 0, mu2 = 0,
                                       sd1 = 1, sd2 = 1,
                                       kurtosis1 = kurtosis1, kurtosis2 = kurtosis2,
                                       skew1 = skew1, skew2 = skew2)
            }
        }
        # Als Test (für den Output) skew und kurtosis einer Variablen
        tempkurtosis <- kurtosis(dat$x)
        tempskew <- skewness(dat$x)
        #
        n_full <- nrow(dat)
        n_subset <- round(n_full * fraction)
        if (is.na(rho_sub)) {
            correlations <- sim_subcor(dat = dat, fraction = fraction)
        } else {
            correlations <- list()
            correlations$cor_full <- cor(dat[(n_subset + 1):n_full, ])[2, 1]
            correlations$cor_subset <- cor(dat[1:n_subset, ])[2, 1]
            correlations$inSubset <- 1:n_subset
        }
        actual_cor <- correlations$cor_full
        pvals_cis_perm <- inference_cor_perm(x = dat$x, y = dat$y, r = correlations$cor_subset,
                                             fraction = fraction, type = permuttype,
                                             n_permut = n_permut,
                                             inSubset = correlations$inSubset)
        res <- c(n_full = n_full,
                 n_subset = n_subset,
                 rho = rho, rho_sub = rho_sub,
                 cor_subset = correlations$cor_subset,
                 cor_full = correlations$cor_full,
                 fraction = fraction,
                 ci95_upper = pvals_cis_perm$ci95_upper,
                 ci95_lower = pvals_cis_perm$ci95_lower,
                 ci99_upper = pvals_cis_perm$ci99_upper,
                 ci99_lower = pvals_cis_perm$ci99_lower,
                 p_less = pvals_cis_perm$p_less,
                 p_greater = pvals_cis_perm$p_greater,
                 p_twosided = pvals_cis_perm$p_twosided,
                 skew = tempskew, kurtosis = tempkurtosis)
    }
    if (progress_bar) close(pb)
    stats <- do.call(rbind, stats)
    stats <- data.frame(stats)
    if (summarize) {
        # Summarize results
        cor_diffs <- stats$cor_subset - stats$cor_full
        inCI95 = cor_diffs >= stats$ci95_lower &
            cor_diffs <= stats$ci95_upper
        inCI99 = cor_diffs >= stats$ci99_lower &
            cor_diffs <= stats$ci99_upper
        return(data.frame(inCI95 = mean(inCI95),
                          inCI99 = mean(inCI99),
                          pgreater_frac_significant = mean(stats$p_greater <= 0.05),
                          pless_frac_significant = mean(stats$p_less <= 0.05),
                          ptwosided_frac_significant = mean(stats$p_twosided <= 0.05),
                          skew = stats$skew, kurtosis = stats$kurtosis,
                          cor_full = stats$cor_full, cor_subset = stats$cor_subset,
                          n_full = stats$n_full, n_subset = stats$n_subset,
                          rho = stats$rho, rho_sub = stats$rho_sub,
                          fraction = stats$fraction))
    } else {
        return(data.frame(
            CI95_lower = stats$ci95_lower,
            CI95_upper = stats$ci95_upper,
            CI99_lower = stats$ci99_lower,
            CI99_upper = stats$ci99_upper,
            pgreater = stats$p_greater,
            pless = stats$p_less,
            ptwosided = stats$p_twosided,
            skew = stats$skew, kurtosis = stats$kurtosis,
            cor_full = stats$cor_full, cor_subset = stats$cor_subset,
            n_full = stats$n_full, n_subset = stats$n_subset,
            rho = stats$rho, rho_sub = stats$rho_sub, fraction = stats$fraction
        ))
    }
}
