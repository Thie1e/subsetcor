#' Density function of sample correlation coefficient based on the closed form
#' solution in Provost (2015)
#'
#' Uses the hypergeometric function 2F1 from the gsl package.
#' @param r Value of the correlation coefficient
#' @param n Sample size
#' @param rho True correlation in the population
dcorr <- function(r, n, rho) {
    require(gsl)
    t1 <- (n-2) * beta((n-1)/2, n/2)**2 * (1-rho**2) ** ((n-1)/2)
    t2 <- 2**(n+1) * beta(n-1, n) * pi
    t3 <- (1-r**2) ** ((n/2) - 2)
    t4 <- hyperg_2F1(n-1, n-1, n-1/2, (1+rho*r)/2)
    t1 / t2 * t3 * t4
}

#' Distribution function of sample correlation coefficient based on the closed form
#' solution in Provost (2015)
#'
#' Uses the hypergeometric function 2F1 from the gsl package.
#' Calculates the distribution directly as the integral of dcorr().
#' @param q Value of the correlation coefficient
#' @param n Sample size
#' @param rho True correlation in the population
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
pcorr <- function(q, n, rho, lower.tail = TRUE) {
    if (lower.tail) {
        lim_lower <- -1
        lim_upper <- q
    } else if (!lower.tail) {
        lim_lower <- q
        lim_upper <- 1
    }
    dcorr2 <- function(r) dcorr(r, n = n, rho = rho)
    integrate(dcorr2, lim_lower, lim_upper)
}

#' Moment functions for the sample correlation coefficient based on the solutions
#' of Anderson (via Provost, 2015)
#'
#' Stops the calculation of the sum term as soon as an infinite value is
#' returned, which may happen very soon or immediately for large n. In that
#' case the results are not reliable.
#'
#' @param rho True correlation coefficient in the population
#' @param n Sample size
#' @param k Moment number
#' @maxterm Maximum number of terms of the infinite sum that the function
#' will try to calculate
moments_R_odd <- function(rho, n, k, maxterms = 1000) {
    if (k %% 2 != 1) stop("k has to be odd")
    i <- 0
    s <- rep(NA, maxterms)
    flag = T
    while(flag == T & i < maxterms) {
        s1 <- suppressWarnings(((2*rho) ** (2 * i + 1)) / factorial(2 * i + 1))
        s2 <- suppressWarnings((gamma(3/2 + (k-1)/2 + i) * gamma(n/2 + i) ** 2) /
                                   gamma((n+1)/2 + (k-1)/2 + i))
        s1s2 <- s1 * s2
        i <- i + 1
        if (is.na(s1s2) | !is.finite(s1s2)) flag = F else s[i] <- s1s2
    }
    s <- na.omit(s)
    terms <- length(s)
    mult <- (1-rho**2) ** ((n-1)/2) / (sqrt(pi) * gamma((n-1)/2))
    moment <- mult * sum(s)
    return(list(value = moment, terms = terms))
}

#' Moment functions for the sample correlation coefficient based on the solutions
#' of Anderson (via Provost, 2015)
#'
#' Stops the calculation of the sum term as soon as an infinite value is
#' returned, which may happen very soon or immediately for large n. In that
#' case the results are not reliable.
#'
#' @param rho True correlation coefficient in the population
#' @param n Sample size
#' @param k Moment number
#' @maxterm Maximum number of terms of the infinite sum that the function
#' will try to calculate
moments_R_even <- function(rho, n, k, maxterms = 1000) {
    if (k %% 2 == 1) stop("k has to be even")
    i <- 0
    s <- rep(NA, maxterms)
    flag = T
    while(flag == T & i < maxterms) {
        s1 <- suppressWarnings(((2*rho) ** (2*i)) / factorial(2*i))
        s2 <- suppressWarnings((gamma(1/2 + k/2 + i) * gamma(n/2 - 1/2 + i) ** 2) /
                                   gamma((n-1)/2 + k/2 + i))
        s1s2 <- s1 * s2
        i <- i + 1
        if (is.na(s1s2) | !is.finite(s1s2)) flag = F else s[i] <- s1s2
    }
    s <- na.omit(s)
    terms <- length(s)
    mult <- (1-rho**2) ** ((n-1)/2) / (sqrt(pi) * gamma((n-1)/2))
    moment <- mult * sum(s)
    return(list(value = moment, terms = terms))
}

#' Expected value of sample correlation coefficient
#'
#' Calculated as the integral of r times the distribution function dcorr in [-1, 1].
#' @param n Sample size
#' @param rho True correlation in the population
E_r_integral <- function(n, rho) {
    e_dcorr <- function(r){r*dcorr(r, n = n, rho = rho)}
    integrate(e_dcorr,-1,1)
}

#' Variance of sample correlation coefficient
#'
#' Calculated as the integral of (r - E(r)) times the distribution function dcorr in [-1, 1].
#' Thus based on the closed form solution in Provost (2015).
#' @param n Sample size
#' @param rho True correlation in the population
var_r_integral <- function(n, rho) {
    mu <- E_r(n, rho)$value
    var_dcorr <- function(r, n. = n, rho. = rho) {
        (r - mu)**2 * dcorr(r, n., rho.)
    }
    integrate(var_dcorr,-1,1)
}

#' p-value of sample correlation coefficient
#'
#' Based on the closed form solution in Provost (2015). May fail for large n (>200).
#' Calculates probability of observing a sample correlation coefficient r given
#' a sample size n and the population correlation coefficient rho. Not suitable
#' for assessing the probability of subset correlation coefficients, i.e. correlation
#' coefficients from finite populations.
#'
#' @param r sample correlation coefficient
#' @param rho population correlation coefficient
#' @param n sample size
#' @param alternative a character string specifying the alternative hypothesis,
#'  must be one of "greater" (default) or "less"
pval_cor <- function(r, rho, n, alternative = "greater") {
    if (!(alternative %in% c("greater", "less"))) stop("alternative != greater or less")
    if (alternative == "greater") {
        return(pcorr(q = r, n = n, rho = rho, lower.tail = F)$value)
    }
    if (alternative == "less") {
        return(pcorr(q = r, n = n, rho = rho, lower.tail = T)$value)
    }
}


#' Confidence interval of sample correlation coefficient
#'
#' Based on the closed form solution in Provost (2015). May fail for large n (>200).
#' Not suitable for constructing confidence intervals for subset correlation coefficients, i.e. correlation
#' coefficients from finite populations.
#'
#' @param alpha confidence level
#' @param rho population correlation coefficient
#' @param n sample size
ci_cor <- function(alpha, rho, n) {
    dist_fun_lower <- function(x) abs(pcorr(q = x, n = n, rho = rho, lower.tail = T)$value - (alpha / 2))
    dist_fun_upper <- function(x) abs(pcorr(q = x, n = n, rho = rho, lower.tail = F)$value - (alpha / 2))
    ci_upper <- optimize(dist_fun_upper, interval = c(-1, 1))
    ci_lower <- optimize(dist_fun_lower, interval = c(-1, 1))
    return(c(cor_ci_lower = ci_lower$minimum, cor_ci_upper = ci_upper$minimum))
}

#' Functions for Fisher's z-transformation
#'
#' Convert a sample correlation coefficient (r) to Fisher's z or the other way
#' around, compute the expected value of z (E_z), compute the standard error of z
#' (se_z), and the standard error of z with a correction for finite samples (se_z_finite) with
#' Npop being the size of the full sample.
#'
#' @param r Sample correlation coefficient
#' @param rho Population correlation coefficient
#' @param n Sample size
#' @param Npop Population size in case of finite populations, e.g. when a subset was
#' extracted from a sample
#' @seealso \code{\link{rtoz}}, \code{\link{ztor}}, \code{\link{E_z}}, \code{\link{se_z}}, \code{\link{se_z_finite}}
rtoz <- function(r) 0.5 * log((1 + r) / (1 - r))

#' Functions for Fisher's z-transformation
#'
#' Convert a sample correlation coefficient (r) to Fisher's z or the other way
#' around, compute the expected value of z (E_z), compute the standard error of z
#' (se_z), and the standard error of z with a correction for finite samples (se_z_finite) with
#' Npop being the size of the full sample.
#'
#' @param r Sample correlation coefficient
#' @param rho Population correlation coefficient
#' @param n Sample size
#' @param Npop Population size in case of finite populations, e.g. when a subset was
#' extracted from a sample
#' @seealso \code{\link{rtoz}}, \code{\link{ztor}}, \code{\link{E_z}}, \code{\link{se_z}}, \code{\link{se_z_finite}}
ztor <- function(z) (exp(2 * z) - 1) / (exp(2 * z) + 1)

#' Functions for Fisher's z-transformation
#'
#' Convert a sample correlation coefficient (r) to Fisher's z or the other way
#' around, compute the expected value of z (E_z), compute the standard error of z
#' (se_z), and the standard error of z with a correction for finite samples (se_z_finite) with
#' Npop being the size of the full sample.
#'
#' @param r Sample correlation coefficient
#' @param rho Population correlation coefficient
#' @param n Sample size
#' @param Npop Population size in case of finite populations, e.g. when a subset was
#' extracted from a sample
#' @seealso \code{\link{rtoz}}, \code{\link{ztor}}, \code{\link{E_z}}, \code{\link{se_z}}, \code{\link{se_z_finite}}
E_z <- function(rho) 0.5 * log((1 + rho) / (1 - rho))

#' Functions for Fisher's z-transformation
#'
#' Convert a sample correlation coefficient (r) to Fisher's z or the other way
#' around, compute the expected value of z (E_z), compute the standard error of z
#' (se_z), and the standard error of z with a correction for finite samples (se_z_finite) with
#' Npop being the size of the full sample.
#'
#' @param r Sample correlation coefficient
#' @param rho Population correlation coefficient
#' @param n Sample size
#' @param Npop Population size in case of finite populations, e.g. when a subset was
#' extracted from a sample
#' @seealso \code{\link{rtoz}}, \code{\link{ztor}}, \code{\link{E_z}}, \code{\link{se_z}}, \code{\link{se_z_finite}}
se_z <- function(n) 1 / (sqrt(n - 3))

#' Functions for Fisher's z-transformation
#'
#' Convert a sample correlation coefficient (r) to Fisher's z or the other way
#' around, compute the expected value of z (E_z), compute the standard error of z
#' (se_z), and the standard error of z with a correction for finite samples (se_z_finite) with
#' Npop being the size of the full sample.
#'
#' @param r Sample correlation coefficient
#' @param rho Population correlation coefficient
#' @param n Sample size
#' @param Npop Population size in case of finite populations, e.g. when a subset was
#' extracted from a sample
#' @seealso \code{\link{rtoz}}, \code{\link{ztor}}, \code{\link{E_z}}, \code{\link{se_z}}, \code{\link{se_z_finite}}
se_z_finite <- function(Npop, n) se_z(n) * sqrt((Npop - n) / (Npop - 1))

#' Confidence interval of sample correlation coefficient (also in finite populations)
#'
#' Based on Fisher's z-transformation and a correction of the standard deviation for
#' finite samples: (Npop - n) / (Npop - 1)
#'
#' @param rho population correlation coefficient
#' @param n sample size
#' @param Npop Population size
#' @param alpha confidence level
ci_cor_fisher <- function(rho, n, Npop, alpha = 0.05) {
    z <- rtoz(rho)
    if (Npop == Inf) {
        sd_z <- se_z(n = n)
    } else {
        sd_z <- se_z_finite(Npop = Npop, n = n)
    }
    # ci_upper <- z + qt(1 - alpha / 2, df = n - 1) * sd_z
    ci_upper <- z + qnorm(1 - alpha / 2) * sd_z
    ci_upper <- ztor(ci_upper)
    # ci_lower <- z + qt(alpha / 2, df = n - 1) * sd_z
    ci_lower <- z + qnorm(alpha / 2) * sd_z
    ci_lower <- ztor(ci_lower)
    return(c(cor_ci_lower = ci_lower, cor_ci_upper = ci_upper))
}

#' p-value of sample correlation coefficient (also in finite populations)
#'
#' Based on Fisher's z-transformation and a correction of the standard deviation for
#' finite samples: (Npop - n) / (Npop - 1)
#' Calculates probability of observing a sample correlation coefficient r given
#' a sample size n, a population size Npop and the population correlation coefficient rho.
#'
#' @param r sample correlation coefficient
#' @param rho population correlation coefficient
#' @param n sample size
#' @param Npop population size
#' @param alternative a character string specifying the alternative hypothesis,
#'  must be one of "two.sided" (default), "greater" or "less"
pval_cor_fisher <- function(r, rho, n, Npop, alternative = "two.sided") {
    if (!(alternative %in% c("greater", "less", "two.sided"))) stop("alternative != greater, less or two.sided")
    r_z <- rtoz(r)
    sd_z <- se_z_finite(Npop, n)
    if (alternative == "greater") {
        # p <- pnorm(q = r_z, mean = E_z(rho), sd = sd_z, lower.tail = F)
        # Das gleiche schöner:
        p <- pnorm(q = (r_z - E_z(rho)) / sd_z, lower.tail = F)
        return(p)
    } else if (alternative == "less") {
        p <- pnorm(q = r_z, mean = E_z(rho), sd = sd_z, lower.tail = T)
        return(p)
    } else if (alternative == "two.sided") {
        pgreater <- pnorm(q = r_z, mean = E_z(rho), sd = sd_z, lower.tail = F)
        pless <- pnorm(q = r_z, mean = E_z(rho), sd = sd_z, lower.tail = T)
        return(2 * min(pgreater, pless))
    }
}

#' Draw a permutation from the full sample with size n_subset and save empirical distribution of cor_subset
#'
#' Confidence intervals can be derived using quantile()
#' p-values of observed cor_subset can be derived using ecdf()
#' Either fraction or n_subset should be specified.
#' @param x Variable x, same length as y
#' @param y Variable y, same length as x
#' @param fraction Size of the subset in relation to the full set in percent
#' @param n_subset Size of the subset
#' @param n_permut number of permutations
permute_ecdf <- function(x, y, fraction = NULL, n_subset = NULL, n_permut = 1000) {
    if (!is.null(fraction) & !is.null(n_subset)) stop("Only enter fraction or n_subset")
    if (is.null(n_subset)) n_subset <- round(length(x) * fraction)
    stopifnot(length(x) == length(y))
    cor_permute <- sapply(1:n_permut, function(r) {
        rand_subset <- sample(1:length(x), size = n_subset, replace = F)
        cor_subset <- cor(x[rand_subset], y[rand_subset])
        cor_full <- cor(x, y)
        return(cor_subset - cor_full)
    })
    return(ecdf(cor_permute))
}

#' Calculate confidence intervals and p-value for a given subsample correlation
#'
#' Given a data set (x and y) and a vector that contains the row numbers of the
#' observations that belong to the subset (inSubset) the functions runs n_permut permutations
#' for determining 95% and 99% confidence intervals and the significance of the
#' given subset correlation coefficient (r)
#'
#' @param x x variable (numeric vector)
#' @param y y variable (numeric vector)
#' @param type "simple" or "studentized" permutation test
#' @param n_permut number of permutations for the permutation test
#' @param inSubset integer vector of the row numbers of all observations that belong
#' to the subset
inference_cor_perm <- function(x, y, type = "simple", n_permut = 1000, inSubset = NULL) {
    stopifnot(all(table(inSubset)) == 1)
    stopifnot(type %in% c("simple", "studentized"))
    stopifnot(length(x) == length(y))
    full_data <- data.frame(x = x, y = y)
    if (type == "simple") {
        fraction <- length(inSubset) / length(x)
        permutation_ecdf <- permute_ecdf(x = x, y = y, fraction = fraction, n_permut = 1000)
        r <- cor(x[inSubset], y[inSubset])
        cor_diff = r - cor(full_data$x, full_data$y)
        p_greater = 1 - permutation_ecdf(cor_diff)
        p_less = permutation_ecdf(cor_diff)
        p_twosided = 1 - permutation_ecdf(abs(cor_diff)) + permutation_ecdf(-abs(cor_diff))
        # confidence intervals
        ci99_lower = as.numeric(quantile(permutation_ecdf, 0.005))
        ci99_upper = as.numeric(quantile(permutation_ecdf, 0.995))
        ci95_lower = as.numeric(quantile(permutation_ecdf, 0.025))
        ci95_upper = as.numeric(quantile(permutation_ecdf, 0.975))
        res <- list(p_twosided = p_twosided, p_greater = p_greater,
                    p_less = p_less, ci99_lower = ci99_lower,
                    ci99_upper = ci99_upper, ci95_lower = ci95_lower,
                    ci95_upper = ci95_upper)
    }
    if (type == "studentized") {
        stopifnot(!is.null(inSubset))
        subset_data <- full_data[inSubset, ]
        res <- test_cor_perm_student(full_data = full_data, subset_data = subset_data, n_permut = n_permut)
    }
    # return(list(p_twosided = p_twosided, p_less = p_less, p_greater = p_greater,
    #             ci95lower = ci95lower, ci95upper = ci95upper,
    #             ci99lower = ci99lower, ci99upper = ci99upper))
    return(res)
}



test_cor_perm_student <- function(full_data, subset_data, n_permut = 1000) {
    #   Function for standardization
    stand <- function(x){
        mean.x <- mean(x);
        res <- (x - mean.x)/sqrt(mean((x-mean.x)^2));
        return(res);
    }

    # #   Sample sizes
    n1 <- nrow(full_data)
    n2 <- nrow(subset_data)
    # n <- n1 + n2 ## Wenn z <- rbind(full_data, subset_data)
    # n <- n1

    # #     Standardization
    full_data <- apply(full_data, 2, stand)
    subset_data <- apply(subset_data, 2, stand)

    # # Sample correlation coefficients (data is already standardized)
    r1 <- mean(full_data[, 1] * full_data[, 2])
    r2 <- mean(subset_data[, 1] * subset_data[, 2])

    # #     Pooled sample correlation coefficient
    # r <- n1/n * r1 + n2/n * r2
    # r <- n1/(n1+n2) * r1 + n2/(n1+n2) * r2
    r <- r1

    # #     Asymptotically linearized observations (Omelka/Pauly Gleichung 9)
    z1 <- full_data[,1]*full_data[,2] - r/2*(full_data[,1]^2 + full_data[,2]^2)
    # z2 <- subset_data[,1]*subset_data[,2] - r/2*(subset_data[,1]^2 + subset_data[,2]^2)

    # #  Estimate of the asymptotic variance
    # Vn <- var(z1)/n1 + var(z2)/n2
    Vn <- var(z1)/n1

    # #     joint sample
    # z <- rbind(full_data, subset_data) ## muss das im Fall von full sample vs. subset sein?
    z <- full_data

    # # Test statistic of the asymptotic test
    stat.as <- abs((r1 - r2)/sqrt(Vn))
    pval.as <- 2*(1 - pnorm(stat.as))

    # #
    # stat.perm <- numeric(n.perm); ## ???
    stat.as <- ((r1 - r2)/sqrt(Vn)) # ohne abs()

    # #     Resampling
    stat.perm <- sapply(1:n_permut, function(r){
        # z <- z[sample(1:n),];
        # x <- z[1:n1,];
        # y <- z[(n1+1):n,];

        # Generate randomly permuted subset
        y <- z
        inSubset <- sample(1:n1, size = n2)
        x <- z[inSubset, ]

        # Next try (zu konservativ, mit z <- rbind(full, subset))
        # z <- z[sample(1:n),]; # shuffling joint sample
        # x <- z[1:n1,]; # new full_data
        # y <- z[(n1+1):n,]; # new subset_data

        r1 <- cor(x)[1,2]
        r2 <- cor(y)[1,2]
        # r <- n1/n * r1 + n2/n * r2
        r <- n1/(n1+n2) * r1 + n2/(n1+n2) * r2

        x <- apply(x, 2, stand)
        y <- apply(y, 2, stand)

        z1 <- x[,1]*x[,2] - r/2*(x[,1]^2 + x[,2]^2)
        # z2 <- y[,1]*y[,2] - r/2*(y[,1]^2 + y[,2]^2)

        # Vn <- var(z1)/n1 + var(z2)/n2
        Vn <- var(z1)/n1

        (r1 - r2)/sqrt(Vn)
    })
    # stat.perm <- abs(stat.perm);
    # #   p-value of the permuation test
    # pval.perm <- (sum(stat.perm >= stat.as)+1)/(n_permut + 1); ### 2-tailed, da stat.as und stat.perm Beträge sind
    ### müsste sich so wie im einfachen Test auf > < erweitern lassen
    # #   Printing results
    # res <- c(stat.as, pval.as, pval.perm);
    # names(res) <- c("test statistic tildeT_n", "p-value of the as. test", "p-value based on perm. test");
    # return(res);

    cor_diff = stat.as
    stat.perm_ecdf <- ecdf(stat.perm)
    p_greater = 1 - stat.perm_ecdf(cor_diff)
    p_less = stat.perm_ecdf(cor_diff)
    p_twosided = 1 - stat.perm_ecdf(abs(cor_diff)) + stat.perm_ecdf(-abs(cor_diff))
    p_as_twosided  = pval.as
    # Mit *sqrt(Vn) um die Differenz der Korrelationen zurückzutransformieren
    ci99_lower = quantile(stat.perm_ecdf, 0.005) * sqrt(Vn)
    ci99_upper = quantile(stat.perm_ecdf, 0.995) * sqrt(Vn)
    ci95_lower = quantile(stat.perm_ecdf, 0.025) * sqrt(Vn)
    ci95_upper = quantile(stat.perm_ecdf, 0.975) * sqrt(Vn)
    res <- list(p_greater = p_greater, p_less = p_less, p_twosided = p_twosided,
                p_as_twosided = p_as_twosided,
                p_asymptotic = pval.as, ci99_lower = ci99_lower, ci99_upper = ci99_upper,
                ci95_lower = ci95_lower, ci95_upper = ci95_upper)
    return(res)
}



# Als weiterer Vergleich zum studentisierten und simplen Permutationstest.
# Diese Funktion permutiert das standardisierte Sample mit
# an * (r2 - r1) als Teststatistik anstelle von (r2 - r1).
test_cor_perm_sakaori <- function(full_data, subset_data, n_permut = 1000) {
    #   Function for standardization
    stand <- function(x){
        mean.x <- mean(x);
        res <- (x - mean.x)/sqrt(mean((x-mean.x)^2));
        return(res);
    }

    # #   Sample sizes
    n1 <- nrow(full_data)
    n2 <- nrow(subset_data)

    an <- sqrt((n1 * n2) / (n1 + n2))

    # #     Standardization
    full_data <- apply(full_data, 2, stand)
    subset_data <- apply(subset_data, 2, stand)

    # #     Resampling
    stat.perm <- sapply(1:n_permut, function(r){
        # Generate randomly permuted subset
        inSubset <- sample(1:n1, size = n2)
        temp_subset <- full_data[inSubset, ]
        r1 <- cor(full_data)[1,2]
        r2 <- cor(temp_subset)[1,2]
        an * (r2 - r1)
    })

    # # Sample correlation coefficients (data is already standardized)
    r1 <- mean(full_data[, 1] * full_data[, 2])
    r2 <- mean(subset_data[, 1] * subset_data[, 2])

    cor_diff = an * (r2 - r1)
    stat.perm_ecdf <- ecdf(stat.perm)
    p_greater = 1 - stat.perm_ecdf(cor_diff)
    p_less = stat.perm_ecdf(cor_diff)
    p_twosided = 1 - stat.perm_ecdf(abs(cor_diff)) + stat.perm_ecdf(-abs(cor_diff))
    ci99_lower = quantile(stat.perm_ecdf, 0.005) / an
    ci99_upper = quantile(stat.perm_ecdf, 0.995) / an
    ci95_lower = quantile(stat.perm_ecdf, 0.025) / an
    ci95_upper = quantile(stat.perm_ecdf, 0.975) / an
    res <- list(p_greater = p_greater, p_less = p_less, p_twosided = p_twosided,
                ci99_lower = ci99_lower, ci99_upper = ci99_upper,
                ci95_lower = ci95_lower, ci95_upper = ci95_upper)
    return(res)
}
