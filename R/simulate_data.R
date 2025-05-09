#' @title Set simulation parameters of a study group.
#'
#' @description This function provides input arguments for each study group needed to
#' simulate data with [simulate_data()]. [simulate_data()] generates data for a two-arms
#' clinical trial with longitudinal continuous outcomes and two intercurrent events (ICEs).
#' ICE1 may be thought of as a discontinuation from study treatment due to study drug or
#' condition related (SDCR) reasons. ICE2 may be thought of as discontinuation from study
#' treatment due to uninformative study drop-out, i.e. due to not study drug or
#' condition related (NSDRC) reasons and outcome data after ICE2 is always missing.
#'
#' @param mu Numeric vector describing the mean outcome trajectory at each visit (including
#' baseline) assuming no ICEs.
#' @param sigma Covariance matrix of the outcome trajectory assuming no ICEs.
#' @param n Number of subjects belonging to the group.
#' @param prob_ice1 Numeric vector that specifies the probability of experiencing ICE1
#' (discontinuation from study treatment due to SDCR reasons) after each visit for a subject
#' with observed outcome at that visit equal to the mean at baseline (`mu[1]`).
#' If a single numeric is provided, then the same probability is applied to each visit.
#' @param or_outcome_ice1 Numeric value that specifies the odds ratio of experiencing ICE1 after
#' each visit corresponding to a +1 higher value of the observed outcome at that visit.
#' @param prob_post_ice1_dropout Numeric value that specifies the probability of study
#' drop-out following ICE1. If a subject is simulated to drop-out after ICE1, all outcomes after
#' ICE1 are set to missing.
#' @param prob_ice2 Numeric that specifies an additional probability that a post-baseline
#' visit is affected by study drop-out. Outcome data at the subject's first simulated visit
#' affected by study drop-out and all subsequent visits are set to missing. This generates
#' a second intercurrent event ICE2, which may be thought as treatment discontinuation due to
#' NSDRC reasons with subsequent drop-out.
#' If for a subject, both ICE1 and ICE2 are simulated to occur,
#' then it is assumed that only the earlier of them counts.
#' In case both ICEs are simulated to occur at the same time, it is assumed that ICE1 counts.
#' This means that a single subject can experience either ICE1 or ICE2, but not both of them.
#' @param prob_miss Numeric value that specifies an additional probability for a given
#' post-baseline observation to be missing. This can be used to produce
#' "intermittent" missing values which are not associated with any ICE.
#'
#' @details For the details, please see [simulate_data()].
#'
#' @return A `simul_pars` object which is a named list containing the simulation parameters.
#'
#' @seealso [simulate_data()]
#'
#' @export
set_simul_pars <- function(
    mu, sigma, n,
    prob_ice1 = 0,
    or_outcome_ice1 = 1,
    prob_post_ice1_dropout = 0,
    prob_ice2 = 0,
    prob_miss = 0
) {

    x <- list(
        mu = mu,
        sigma = sigma,
        n = n,
        prob_ice1 = prob_ice1,
        or_outcome_ice1 = or_outcome_ice1,
        prob_post_ice1_dropout = prob_post_ice1_dropout,
        prob_ice2 = prob_ice2,
        prob_miss = prob_miss
    )

    class(x) <- "simul_pars"

    return(x)
}

#' @title Validate a `simul_pars` object
#'
#' @param x An `simul_pars` object as generated by [set_simul_pars()].
#' @param ... Not used.
#'
#' @export
validate.simul_pars <- function(x, ...) {

    expected_names <- c(
        "mu", "sigma", "n", "prob_ice1", "or_outcome_ice1",
        "prob_post_ice1_dropout", "prob_ice2", "prob_miss"
    )

    assert_that(
        length(x) == 8,
        all(names(x) %in% expected_names),
        msg = "`x` must be a named list of length 8"
    )
    assert_that(
        all(!sapply(x, function(y) any(is.na(y)))),
        msg = "some element of `x` contains missing values"
    )
    assert_that(
        is.vector(x$mu) &&
            is.matrix(x$sigma) &&
            all(length(x$mu) == dim(x$sigma)),
        msg = "`mu` must be a numeric vector with length equal to the dimension of `sigma`"
    )
    assert_that(
        is.numeric(x$n) && x$n > 0,
        msg = "`n` must be a positive number"
    )
    assert_that(
        is.numeric(x$prob_ice1) && all(x$prob_ice1 >= 0 & x$prob_ice1 <= 1),
        msg = "each element of `prob_ice1` must be a number in [0,1]"
    )
    assert_that(
        length(x$prob_ice1) == 1 || length(x$prob_ice1) == length(x$mu) - 1,
        msg = "`prob_ice1` must have either length 1 or equal to `length(mu) - 1`"
    )
    assert_that(
        is.numeric(x$or_outcome_ice1),
        x$or_outcome_ice1 > 0,
        msg = "`or_outcome_ice1` must be a positive number"
    )
    assert_that(
        is.numeric(x$prob_post_ice1_dropout),
        (x$prob_post_ice1_dropout >= 0 && x$prob_post_ice1_dropout <= 1),
        msg = "`prob_post_ice1_dropout` must be a number in [0,1]"
    )
    assert_that(
        is.numeric(x$prob_ice2),
        (x$prob_ice2 >= 0 && x$prob_ice2 <= 1),
        msg = "`prob_ice2` must be a number in [0,1]"
    )
    assert_that(
        is.numeric(x$prob_miss),
        (x$prob_miss >= 0 && x$prob_miss <= 1),
        msg = "`prob_miss` must be a number in [0,1]"
    )
}


#' @title Generate data
#'
#' @description Generate data for a two-arms clinical trial with longitudinal continuous
#' outcome and two intercurrent events (ICEs).
#' ICE1 may be thought of as a discontinuation from study treatment due to study drug or
#' condition related (SDCR) reasons.
#' ICE2 may be thought of as discontinuation from study treatment due to uninformative
#' study drop-out, i.e. due to not study drug or
#' condition related (NSDRC) reasons and outcome data after ICE2 is always missing.
#'
#' @param pars_c A `simul_pars` object as generated by [set_simul_pars()]. It specifies
#' the simulation parameters of the control arm.
#' @param pars_t A `simul_pars` object as generated by [set_simul_pars()]. It specifies
#' the simulation parameters of the treatment arm.
#' @param post_ice1_traj A string which specifies how observed outcomes occurring after
#' ICE1 are simulated.
#' Must target a function included in `strategies`. Possible choices are: Missing At
#' Random `"MAR"`, Jump to Reference `"JR"`,
#' Copy Reference `"CR"`, Copy Increments in Reference `"CIR"`, Last Mean Carried
#' Forward `"LMCF"`. User-defined strategies
#' could also be added. See [getStrategies()] for details.
#' @param strategies A named list of functions. Default equal to [getStrategies()].
#' See [getStrategies()] for details.
#'
#'
#' @details
#' The data generation works as follows:
#'
#' - Generate outcome data for all visits (including baseline) from a multivariate
#' normal distribution with parameters `pars_c$mu` and `pars_c$sigma`
#' for the control arm and parameters `pars_t$mu` and `pars_t$sigma` for the treatment
#' arm, respectively.
#' Note that for a randomized trial, outcomes have the same distribution at baseline
#' in both treatment groups, i.e. one should set
#' `pars_c$mu[1]=pars_t$mu[1]` and `pars_c$sigma[1,1]=pars_t$sigma[1,1]`.
#' - Simulate whether ICE1 (study treatment discontinuation due to SDCR reasons) occurs
#' after each visit according to parameters `pars_c$prob_ice1` and `pars_c$or_outcome_ice1`
#' for the control arm and `pars_t$prob_ice1` and `pars_t$or_outcome_ice1` for the
#' treatment arm, respectively.
#' - Simulate drop-out following ICE1 according to `pars_c$prob_post_ice1_dropout` and
#' `pars_t$prob_post_ice1_dropout`.
#' - Simulate an additional uninformative study drop-out with probabilities `pars_c$prob_ice2`
#' and `pars_t$prob_ice2` at each visit. This generates a second intercurrent event ICE2, which
#' may be thought as treatment discontinuation due to NSDRC reasons with subsequent drop-out.
#' The simulated time of drop-out is the subject's first visit which is affected by
#' drop-out and data from this visit and all subsequent visits are consequently set to missing.
#' If for a subject, both ICE1 and ICE2 are simulated to occur,
#' then it is assumed that only the earlier of them counts.
#' In case both ICEs are simulated to occur at the same time, it is assumed that ICE1 counts.
#' This means that a single subject can experience either ICE1 or ICE2, but not both of them.
#' - Adjust trajectories after ICE1 according to the given assumption expressed with
#' the `post_ice1_traj` argument. Note that only post-ICE1 outcomes in the intervention arm can be
#' adjusted. Post-ICE1 outcomes from the control arm  are not adjusted.
#' - Simulate additional intermittent missing outcome data as per arguments `pars_c$prob_miss`
#' and `pars_t$prob_miss`.
#'
#' The probability of the ICE after each visit is modeled according to the following
#' logistic regression model:
#' `~ 1 + I(visit == 0) + ... + I(visit == n_visits-1) + I((x-alpha))` where:
#' - `n_visits` is the number of visits (including baseline).
#' - `alpha` is the baseline outcome mean.
#' The term `I((x-alpha))` specifies the dependency of the probability of the ICE on
#' the current outcome value.
#' The corresponding regression coefficients of the logistic model are defined as follows:
#' The intercept is set to 0, the coefficients corresponding to discontinuation after
#' each visit for a subject with outcome equal to
#' the mean at baseline are set according to parameters `pars_c$prob_ice1` (`pars_t$prob_ice1`),
#' and the regression coefficient associated with the covariate `I((x-alpha))` is set
#' to `log(pars_c$or_outcome_ice1)` (`log(pars_t$or_outcome_ice1)`).
#'
#' Please note that the baseline outcome cannot be missing nor be affected by any ICEs.
#'
#' @returns A `data.frame` containing the simulated data. It includes the following variables:
#' - `id`: Factor variable that specifies the id of each subject.
#' - `visit`: Factor variable that specifies the visit of each assessment. Visit `0` denotes
#' the baseline visit.
#' - `group`: Factor variable that specifies which treatment group each subject belongs to.
#' - `outcome_bl`: Numeric variable that specifies the baseline outcome.
#' - `outcome_noICE`: Numeric variable that specifies the longitudinal outcome assuming
#' no ICEs.
#' - `ind_ice1`: Binary variable that takes value `1` if the corresponding visit is
#' affected by ICE1 and `0` otherwise.
#' - `dropout_ice1`: Binary variable that takes value `1` if the corresponding visit is
#' affected by the drop-out following ICE1 and `0` otherwise.
#' - `ind_ice2`: Binary variable that takes value `1` if the corresponding visit is affected
#' by ICE2.
#' - `outcome`: Numeric variable that specifies the longitudinal outcome including ICE1, ICE2
#' and the intermittent missing values.
#'
#' @export
simulate_data <- function(pars_c, pars_t, post_ice1_traj, strategies = getStrategies()) {

    assert_that(
        has_class(pars_c, "simul_pars"),
        has_class(pars_t, "simul_pars"),
        msg = "`pars_c` and `pars_t` must be of class `simul_pars`"
    )
    validate(pars_c)
    validate(pars_t)
    assert_that(
        length(pars_c$mu) == length(pars_t$mu),
        msg = "The two groups must have the same number of visits"
    )

    match.arg(
        arg = post_ice1_traj,
        choices = names(strategies)
    )
    strategy_fun <- strategies[[post_ice1_traj]]

    data_c <- generate_data_single(pars_c)
    data_t <- generate_data_single(
        pars_t,
        strategy_fun,
        distr_pars_ref = list(mu = pars_c$mu, sigma = pars_c$sigma)
    )

    data <- rbind(data_c, data_t)

    # overwrite ids to have unique ids
    n_visits <- length(pars_c$mu)
    unique_ids <- paste0("id_", seq.int(pars_c$n + pars_t$n))
    data$id <- factor(
        rep(unique_ids, each = n_visits),
        levels = unique_ids
    )

    # add group variable
    data$group <- factor(
        rep(
            c(
                rep("Control", pars_c$n),
                rep("Intervention", pars_t$n)
            ),
            each = n_visits
        ),
        levels = c("Control", "Intervention")
    )

    return(data)
}


#' Generate data for a single group
#'
#' @param pars_group A `simul_pars` object as generated by [set_simul_pars()]. It specifies
#' the simulation parameters of the given group.
#' @param strategy_fun Function implementing trajectories after the intercurrent event (ICE).
#' Must be one of [getStrategies()]. See [getStrategies()] for details. If `NULL` then post-ICE
#' outcomes are untouched.
#' @param distr_pars_ref Optional. Named list containing the simulation parameters of the
#' reference arm. It contains the following elements:
#' - `mu`: Numeric vector indicating the mean outcome trajectory assuming no ICEs. It should
#' include the outcome at baseline.
#' - `sigma` Covariance matrix of the outcome trajectory assuming no ICEs.
#' If `NULL`, then these parameters are inherited from `pars_group`.
#'
#' @inherit simulate_data return
#'
#' @seealso [simulate_data()].
generate_data_single <- function(pars_group, strategy_fun = NULL, distr_pars_ref = NULL) {

    n_visits <- length(pars_group$mu)

    data <- data.frame(
        id = as.factor(rep(paste0("id_", seq.int(pars_group$n)), each = n_visits)),
        visit = as.factor(rep(seq.int(n_visits) - 1, pars_group$n)),
        group = as.factor(rep("g", pars_group$n)),
        outcome_bl = NA,
        outcome_noICE = c(replicate(pars_group$n, sample_mvnorm(pars_group$mu, pars_group$sigma)))
    )
    data$outcome_bl <- rep(data$outcome_noICE[data$visit == "0"], each = n_visits)

    data$id <- factor(data$id, levels = unique(data$id))

    data$ind_ice1 <- simulate_ice(
        outcome = data$outcome_noICE,
        visits = data$visit,
        ids = data$id,
        prob_ice = pars_group$prob_ice1,
        or_outcome_ice = pars_group$or_outcome_ice1,
        baseline_mean = pars_group$mu[1]
    )

    # ICE2 happens if a subject is on treatment and drops-out
    data$ind_ice2 <- simulate_dropout(
        prob_dropout = pars_group$prob_ice2,
        ids = data$id,
        subset = ifelse(data$ind_ice1 == 0, 1, 0)
    )

    # remove ICE1 if it happened after ICE2
    data$ind_ice1 <- unlist(
        lapply(
            split(
                data[, c("ind_ice1", "ind_ice2")],
                data$id
            ),
            function(x) {
                if (sum(x$ind_ice1) < sum(x$ind_ice2)) {
                    x$ind_ice1 <- 0
                }
                return(x$ind_ice1)
            }
        )
    )

    # return binary vector with "1" where first ICE1 visit
    first_ice1_visit <- unlist(
        tapply(
            data$ind_ice1,
            data$id,
            function(x) {
                subset <- rep(0, n_visits)
                subset[which(x == 1)[1]] <- 1
                return(subset)
            }
        ),
        use.names = FALSE
    )

    data$dropout_ice1 <- simulate_dropout(
        prob_dropout = pars_group$prob_post_ice1_dropout,
        ids = data$id,
        subset = first_ice1_visit
    )

    if (!is.null(strategy_fun)) {
        data$outcome <- adjust_trajectories(
            distr_pars_group = list(mu = pars_group$mu, sigma = pars_group$sigma),
            outcome = data$outcome_noICE,
            ids = data$id,
            ind_ice = data$ind_ice1,
            strategy_fun = strategy_fun,
            distr_pars_ref = distr_pars_ref
        )
    } else {
        data$outcome <- data$outcome_noICE
    }

    data$outcome[data$dropout_ice1 == 1 | data$ind_ice2 == 1] <- NA

    rand_miss <- rbinom(n = pars_group$n * (n_visits - 1), size = 1, prob = pars_group$prob_miss)
    data$outcome[data$visit != "0"][rand_miss == 1] <- NA

    return(data)
}


#' Simulate intercurrent event
#'
#' @param outcome Numeric variable that specifies the longitudinal outcome for a single group.
#' @param visits Factor variable that specifies the visit of each assessment.
#' @param ids Factor variable that specifies the id of each subject.
#' @param prob_ice Numeric vector that specifies for each visit the probability of experiencing
#' the ICE after the current visit for a subject with outcome equal to the mean at baseline.
#' If a single numeric is provided, then the same probability is applied to each visit.
#' @param or_outcome_ice Numeric value that specifies the odds ratio of the ICE corresponding to
#' a +1 higher value of the outcome at the visit.
#' @param baseline_mean Mean outcome value at baseline.
#'
#' @details The probability of the ICE after each visit is modeled according to the following
#' logistic regression model:
#' `~ 1 + I(visit == 0) + ... + I(visit == n_visits-1) + I((x-alpha))` where:
#' - `n_visits` is the number of visits (including baseline).
#' - `alpha` is the baseline outcome mean set via argument `baseline_mean`.
#' The term `I((x-alpha))` specifies the dependency of the probability of the ICE on the current
#' outcome value.
#' The corresponding regression coefficients of the logistic model are defined as follows:
#' The intercept is set to 0, the coefficients corresponding to discontinuation after each visit
#' for a subject with outcome equal to
#' the mean at baseline are set according to parameter `or_outcome_ice`,
#' and the regression coefficient associated with the covariate `I((x-alpha))` is set to
#' `log(or_outcome_ice)`.
#'
#' @return A binary variable that takes value `1` if the corresponding outcome is affected
#' by the ICE and `0` otherwise.
#'
#' @importFrom stats rbinom model.matrix binomial
simulate_ice <- function(outcome, visits, ids, prob_ice, or_outcome_ice, baseline_mean) {

    assert_that(
        is.factor(visits) && length(visits) == length(outcome),
        msg = "`visits` must be a factor of length equal to the length of `outcome`"
    )

    assert_that(
        is.factor(ids) && length(ids) == length(outcome),
        msg = "`ids` must be a factor of length equal to the length of `outcome`"
    )

    n_visits <- nlevels(visits)
    if (length(prob_ice) == 1) {
        prob_ice <- rep(prob_ice, n_visits - 1)
    }
    prob_ice[prob_ice == 0] <- 1e-20
    prob_ice[prob_ice == 1] <- 1 - 1e-15

    model_coef <- c(log(prob_ice / (1 - prob_ice)), log(or_outcome_ice))

    dat <- data.frame(
        visits = visits,
        x = outcome - baseline_mean
    )
    no_lastvisit <- dat$visits != levels(dat$visits)[n_visits]
    dat <- dat[no_lastvisit, ]
    dat$visits <- factor(dat$visits, levels = unique(dat$visits))

    mod <- model.matrix(~ 0 + visits + x, dat)
    lp <- mod %*% model_coef
    probs_ice <- binomial(link = "logit")$linkinv(lp)
    ind_ice <- rbinom(n = length(probs_ice), size = 1, prob = probs_ice)
    ind_ice <- unlist(
        tapply(
            ind_ice, ids[no_lastvisit],
            function(x) c(0, cummax(x))
        ),
        use.names = FALSE
    )

    return(ind_ice)
}


#' Simulate drop-out
#'
#' @param prob_dropout Numeric that specifies the probability that a post-baseline visit is
#' affected by study drop-out.
#' @param ids Factor variable that specifies the id of each subject.
#' @param subset Binary variable that specifies the subset that could be affected by drop-out.
#' I.e. `subset` is a binary vector
#' of length equal to the length of `ids` that takes value `1` if the corresponding visit could
#' be affected by drop-out and `0` otherwise.
#'
#' @return A binary vector of length equal to the length of `ids` that takes value `1` if the
#' corresponding outcome is
#' affected by study drop-out.
#'
#' @details `subset` can be used to specify outcome values that cannot be affected by the
#' drop-out. By default
#' `subset` will be set to `1` for all the values except the values corresponding to the
#' baseline outcome, since baseline is supposed to not be affected by drop-out.
#' Even if `subset` is specified by the user, the values corresponding to the baseline
#' outcome are still hard-coded to be `0`.
#'
#' @importFrom stats rbinom
simulate_dropout <- function(prob_dropout, ids, subset = rep(1, length(ids))) {


    # baseline values cannot be missing
    subset <- unlist(tapply(subset, ids, function(x) {
        x[1] <- 0
        return(x)
    }))

    dropout <- rep(0, length(ids))
    dropout[subset == 1] <- rbinom(n = sum(subset), size = 1, prob = prob_dropout)
    dropout <- unlist(
        tapply(
            dropout,
            ids,
            function(x) pmin(cumsum(x), 1)
        ),
        use.names = FALSE
    )
    return(dropout)
}


#' Adjust trajectories due to the intercurrent event (ICE)
#'
#' @param distr_pars_group Named list containing the simulation parameters of the multivariate
#' normal distribution assumed for the given treatment group. It contains the following elements:
#' - `mu`: Numeric vector indicating the mean outcome trajectory. It should include the outcome
#' at baseline.
#' - `sigma` Covariance matrix of the outcome trajectory.
#' @param outcome Numeric variable that specifies the longitudinal outcome.
#' @param ids Factor variable that specifies the id of each subject.
#' @param ind_ice A binary variable that takes value `1` if the corresponding outcome is affected
#' by the ICE and `0` otherwise.
#' @param strategy_fun Function implementing trajectories after the intercurrent event (ICE). Must
#' be one of [getStrategies()]. See [getStrategies()] for details.
#' @param distr_pars_ref Optional. Named list containing the simulation parameters of the
#' reference arm. It contains the following elements:
#' - `mu`: Numeric vector indicating the mean outcome trajectory assuming no ICEs. It should
#' include the outcome at baseline.
#' - `sigma` Covariance matrix of the outcome trajectory assuming no ICEs.
#'
#' @return A numeric vector containing the adjusted trajectories.
#'
#' @seealso [adjust_trajectories_single()].
#'
adjust_trajectories <- function(
    distr_pars_group,
    outcome,
    ids,
    ind_ice,
    strategy_fun,
    distr_pars_ref = NULL
) {

    assert_that(
        all(!is.na(outcome)),
        msg = "`outcome` contains missing values"
    )

    outcome[ind_ice == 1] <- NA
    outcome <- unlist(
        tapply(
            outcome,
            ids,
            function(x) {
                adjust_trajectories_single(
                    outcome = x,
                    distr_pars_group = distr_pars_group,
                    strategy_fun = strategy_fun,
                    distr_pars_ref = distr_pars_ref
                )
            }
        ),
        use.names = FALSE
    )

    return(outcome)
}


#' Adjust trajectory of a subject's outcome due to the intercurrent event (ICE)
#'
#' @param distr_pars_group Named list containing the simulation parameters of the multivariate
#' normal distribution assumed for the given treatment group. It contains the following elements:
#' - `mu`: Numeric vector indicating the mean outcome trajectory. It should include the
#' outcome at baseline.
#' - `sigma` Covariance matrix of the outcome trajectory.
#' @param outcome Numeric variable that specifies the longitudinal outcome.
#' @param strategy_fun Function implementing trajectories after the intercurrent event (ICE).
#' Must be one of [getStrategies()]. See [getStrategies()] for details.
#' @param distr_pars_ref Optional. Named list containing the simulation parameters of the
#' reference arm. It contains the following elements:
#' - `mu`: Numeric vector indicating the mean outcome trajectory assuming no ICEs. It should
#' include the outcome at baseline.
#' - `sigma` Covariance matrix of the outcome trajectory assuming no ICEs.
#'
#' @return A numeric vector containing the adjusted trajectory for a single subject.
#'
#' @details `outcome` should be specified such that all-and-only the post-ICE observations
#' (i.e. the
#' observations to be adjusted) are set to `NA`.
#'
adjust_trajectories_single <- function(
    distr_pars_group,
    outcome,
    strategy_fun,
    distr_pars_ref = NULL
) {

    is_post_ice <- is.na(outcome)
    if (all(!is_post_ice)) {
        return(outcome)
    }

    if (is.null(distr_pars_ref)) {
        distr_pars_ref <- distr_pars_group
    }

    pars <- strategy_fun(
        pars_group = distr_pars_group,
        pars_ref = distr_pars_ref,
        index_mar = !is_post_ice
    )

    conditional_parameters <- get_conditional_parameters(pars, outcome)
    outcome[is_post_ice] <- impute_outcome(conditional_parameters)[[1]]
    return(outcome)
}


#' Simulate a realistic example dataset
#'
#' @description
#' Simulate a realistic example dataset using [simulate_data()] with hard-coded
#' values of all the input arguments.
#'
#' @details
#' [get_example_data()] simulates a 1:1 randomized trial of
#' an active drug (intervention) versus placebo (control) with 100 subjects per
#' group and 6 post-baseline assessments (bi-monthly visits until 12 months).
#' One intercurrent event corresponding to treatment discontinuation is also simulated.
#' Specifically, data are simulated under the following assumptions:
#'
#'- The mean outcome trajectory in the placebo group increases linearly from
#' 50 at baseline (visit 0) to 60 at visit 6, i.e. the slope is 10 points/year.
#'- The mean outcome trajectory in the intervention group is identical to the
#' placebo group up to visit 2. From visit 2 onward, the slope decreases by 50% to 5 points/year.
#'- The covariance structure of the baseline and follow-up values in both groups
#' is implied by a random intercept and slope model with a standard deviation of 5
#' for both the intercept and the slope, and a correlation of 0.25.
#' In addition, an independent residual error with standard deviation 2.5 is added
#' to each assessment.
#'- The probability of study drug discontinuation after each visit is calculated
#' according to a logistic model which depends on the observed outcome at that visit.
#' Specifically, a visit-wise discontinuation probability of 2% and 3% in the control
#' and intervention group, respectively, is specified in case the observed outcome is
#' equal to 50 (the mean value at baseline). The odds of a discontinuation is simulated
#' to increase by +10% for each +1 point increase of the observed outcome.
#'- Study drug discontinuation is simulated to have no effect on the mean trajectory in
#' the placebo group. In the intervention group, subjects who discontinue follow
#' the slope of the mean trajectory from the placebo group from that time point onward.
#' This is compatible with a copy increments in reference (CIR) assumption.
#'- Study drop-out at the study drug discontinuation visit occurs with a probability
#' of 50% leading to missing outcome data from that time point onward.
#'
#' @seealso [simulate_data()], [set_simul_pars()]
#'
#' @export
get_example_data <- function() {

    n <- 100
    time <- c(0, 2, 4, 6, 8, 10, 12)

    # Mean trajectory control
    muC <- c(50.0, 51.66667, 53.33333, 55.0, 56.66667, 58.33333, 60.0)

    # Mean trajectory intervention
    muT <- c(50.0, 51.66667, 53.33333, 54.16667, 55.0, 55.83333, 56.66667)

    # Create Sigma
    sd_error <- 2.5
    covRE <- rbind(
        c(25.0, 6.25),
        c(6.25, 25.0)
    )

    Sigma <- cbind(1, time / 12) %*% covRE %*% rbind(1, time / 12) + diag(sd_error^2, nrow = length(time))

    # Set probability of discontinuation
    probDisc_C <- 0.02
    probDisc_T <- 0.03
    or_outcome <- 1.10 # +1 point increase => +10% odds of discontinuation

    # Set drop-out rate following discontinuation
    prob_dropout <- 0.5

    # Set simulation parameters of the control group
    parsC <- set_simul_pars(
        mu = muC,
        sigma = Sigma,
        n = n,
        prob_ice1 = probDisc_C,
        or_outcome_ice1 = or_outcome,
        prob_post_ice1_dropout = prob_dropout
    )

    # Set simulation parameters of the intervention group
    parsT <- parsC
    parsT$mu <- muT
    parsT$prob_ice1 <- probDisc_T

    # Set assumption about post-ice trajectory
    post_ice_traj <- "CIR"

    # Simulate data
    data <- simulate_data(
        pars_c = parsC,
        pars_t = parsT,
        post_ice1_traj = post_ice_traj
    )

    return(data)
}
