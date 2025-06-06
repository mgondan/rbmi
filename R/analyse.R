
#' Analyse Multiple Imputed Datasets
#'
#' @description
#' This function takes multiple imputed datasets (as generated by
#' the [impute()] function) and runs an analysis function on
#' each of them.
#'
#' @details
#' This function works by performing the following steps:
#'
#' 1. Extract a dataset from the `imputations` object.
#' 2. Apply any delta adjustments as specified by the `delta` argument.
#' 3. Run the analysis function `fun` on the dataset.
#' 4. Repeat steps 1-3 across all of the datasets inside the `imputations`
#' object.
#' 5. Collect and return all of the analysis results.
#'
#' The analysis function `fun` must take a `data.frame` as its first
#' argument. All other options to [analyse()] are passed onto `fun`
#' via `...`.
#' `fun` must return a named list with each element itself being a
#' list containing a single
#' numeric element called `est` (or additionally `se` and `df` if
#' you had originally specified [method_bayes()] or [method_approxbayes()])
#' i.e.:
#' ```
#' myfun <- function(dat, ...) {
#'     mod_1 <- lm(data = dat, outcome ~ group)
#'     mod_2 <- lm(data = dat, outcome ~ group + covar)
#'     x <- list(
#'         trt_1 = list(
#'             est = coef(mod_1)[[group]],
#'             se = sqrt(vcov(mod_1)[group, group]),
#'             df = df.residual(mod_1)
#'         ),
#'         trt_2 = list(
#'             est = coef(mod_2)[[group]],
#'             se = sqrt(vcov(mod_2)[group, group]),
#'             df = df.residual(mod_2)
#'         )
#'      )
#'      return(x)
#'  }
#' ```
#'
#' Please note that the `vars$subjid` column (as defined in the original call to
#' [draws()]) will be scrambled in the data.frames that are provided to `fun`.
#' This is to say they will not contain the original subject values and as such
#' any hard coding of subject ids is strictly to be avoided.
#'
#' By default `fun` is the [ancova()] function.
#' Please note that this function
#' requires that a `vars` object, as created by [set_vars()], is provided via
#' the `vars` argument e.g. `analyse(imputeObj, vars = set_vars(...))`. Please
#' see the documentation for [ancova()] for full details.
#' Please also note that the theoretical justification for the conditional mean imputation
#' method (`method = method_condmean()` in [draws()]) relies on the fact that ANCOVA is
#' a linear transformation of the outcomes.
#' Thus care is required when applying alternative analysis functions in this setting.
#'
#' The `delta` argument can be used to specify offsets to be applied
#' to the outcome variable in the imputed datasets prior to the analysis.
#' This is typically used for sensitivity or tipping point analyses. The
#' delta dataset must contain columns `vars$subjid`, `vars$visit` (as specified
#' in the original call to [draws()]) and `delta`. Essentially this `data.frame`
#' is merged onto the imputed dataset by `vars$subjid` and `vars$visit` and then
#' the outcome variable is modified by:
#'
#' ```
#' imputed_data[[vars$outcome]] <- imputed_data[[vars$outcome]] + imputed_data[["delta"]]
#' ```
#'
#' Please note that in order to provide maximum flexibility, the `delta` argument
#' can be used to modify any/all outcome values including those that were not
#' imputed. Care must be taken when defining offsets. It is recommend that you
#' use the helper function [delta_template()] to define the delta datasets as
#' this provides utility variables such as `is_missing` which can be used to identify
#' exactly which visits have been imputed.
#'
#' @seealso [extract_imputed_dfs()] for manually extracting imputed
#' datasets.
#' @seealso [delta_template()] for creating delta data.frames.
#' @seealso [ancova()] for the default analysis function.
#'
#' @param imputations An `imputations` object as created by [impute()].
#' @param fun An analysis function to be applied to each imputed dataset. See details.
#' @param delta A `data.frame` containing the delta transformation to be applied to the imputed
#' datasets prior to running `fun`. See details.
#' @param ... Additional arguments passed onto `fun`.
#' @param ncores The number of parallel processes to use when running this function. Can also be a
#' cluster object created by [`make_rbmi_cluster()`]. See the parallisation section below.
#' @param .validate Should `inputations` be checked to ensure it conforms to the required format
#' (default = `TRUE`) ? Can gain a small performance increase if this is set to `FALSE` when
#' analysing a large number of samples.
#'
#' @section Parallelisation:
#' To speed up the evaluation of `analyse()` you can use the `ncores` argument to enable parallelisation.
#' Simply providing an integer will get `rbmi` to automatically spawn that many background processes
#' to parallelise across. If you are using a custom analysis function then you need to ensure
#' that any libraries or global objects required by your function are available in the
#' sub-processes. To do this you need to use the [`make_rbmi_cluster()`] function for example:
#' ```
#' my_custom_fun <- function(...) <some analysis code>
#' cl <- make_rbmi_cluster(
#'     4,
#'     objects = list("my_custom_fun" = my_custom_fun),
#'     packages = c("dplyr", "nlme")
#' )
#' analyse(
#'     imputations = imputeObj,
#'     fun = my_custom_fun,
#'     ncores = cl
#' )
#' parallel::stopCluster(cl)
#' ```
#'
#' Note that there is significant overhead both with setting up the sub-processes and with
#' transferring data back-and-forth between the main process and the sub-processes. As such
#' parallelisation of the `analyse()` function tends to only be worth it when you have
#' `> 2000` samples generated by [`draws()`]. Conversely using parallelisation if your samples
#' are smaller than this may lead to longer run times than just running it sequentially.
#'
#' It is important to note that the implementation of parallel processing within [`analyse()`] has
#' been optimised around the assumption that the parallel processes will be spawned on the same
#' machine and not a remote cluster. One such optimisation is that the required data is saved to
#' a temporary file on the local disk from which it is then read into each sub-process. This is
#' done to avoid the overhead of transferring the data over the network. Our assumption is that
#' if you are at the stage where you need to be parallelising your analysis over a remote cluster
#' then you would likely be better off parallelising across multiple `rbmi` runs rather than within
#' a single `rbmi` run.
#'
#' Finally, if you are doing a tipping point analysis you can get a reasonable performance
#' improvement by re-using the cluster between each call to `analyse()` e.g.
#' ```
#' cl <- make_rbmi_cluster(4)
#' ana_1 <- analyse(
#'     imputations = imputeObj,
#'     delta = delta_plan_1,
#'     ncores = cl
#' )
#' ana_2 <- analyse(
#'     imputations = imputeObj,
#'     delta = delta_plan_2,
#'     ncores = cl
#' )
#' ana_3 <- analyse(
#'     imputations = imputeObj,
#'     delta = delta_plan_3,
#'     ncores = cl
#' )
#' parallel::clusterStop(cl)
#' ```
#'
#' @examples
#' \dontrun{
#' vars <- set_vars(
#'     subjid = "subjid",
#'     visit = "visit",
#'     outcome = "outcome",
#'     group = "group",
#'     covariates = c("sex", "age", "sex*age")
#' )
#'
#' analyse(
#'     imputations = imputeObj,
#'     vars = vars
#' )
#'
#' deltadf <- data.frame(
#'     subjid = c("Pt1", "Pt1", "Pt2"),
#'     visit = c("Visit_1", "Visit_2", "Visit_2"),
#'     delta = c( 5, 9, -10)
#' )
#'
#' analyse(
#'     imputations = imputeObj,
#'     delta = deltadf,
#'     vars = vars
#' )
#' }
#' @export
analyse <- function(
    imputations,
    fun = ancova,
    delta = NULL,
    ...,
    ncores = 1,
    .validate = TRUE
) {

    if (.validate) validate(imputations)

    assert_that(
        is.function(fun),
        msg = "`fun` must be a function"
    )

    assert_that(
        is.null(delta) | is.data.frame(delta),
        msg = "`delta` must be NULL or a data.frame"
    )

    vars <- imputations$data$vars

    if (.validate) devnull <- lapply(imputations$imputations, function(x) validate(x))

    if (!is.null(delta)) {
        expected_vars <- c(
            vars$subjid,
            vars$visit,
            "delta"
        )
        assert_that(
            all(expected_vars %in% names(delta)),
            msg = sprintf(
                "The following variables must exist witin `delta`: `%s`",
                paste0(expected_vars, collapse = "`, `")
            )
        )
    }

    # Mangle name to avoid any conflicts with user defined objects if running in a cluster
    ..rbmi..analysis..imputations <- imputations
    ..rbmi..analysis..delta <- delta
    ..rbmi..analysis..fun <- fun
    objects <- list(
        ..rbmi..analysis..imputations = ..rbmi..analysis..imputations,
        ..rbmi..analysis..delta = ..rbmi..analysis..delta,
        ..rbmi..analysis..fun = ..rbmi..analysis..fun
    )

    cl <- make_rbmi_cluster(ncores)

    if (is(cl, "cluster")) {
        ..rbmi..analysis..data..path <- tempfile()
        saveRDS(objects, file = ..rbmi..analysis..data..path, compress = FALSE)
        devnull <- parallel::clusterExport(cl, "..rbmi..analysis..data..path", environment())
        devnull <- parallel::clusterEvalQ(
            cl,
            {
                ..rbmi..analysis..objects <- readRDS(..rbmi..analysis..data..path)
                list2env(..rbmi..analysis..objects, envir = environment())
            }
        )
    }

    # If the user provided the clusters object directly then do not close it on completion
    if (!is(ncores, "cluster")) {
        on.exit(
            { if (!is.null(cl)) parallel::stopCluster(cl) },
            add = TRUE,
            after = FALSE
        )
    }

    # Chunk up requests for significant speed improvement when running in parallel
    number_of_cores <- ifelse(is.null(cl), 1, length(cl))
    indexes <- seq_along(imputations$imputations)
    indexes_split <- split(indexes, (indexes %% number_of_cores) + 1)

    results_list <- par_lapply(
        cl,
        function(indicies, ...) {
            inner_fun <- function(idx, ...) {
                dat2 <- extract_imputed_df(
                    ..rbmi..analysis..imputations$imputations[[idx]],
                    ..rbmi..analysis..imputations$data,
                    ..rbmi..analysis..delta
                )
                ..rbmi..analysis..fun(dat2, ...)
            }
            lapply(indicies, inner_fun, ...)
        },
        indexes_split,
        ...
    )
    results <- unlist(results_list, recursive = FALSE, use.names = FALSE)

    # Re-order to ensure results are returned in same order as imputations
    results <- results[order(unlist(indexes_split, use.names = FALSE))]
    names(results) <- NULL

    fun_name <- deparse(substitute(fun))
    if (length(fun_name) > 1) {
        fun_name <- "<Anonymous Function>"
    } else if (is.null(fun_name)) {
        fun_name <- "<NULL>"
    }

    ret <- as_analysis(
        results = results,
        fun_name = fun_name,
        delta = delta,
        fun = fun,
        method = imputations$method
    )
    validate(ret)
    return(ret)
}



#' Extract imputed datasets
#'
#' @description
#' Extracts the imputed datasets contained within an `imputations` object generated
#' by [impute()].
#'
#' @param imputations An `imputations` object as created by [impute()].
#' @param index The indexes of the imputed datasets to return. By default,
#' all datasets within the `imputations` object will be returned.
#' @param delta A `data.frame` containing the delta transformation to be
#' applied to the imputed dataset. See [analyse()] for details on the
#' format and specification of this `data.frame`.
#' @param idmap Logical. The subject IDs in the imputed `data.frame`'s are
#' replaced with new IDs to ensure they are unique. Setting this argument to
#' `TRUE` attaches an attribute, called `idmap`, to the returned `data.frame`'s
#' that will provide a map from the new subject IDs to the old subject IDs.
#'
#' @examples
#' \dontrun{
#' extract_imputed_dfs(imputeObj)
#' extract_imputed_dfs(imputeObj, c(1:3))
#' }
#' @returns
#' A list of data.frames equal in length to the `index` argument.
#'
#' @seealso [delta_template()] for creating delta data.frames.
#' @seealso [analyse()].
#' @export
extract_imputed_dfs <- function(
    imputations,
    index = seq_along(imputations$imputations),
    delta = NULL,
    idmap = FALSE
) {
    x <- imputations$imputations[index]
    lapply(
        x,
        function(x) extract_imputed_df(x, imputations$data, delta, idmap)
    )
}


#' Extract imputed dataset
#'
#' @description
#' Takes an imputation object as generated by [imputation_df()] and uses
#' this to extract a completed dataset from a `longdata` object as created
#' by [longDataConstructor()]. Also applies a delta transformation
#' if a `data.frame` is provided to the `delta` argument. See [analyse()] for
#' details on the structure of this `data.frame`.
#'
#' Subject IDs in the returned `data.frame` are scrambled i.e. are not the original
#' values.
#'
#' @param imputation An imputation object as generated by [imputation_df()].
#' @param ld A `longdata` object as generated by [longDataConstructor()].
#' @param delta Either `NULL` or a `data.frame`. Is used to offset outcome values in the imputed dataset.
#' @param idmap Logical. If `TRUE` an attribute called "idmap" is attached to
#' the return object which contains a `list` that maps the old subject ids
#' the new subject ids.
#' @returns A `data.frame`.
extract_imputed_df <- function(imputation, ld, delta = NULL, idmap = FALSE) {

    vars <- ld$vars
    dat <- ld$get_data(imputation, idmap = TRUE)
    id_map <- attr(dat, "idmap")

    if (!is.null(delta)) {
        # We are injecting a variable into the dataset so are using a obscured variable
        # name to remove the chance of a clash
        oldvar <- "old_subject_variable_zkfed1fgkadwni6g4oajd2aw"
        dat[[oldvar]] <- id_map[dat[[vars$subjid]]]
        delta[[oldvar]] <- delta[[vars$subjid]]
        dat2 <- apply_delta(
            dat,
            delta,
            group = c(oldvar, vars$visit),
            outcome = vars$outcome
        )
        dat2[[oldvar]] <- NULL
    } else {
        dat2 <- dat
    }

    dat2 <- as_dataframe(dat2)

    if (idmap) {
        attr(dat2, "idmap") <- id_map
    } else {
        attr(dat2, "idmap") <- NULL
    }

    return(dat2)
}



#' Construct an `analysis` object
#'
#' @description
#' Creates an analysis object ensuring that all components are
#' correctly defined.
#'
#' @param results A list of lists contain the analysis results for each imputation
#' See [analyse()] for details on what this object should look like.
#' @param method The method object as specified in [draws()].
#' @param delta The delta dataset used. See [analyse()] for details on how this
#' should be specified.
#' @param fun The analysis function that was used.
#' @param fun_name The character name of the analysis function (used for printing)
#' purposes.
as_analysis <- function(results, method, delta = NULL, fun = NULL, fun_name = NULL) {

    next_class <- switch(class(method)[[2]],
        bayes = "rubin",
        approxbayes = "rubin",
        condmean = ifelse(
            method$type == "jackknife",
            "jackknife",
            "bootstrap"
        ),
        bmlmi = "bmlmi"
    )

    assert_that(
        is.list(results),
        length(next_class) == 1,
        is.character(next_class),
        next_class %in% c("jackknife", "bootstrap", "rubin", "bmlmi")
    )

    x <- list(
        results = as_class(results, c(next_class, "list")),
        delta = delta,
        fun = fun,
        fun_name = fun_name,
        method = method
    )
    class(x) <- c("analysis", "list")
    validate(x)
    return(x)
}



#' Print `analysis` object
#'
#' @param x An `analysis` object generated by [analyse()].
#' @param ... Not used.
#' @importFrom utils capture.output
#' @export
print.analysis <- function(x, ...) {

    n_samp <- length(x$results)
    n_samp_string <- ife(
        has_class(x$results, "bootstrap") | has_class(x$results, "jackknife"),
        sprintf("1 + %s", n_samp - 1),
        n_samp
    )

    string <- c(
        "",
        "Analysis Object",
        "---------------",
        sprintf("Number of Results: %s", n_samp_string),
        sprintf("Analysis Function: %s", x$fun_name),
        sprintf("Delta Applied: %s", !is.null(x$delta)),
        "Analysis Estimates:",
        sprintf("    %s", names(x$results[[1]])),
        ""
    )

    cat(string, sep = "\n")
    return(invisible(x))
}



#' Validate `analysis` objects
#'
#' Validates the return object of the [analyse()] function.
#'
#' @param x An `analysis` results object (of class `"jackknife"`, `"bootstrap"`, `"rubin"`).
#' @param ... Not used.
#' @export
validate.analysis <- function(x, ...) {

    next_class <- class(x$results)[[1]]

    assert_that(
        next_class %in% c("jackknife", "bootstrap", "rubin", "bmlmi"),
        msg = "`results` must be of class 'jackknife', 'bootstrap', 'rubin' or 'bmlmi'"
    )

    if (next_class %in% c("bootstrap", "rubin")) {
        nsamp <- ife(
            next_class %in% c("bootstrap"),
            x$method$n_samples + 1,
            x$method$n_samples
        )
        assert_that(
            length(x$results) == nsamp
        )
    }

    assert_that(
        is.list(x$results),
        is.null(x$delta) | is.data.frame(x$delta),
        is.null(x$fun) | is.function(x$fun),
        is.null(x$fun_name) | is.character(x$fun_name)
    )

    validate(x$results)
}


#' @export
validate.jackknife <- function(x, ...) {
    validate_analyse_pars(x, get_pool_components("jackknife"))
}


#' @export
validate.bootstrap <- function(x, ...) {
    validate_analyse_pars(x, get_pool_components("bootstrap"))
}

#' @export
validate.rubin <- function(x, ...) {
    validate_analyse_pars(x, get_pool_components("rubin"))
}

#' @export
validate.bmlmi <- function(x, ...) {
    validate_analyse_pars(x, get_pool_components("bmlmi"))
}


#' Validate analysis results
#'
#' Validates analysis results generated by [analyse()].
#'
#'
#' @param results A list of results generated by the analysis `fun`
#' used in [analyse()].
#' @param pars A list of expected parameters in each of the analysis.
#' lists i.e. `c("est", "se", "df")`.
#'
#' @export
validate_analyse_pars <- function(results, pars) {

    assert_that(
        length(results) != 0,
        is.list(results),
        all(vapply(results, is.list, logical(1))),
        msg = "Analysis results must be a list of lists"
    )

    assert_that(
        length(names(results[[1]])) != 0,
        all(vapply(results, function(x) !is.null(names(x)) & all(names(x) != ""), logical(1))),
        msg = "Individual analysis results must be named lists"
    )

    results_names <- lapply(results, function(x) unique(names(x)))
    results_names_flat <- unlist(results_names, use.names = FALSE)
    results_names_count <- table(results_names_flat)

    assert_that(
        all(results_names_count == length(results)),
        msg = "Each individual analysis result must contain identically named elements"
    )

    results_unnested <- unlist(results, recursive = FALSE, use.names = FALSE)

    devnull <- lapply(
        results_unnested,
        function(x) {
            assert_that(
                is.list(x),
                all(pars %in% names(x)),
                msg = sprintf(
                    "Each individual analysis result element must be a list with elements `%s`",
                    paste0(pars, collapse = "`, `")
                )
            )
        }
    )

    for (par in pars) {
        if (!par %in% c("df", "se")) {
            assert_that(
                all(!is.na(vapply(results_unnested, function(x) x[[par]], numeric(1)))),
                msg = sprintf("Parameter `%s` contains missing values", par)
            )
        } else {
            assert_that(
                all(!is.na(vapply(results_unnested, function(x) x[[par]], numeric(1)))) ||
                all(is.na(vapply(results_unnested, function(x) x[[par]], numeric(1)))) ,
                msg = sprintf("Parameter `%s` contains both missing and observed values", par)
            )
        }
    }

    return(invisible(TRUE))
}
