library(rbmi)

set.seed(4711)
d = antidepressant_data

# Add Therapist (different for DRUG and PLACEBO)
d$CLUSTER = numeric(nrow(d))
d$CLUSTER[d$THERAPY == "DRUG"] = (1:10)[1 + as.numeric(d$PATIENT[d$THERAPY == "DRUG"]) %% 10]
d$CLUSTER[d$THERAPY == "PLACEBO"] = (11:20)[1 + as.numeric(d$PATIENT[d$THERAPY == "PLACEBO"]) %% 10]

# Add cluster effects
clu_eff = rnorm(20, mean=0, sd=1)
d$CHANGE = d$CHANGE + clu_eff[d$CLUSTER]

# Missing data = missing (e.g., 1513, 1514, 1804)
d[c(d$PATIENT %in% c(1511, 1513, 1514, 1804)), ]

# Replace by lines with NA
d = expand_locf(d, PATIENT=levels(d$PATIENT),
                VISIT=levels(d$VISIT),
                vars=c("BASVAL", "THERAPY", "CLUSTER"),
                group=c("PATIENT"),
                order=c("PATIENT", "VISIT"))

# See result
d[c(d$PATIENT %in% c(1513, 1514, 1804)), ]

# Complicated code that searches for the first NA in each patient
# and adds an entry "JR" to initiate a jump to the reference.
library(dplyr)
ice = d %>%
    arrange(PATIENT, VISIT) %>%
    filter(is.na(CHANGE)) %>%
    group_by(PATIENT) %>%
    slice(1) %>% ungroup() %>%
    select(PATIENT, VISIT) %>%
    mutate(strategy="JR")

# In 1513, 1514, 1517, JR at Visit 5. In 1804, JR at Visit 7.
head(ice)

# Roles of the variables in the imputation model. Note that we can have
# therapy x time interactions in the imputation model.
vars = set_vars(outcome="CHANGE",
                visit="VISIT", subjid="PATIENT", group="THERAPY",
                strata="THERAPY",
                cluster="CLUSTER",
                covariates=c("BASVAL*VISIT", "THERAPY*VISIT"))

# Estimate variance using "stratified" bootstrap
method = method_condmean(type="bootstrap", n_samples=999)
drw = draws(data=d, data_ice=ice, vars=vars, method=method, quiet=TRUE)

# Imputation
imp = impute(drw, references=c(DRUG="PLACEBO", PLACEBO="PLACEBO"))

# Roles of the variables in the primary analysis
vars = set_vars(subjid="PATIENT", outcome="CHANGE", visit="VISIT",
                # cluster="CLUSTER",
                group="THERAPY", covariates="BASVAL")
ana = analyse(imp, ancova, vars=vars)

# Pool results of imputations. The primary outcome is the last
# visit (trt_7)
pool(ana, conf.level=0.95, alternative="two.sided")





# jregress(formula, data, clusters, fe = FALSE, df.adjust = TRUE, df.display = FALSE, level = .95, tolerance = 1e-8, print = TRUE, digits = 4L, signif.stars = TRUE)
#
# jreg(X, Y, clusterid, fe = FALSE, df.adjust = TRUE, df.display = FALSE, level = .95, tolerance = 1e-8, print = TRUE, digits = 4L, signif.stars = TRUE)
#
# jregress arguments
#	formula		regression formula, similar to lm function
#	data		data file
#	clusters	cluster variable name (if omitted, regression is treated as non-clustered)
#
# jreg arguments
#	X		nxk regressor matrix
#			if fe=FALSE, X should (generally) contain an intercept column
#			if fe=TRUE,  X should not contain an intercept column
#	Y		nx1 dependent variable
#	clusterid	nx1 cluster variable
#
# common arguments
#	fe		TRUE to include cluster-level fixed effects (default is FALSE)
#	df.adjust	TRUE for degree-of-freedom adjustment (default)
#			FALSE for t(G-1) or t(n-k) distribution
#	df.display	TRUE to display adjustment coefficients (K and a)
#			FALSE to not display adjustment coefficients (default)
#	level		confidence level (default = .95)
#	tolerance	tolerance for singularity of X'X (default = 1e-8)
#	print		TRUE to print results (default)
#	digits		number of digits for display (default = 4L)
#	signif.stars	TRUE to display significance stars (default)
#
# Outputs
#	coefficients	vector of least squares coefficient estimates
#	se		vector of jackknife standard errors
#	pvalue		vector of p-values for coefficients
#	L		vector of lower bounds of confidence intervals for coefficients
#	U		vector of upper bounds of confidence intervals for coefficients
#	cov.mat		estimated covariance matrix of coefficient estimates
#	residuals	vector of least squares residuals
#	loo.residuals	vector of leave-cluster-out residuals
#	K		vector of adjustment degree-of-freedom coefficients
#	a		vector of scale adjustment coefficients
#	CV		leave-cluster-out cross-validation criterion
#---------------------------------------------------------------------------------#

jregress <- function(formula, data, clusters, fe = FALSE, df.adjust = TRUE, df.display = FALSE, level = .95, tolerance = 1e-8, print = TRUE, digits = 4L, signif.stars = TRUE){

    ## --- Parsing Arguments --- ## (add method for dealing w/ NA's)
    cl <- match.call()
    m <- match(c("formula","data","clusters"), names(cl), 0L)
    mf <- cl[c(1L,m)]
    mf[[1L]] <- quote(stats::model.frame)

    # create model frame and terms objects
    mf <- eval(mf, envir = parent.frame())
    mt <- attr(mf, "terms")

    attr(mt,"clustervar") <- cl$clusters  # add cluster id (if specified) to terms object
    clustervar <- mf$`(clusters)`         # extract cluster variable (if specified)

    # Error checking
    if (!is.null(clustervar) & !is.numeric(clustervar)){ # non-numeric cluster variable
        stop("Non-numeric cluster variable")
    }

    if (fe == TRUE){
        if (is.null(clustervar)){
            stop("Fixed effects option specified without a cluster variable")
        }

        cluster_id <- as.character(attr(mt,"clustervar"))
        varlist <- all.vars(attr(mt,"variables"))
        if (match(cluster_id,varlist,0L) > 0L){
            stop("Cluster variable cannot be included as an (in)dependent variable when fixed effects option is specified.")
        }
    }

    # generate response variable and design matrix
    Y <- model.response(mf, type = "numeric")
    X <- model.matrix(mt,mf)  #Q: see if we can account for collinearity w/ fe here...
    #A: try the contrasts.arg option, see link: https://www.r-bloggers.com/2022/02/how-to-include-all-levels-of-a-factor-variable-in-a-model-matrix-in-r/

    # omit intercept from model if fe option specified
    if (fe == TRUE){
        attr(mt,"intercept") <- 0	# update terms object; Q: necessary?
        X <- X[,-1]
    }

    ## --- Estimation and Output --- ##
    if (print == TRUE)   cat("\nCall:\n",paste(deparse(cl),sep="\n",collapse="\n"),"\n",sep="")

    out <- jreg(X,Y,clustervar,fe,df.adjust,df.display,level,tolerance,print,digits,signif.stars)
    out$call <- cl  # formula call
    out$terms <- mt # terms object
    return(out)
}

#---------------------------------------------------------------------------------#
### Ginv: Moore-Penrose pseudoinverse for symmetric matrix A ###
# 	Ginv(A,tolerance)
#
#	Argument:
#	A		symmetric matrix
#	tolerance	default = 1e-15
#
# 	Outputs:
#	inverse		pseudoinverse
#	rcond		condition number
#---------------------------------------------------------------------------------#

Ginv <- function(A,tolerance=1e-15){
    tol <- tolerance * max(A) * nrow(A)
    ei <- eigen(A,symmetric=TRUE)
    d <- ei$values
    H <- ei$vectors
    n <- (d < tol)
    di <- 1/(d+n)
    di[n] <- 0
    B <- H %*% (di*t(H))
    c <- min(d)/max(d)
    return(list(inverse=B,rcond=c))
}

#---------------------------------------------------------------------------------#
# jreg
#---------------------------------------------------------------------------------#

jreg <- function(X,Y,clusterid,fe = FALSE,df.adjust = TRUE,df.display = FALSE,level=.95,tolerance=1e-15,print = TRUE,digits = 4L,signif.stars = TRUE){
    n <- length(Y)
    k <- kF <- ncol(X)
    if (!is.null(clusterid)){
        c_id <- unique(clusterid)  # unique clusters
        G <- length(c_id)          # no. of clusters
        G_tab <- table(clusterid)   # no. obs. per cluster
    }
    ## within transformation ##
    if (fe == TRUE){
        YX <- cbind(Y,X)
        YX <- do.call(rbind, by(YX,clusterid,function(x){sweep(x,2,colMeans(x),FUN = "-")}))
        Y <- YX[,1]
        X <- as.matrix(YX[,-1])
        clusterid <- sort(clusterid) # sort clusterid to match the above lines
        kF <- kF + G
    }
    ## OLS estimates ##
    XX <- crossprod(X)
    if (rcond(XX) < tolerance) stop("Singular design matrix. (Condition number below threshold.) Recommend re-specify regression model.")
    Q <- solve(XX)
    XY <- crossprod(X,Y)
    betahat <- Q%*%XY
    ehat <- Y - X%*%betahat
    etilde <- ehat
    r.squared <- 1 - var(ehat)/var(Y)

    ## VCE estimate ##
    if (!is.null(clusterid)){ # clustered case

        beta_loo <- matrix(NA,G,k)
        rcond_loo <- matrix(NA,G,1)
        if (df.adjust == TRUE) {
            U <- V <- array(0,c(G,k,k))
            S <- matrix(0,G,k)
        }
        for (g in 1:G) {
            gi <- (clusterid==c_id[g])
            ni <- sum(gi)
            Yg <- matrix(Y[gi],nrow=ni)
            Xg <- matrix(X[gi,],nrow=ni,ncol=k)
            XXg <- crossprod(Xg)
            Ig <- Ginv(XX - XXg)
            Qg <- Ig$inverse
            rcond_loo[g] <- Ig$rcond
            beta_g <- Qg%*%(XY - crossprod(Xg,Yg))
            beta_loo[g,] <- beta_g - betahat
            etilde[gi] <- Yg - Xg%*%beta_g
            if (df.adjust == TRUE) {
                U[g,,] <- Ug <- Qg%*%XXg%*%Q
                V[g,,] <- Vg <- XXg%*%(Q + Ug)
                S[g,]  <- colSums((Q+Ug)*Vg)      }
        }

    } else { # non-clustered case

        beta_loo <- matrix(NA,n,k)
        rcond_loo <- matrix(NA,n,1)
        if (df.adjust == TRUE) {
            U <- V <- array(0,c(n,k,k))
            S <- matrix(0,n,k)
        }
        for (g in 1:n) {
            Xg <- matrix(X[g,],nrow=1,ncol=k)
            Yg <- Y[g]
            XXg <- crossprod(Xg)
            Ig <- Ginv(XX - XXg)
            Qg <- Ig$inverse
            rcond_loo[g] <- Ig$rcond
            beta_g <- Qg%*%(XY - crossprod(Xg,Yg))
            beta_loo[g,] <- beta_g - betahat
            etilde[g] <- Yg - Xg%*%beta_g
            if (df.adjust == TRUE) {
                U[g,,] <- Ug <- Qg%*%XXg%*%Q
                V[g,,] <- Vg <- XXg%*%(Q + Ug)
                S[g,]  <- colSums((Q+Ug)*Vg)
            }
        }
    }

    Vhat <- crossprod(beta_loo)
    se <- matrix(sqrt(diag(Vhat)))
    tstat <- betahat/se
    CV <- sum(etilde^2)
    noninvertible <- sum(rcond_loo < tolerance)

    if (df.adjust == FALSE) {
        if (!is.null(clusterid)) K <- G-1 else K <- n-kF
        a <- 1
    } else {
        a <- rep(0,k)
        K <- rep(0,k)
        for (j in 1:k) {
            Sj <- S[,j]
            Uj <- U[,,j]
            Vj <- V[,,j]
            Wj <- sweep(Uj,1,Sj,FUN="*")
            UU <- crossprod(Uj)
            UV <- crossprod(Uj,Vj)
            UW <- crossprod(Uj,Wj)
            VV <- crossprod(Vj)
            UUXX <- UU%*%XX
            trL <- sum(Sj) + sum(diag(UUXX)) - 2*sum(diag(UV))
            # The following uses the matrix fact:   sum(diag(A%*%B)) = sum(t(A)*B))
            trLL <- sum(Sj^2) + sum(t(UUXX)*UUXX) + 2*sum(t(UV)*UV) + 2*sum(XX*UW) - 4*sum(Wj*Vj) - 4*sum(t(UUXX)*UV) + 2*sum(UU*VV)
            K[j] <- (trL^2)/trLL
            a[j] <- sqrt(trL/Q[j,j])
        }
    }
    c <- qt((1+level)/2,K)/a
    Lb <- betahat - c*se
    Ub <- betahat + c*se
    pvalue <- 1 - pf((a*tstat)^2,1,K)

    # print output
    if (print == TRUE){
        cat("\nNumber of observations:", n,"\n")
        cat("Number of regressors:", k,"\n")
        if (!is.null(clusterid)){
            G.tab = G_tab
            G.summary = c(min(G.tab),mean(G.tab),median(G.tab),max(G.tab))
            G.info = structure(G.summary, names = c("Min","Mean","Median","Max"))
            cat("Number of clusters:", G,"\n")
            cat("\nNumber of observations per cluster:\n")
            print(G.info, digits = digits)
        }
        coef.table <- cbind(betahat, se, tstat, Lb, Ub, pvalue)
        if (df.adjust == TRUE & df.display == TRUE) {
            coef.table <- cbind(coef.table[,1:5], K, a, coef.table[,6])
            colnames(coef.table) <- c("Estimate", "Std. Error", "t-ratio","Lower","Upper","df","scale", "Pr(>|t|)")
        } else colnames(coef.table) <- c("Estimate", "Std. Error", "t-ratio","Lower","Upper","Pr(>|t|)")
        cat("\nCoefficients:\n")
        printCoefmat(coef.table, digits = digits, signif.stars = signif.stars)

        if (!is.null(clusterid)) cat("\nStandard errors calculated by the delete-cluster jackknife.")
        else cat("\nStandard errors calculated by the jackknife.")

        if (df.adjust == TRUE) cat("\nConfidence intervals and p-values calculated using adjusted student t distribution.\n")
        else {
            if (!is.null(clusterid)) cat("\nConfidence intervals and p-values calculated using t(G-1) distribution.\n")
            else cat("\nConfidence intervals and p-values calculated using t(n-k) distribution.\n")
        }

        cat("\nR-squared: ", formatC(r.squared, digits = digits))
        if (!is.null(clusterid)) cat("\nCluster Cross-Validation Criterion: ", formatC(CV, digits = digits, format = "f"))
        else                     cat("\nCross-Validation Criterion: ", formatC(CV, digits = digits, format = "f"))
        cat("\n")
        if (noninvertible > 0) {
            cat("\nWarning: Jackknife coefficients not identified in", formatC(noninvertible,digits=0,format = "f"))
            cat("/")
            if (!is.null(clusterid)) cat(formatC(G,digits=0,format = "f"),"clusters.")
            else                     cat(formatC(n,digits=0,format = "f"),"iterations.")
            cat("\nJackknife standard errors may be conservative. Recommend respecification of regression model.")
            cat("\n")
        }
    }
    out <- list(coefficients=betahat,se=se,pvalue=pvalue,cov.mat=Vhat,residuals=ehat,loo.residuals=etilde,K=K,a=a,CV=CV,L=Lb,U=Ub)
    return(out)
}

jregress(CHANGE ~ BASVAL + THERAPY, data=d[d$VISIT==7, ], clusters=NULL)
jregress(CHANGE ~ BASVAL + THERAPY, data=d[d$VISIT==7, ], clusters=CLUSTER)

