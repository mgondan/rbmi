library(rbmi)
library(mvtnorm)

sd_r = c(10, 15, 17, 18)
rho_r = 0.8
cor_r = cbind(
    c(1,         rho_r, rho_r^2, rho_r^3),
    c(rho_r,         1,   rho_r, rho_r^2),
    c(rho_r^2,   rho_r,       1,   rho_r),
    c(rho_r^3, rho_r^2,   rho_r,       1))

Mu_r = c(T0=100, T1=110, T2=112, T3=114)
Sigma_r = diag(sd_r) %*% cor_r %*% diag(sd_r)

n_r = 5000
d_r = cbind(X="R", as.data.frame(rmvnorm(n_r, Mu_r, Sigma_r)))

sd_a = c(10, 20, 23, 28)
rho_a = 0.6
cor_a = cbind(
    c(1,         rho_a, rho_a^2, rho_a^3),
    c(rho_a,         1,   rho_a, rho_a^2),
    c(rho_a^2,   rho_a,       1,   rho_a),
    c(rho_a^3, rho_a^2,   rho_a,       1))

Mu_a = c(T0=100, T1=120, T2=122, T3=124)
Sigma_a = diag(sd_a) %*% cor_a %*% diag(sd_a)

n_a = 5000
d_a = cbind(X="A", as.data.frame(rmvnorm(n_a, Mu_a, Sigma_a)))

N = n_r + n_a
d = rbind(d_r, d_a)
d$Id = factor(1:N)
rownames(d) = d$Id

# Covariance
round(Sigma_r)
round(cov(d[d$X == "R", c("T0", "T1", "T2", "T3")]))

round(Sigma_a)
round(cov(d[d$X == "A", c("T0", "T1", "T2", "T3")]))

# Ampute
i = as.character(sample(d$Id, size=N/3))
d_amp = d
d_amp[i, c("T2", "T3")] = NA

# Impute
long = reshape(d_amp, idvar="Id",
    times=c("T0", "T1", "T2", "T3"), timevar="Time",
    varying=list(c("T0", "T1", "T2", "T3")), v.names="Outcome",
    direction="long")

long$Id = factor(long$Id)
long$Time = factor(long$Time)
long$X = factor(long$X)

library(dplyr)
ice = long %>%
    arrange(Id, Time) %>%
    filter(is.na(Outcome)) %>%
    group_by(Id) %>%
    slice(1) %>%
    ungroup() %>%
    select(Id, Time) %>%
    mutate(strategy="MAR")

# Roles of the variables in the imputation model. Note that we can have
# therapy x time interactions in the imputation model.
vars = set_vars(outcome="Outcome",
    visit="Time", subjid="Id",
    group="X",
    strata="X",
    covariates="X*Time")

# Conditional mean imputation
B = 999
method = method_bayes(n_samples=B)
drw = draws(data=long, data_ice=ice, vars=vars, method=method, quiet=TRUE)
imp = impute(drw)

d_imp = list()
for(b in 1:B)
{
    ex = extract_imputed_dfs(imp, b)[[1]]
    d_imp[[b]] = reshape(ex, idvar="Id",
        times=c("T0", "T1", "T2", "T3"), timevar="Time",
        varying=list(c("T0", "T1", "T2", "T3")), v.names="Outcome",
        direction="wide")
}

# Estimate covariance for each imputation
cov_r = cov_a = array(NA, dim=c(4, 4, B))
for(b in 1:B)
{
    cov_r[, , b] = cov(d_imp[[b]][d_imp[[b]]$X == "R", c("T0", "T1", "T2", "T3")])
    cov_a[, , b] = cov(d_imp[[b]][d_imp[[b]]$X == "A", c("T0", "T1", "T2", "T3")])
}

# "Can you average them?" – "Sure. I can also..."
round(Sigma_r)
round(apply(cov_r, MARGIN=c(1, 2), mean))

round(Sigma_a)
round(apply(cov_a, MARGIN=c(1, 2), mean))
