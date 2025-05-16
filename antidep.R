library(rbmi)
d = antidepressant_data

# Missing data = missing (e.g., 1513, 1514, 1804)
d[c(d$PATIENT %in% c(1511, 1513, 1514, 1804)), ]

# Replace missing lines by NA
d = expand_locf(d, PATIENT=levels(d$PATIENT),
  VISIT=levels(d$VISIT),
  vars=c("BASVAL", "THERAPY"),
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
  mutate(strategy="MAR")

# In 1513, 1514, 1517, JR at Visit 5. In 1804, JR at Visit 7.
head(ice)

# Roles of the variables in the imputation model. Note that we can have
# therapy x time interactions in the imputation model.
vars = set_vars(outcome="CHANGE",
  visit="VISIT", subjid="PATIENT", group="THERAPY",
  strata=c("THERAPY"),
  covariates=c("BASVAL*VISIT", "THERAPY*VISIT"))

# Estimate variance using jackknife
method = method_condmean(type="jackknife")
drw = draws(data=d, data_ice=ice, vars=vars, method=method, quiet=TRUE)

# Jump to reference imputation
imp = impute(drw, references=c(DRUG="PLACEBO", PLACEBO="PLACEBO"))

# Roles of the variables in the primary analysis
vars = set_vars(subjid="PATIENT", outcome="CHANGE", visit="VISIT",
  group="THERAPY", covariates="BASVAL")
ana = analyse(imp, ancova, vars=vars)

# Pool results of imputations. The primary outcome is the last
# visit (trt_7)
pool(ana, conf.level=0.95, alternative="two.sided")
