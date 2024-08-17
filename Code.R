# Load necessary libraries
require(meta)
require(metafor)

# Conduct meta-analysis
meta <- metagen(
  TE = logES,
  seTE = SE,
  studlab = studlab,
  data = data,
  sm = "ES", 
  fixed = TRUE,
  random = TRUE,
  comb.random = TRUE
)

# Perform Egger's test
meta.eggers <- eggers.test(meta)

# Excess significance test
m.rma <- rma(
  yi = meta$TE,
  sei = meta$seTE,
  method = meta$method.tau,
  control = list(maxiter = 100000)
)
es <- tes(m.rma)

# Find outliers
outliers <- find.outliers(
  meta, 
  control = list(maxiter = 10000, stepadj = 0.5)
)

# Calculate p-values under different credibility ceilings
credibility_ceilings <- c(0.05, 0.10, 0.15, 0.20)
p_values <- numeric(length(credibility_ceilings))

for (i in seq_along(credibility_ceilings)) {
  cc <- credibility_ceilings[i]
  vi <- meta$seTE^2
  wi <- 1 / vi
  cc_weights <- (1 - cc) / (cc * (length(meta$TE) - 1))
  wi[wi > cc_weights] <- cc_weights
  res <- rma(
    yi = meta$TE, 
    sei = meta$seTE, 
    weights = wi, 
    method = meta$method.tau
  )
  p_values[i] <- coef(summary(res))[1, 4]
}