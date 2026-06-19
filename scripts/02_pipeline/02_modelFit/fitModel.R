# =============================================================================
# Multispecies Occupancy Model — Point-level Spatial NNGP
# =============================================================================
#
# Hierarchy redefinition
# ----------------------
#   Sites (J)      = unique survey POINTS  — NNGP spatial GP lives here
#   Replicates (R) = YEARS surveyed at each point
#   Route          = grouping random intercept in occ.formula
#
# This framing treats each point as the primary spatial unit, letting the
# NNGP capture fine-scale spatial structure in occupancy (e.g. habitat
# patchiness at the 100–500 m scale) rather than coarser route-level patterns.
#
# Submodel structure
# ------------------
#   occ.formula:  ~ route_mean_DEM + (1 | route_id)
#       ψ varies with broad-scale elevation (route mean) plus a route-level
#       random intercept that absorbs unmeasured between-route heterogeneity.
#       The NNGP spatial random effect captures residual point-level spatial
#       structure after route and DEM are accounted for.
#
#   det.formula:  ~ point_DEM
#       Detection p varies with point-level elevation — separates local
#       topographic influence on detectability from occupancy-level effects.
#
# spOccupancy function
# --------------------
#   spMsPGOcc()  with NNGP = TRUE, coords = point coordinates
#
# Array dimensions
# ----------------
#   y[N_species, J_points, R_years]   — 3-D (no separate subunit axis)
#   NAs where a point was not surveyed in a given year
# =============================================================================


# ── 0. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("spOccupancy", "dplyr", "tidyr", "coda")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))


# =============================================================================
# 1. LOAD & PARSE DATA
# =============================================================================

dataForModelling <- ecosystemBreedingBird

if(ecosystemType %in% "Forest"){
  species_guild_groupings <- readxl::read_excel(file.path(path_finder("expertData"),
                                                          "processedData",
                                                          "speciesGroupings",
                                                          "forestBirdSpecies.xlsx"))
} else if(ecosystemType %in% "Mountain"){
  species_guild_groupings <- readxl::read_excel(file.path(path_finder("expertData"),
                                                          "processedData",
                                                          "speciesGroupings",
                                                          "mountainBirdSpecies.xlsx"))
}

group <- unique(species_guild_groupings$broad_guild)[2]
speciesInGroup <- species_guild_groupings$Species[species_guild_groupings$broad_guild %in% group]

dat <- dataForModelling %>%
  dplyr::filter(Species %in% speciesInGroup) %>%
  mutate(
    route_id  = sub("_.*", "", Route_point_id),
    point_id  = sub(".*_", "", Route_point_id),
    year_chr  = as.character(Year),
    detect    = as.integer(TotalCount > 0)
  )

#Remove routes with missing NAs in lon and lat
tt <- dat %>%
  dplyr::filter(is.na(Longitude) | is.na(Latitude)) %>%
  dplyr::select(-Longitude, -Latitude)

routes_with_missing_point_coords <- unique(tt$route_id)
cat("Number of routes taken from modelling:", length(routes_with_missing_point_coords),
    "out of ", length(unique(dat$route_id)), "due to missing coordinates")


dat <- dat %>%
  filter(!route_id %in% routes_with_missing_point_coords)

# Index vectors — points are the primary spatial unit
points  <- sort(unique(dat$Route_point_id))   # full "route_point" label as site ID
routes  <- sort(unique(dat$route_id))
years   <- sort(unique(dat$year_chr))
species <- sort(unique(dat$Species))

J <- length(points)    # number of unique survey points
R <- length(years)     # number of replicate years
N <- length(species)

cat(sprintf("Points (J): %d  |  Years (R): %d  |  Species (N): %d  |  Routes: %d\n",
            J, R, N, length(routes)))


# =============================================================================
# 2. BUILD THE 3-D DETECTION ARRAY  y[N, J, R]
# =============================================================================
# Dimensions: species × point × year
# If point-year visited by at least one species → missing species get 0
# If point-year not visited at all → all species remain NA
# =============================================================================

y <- array(
  NA_integer_,
  dim      = c(N, J, R),
  dimnames = list(species, points, years)
)

# Fill observed detections
for (i in seq_len(nrow(dat))) {
  sp  <- dat$Species[i]
  pt  <- dat$Route_point_id[i]
  yr  <- dat$year_chr[i]
  val <- dat$detect[i]

  cur <- y[sp, pt, yr]

  y[sp, pt, yr] <- if (is.na(cur)) {
    val
  } else {
    max(cur, val, na.rm = TRUE)
  }
}

# Identify point-year combinations that were visited
visited_point_year <- dat %>%
  dplyr::distinct(Route_point_id, year_chr) %>%
  dplyr::mutate(visited = TRUE)

visited_mat <- matrix(
  FALSE,
  nrow = J,
  ncol = R,
  dimnames = list(points, years)
)

visited_mat[
  as.matrix(
    visited_point_year %>%
      dplyr::transmute(
        point_i = match(Route_point_id, points),
        year_i  = match(year_chr, years)
      )
  )
] <- TRUE

# For visited point-years, replace missing species records with 0
for (pt in points) {
  for (yr in years) {
    if (visited_mat[pt, yr]) {
      y[is.na(y[, pt, yr]), pt, yr] <- 0L
    }
  }
}

cat("y dimensions      :", paste(dim(y), collapse = " × "), "\n")
cat("Naive detection   :", round(mean(y, na.rm = TRUE), 3), "\n")
cat("Proportion of NAs :", round(mean(is.na(y)), 3), "\n")


# =============================================================================
# 3. SITE-LEVEL (POINT-LEVEL) COVARIATES
# =============================================================================
# Each row in occ.covs / det.covs corresponds to one POINT (J rows total).
#
# occ.covs
#   route_mean_DEM : mean DEM across all points in the same route, standardised
#                    Represents broad-scale elevation for occupancy ψ
#   route_id       : factor — used as grouping variable for random intercept
#
# det.covs
#   point_DEM      : DEM at this specific point, standardised
#                    Represents local elevation effect on detection p
# =============================================================================

# ── Point-level summary table (one row per unique point) ─────────────────────
point_info <- dat %>%
  group_by(Route_point_id, route_id, point_id) %>%
  summarise(
    point_DEM = mean(DEM,       na.rm = TRUE),
    Lat       = mean(Latitude,  na.rm = TRUE),
    Lon       = mean(Longitude, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(match(Route_point_id, points))   # align to 'points' index

# ── Route-mean DEM: average the point DEMs within each route ─────────────────
route_mean_dem <- dat %>%
  group_by(route_id) %>%
  summarise(route_mean_DEM = mean(DEM, na.rm = TRUE), .groups = "drop")

point_info <- point_info %>%
  left_join(route_mean_dem, by = "route_id")

# ── Standardise both DEM variables (important: use separate scale() calls) ───
point_info <- point_info %>%
  mutate(
    route_mean_DEM_sc = as.numeric(scale(route_mean_DEM)),
    point_DEM_sc      = as.numeric(scale(point_DEM))
  )

# ── Verify alignment with points index ───────────────────────────────────────
stopifnot(all(point_info$Route_point_id == points))

# ── occ.covs : J rows, columns = occurrence covariates + grouping factor ──────
occ.covs <- data.frame(
  route_mean_DEM = point_info$route_mean_DEM_sc,
  route_id       = as.numeric(factor(point_info$route_id)),#, levels = routes)),
  row.names      = points
)



cat("\nOccurrence covariates summary:\n")
print(summary(occ.covs))

# ── det.covs : named list; each element is a [J × R] matrix ──────────────────
# point_DEM is constant across years for a given point → replicate the column
point_dem_mat <- matrix(
  rep(point_info$point_DEM_sc, times = R),
  nrow = J, ncol = R,
  dimnames = list(points, years)
)

det.covs <- list(point_DEM = point_dem_mat)

cat("\nDetection covariate (point DEM) summary:\n")
print(summary(as.vector(point_dem_mat)))


# =============================================================================
# 4. COORDINATES  (point-level — for NNGP)
# =============================================================================

coords <- as.matrix(point_info[, c("Lon", "Lat")])
rownames(coords) <- points

cat(sprintf("\nCoordinate range — Lon: [%.4f, %.4f]  Lat: [%.4f, %.4f]\n",
            min(coords[,1]), max(coords[,1]),
            min(coords[,2]), max(coords[,2])))

# Compute pairwise distance range for setting phi prior
dists  <- dist(coords)
max_d  <- max(dists)
min_d  <- quantile(dists, 0.05)   # 5th percentile avoids near-zero pairs
phi_lo <- 3 / max_d               # effective range = maximum extent
phi_hi <- 3 / min_d               # effective range = near-neighbour distance

cat(sprintf("Spatial decay prior: Uniform(%.5f, %.5f)\n", phi_lo, phi_hi))
cat(sprintf("Effective range covers %.1f to %.1f coordinate units\n",
            3 / phi_hi, 3 / phi_lo))


# =============================================================================
# 5. DATA LIST
# =============================================================================

data.list <- list(
  y        = y,
  occ.covs = occ.covs,
  det.covs = det.covs,
  coords   = coords
)


# y dimensions: species x site x year
visit_na <- is.na(data.list$y[1, , ])

# Apply the same missingness mask to every detection covariate
for (nm in names(data.list$det.covs)) {
  data.list$det.covs[[nm]][visit_na] <- NA
}

dim(data.list$y)
lapply(data.list$det.covs, dim)

stopifnot(all(
  sapply(data.list$det.covs, function(x) {
    identical(dim(x), dim(data.list$y)[2:3])
  })
))

#####
# Missing values in covariate for detection and occupancy
# Else there will be error in the model fit
#####
# Identify observed y values
y_observed <- !is.na(data.list$y[1, , ])

# Check problem cells
bad <- is.na(data.list$det.covs$point_DEM) & y_observed

sum(bad)

# Mean-impute point_DEM where y is observed
mean_point_DEM <- mean(data.list$det.covs$point_DEM, na.rm = TRUE)

data.list$det.covs$point_DEM[bad] <- mean_point_DEM
# =============================================================================
# 6. PRIORS
# =============================================================================
# beta  : occurrence coefficients  (intercept + route_mean_DEM)
# alpha : detection coefficients   (intercept + point_DEM)
# sigma.sq.psi : variance of route_id random intercept
# sigma.sq     : spatial variance (NNGP, one per species)
# phi          : spatial decay    (NNGP, one per species)
# =============================================================================

priors <- list(
  # Community-level occurrence (hyperpriors over species-level betas)
  beta.comm.normal  = list(mean = 0, var = 2.72),
  tau.sq.beta.ig    = list(a = 0.1, b = 0.1),

  # Community-level detection (hyperpriors over species-level alphas)
  alpha.comm.normal = list(mean = 0, var = 2.72),
  tau.sq.alpha.ig   = list(a = 0.1, b = 0.1),

  # Route random intercept variance: IG(2, 0.5) — moderately regularising
  sigma.sq.psi.ig   = list(a = 2, b = 0.5),

  # Spatial parameters
  sigma.sq.ig       = list(a = 2, b = 1),
  phi.unif          = list(a = phi_lo, b = phi_hi)
)


# =============================================================================
# 7. INITIAL VALUES
# =============================================================================

# Naive z: 1 if species detected at this point in any year
z.init <- apply(y, c(1, 2), function(x) as.integer(any(x == 1, na.rm = TRUE)))

inits <- list(
  # Community means: intercept + route_mean_DEM
  beta.comm    = c(0, 0),
  # Community means: intercept + point_DEM
  alpha.comm   = c(0, 0),

  tau.sq.beta  = c(1, 1),
  tau.sq.alpha = c(1, 1),

  # Route random intercept variance
  sigma.sq.psi = 1,

  # Spatial: one value per species
  sigma.sq     = rep(1,             N),
  phi          = rep(3 / (max_d / 2), N),   # start at mid-range decay

  z            = z.init
)


# =============================================================================
# 8. FIT MODEL  —  spMsPGOcc  (spatial multispecies, point-level NNGP)
# =============================================================================
#
# occ.formula:
#   ~ route_mean_DEM   — broad-scale elevation effect on occupancy
#   + (1 | route_id)   — route-level random intercept (lme4 syntax)
#   [NNGP spatial GP is added automatically by spMsPGOcc via coords]
#
# det.formula:
#   ~ point_DEM        — local elevation effect on detection
#
# NNGP settings:
#   n.neighbors = 10  — reasonable for dense point networks; increase to 15
#                       for better GP approximation if J is large
#   cov.model = "exponential" — standard choice; try "matern" if residuals
#                               show short-range structure
#
# MCMC settings:
#   Increase n.samples to >= 20000 and n.burn to >= 5000 for production runs.
#   n.chains >= 3 required for reliable R-hat convergence diagnostics.
# =============================================================================

cat("\n=== Fitting spMsPGOcc: point-level NNGP + route random intercept ===\n")

out <- spMsPGOcc(
  occ.formula   = ~ route_mean_DEM + (1 | route_id),
  det.formula   = ~ point_DEM,
  data          = data.list,
  inits         = inits,
  priors        = priors,
  cov.model     = "exponential",
  NNGP          = TRUE,
  n.neighbors   = 3,
  n.batch = 100, #change to 800
  batch.length = 5, #change to 25
  n.burn = 100, #change to 10000
  n.thin = 2, #change to 10
  n.chains      = 2,
  n.omp.threads = 2,
  verbose       = TRUE,
  n.report      = 1000
)

cat("\n=== Model summary ===\n")
summary(out)


# =============================================================================
# 9. CONVERGENCE
# =============================================================================

n.post <- nrow(out$beta.comm.samples)   # total post-burn posterior samples

check_rhat <- function(samples, parnames, n.chains) {
  chains <- lapply(seq_len(n.chains), function(ch) {
    idx <- seq(ch, nrow(samples), by = n.chains)
    as.mcmc(samples[idx, , drop = FALSE])
  })
  rhat <- gelman.diag(as.mcmc.list(chains))$psrf
  rownames(rhat) <- parnames
  rhat
}

cat("\n--- R-hat: community occurrence coefficients ---\n")
print(round(check_rhat(out$beta.comm.samples,
                       c("Intercept", "route_mean_DEM"), 2), 3))

cat("\n--- R-hat: community detection coefficients ---\n")
print(round(check_rhat(out$alpha.comm.samples,
                       c("Intercept", "point_DEM"), 2), 3))

cat("\n--- R-hat: route random intercept variance ---\n")
# sigma.sq.psi.samples: [n.post × 1]
sigma_samples <- matrix(out$sigma.sq.psi.samples, ncol = 1)
print(round(check_rhat(sigma_samples, "sigma.sq.psi (route)", 2), 3))


# =============================================================================
# 10. POSTERIOR PREDICTIVE CHECK  &  WAIC
# =============================================================================

ppc <- ppcOcc(out, fit.stat = "freeman-tukey", group = 1)
bp_species <- colMeans(ppc$fit.y.rep > ppc$fit.y, na.rm = TRUE)
bp_speciesForSpecies <- data.frame(
  Species = ppc$sp.names,
  Bayesian_p = bp_species
)

bp_overall <- mean(ppc$fit.y.rep > ppc$fit.y, na.rm = TRUE)

cat(sprintf(
  "\nOverall Bayesian p-value: %.3f  (0.3–0.7 = good fit)\n",
  bp_overall
))

waic <- waicOcc(out)
cat(sprintf("WAIC            : %.1f\n", waic[3]))


# =============================================================================
# 11. POSTERIOR SUMMARIES
# =============================================================================

post_summary <- function(samples, parnames) {
  tab <- t(apply(samples, 2, function(x)
    c(Mean    = round(mean(x), 3),
      SD      = round(sd(x), 3),
      `2.5%`  = round(quantile(x, 0.025), 3),
      `97.5%` = round(quantile(x, 0.975), 3),
      `Pr(>0)`= round(mean(x > 0), 3))))
  rownames(tab) <- parnames
  tab
}

cat("\n--- Community occurrence coefficients (logit scale) ---\n")
print(post_summary(out$beta.comm.samples,
                   c("Intercept", "route_mean_DEM")))

cat("\n--- Community detection coefficients (logit scale) ---\n")
print(post_summary(out$alpha.comm.samples,
                   c("Intercept", "point_DEM")))

cat("\n--- Route random intercept variance ---\n")
cat(sprintf("  Mean: %.3f  SD: %.3f  95%% CI: [%.3f, %.3f]\n",
            mean(out$sigma.sq.psi.samples),
            sd(out$sigma.sq.psi.samples),
            quantile(out$sigma.sq.psi.samples, 0.025),
            quantile(out$sigma.sq.psi.samples, 0.975)))

# ── Species-specific DEM effects on occupancy ─────────────────────────────────
# beta.samples: [n.post × N × p.occ]  — p.occ = 2 (intercept, route_mean_DEM)
cat("\n--- Species-specific route_mean_DEM effect on occupancy ---\n")

dem_cols <- grep("^route_mean_DEM-", colnames(out$beta.samples), value = TRUE)

dem_occ <- t(apply(out$beta.samples[, dem_cols, drop = FALSE], 2, function(x) {
  c(
    Mean     = round(mean(x), 3),
    `2.5%`   = round(quantile(x, 0.025), 3),
    `97.5%`  = round(quantile(x, 0.975), 3),
    `Pr(>0)` = round(mean(x > 0), 3)
  )
}))

rownames(dem_occ) <- sub("^route_mean_DEM-", "", dem_cols)

print(dem_occ)

# ── Species-specific point DEM effects on detection ──────────────────────────
# alpha.samples: [n.post × N × p.det]  — p.det = 2 (intercept, point_DEM)
cat("\n--- Species-specific point_DEM effect on detection ---\n")

dem_cols <- grep("^point_DEM-", colnames(out$alpha.samples), value = TRUE)

dem_det <- t(
  apply(
    out$alpha.samples[, dem_cols, drop = FALSE],
    2,
    function(x) {
      c(
        Mean     = round(mean(x), 3),
        `2.5%`   = round(quantile(x, 0.025), 3),
        `97.5%`  = round(quantile(x, 0.975), 3),
        `Pr(>0)` = round(mean(x > 0), 3)
      )
    }
  )
)

rownames(dem_det) <- sub("^point_DEM-", "", dem_cols)

print(dem_det)

# ── Spatial parameters (sigma.sq, phi) per species ───────────────────────────
# sigma.sq.samples: [n.post × N]
# phi.samples:      [n.post × N]
cat("\n--- Spatial parameters per species ---\n")

sigma_cols <- grep("^sigma.sq-", colnames(out$theta.samples), value = TRUE)
phi_cols   <- grep("^phi-", colnames(out$theta.samples), value = TRUE)

sigma_mean <- colMeans(out$theta.samples[, sigma_cols, drop = FALSE])
phi_mean   <- colMeans(out$theta.samples[, phi_cols, drop = FALSE])

sp_spatial <- data.frame(
  species = sub("^sigma.sq-", "", sigma_cols),
  sigma_sq = round(sigma_mean, 3),
  phi = round(phi_mean, 3),
  eff_range = round(3 / phi_mean, 4)
)

print(sp_spatial)


# =============================================================================
# 12. POINT-LEVEL OCCUPANCY  (posterior mean ψ per point)
# =============================================================================

# psi.samples: [n.post × N × J]
psi.mean <- apply(out$psi.samples, c(2, 3), mean)   # [N × J]
psi.ci   <- apply(out$psi.samples, c(2, 3),
                  function(x) quantile(x, c(0.025, 0.975)))
dimnames(psi.mean) <- list(species, points)

cat("\nMean occupancy probability (across points) per species:\n")
print(round(sort(rowMeans(psi.mean), decreasing = TRUE), 3))


# =============================================================================
# 13. SPATIAL RANDOM EFFECTS  (point-level)
# =============================================================================
# w.samples: [n.post × N × J]
# These are the NNGP spatial random effects on the logit(ψ) scale,
# capturing spatial dependence in occupancy after DEM and route are removed.
# =============================================================================

w.mean <- apply(out$w.samples, c(2, 3), mean)   # [N × J]
dimnames(w.mean) <- list(species, points)

cat("\nSpatial random effect range per species (min → max):\n")
for (n in seq_len(N)) {
  cat(sprintf("  %-20s  [%+.3f, %+.3f]\n",
              species[n], min(w.mean[n,]), max(w.mean[n,])))
}


# =============================================================================
# 14. ROUTE-LEVEL RANDOM INTERCEPTS
# =============================================================================
# beta.star.samples: [n.post × J]  — realised random effects (one per point,
# all points in the same route share the same drawn value from N(0, sigma.sq.psi))
# Marginalise to route level by averaging over points within each route.
# =============================================================================

# beta.star.samples shape: [n.post × J] (lme4-style, one entry per site)
route_re_mean <- sapply(as.numeric(as.factor(routes)), function(rt) {
  #pt_idx <- which(point_info$route_id == rt)
  rowMeans(out$beta.star.samples[, rt, drop = FALSE])
})  # [n.post × n.routes]

cat("\nRoute-level random intercepts (posterior mean, sorted):\n")
route_re <- sort(round(colMeans(route_re_mean), 3), decreasing = TRUE)
print(head(route_re, 10))


# =============================================================================
# 15. POINT-LEVEL SPECIES RICHNESS
# =============================================================================

# z.samples: [n.post × N × J]
richness.samples <- apply(out$z.samples, c(1, 3), sum)   # [n.post × J]
richness.mean    <- colMeans(richness.samples)
richness.ci      <- apply(richness.samples, 2,
                          function(x) quantile(x, c(0.025, 0.975)))
names(richness.mean) <- points

cat("\nTop 10 points by estimated species richness:\n")
top10 <- order(richness.mean, decreasing = TRUE)[1:min(10, J)]
print(data.frame(
  Point   = points[top10],
  Route   = point_info$route_id[top10],
  Mean    = round(richness.mean[top10], 1),
  CI_lo   = round(richness.ci[1, top10], 1),
  CI_hi   = round(richness.ci[2, top10], 1),
  DEM     = round(point_info$point_DEM[top10], 0)
), row.names = FALSE)


# =============================================================================
# 16. PLOTS
# =============================================================================

op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# ── 1. Spatial map of richness ────────────────────────────────────────────────
pal <- colorRampPalette(c("#E6F1FB", "#378ADD", "#042C53"))(100)
rich_col <- pal[cut(richness.mean, 100, labels = FALSE)]
plot(coords[,1], coords[,2], col = rich_col, pch = 16, cex = 1.2,
     xlab = "Longitude", ylab = "Latitude",
     main = "Estimated species richness\n(point level)")
legend("bottomleft",
       legend = c("Low", "High"),
       col = c(pal[1], pal[100]), pch = 16, bty = "n", cex = 0.8)

# ── 2. Spatial random effect — species 1 ─────────────────────────────────────
w1  <- w.mean[1,]
w_breaks <- seq(min(w1), max(w1), length.out = 101)
w_col <- colorRampPalette(c("#185FA5", "#F1EFE8", "#A32D2D"))(100)
w_idx <- findInterval(w1, w_breaks, all.inside = TRUE)
plot(coords[,1], coords[,2], col = w_col[w_idx], pch = 16, cex = 1.2,
     xlab = "Longitude", ylab = "Latitude",
     main = paste0("Spatial random effect\n(", species[1], ")"))
legend("bottomleft", legend = c("Negative", "Positive"),
       col = c(w_col[1], w_col[100]), pch = 16, bty = "n", cex = 0.8)

# ── 3. Forest plot — route_mean_DEM effect on occupancy ──────────────────────
sp_ord  <- order(dem_occ[, "Mean"])
n_sp    <- nrow(dem_occ)
plot(dem_occ[sp_ord, "Mean"], seq_len(n_sp),
     xlim = range(c(dem_occ[, "2.5%"], dem_occ[, "97.5%"])) * 1.15,
     yaxt = "n", ylab = "", xlab = "Coefficient (logit)",
     main = "Route DEM effect\non occupancy (ψ)",
     pch = 16, col = "#0C447C")
abline(v = 0, lty = 2, col = "grey60")
segments(dem_occ[sp_ord, "2.5%"], seq_len(n_sp),
         dem_occ[sp_ord, "97.5%"], seq_len(n_sp),
         col = "#0C447C", lwd = 1.2)
axis(2, at = seq_len(n_sp),
     labels = abbreviate(rownames(dem_occ)[sp_ord], 8),
     las = 2, cex.axis = 0.65)

# ── 4. Forest plot — point_DEM effect on detection ───────────────────────────
sp_ord2 <- order(dem_det[, "Mean"])
plot(dem_det[sp_ord2, "Mean"], seq_len(n_sp),
     xlim = range(c(dem_det[, "2.5%"], dem_det[, "97.5%"])) * 1.15,
     yaxt = "n", ylab = "", xlab = "Coefficient (logit)",
     main = "Point DEM effect\non detection (p)",
     pch = 16, col = "#3B6D11")
abline(v = 0, lty = 2, col = "grey60")
segments(dem_det[sp_ord2, "2.5%"], seq_len(n_sp),
         dem_det[sp_ord2, "97.5%"], seq_len(n_sp),
         col = "#3B6D11", lwd = 1.2)
axis(2, at = seq_len(n_sp),
     labels = abbreviate(rownames(dem_det)[sp_ord2], 8),
     las = 2, cex.axis = 0.65)

# ── 5. Effective spatial range per species ────────────────────────────────────
eff_range <- 3 / colMeans(out$phi.samples)
barplot(sort(eff_range, decreasing = TRUE),
        names.arg = abbreviate(species[order(eff_range, decreasing = TRUE)], 6),
        las = 2, cex.names = 0.7,
        xlab = "Species", ylab = "Effective range (coord. units)",
        main = "NNGP effective spatial range\nper species",
        col = "#534AB7")

# ── 6. MCMC trace — community route_mean_DEM effect ──────────────────────────
plot(out$beta.comm.samples[, 2], type = "l", col = "#185FA5",
     xlab = "MCMC iteration", ylab = "Coefficient",
     main = "Trace: community route_mean_DEM\n(occurrence)")
abline(h = 0, lty = 2, col = "grey60")

par(op)


# =============================================================================
# 17. SAVE RESULTS
# =============================================================================

results <- list(
  model         = out,
  waic          = waic,
  psi_mean      = psi.mean,          # [N × J] point-level occupancy
  w_mean        = w.mean,            # [N × J] spatial random effects
  richness_mean = richness.mean,     # [J] point-level richness
  richness_ci   = richness.ci,
  dem_occ       = dem_occ,           # species DEM effects on occupancy
  dem_det       = dem_det,           # species DEM effects on detection
  sp_spatial    = sp_spatial,        # spatial parameters per species
  point_info    = point_info,        # covariate + coordinate table
  dims          = list(N=N, J=J, R=R, species=species,
                       points=points, routes=routes, years=years)
)

saveRDS(results, "ms_point_spatial_results.rds")
cat("\nResults saved to ms_point_spatial_results.rds\n")
cat("✓ Script complete.\n")


# =============================================================================
# DESIGN NOTES
# =============================================================================
#
# Why points as sites?
#   The NNGP spatial process is defined over the site coordinates. Placing
#   it at the point level (finer scale) allows it to capture habitat
#   patchiness within routes — variation that a route-centroid model averages
#   away. The route random intercept (1 | route_id) in occ.formula absorbs
#   the residual between-route clustering that the NNGP doesn't explain.
#
# Why separate DEM variables?
#   route_mean_DEM in occ.formula captures macroclimatic elevation effects
#   on species range (e.g. montane vs lowland communities).
#   point_DEM in det.formula captures microtopographic effects on detection
#   (e.g. valley shadows reducing acoustic detectability). Using the same
#   variable at both levels would confound these two mechanisms.
#
# Why not point-level spatial in det.formula?
#   spMsPGOcc() places the spatial GP in the occupancy submodel only. Point-
#   level detection heterogeneity beyond point_DEM is handled by point_DEM
#   itself plus the species-specific detection intercept.
#
# Scaling up (J > 1000 points)
#   Increase n.neighbors to 15. Consider sfMsPGOcc() (shared spatial factors)
#   which scales as O(J × n.factors) rather than O(J × N).
# =============================================================================
