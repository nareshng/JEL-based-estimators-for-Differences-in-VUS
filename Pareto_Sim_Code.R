## ============================================================================
##
## FGM-copula-based bivariate Pareto model (FGMPM)
##
## Methods:
##   1. Jackknife empirical likelihood (JEL)
##   2. Normal approximation with stratified jackknife variance
##   3. Kernel-smoothed bootstrap percentile confidence interval
##
##
## ============================================================================



OUTPUT_DIR <- getwd()

MASTER_SEED <- 2026
ITERATIONS <- 1000L
BOOT <- 1000L
CONFIDENCE_LEVELS <- c(0.90, 0.95)
GRID_LEN_JEL <- 401L

## Set TRUE for a short diagnostic run before starting the full simulation.
QUICK_TEST <- FALSE


KEEP_REPLICATION_RESULTS <- TRUE


FAIL_FAST <- QUICK_TEST

if (QUICK_TEST) {
  ITERATIONS <- 10L
  BOOT <- 100L
  GRID_LEN_JEL <- 201L
}

LEVEL_TAG <- paste0(
  sprintf("%02d", round(100 * CONFIDENCE_LEVELS)),
  collapse = "_"
)

RUN_TAG <- paste0(
  if (QUICK_TEST) "quick_test" else "full",
  "_L",
  LEVEL_TAG
)

OUTFILE <- file.path(
  OUTPUT_DIR,
  paste0(
    "VUS_FGMPM_Pareto_simulation_corrected_",
    RUN_TAG,
    ".xlsx"
  )
)

CHECKPOINT_DIR <- file.path(
  OUTPUT_DIR,
  paste0("VUS_FGMPM_Pareto_checkpoints_", RUN_TAG)
)

## ----------------------------- packages ------------------------------------

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package could not be installed: ", pkg)
  }
}

need_pkg("openxlsx")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(OUTPUT_DIR) || file.access(OUTPUT_DIR, mode = 2) != 0) {
  stop("OUTPUT_DIR does not exist or is not writable: ", OUTPUT_DIR)
}

if (!dir.exists(CHECKPOINT_DIR) ||
    file.access(CHECKPOINT_DIR, mode = 2) != 0) {
  stop(
    "CHECKPOINT_DIR does not exist or is not writable: ",
    CHECKPOINT_DIR
  )
}

RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

## ----------------------------- validation --------------------------

assert_scalar_integer <- function(x, name, lower = 1L) {
  if (length(x) != 1L ||
      !is.finite(x) ||
      x != floor(x) ||
      x < lower) {
    stop(name, " must be a single integer greater than or equal to ", lower, ".")
  }
  invisible(TRUE)
}

assert_probability <- function(x, name, open = FALSE) {
  ok <- length(x) == 1L && is.finite(x)
  
  if (open) {
    ok <- ok && x > 0 && x < 1
  } else {
    ok <- ok && x >= 0 && x <= 1
  }
  
  if (!ok) {
    stop(name, if (open) " must lie in (0,1)." else " must lie in [0,1].")
  }
  invisible(TRUE)
}

safe_mean <- function(x) {
  keep <- is.finite(x)
  if (!any(keep)) return(NA_real_)
  mean(x[keep])
}

safe_mcse_mean <- function(x) {
  keep <- is.finite(x)
  if (sum(keep) <= 1L) return(NA_real_)
  stats::sd(x[keep]) / sqrt(sum(keep))
}

safe_binary_mcse <- function(x) {
  keep <- is.finite(x)
  if (!any(keep)) return(NA_real_)
  p <- mean(x[keep])
  sqrt(p * (1 - p) / sum(keep))
}

STATUS_CODEBOOK <- data.frame(
  Method = c(
    "JEL", "JEL", "JEL", "JEL", "JEL",
    "Normal", "Normal", "Normal",
    "Kernel", "Kernel", "Kernel", "Kernel"
  ),
  Code = c(
    0L, 1L, 2L, 3L, 99L,
    0L, 1L, 99L,
    0L, 1L, 2L, 99L
  ),
  Status = c(
    "ok",
    "degenerate_estimating_equations",
    "likelihood_not_zero_at_U0",
    "central_component_not_unique",
    "unhandled_replication_error",
    "ok",
    "invalid_variance",
    "unhandled_replication_error",
    "ok",
    "nonfinite_bootstrap_statistics",
    "invalid_quantile",
    "unhandled_replication_error"
  ),
  stringsAsFactors = FALSE
)

method_status_code <- function(method, status) {
  selected <- STATUS_CODEBOOK$Method == method &
    STATUS_CODEBOOK$Status == status
  
  if (sum(selected) == 1L) {
    return(STATUS_CODEBOOK$Code[selected])
  }
  
  99L
}

## Independent seeds make the simulated data invariant to BOOT and method order.
make_replication_seed <- function(setting_id, design_id, replication_id, stream) {
  assert_scalar_integer(setting_id, "setting_id")
  assert_scalar_integer(design_id, "design_id")
  assert_scalar_integer(replication_id, "replication_id")
  assert_scalar_integer(stream, "stream")
  
  if (design_id > 99L) {
    stop("The documented seed scheme supports at most 99 designs.")
  }
  
  if (replication_id > 9999L) {
    stop("The documented seed scheme supports at most 9999 replications.")
  }
  
  if (stream > 9L) {
    stop("stream must lie between 1 and 9 under the documented seed scheme.")
  }
  
  seed_value <- as.double(MASTER_SEED) +
    10000000 * setting_id +
    100000 * design_id +
    10 * replication_id +
    stream
  
  if (seed_value > .Machine$integer.max) {
    stop("Constructed random seed exceeds the R integer limit.")
  }
  
  as.integer(seed_value)
}

## ----------------------------- U-statistic  --------------------------

trip_count <- function(x, y, z) {
  if (anyNA(x) || anyNA(y) || anyNA(z) ||
      any(!is.finite(x)) || any(!is.finite(y)) || any(!is.finite(z))) {
    stop("trip_count requires finite observations without missing values.")
  }
  
  xs <- sort(x)
  zs <- sort(z)
  
  nx <- findInterval(y, xs, left.open = TRUE)  # number of x values satisfying x < y
  nz <- length(zs) - findInterval(y, zs)       # number of z values satisfying z > y
  
  sum(nx * nz)
}

u_stat <- function(x, y, z) {
  if (!is.matrix(x) || !is.matrix(y) || !is.matrix(z) ||
      ncol(x) < 2L || ncol(y) < 2L || ncol(z) < 2L) {
    stop("x, y, and z must be matrices with at least two columns.")
  }
  
  n1 <- nrow(x)
  n2 <- nrow(y)
  n3 <- nrow(z)
  
  if (any(c(n1, n2, n3) < 1L)) {
    stop("All three samples must be nonempty.")
  }
  
  c1 <- trip_count(x[, 1], y[, 1], z[, 1])
  c2 <- trip_count(x[, 2], y[, 2], z[, 2])
  
  (c1 - c2) / (n1 * n2 * n3)
}

## ----------------------------- kernel estimator ----------------------------

bw <- function(v, n) {
  assert_scalar_integer(n, "bandwidth sample size", lower = 2L)
  
  s <- min(
    stats::sd(v),
    stats::IQR(v) / 1.34,
    na.rm = TRUE
  )
  
  if (!is.finite(s) || s <= 0) {
    s <- stats::sd(v, na.rm = TRUE)
  }
  
  if (!is.finite(s) || s <= 0) {
    s <- 1e-8
  }
  
  0.9 * s * n^(-0.2)
}

kernel_trip_sum <- function(x, y, z, s12, s23) {
  if (!is.finite(s12) || !is.finite(s23) || s12 <= 0 || s23 <= 0) {
    stop("Kernel bandwidth combinations must be finite and positive.")
  }
  
  left_part <- rowSums(stats::pnorm(outer(y, x, "-") / s12))
  right_part <- colSums(stats::pnorm(outer(z, y, "-") / s23))
  
  sum(left_part * right_part)
}

kernel_stat <- function(x, y, z) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  n3 <- nrow(z)
  
  h1 <- bw(x[, 1], n1)
  h2 <- bw(y[, 1], n2)
  h3 <- bw(z[, 1], n3)
  h4 <- bw(x[, 2], n1)
  h5 <- bw(y[, 2], n2)
  h6 <- bw(z[, 2], n3)
  
  s12_1 <- sqrt(h1^2 + h2^2)
  s23_1 <- sqrt(h2^2 + h3^2)
  s12_2 <- sqrt(h4^2 + h5^2)
  s23_2 <- sqrt(h5^2 + h6^2)
  
  marker1 <- kernel_trip_sum(
    x[, 1], y[, 1], z[, 1], s12_1, s23_1
  )
  
  marker2 <- kernel_trip_sum(
    x[, 2], y[, 2], z[, 2], s12_2, s23_2
  )
  
  (marker1 - marker2) / (n1 * n2 * n3)
}

kernel_ci <- function(x0, y0, z0, boot = 1000L, level = 0.95, seed = NULL) {
  assert_scalar_integer(boot, "boot", lower = 20L)
  assert_probability(level, "level", open = TRUE)
  
  n1 <- nrow(x0)
  n2 <- nrow(y0)
  n3 <- nrow(z0)
  
  if (!is.null(seed)) {
    assert_scalar_integer(seed, "bootstrap seed")
    set.seed(seed)
  }
  
  kstat <- numeric(boot)
  
  for (b in seq_len(boot)) {
    x <- x0[
      sample.int(n1, size = n1, replace = TRUE),
      ,
      drop = FALSE
    ]
    
    y <- y0[
      sample.int(n2, size = n2, replace = TRUE),
      ,
      drop = FALSE
    ]
    
    z <- z0[
      sample.int(n3, size = n3, replace = TRUE),
      ,
      drop = FALSE
    ]
    
    kstat[b] <- kernel_stat(x, y, z)
  }
  
  finite_boot <- is.finite(kstat)
  number_valid <- sum(finite_boot)
  
  if (number_valid != boot) {
    return(list(
      ci = c(NA_real_, NA_real_),
      valid = FALSE,
      boot_requested = boot,
      boot_valid = number_valid,
      status = "nonfinite_bootstrap_statistics"
    ))
  }
  
  alpha <- 1 - level
  
  ci <- as.numeric(stats::quantile(
    kstat,
    probs = c(alpha / 2, 1 - alpha / 2),
    type = 1,
    names = FALSE,
    na.rm = FALSE
  ))
  
  ci_valid <-
    length(ci) == 2L && all(is.finite(ci)) && ci[1] <= ci[2]
  
  list(
    ci = ci,
    valid = ci_valid,
    boot_requested = boot,
    boot_valid = number_valid,
    status = if (ci_valid) "ok" else "invalid_quantile"
  )
}

## ----------------------------- stratified jackknife ------------------------

jackknife_values_stratified <- function(x, y, z) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  n3 <- nrow(z)
  
  if (any(c(n1, n2, n3) <= 1L)) {
    stop("Each group must contain at least two observations for jackknifing.")
  }
  
  U0 <- u_stat(x, y, z)
  
  vx <- numeric(n1)
  vy <- numeric(n2)
  vz <- numeric(n3)
  
  for (i in seq_len(n1)) {
    x_leave <- x[-i, , drop = FALSE]
    U_leave <- u_stat(x_leave, y, z)
    vx[i] <- n1 * U0 - (n1 - 1) * U_leave
  }
  
  for (j in seq_len(n2)) {
    y_leave <- y[-j, , drop = FALSE]
    U_leave <- u_stat(x, y_leave, z)
    vy[j] <- n2 * U0 - (n2 - 1) * U_leave
  }
  
  for (k in seq_len(n3)) {
    z_leave <- z[-k, , drop = FALSE]
    U_leave <- u_stat(x, y, z_leave)
    vz[k] <- n3 * U0 - (n3 - 1) * U_leave
  }
  
  tolerance <- 1e-8
  
  if (abs(mean(vx) - U0) > tolerance ||
      abs(mean(vy) - U0) > tolerance ||
      abs(mean(vz) - U0) > tolerance) {
    stop("Stratified pseudo-values do not average to U0.")
  }
  
  list(
    U0 = U0,
    vx = vx,
    vy = vy,
    vz = vz
  )
}

normal_ci_from_jackknife <- function(jk, level = 0.95) {
  assert_probability(level, "level", open = TRUE)
  
  n1 <- length(jk$vx)
  n2 <- length(jk$vy)
  n3 <- length(jk$vz)
  
  group_variances <- c(
    stats::var(jk$vx),
    stats::var(jk$vy),
    stats::var(jk$vz)
  )
  
  variance_estimate <-
    group_variances[1] / n1 +
    group_variances[2] / n2 +
    group_variances[3] / n3
  
  ## A zero estimated variance gives a degenerate point interval and does not
  ## satisfy the nondegeneracy conditions used for the normal approximation.
  
  valid <- is.finite(variance_estimate) && variance_estimate > 1e-14
  
  if (!valid) {
    return(list(
      ci_raw = c(NA_real_, NA_real_),
      ci_truncated = c(NA_real_, NA_real_),
      valid = FALSE,
      variance = variance_estimate,
      near_degenerate = TRUE,
      status = "invalid_variance"
    ))
  }
  
  zcrit <- stats::qnorm(1 - (1 - level) / 2)
  standard_error <- sqrt(variance_estimate)
  
  ci_raw <- c(
    jk$U0 - zcrit * standard_error,
    jk$U0 + zcrit * standard_error
  )
  
  ci_truncated <- c(
    max(-1, ci_raw[1]),
    min(1, ci_raw[2])
  )
  
  list(
    ci_raw = ci_raw,
    ci_truncated = ci_truncated,
    valid = all(is.finite(ci_raw)) && ci_raw[1] <= ci_raw[2],
    variance = variance_estimate,
    near_degenerate = any(group_variances < 1e-12),
    status = "ok"
  )
}

## ----------------------------- pooled JEL ----------------------------------

pooled_jackknife_from_stratified <- function(jk) {
  n1 <- length(jk$vx)
  n2 <- length(jk$vy)
  n3 <- length(jk$vz)
  n <- n1 + n2 + n3
  
  if (n <= 3L) {
    stop("The pooled sample size must exceed three.")
  }
  
  ## Exact algebraic transformation 
  
  pooled_x <- (n / (n - 3)) * (
    ((n - 1) / n1) * jk$vx - 2 * jk$U0
  )
  
  pooled_y <- (n / (n - 3)) * (
    ((n - 1) / n2) * jk$vy - 2 * jk$U0
  )
  
  pooled_z <- (n / (n - 3)) * (
    ((n - 1) / n3) * jk$vz - 2 * jk$U0
  )
  
  v <- c(pooled_x, pooled_y, pooled_z)
  
  coefficient <- function(group_size) {
    (n / (n - 3)) * (
      n - 3 - ((n - 1) * (group_size - 1) / group_size)
    )
  }
  
  a <- c(
    rep(coefficient(n1), n1),
    rep(coefficient(n2), n2),
    rep(coefficient(n3), n3)
  )
  
  if (abs(mean(v) - jk$U0) > 1e-8) {
    stop("Pooled pseudo-values do not average to U0.")
  }
  
  if (abs(mean(a) - 1) > 1e-10) {
    stop("The pooled centering coefficients do not average to one.")
  }
  
  list(
    U0 = jk$U0,
    v = v,
    a = a
  )
}

el_neg2log <- function(v, a, theta) {
  if (length(v) != length(a) ||
      any(!is.finite(v)) ||
      any(!is.finite(a)) ||
      !is.finite(theta)) {
    return(Inf)
  }
  
  g <- v - a * theta
  
  scale_g <- max(
    1,
    max(abs(v)),
    max(abs(a * theta))
  )
  
  ## If every estimating equation is zero, uniform empirical weights satisfy
  ## the constraint and -2 log R(theta) is exactly zero.
  
  if (max(abs(g)) <= 1e-12 * scale_g) {
    return(0)
  }
  
  ## Zero must lie in the interior of the scalar convex hull.
  if (!(min(g) < 0 && max(g) > 0)) {
    return(Inf)
  }
  
  lower_boundary <- max(-1 / g[g > 0])
  upper_boundary <- min(-1 / g[g < 0])
  width <- upper_boundary - lower_boundary
  
  if (!is.finite(width) || width <= 0) {
    return(Inf)
  }
  
  movement <- max(
    width * 1e-12,
    .Machine$double.eps * max(
      1,
      abs(lower_boundary),
      abs(upper_boundary)
    )
  )
  movement <- min(movement, width / 4)
  
  lower <- lower_boundary + movement
  upper <- upper_boundary - movement
  
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    return(Inf)
  }
  
  score <- function(lambda) {
    denominator <- 1 + lambda * g
    
    if (any(denominator <= 0) || any(!is.finite(denominator))) {
      return(NA_real_)
    }
    
    sum(g / denominator)
  }
  
  score_at_zero <- sum(g)
  score_scale <- max(1, sum(abs(g)))
  
  if (abs(score_at_zero) <= 1e-12 * score_scale) {
    lambda <- 0
  } else {
    score_lower <- score(lower)
    score_upper <- score(upper)
    
    if (!is.finite(score_lower) ||
        !is.finite(score_upper) ||
        score_lower * score_upper > 0) {
      return(Inf)
    }
    
    lambda <- tryCatch(
      stats::uniroot(
        score,
        lower = lower,
        upper = upper,
        tol = 1e-12
      )$root,
      error = function(e) NA_real_
    )
  }
  
  if (!is.finite(lambda)) {
    return(Inf)
  }
  
  log_arguments <- lambda * g
  
  if (any(1 + log_arguments <= 0) || any(!is.finite(log_arguments))) {
    return(Inf)
  }
  
  value <- 2 * sum(log1p(log_arguments))
  
  if (!is.finite(value)) {
    return(Inf)
  }
  
  if (value < -1e-8) {
    stop("The empirical-likelihood statistic is materially negative.")
  }
  
  max(value, 0)
}

refine_inside_outside_boundary <- function(
    outside_point,
    inside_point,
    is_inside,
    tolerance = 1e-10,
    max_iterations = 100L
) {
  outside <- outside_point
  inside <- inside_point
  
  if (is_inside(outside) || !is_inside(inside)) {
    stop("Invalid inside/outside bracket supplied for boundary refinement.")
  }
  
  for (iteration in seq_len(max_iterations)) {
    midpoint <- (outside + inside) / 2
    
    if (is_inside(midpoint)) {
      inside <- midpoint
    } else {
      outside <- midpoint
    }
    
    if (abs(inside - outside) <= tolerance) {
      break
    }
  }
  
  (inside + outside) / 2
}

jel_confidence_set_from_jackknife <- function(
    jk,
    level = 0.95,
    lower = -1,
    upper = 1,
    grid_len = 401L
) {
  assert_probability(level, "level", open = TRUE)
  assert_scalar_integer(grid_len, "grid_len", lower = 101L)
  
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop("Invalid JEL parameter bounds.")
  }
  
  pooled <- pooled_jackknife_from_stratified(jk)
  v <- pooled$v
  a <- pooled$a
  U0 <- pooled$U0
  critical_value <- stats::qchisq(level, df = 1)
  
  g_at_estimate <- v - a * U0
  scale_at_estimate <- max(1, max(abs(v)), max(abs(a * U0)))
  degenerate_at_estimate <-
    max(abs(g_at_estimate)) <= 1e-12 * scale_at_estimate
  
  if (degenerate_at_estimate) {
    return(list(
      intervals = matrix(numeric(0), ncol = 2L),
      central_interval = c(NA_real_, NA_real_),
      hull = c(NA_real_, NA_real_),
      total_length = NA_real_,
      valid = FALSE,
      components = 0L,
      degenerate = TRUE,
      lr_at_estimate = 0,
      U0 = U0,
      status = "degenerate_estimating_equations"
    ))
  }
  
  lr_value <- function(theta) {
    el_neg2log(v, a, theta)
  }
  
  is_inside <- function(theta) {
    value <- lr_value(theta)
    is.finite(value) && value <= critical_value + 1e-10
  }
  
  grid <- sort(unique(c(
    seq(lower, upper, length.out = grid_len),
    min(max(U0, lower), upper)
  )))
  
  likelihood_values <- vapply(grid, lr_value, numeric(1))
  inside <- is.finite(likelihood_values) &
    likelihood_values <= critical_value + 1e-10
  
  center_index <- which.min(abs(grid - U0))
  lr_at_estimate <- likelihood_values[center_index]
  
  if (!inside[center_index] ||
      !is.finite(lr_at_estimate) ||
      lr_at_estimate > 1e-7) {
    return(list(
      intervals = matrix(numeric(0), ncol = 2L),
      central_interval = c(NA_real_, NA_real_),
      hull = c(NA_real_, NA_real_),
      total_length = NA_real_,
      valid = FALSE,
      components = 0L,
      degenerate = degenerate_at_estimate,
      lr_at_estimate = lr_at_estimate,
      U0 = U0,
      status = "likelihood_not_zero_at_U0"
    ))
  }
  
  starts <- which(inside & !c(FALSE, head(inside, -1L)))
  ends <- which(inside & !c(tail(inside, -1L), FALSE))
  
  number_components <- length(starts)
  intervals <- matrix(
    NA_real_,
    nrow = number_components,
    ncol = 2L,
    dimnames = list(NULL, c("lower", "upper"))
  )
  
  for (component in seq_len(number_components)) {
    start_index <- starts[component]
    end_index <- ends[component]
    
    if (start_index == 1L) {
      lower_endpoint <- lower
    } else {
      lower_endpoint <- refine_inside_outside_boundary(
        outside_point = grid[start_index - 1L],
        inside_point = grid[start_index],
        is_inside = is_inside
      )
    }
    
    if (end_index == length(grid)) {
      upper_endpoint <- upper
    } else {
      upper_endpoint <- refine_inside_outside_boundary(
        outside_point = grid[end_index + 1L],
        inside_point = grid[end_index],
        is_inside = is_inside
      )
    }
    
    intervals[component, ] <- c(lower_endpoint, upper_endpoint)
  }
  
  central_component <- which(
    starts <= center_index & ends >= center_index
  )
  
  if (length(central_component) != 1L) {
    return(list(
      intervals = intervals,
      central_interval = c(NA_real_, NA_real_),
      hull = c(intervals[1, 1], intervals[nrow(intervals), 2]),
      total_length = sum(intervals[, 2] - intervals[, 1]),
      valid = FALSE,
      components = number_components,
      degenerate = degenerate_at_estimate,
      lr_at_estimate = lr_at_estimate,
      U0 = U0,
      status = "central_component_not_unique"
    ))
  }
  
  central_interval <- intervals[central_component, ]
  
  list(
    intervals = intervals,
    central_interval = as.numeric(central_interval),
    hull = c(intervals[1, 1], intervals[nrow(intervals), 2]),
    total_length = sum(intervals[, 2] - intervals[, 1]),
    valid = all(is.finite(intervals)) &&
      all(intervals[, 1] <= intervals[, 2]),
    components = number_components,
    degenerate = degenerate_at_estimate,
    lr_at_estimate = lr_at_estimate,
    U0 = U0,
    status = "ok"
  )
}

value_in_intervals <- function(value, intervals, tolerance = 1e-12) {
  if (!is.finite(value) || nrow(intervals) == 0L) {
    return(FALSE)
  }
  
  any(
    value >= intervals[, 1] - tolerance &
      value <= intervals[, 2] + tolerance
  )
}

## ----------------------------- FGM Pareto model ----------------------------

fgm_sample_uv <- function(n, theta) {
  assert_scalar_integer(n, "n")
  
  if (length(theta) != 1L || !is.finite(theta) || abs(theta) > 1) {
    stop("The FGM copula parameter theta must lie in [-1,1].")
  }
  
  u <- stats::runif(n)
  w <- stats::runif(n)
  
  a <- theta * (1 - 2 * u)
  v <- numeric(n)
  
  near_zero <- abs(a) < 1e-12
  v[near_zero] <- w[near_zero]
  
  nonzero <- !near_zero
  
  if (any(nonzero)) {
    discriminant <-
      (1 + a[nonzero])^2 - 4 * a[nonzero] * w[nonzero]
    
    if (any(discriminant < -1e-12)) {
      stop("Negative discriminant in the FGM conditional inversion.")
    }
    
    discriminant <- pmax(discriminant, 0)
    
    ## Rationalized root avoids subtractive cancellation near a = 0.
    denominator <-
      (1 + a[nonzero]) + sqrt(discriminant)
    
    zero_over_zero <- denominator == 0 & w[nonzero] == 0
    
    if (any(!is.finite(denominator)) ||
        any(denominator[!zero_over_zero] <= 0)) {
      stop("Invalid denominator in the FGM conditional inversion.")
    }
    
    conditional_v <- numeric(sum(nonzero))
    conditional_v[zero_over_zero] <- 0
    conditional_v[!zero_over_zero] <-
      2 * w[nonzero][!zero_over_zero] /
      denominator[!zero_over_zero]
    v[nonzero] <- conditional_v
  }
  
  if (any(!is.finite(v)) ||
      any(v < -1e-10) ||
      any(v > 1 + 1e-10)) {
    stop("Invalid FGM uniform variates were generated.")
  }
  
  v <- pmin(pmax(v, 0), 1)
  
  cbind(u, v)
}

fgmpm_sample <- function(n, lambda1, lambda2, alpha1, alpha2, theta) {
  assert_scalar_integer(n, "n")
  
  parameters <- c(lambda1, lambda2, alpha1, alpha2)
  
  if (any(!is.finite(parameters)) ||
      any(c(lambda1, lambda2) <= 0) ||
      any(c(alpha1, alpha2) <= 0)) {
    stop("Pareto scale and shape parameters must be finite and positive.")
  }
  
  uniforms <- fgm_sample_uv(n, theta)
  
  ## Stable inverse Pareto transformation.
  x1 <- exp(log(lambda1) - log1p(-uniforms[, 1]) / alpha1)
  x2 <- exp(log(lambda2) - log1p(-uniforms[, 2]) / alpha2)
  
  if (any(!is.finite(x1)) || any(!is.finite(x2))) {
    stop("Nonfinite Pareto observations were generated.")
  }
  
  cbind(x1, x2)
}

## Exact expression for P(X < Y < Z) under independent Pareto marginals.
vus_pareto_exact <- function(lx, ax, ly, ay, lz, az) {
  parameters <- c(lx, ax, ly, ay, lz, az)
  
  if (any(!is.finite(parameters)) ||
      any(c(lx, ly, lz) <= 0) ||
      any(c(ax, ay, az) <= 0)) {
    stop("Pareto scale and shape parameters must be finite and positive.")
  }
  
  lower_limit <- max(lx, ly)
  result <- 0
  
  ## For lower_limit <= t < lz, S_Z(t) = 1.
  if (lower_limit < lz) {
    result <- result + ay * ly^ay * (
      (lower_limit^(-ay) - lz^(-ay)) / ay -
        lx^ax * (
          lower_limit^(-(ay + ax)) -
            lz^(-(ay + ax))
        ) / (ay + ax)
    )
  }
  
  
  tail_limit <- max(lower_limit, lz)
  
  result <- result + ay * ly^ay * lz^az * (
    tail_limit^(-(ay + az)) / (ay + az) -
      lx^ax *
      tail_limit^(-(ay + az + ax)) /
      (ay + az + ax)
  )
  
  if (!is.finite(result) || result < -1e-12 || result > 1 + 1e-12) {
    stop("The exact Pareto VUS calculation produced an invalid probability.")
  }
  
  min(max(result, 0), 1)
}

## Copula parameters do not enter the true marginal VUS values. 
true_fgmpm_values <- function(
    lx1, lx2, ax1, ax2,
    ly1, ly2, ay1, ay2,
    lz1, lz2, az1, az2
) {
  vus1 <- vus_pareto_exact(lx1, ax1, ly1, ay1, lz1, az1)
  vus2 <- vus_pareto_exact(lx2, ax2, ly2, ay2, lz2, az2)
  
  list(
    vus1 = vus1,
    vus2 = vus2,
    theta = vus1 - vus2
  )
}

## ----------------------------- simulation settings -------------------------

pareto_settings <- data.frame(
  setting_id = 1:5,
  scenario_class = c(
    "regular_benchmark",
    "regular_benchmark",
    "stress",
    "stress",
    "extreme_stress"
  ),
  selection_basis = c(
    "Prespecified symmetric regular case",
    "Prespecified separated-target regular case",
    "Prespecified heterogeneous heavy-tail stress case",
    "Prespecified heterogeneous tail/dependence stress case",
    "Prespecified extreme heavy-tail/dependence stress case"
  ),
  label = c(
    "(1,1,1,1,0.5;1,1,1,1,0.5;1,1,1,1,-0.5)",
    "(1,2,1,2,0.5;2,1,2,1,0.5;2,2,1,1,-0.5)",
    "(5,1,1,5,0.2;5,1,2,2,0.5;1,5,5,1,0.9)",
    "(0.5,1,0.5,1,-0.5;1,1,0.5,0.5,-0.5;1,1,5,1,-0.2)",
    "(10,1,0.5,10,-0.9;1,1,15,0.2,-0.1;0.2,1,5,0.1,0.9)"
  ),
  lx1 = c(1, 1, 5, 0.5),
  lx2 = c(1, 2, 1, 1),
  ax1 = c(1, 1, 1, 0.5),
  ax2 = c(1, 2, 5, 1),
  tx = c(0.5, 0.5, 0.2, -0.5),
  ly1 = c(1, 2, 5, 1),
  ly2 = c(1, 1, 1, 1),
  ay1 = c(1, 2, 2, 0.5),
  ay2 = c(1, 1, 2, 0.5),
  ty = c(0.5, 0.5, 0.5, -0.5),
  lz1 = c(1, 2, 1, 1),
  lz2 = c(1, 2, 5, 1),
  az1 = c(1, 1, 5, 5),
  az2 = c(1, 1, 1, 1),
  tz = c(-0.5, -0.5, 0.9, -0.2),
  stringsAsFactors = FALSE
)

sample_sizes <- data.frame(
  design_id = 1:6,
  n1 = c(10, 30, 30, 60, 50, 30),
  n2 = c(20, 20, 30, 70, 80, 120),
  n3 = c(30, 10, 30, 40, 100, 50)
)

validate_simulation_settings <- function() {
  assert_scalar_integer(MASTER_SEED, "MASTER_SEED", lower = 0L)
  assert_scalar_integer(ITERATIONS, "ITERATIONS")
  assert_scalar_integer(BOOT, "BOOT", lower = 20L)
  if (length(CONFIDENCE_LEVELS) < 1L ||
      anyDuplicated(CONFIDENCE_LEVELS)) {
    stop("CONFIDENCE_LEVELS must contain unique confidence levels.")
  }
  invisible(lapply(
    CONFIDENCE_LEVELS,
    assert_probability,
    name = "each CONFIDENCE_LEVELS entry",
    open = TRUE
  ))
  assert_scalar_integer(GRID_LEN_JEL, "GRID_LEN_JEL", lower = 101L)
  
  if (ITERATIONS > 9999L) {
    stop(
      "ITERATIONS must not exceed 9999 under the documented seed scheme."
    )
  }
  
  size_matrix <- as.matrix(sample_sizes[, c("n1", "n2", "n3")])
  
  if (any(!is.finite(size_matrix)) ||
      any(size_matrix <= 1) ||
      any(size_matrix != floor(size_matrix))) {
    stop("Every sample size must be an integer greater than one.")
  }
  
  copula_parameters <- unlist(
    pareto_settings[, c("tx", "ty", "tz")],
    use.names = FALSE
  )
  
  if (any(!is.finite(copula_parameters)) ||
      any(abs(copula_parameters) > 1)) {
    stop("Every FGM copula parameter must lie in [-1,1].")
  }
  
  positive_columns <- c(
    "lx1", "lx2", "ax1", "ax2",
    "ly1", "ly2", "ay1", "ay2",
    "lz1", "lz2", "az1", "az2"
  )
  
  positive_values <- unlist(
    pareto_settings[, positive_columns],
    use.names = FALSE
  )
  
  if (any(!is.finite(positive_values)) || any(positive_values <= 0)) {
    stop("Every Pareto scale and shape parameter must be positive.")
  }
  
  if (anyDuplicated(pareto_settings$setting_id) ||
      anyDuplicated(sample_sizes$design_id)) {
    stop("setting_id and design_id values must be unique.")
  }
  
  invisible(lapply(
    pareto_settings$setting_id,
    assert_scalar_integer,
    name = "setting_id"
  ))
  
  invisible(lapply(
    sample_sizes$design_id,
    assert_scalar_integer,
    name = "design_id"
  ))
  
  invisible(TRUE)
}

## ----------------------------- one replication -----------------------------

run_one_dataset <- function(
    x,
    y,
    z,
    theta,
    level,
    boot,
    grid_len_jel,
    bootstrap_seed
) {
  
  method_errors <- character(0)
  
  run_timed_method <- function(method, expression, fallback) {
    method_start <- proc.time()[["elapsed"]]
    
    value <- tryCatch(
      force(expression),
      error = function(error_condition) {
        if (FAIL_FAST) {
          stop(error_condition)
        }
        
        method_errors <<- c(
          method_errors,
          paste0(method, ": ", conditionMessage(error_condition))
        )
        
        fallback
      }
    )
    
    list(
      value = value,
      elapsed = proc.time()[["elapsed"]] - method_start
    )
  }
  
  jel_fallback <- list(
    intervals = matrix(numeric(0), ncol = 2L),
    total_length = NA_real_,
    valid = FALSE,
    components = 0L,
    degenerate = FALSE,
    status = "unhandled_replication_error"
  )
  
  normal_fallback <- list(
    ci_raw = c(NA_real_, NA_real_),
    ci_truncated = c(NA_real_, NA_real_),
    valid = FALSE,
    near_degenerate = FALSE,
    status = "unhandled_replication_error"
  )
  
  kernel_fallback <- list(
    ci = c(NA_real_, NA_real_),
    valid = FALSE,
    boot_valid = 0L,
    status = "unhandled_replication_error"
  )
  
  ## The stratified pseudo-values are shared by JEL and Normal, but Kernel does
  ## not depend on them. A shared-jackknife failure therefore invalidates only
  ## JEL and Normal.
  jackknife_start <- proc.time()[["elapsed"]]
  jk <- tryCatch(
    jackknife_values_stratified(x, y, z),
    error = function(error_condition) {
      if (FAIL_FAST) {
        stop(error_condition)
      }
      
      method_errors <<- c(
        method_errors,
        paste0(
          "Shared Jackknife: ",
          conditionMessage(error_condition)
        )
      )
      NULL
    }
  )
  jackknife_time <- proc.time()[["elapsed"]] - jackknife_start
  
  jel_run <- if (is.null(jk)) {
    list(value = jel_fallback, elapsed = 0)
  } else {
    run_timed_method(
      method = "JEL",
      expression = jel_confidence_set_from_jackknife(
        jk = jk,
        level = level,
        lower = -1,
        upper = 1,
        grid_len = grid_len_jel
      ),
      fallback = jel_fallback
    )
  }
  jel <- jel_run$value
  jel_incremental_time <- jel_run$elapsed
  jel_time <- jackknife_time + jel_incremental_time
  
  normal_run <- if (is.null(jk)) {
    list(value = normal_fallback, elapsed = 0)
  } else {
    run_timed_method(
      method = "Normal",
      expression = normal_ci_from_jackknife(jk, level = level),
      fallback = normal_fallback
    )
  }
  normal <- normal_run$value
  normal_incremental_time <- normal_run$elapsed
  normal_time <- jackknife_time + normal_incremental_time
  
  kernel_run <- run_timed_method(
    method = "Kernel",
    expression = kernel_ci(
      x0 = x,
      y0 = y,
      z0 = z,
      boot = boot,
      level = level,
      seed = bootstrap_seed
    ),
    fallback = kernel_fallback
  )
  kernel <- kernel_run$value
  kernel_time <- kernel_run$elapsed
  
  jel_valid <- isTRUE(jel$valid)
  normal_valid <- isTRUE(normal$valid)
  kernel_valid <- isTRUE(kernel$valid)
  
  normal_ci_raw <- normal$ci_raw
  normal_ci_truncated <- normal$ci_truncated
  kernel_interval <- kernel$ci
  
  result <- c(
    U0 = if (is.null(jk)) NA_real_ else jk$U0,
    Shared_jackknife_time_seconds = jackknife_time,
    
    JEL_valid = as.numeric(jel_valid),
    JEL_cover = if (jel_valid) {
      as.numeric(value_in_intervals(theta, jel$intervals))
    } else {
      0
    },
    JEL_len = if (jel_valid) jel$total_length else NA_real_,
    JEL_components = if (jel_valid) jel$components else NA_real_,
    JEL_multicomponent = as.numeric(jel_valid && jel$components > 1L),
    JEL_degenerate = as.numeric(isTRUE(jel$degenerate)),
    JEL_status_code = method_status_code("JEL", jel$status),
    JEL_time_seconds = jel_time,
    
    Normal_valid = as.numeric(normal_valid),
    Normal_cover = if (normal_valid) {
      as.numeric(
        theta >= normal_ci_raw[1] && theta <= normal_ci_raw[2]
      )
    } else {
      0
    },
    Normal_len = if (normal_valid) diff(normal_ci_raw) else NA_real_,
    Normal_cover_truncated = if (normal_valid) {
      as.numeric(
        theta >= normal_ci_truncated[1] &&
          theta <= normal_ci_truncated[2]
      )
    } else {
      0
    },
    Normal_len_truncated = if (normal_valid) {
      diff(normal_ci_truncated)
    } else {
      NA_real_
    },
    Normal_near_degenerate = as.numeric(isTRUE(normal$near_degenerate)),
    Normal_status_code = method_status_code("Normal", normal$status),
    Normal_time_seconds = normal_time,
    
    Kernel_valid = as.numeric(kernel_valid),
    Kernel_cover = if (kernel_valid) {
      as.numeric(
        theta >= kernel_interval[1] && theta <= kernel_interval[2]
      )
    } else {
      0
    },
    Kernel_len = if (kernel_valid) diff(kernel_interval) else NA_real_,
    Kernel_boot_valid = kernel$boot_valid,
    Kernel_status_code = method_status_code("Kernel", kernel$status),
    Kernel_time_seconds = kernel_time
  )
  
  attr(result, "method_errors") <- method_errors
  result
}

failed_replication_result <- function() {
  c(
    U0 = NA_real_,
    Shared_jackknife_time_seconds = NA_real_,
    JEL_valid = 0,
    JEL_cover = 0,
    JEL_len = NA_real_,
    JEL_components = NA_real_,
    JEL_multicomponent = 0,
    JEL_degenerate = 0,
    JEL_status_code = 99,
    JEL_time_seconds = NA_real_,
    Normal_valid = 0,
    Normal_cover = 0,
    Normal_len = NA_real_,
    Normal_cover_truncated = 0,
    Normal_len_truncated = NA_real_,
    Normal_near_degenerate = 0,
    Normal_status_code = 99,
    Normal_time_seconds = NA_real_,
    Kernel_valid = 0,
    Kernel_cover = 0,
    Kernel_len = NA_real_,
    Kernel_boot_valid = 0,
    Kernel_status_code = 99,
    Kernel_time_seconds = NA_real_
  )
}

## ----------------------------- scenario summary ----------------------------

summarize_scenario <- function(
    tmp,
    setting,
    design,
    true_values,
    first_data_seed,
    level
) {
  jel_valid <- tmp[, "JEL_valid"] == 1
  normal_valid <- tmp[, "Normal_valid"] == 1
  kernel_valid <- tmp[, "Kernel_valid"] == 1
  
  jel_cover_all <- tmp[, "JEL_cover"]
  normal_cover_all <- tmp[, "Normal_cover"]
  normal_cover_truncated_all <- tmp[, "Normal_cover_truncated"]
  kernel_cover_all <- tmp[, "Kernel_cover"]
  unhandled_error <-
    tmp[, "JEL_status_code"] == 99 |
    tmp[, "Normal_status_code"] == 99 |
    tmp[, "Kernel_status_code"] == 99
  
  data.frame(
    Nominal_Level = level,
    Setting_ID = setting$setting_id,
    Setting = setting$label,
    Scenario_Class = setting$scenario_class,
    Selection_Basis = setting$selection_basis,
    Design_ID = design$design_id,
    n1 = design$n1,
    n2 = design$n2,
    n3 = design$n3,
    True_VUS1 = true_values$vus1,
    True_VUS2 = true_values$vus2,
    True_theta = true_values$theta,
    R_total = nrow(tmp),
    First_data_seed = first_data_seed,
    Unhandled_error_n = sum(unhandled_error),
    Valid_for_reporting = sum(unhandled_error) == 0L,
    
    JEL_valid_n = sum(jel_valid),
    JEL_failure_percent = 100 * mean(!jel_valid),
    JEL_CP = 100 * mean(jel_cover_all),
    JEL_CP_MCSE = 100 * safe_binary_mcse(jel_cover_all),
    JEL_CP_valid = if (any(jel_valid)) {
      100 * mean(jel_cover_all[jel_valid])
    } else {
      NA_real_
    },
    JEL_AL = safe_mean(tmp[, "JEL_len"]),
    JEL_AL_MCSE = safe_mcse_mean(tmp[, "JEL_len"]),
    JEL_multicomponent_percent = if (any(jel_valid)) {
      100 * mean(tmp[jel_valid, "JEL_multicomponent"])
    } else {
      NA_real_
    },
    JEL_degenerate_percent = 100 * mean(tmp[, "JEL_degenerate"]),
    JEL_time_observed_n = sum(is.finite(tmp[, "JEL_time_seconds"])),
    JEL_total_time_seconds = sum(
      tmp[, "JEL_time_seconds"],
      na.rm = TRUE
    ),
    
    Normal_valid_n = sum(normal_valid),
    Normal_failure_percent = 100 * mean(!normal_valid),
    Normal_CP = 100 * mean(normal_cover_all),
    Normal_CP_MCSE = 100 * safe_binary_mcse(normal_cover_all),
    Normal_CP_valid = if (any(normal_valid)) {
      100 * mean(normal_cover_all[normal_valid])
    } else {
      NA_real_
    },
    Normal_AL = safe_mean(tmp[, "Normal_len"]),
    Normal_AL_MCSE = safe_mcse_mean(tmp[, "Normal_len"]),
    Normal_CP_truncated = 100 * mean(
      normal_cover_truncated_all
    ),
    Normal_CP_truncated_MCSE = 100 * safe_binary_mcse(
      normal_cover_truncated_all
    ),
    Normal_AL_truncated = safe_mean(
      tmp[, "Normal_len_truncated"]
    ),
    Normal_near_degenerate_percent = 100 * mean(
      tmp[, "Normal_near_degenerate"]
    ),
    Normal_time_observed_n = sum(is.finite(tmp[, "Normal_time_seconds"])),
    Normal_total_time_seconds = sum(
      tmp[, "Normal_time_seconds"],
      na.rm = TRUE
    ),
    
    Kernel_valid_n = sum(kernel_valid),
    Kernel_failure_percent = 100 * mean(!kernel_valid),
    Kernel_CP = 100 * mean(kernel_cover_all),
    Kernel_CP_MCSE = 100 * safe_binary_mcse(kernel_cover_all),
    Kernel_CP_valid = if (any(kernel_valid)) {
      100 * mean(kernel_cover_all[kernel_valid])
    } else {
      NA_real_
    },
    Kernel_AL = safe_mean(tmp[, "Kernel_len"]),
    Kernel_AL_MCSE = safe_mcse_mean(tmp[, "Kernel_len"]),
    Kernel_boot_valid_min = min(tmp[, "Kernel_boot_valid"]),
    Kernel_time_observed_n = sum(is.finite(tmp[, "Kernel_time_seconds"])),
    Kernel_total_time_seconds = sum(
      tmp[, "Kernel_time_seconds"],
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

## ----------------------------- full simulation -----------------------------

run_pareto_all <- function(level) {
  validate_simulation_settings()
  assert_probability(level, "level", open = TRUE)
  
  number_scenarios <- nrow(pareto_settings) * nrow(sample_sizes)
  scenario_results <- vector("list", number_scenarios)
  replication_results <- if (KEEP_REPLICATION_RESULTS) {
    vector("list", number_scenarios)
  } else {
    NULL
  }
  
  scenario_index <- 0L
  
  for (setting_index in seq_len(nrow(pareto_settings))) {
    setting <- pareto_settings[setting_index, , drop = FALSE]
    setting_id <- as.integer(setting$setting_id)
    
    true_values <- true_fgmpm_values(
      lx1 = setting$lx1,
      lx2 = setting$lx2,
      ax1 = setting$ax1,
      ax2 = setting$ax2,
      ly1 = setting$ly1,
      ly2 = setting$ly2,
      ay1 = setting$ay1,
      ay2 = setting$ay2,
      lz1 = setting$lz1,
      lz2 = setting$lz2,
      az1 = setting$az1,
      az2 = setting$az2
    )
    
    for (design_index in seq_len(nrow(sample_sizes))) {
      design <- sample_sizes[design_index, , drop = FALSE]
      design_id <- as.integer(design$design_id)
      
      n1 <- as.integer(design$n1)
      n2 <- as.integer(design$n2)
      n3 <- as.integer(design$n3)
      
      scenario_index <- scenario_index + 1L
      
      first_data_seed <- make_replication_seed(
        setting_id = setting_id,
        design_id = design_id,
        replication_id = 1L,
        stream = 1L
      )
      
      cat(sprintf(
        paste0(
          "FGM Pareto setting %d/%d, design %d/%d, ",
          "(n1,n2,n3) = (%d,%d,%d)\n"
        ),
        setting_index,
        nrow(pareto_settings),
        design_index,
        nrow(sample_sizes),
        n1,
        n2,
        n3
      ))
      
      tmp <- NULL
      data_seeds <- integer(ITERATIONS)
      bootstrap_seeds <- integer(ITERATIONS)
      replication_errors <- rep("", ITERATIONS)
      
      for (replication in seq_len(ITERATIONS)) {
        data_seed <- make_replication_seed(
          setting_id = setting_id,
          design_id = design_id,
          replication_id = replication,
          stream = 1L
        )
        
        bootstrap_seed <- make_replication_seed(
          setting_id = setting_id,
          design_id = design_id,
          replication_id = replication,
          stream = 2L
        )
        
        data_seeds[replication] <- data_seed
        bootstrap_seeds[replication] <- bootstrap_seed
        
        result <- tryCatch(
          {
            set.seed(data_seed)
            
            x <- fgmpm_sample(
              n = n1,
              lambda1 = setting$lx1,
              lambda2 = setting$lx2,
              alpha1 = setting$ax1,
              alpha2 = setting$ax2,
              theta = setting$tx
            )
            
            y <- fgmpm_sample(
              n = n2,
              lambda1 = setting$ly1,
              lambda2 = setting$ly2,
              alpha1 = setting$ay1,
              alpha2 = setting$ay2,
              theta = setting$ty
            )
            
            z <- fgmpm_sample(
              n = n3,
              lambda1 = setting$lz1,
              lambda2 = setting$lz2,
              alpha1 = setting$az1,
              alpha2 = setting$az2,
              theta = setting$tz
            )
            
            stopifnot(
              nrow(x) == n1,
              nrow(y) == n2,
              nrow(z) == n3
            )
            
            run_one_dataset(
              x = x,
              y = y,
              z = z,
              theta = true_values$theta,
              level = level,
              boot = BOOT,
              grid_len_jel = GRID_LEN_JEL,
              bootstrap_seed = bootstrap_seed
            )
          },
          error = function(error_condition) {
            if (FAIL_FAST) {
              stop(error_condition)
            }
            
            replication_errors[replication] <<-
              conditionMessage(error_condition)
            
            message(sprintf(
              paste0(
                "Recorded failure: setting %d, design %d, ",
                "replication %d: %s"
              ),
              setting_index,
              design_index,
              replication,
              replication_errors[replication]
            ))
            
            failed_replication_result()
          }
        )
        
        method_errors <- attr(result, "method_errors", exact = TRUE)
        
        if (length(method_errors) > 0L) {
          replication_errors[replication] <- paste(
            method_errors,
            collapse = "; "
          )
        }
        
        if (is.null(tmp)) {
          tmp <- matrix(
            NA_real_,
            nrow = ITERATIONS,
            ncol = length(result),
            dimnames = list(NULL, names(result))
          )
        }
        
        tmp[replication, ] <- result
      }
      
      scenario_summary <- summarize_scenario(
        tmp = tmp,
        setting = setting,
        design = design,
        true_values = true_values,
        first_data_seed = first_data_seed,
        level = level
      )
      
      scenario_results[[scenario_index]] <- scenario_summary
      
      replication_frame <- data.frame(
        Nominal_Level = level,
        Setting_ID = setting$setting_id,
        Design_ID = design$design_id,
        n1 = n1,
        n2 = n2,
        n3 = n3,
        Replication = seq_len(ITERATIONS),
        Data_seed = data_seeds,
        Bootstrap_seed = bootstrap_seeds,
        True_theta = true_values$theta,
        Replication_error = replication_errors,
        as.data.frame(tmp),
        check.names = FALSE
      )
      
      if (KEEP_REPLICATION_RESULTS) {
        replication_results[[scenario_index]] <- replication_frame
      }
      
      checkpoint_name <- sprintf(
        "FGMPM_L%02d_setting_%02d_design_%02d.rds",
        round(100 * level),
        setting_id,
        design_id
      )
      
      saveRDS(
        list(
          summary = scenario_summary,
          replications = replication_frame,
          configuration = list(
            MASTER_SEED = MASTER_SEED,
            ITERATIONS = ITERATIONS,
            BOOT = BOOT,
            LEVEL = level,
            GRID_LEN_JEL = GRID_LEN_JEL,
            setting = setting,
            design = design
          )
        ),
        file = file.path(CHECKPOINT_DIR, checkpoint_name)
      )
    }
  }
  
  list(
    summary = do.call(rbind, scenario_results),
    replications = if (KEEP_REPLICATION_RESULTS) {
      do.call(rbind, replication_results)
    } else {
      NULL
    }
  )
}

## ----------------------------- self-tests ----------------------------------

run_self_tests <- function() {
  expected_theta <- c(
    0,
    0.2916666666666667,
    -0.7009414095238095,
    -0.1013498075231214,
    -0.6472491909385113
  )
  
  calculated_theta <- numeric(nrow(pareto_settings))
  
  for (setting_index in seq_len(nrow(pareto_settings))) {
    setting <- pareto_settings[setting_index, , drop = FALSE]
    
    truth <- true_fgmpm_values(
      setting$lx1,
      setting$lx2,
      setting$ax1,
      setting$ax2,
      setting$ly1,
      setting$ly2,
      setting$ay1,
      setting$ay2,
      setting$lz1,
      setting$lz2,
      setting$az1,
      setting$az2
    )
    
    calculated_theta[setting_index] <- truth$theta
  }
  
  if (max(abs(calculated_theta - expected_theta)) > 1e-10) {
    stop("The exact Pareto VUS self-test failed.")
  }
  
  set.seed(123456L)
  x <- fgmpm_sample(7L, 1, 2, 1, 2, 0.5)
  y <- fgmpm_sample(11L, 2, 1, 2, 1, -0.4)
  z <- fgmpm_sample(13L, 1, 1, 1.5, 0.8, 0.7)
  
  jk <- jackknife_values_stratified(x, y, z)
  pooled <- pooled_jackknife_from_stratified(jk)
  lr_at_U0 <- el_neg2log(pooled$v, pooled$a, pooled$U0)
  
  if (!is.finite(lr_at_U0) || abs(lr_at_U0) > 1e-8) {
    stop("The empirical-likelihood self-test failed at U0.")
  }
  
  set.seed(MASTER_SEED)
  message("All internal self-tests passed.")
  invisible(TRUE)
}

## ----------------------------- execute simulation --------------------------

validate_simulation_settings()
run_self_tests()

start_all <- Sys.time()

simulation_outputs <- lapply(CONFIDENCE_LEVELS, run_pareto_all)
pareto_result <- do.call(
  rbind,
  lapply(simulation_outputs, function(output) output$summary)
)
replication_result <- if (KEEP_REPLICATION_RESULTS) {
  do.call(
    rbind,
    lapply(simulation_outputs, function(output) output$replications)
  )
} else {
  NULL
}

if (any(!pareto_result$Valid_for_reporting)) {
  warning(
    paste0(
      "At least one scenario contains an unhandled error (status code 99). ",
      "Do not use that scenario for inference until the error is resolved."
    ),
    call. = FALSE
  )
}

end_all <- Sys.time()
total_time <- end_all - start_all

## ----------------------------- paper-format table --------------------------

make_paper_table <- function(df) {
  data.frame(
    Nominal_Level = df$Nominal_Level,
    Scenario_Class = df$Scenario_Class,
    Parameter_setting = df$Setting,
    Sample_Size = paste0("(", df$n1, ",", df$n2, ",", df$n3, ")"),
    True_theta = round(df$True_theta, 6),
    Valid_for_reporting = df$Valid_for_reporting,
    Unhandled_error_n = df$Unhandled_error_n,
    JEL_CP_percent = round(df$JEL_CP, 2),
    JEL_CP_MCSE = round(df$JEL_CP_MCSE, 2),
    JEL_AL = round(df$JEL_AL, 4),
    JEL_failure_percent = round(df$JEL_failure_percent, 2),
    JEL_multicomponent_percent = round(
      df$JEL_multicomponent_percent,
      2
    ),
    ## The manuscript defines the raw normal interval, so it is primary here.
    Normal_CP_percent = round(df$Normal_CP, 2),
    Normal_CP_MCSE = round(df$Normal_CP_MCSE, 2),
    Normal_AL = round(df$Normal_AL, 4),
    Normal_CP_truncated_percent = round(df$Normal_CP_truncated, 2),
    Normal_AL_truncated = round(df$Normal_AL_truncated, 4),
    Normal_failure_percent = round(df$Normal_failure_percent, 2),
    Kernel_CP_percent = round(df$Kernel_CP, 2),
    Kernel_CP_MCSE = round(df$Kernel_CP_MCSE, 2),
    Kernel_AL = round(df$Kernel_AL, 4),
    Kernel_failure_percent = round(df$Kernel_failure_percent, 2),
    stringsAsFactors = FALSE
  )
}

pareto_table_90 <- make_paper_table(
  pareto_result[pareto_result$Nominal_Level == 0.90, , drop = FALSE]
)
pareto_table_95 <- make_paper_table(
  pareto_result[pareto_result$Nominal_Level == 0.95, , drop = FALSE]
)

make_calibration_comparison <- function(data) {
  nominal_percent <- 100 * data$Nominal_Level
  errors <- cbind(
    JEL = abs(data$JEL_CP - nominal_percent),
    Normal = abs(data$Normal_CP - nominal_percent),
    Kernel = abs(data$Kernel_CP - nominal_percent)
  )
  errors[!data$Valid_for_reporting, ] <- NA_real_
  closest_method <- apply(
    errors,
    1L,
    function(row) {
      if (any(!is.finite(row))) return(NA_character_)
      minimum <- min(row)
      paste(names(row)[abs(row - minimum) <= 1e-12], collapse = "=")
    }
  )

  data.frame(
    Nominal_Level = data$Nominal_Level,
    Scenario_Class = data$Scenario_Class,
    Setting_ID = data$Setting_ID,
    Design_ID = data$Design_ID,
    Sample_Size = paste0("(", data$n1, ",", data$n2, ",", data$n3, ")"),
    JEL_absolute_CP_error = round(errors[, "JEL"], 3),
    Normal_absolute_CP_error = round(errors[, "Normal"], 3),
    Kernel_absolute_CP_error = round(errors[, "Kernel"], 3),
    Closest_to_nominal = closest_method,
    JEL_closer_than_Normal = errors[, "JEL"] < errors[, "Normal"],
    JEL_closer_than_Kernel = errors[, "JEL"] < errors[, "Kernel"],
    Note = paste0(
      "Descriptive point-estimate comparison over every prespecified ",
      "scenario; no rows are selected or suppressed. Use the reported ",
      "MCSEs before interpreting small differences."
    ),
    stringsAsFactors = FALSE
  )
}

calibration_comparison <- make_calibration_comparison(pareto_result)

true_value_table <- do.call(
  rbind,
  lapply(seq_len(nrow(pareto_settings)), function(setting_index) {
    setting <- pareto_settings[setting_index, , drop = FALSE]
    
    truth <- true_fgmpm_values(
      setting$lx1,
      setting$lx2,
      setting$ax1,
      setting$ax2,
      setting$ly1,
      setting$ly2,
      setting$ay1,
      setting$ay2,
      setting$lz1,
      setting$lz2,
      setting$az1,
      setting$az2
    )
    
    data.frame(
      Setting_ID = setting$setting_id,
      Setting = setting$label,
      Scenario_Class = setting$scenario_class,
      Selection_Basis = setting$selection_basis,
      True_VUS1 = truth$vus1,
      True_VUS2 = truth$vus2,
      True_theta = truth$theta,
      stringsAsFactors = FALSE
    )
  })
)

## ----------------------------- Excel output --------------------------------

workbook <- openxlsx::createWorkbook()

openxlsx::addWorksheet(workbook, "README")
openxlsx::addWorksheet(workbook, "Table Pareto 90")
openxlsx::addWorksheet(workbook, "Table Pareto 95")
openxlsx::addWorksheet(workbook, "Calibration Comparison")
openxlsx::addWorksheet(workbook, "Scenario Summary")
openxlsx::addWorksheet(workbook, "Parameter Pareto")
openxlsx::addWorksheet(workbook, "Sample Sizes")
openxlsx::addWorksheet(workbook, "True Values")
openxlsx::addWorksheet(workbook, "Status Codes")

if (KEEP_REPLICATION_RESULTS && !is.null(replication_result)) {
  openxlsx::addWorksheet(workbook, "Replication Results")
}

readme <- data.frame(
  Item = c(
    "Master seed actually used",
    "Data seed formula",
    "Bootstrap seed formula",
    "Monte Carlo iterations",
    "Bootstrap replications",
    "Confidence levels",
    "JEL grid length",
    "Quick-test mode",
    "Coverage convention for every method",
    "Average-length convention for every method",
    "JEL multicomponent convention",
    "JEL grid sensitivity",
    "JEL calibration note",
    "Bootstrap quantile type",
    "Normal interval in main table",
    "JEL and Normal timing convention",
    "Unhandled-error reporting rule",
    "Unequal-size table note",
    "Started",
    "Finished",
    "Total time",
    "R version",
    "RNG kind",
    "openxlsx version",
    "Working directory",
    "Output file",
    "Checkpoint directory"
  ),
  Value = c(
    MASTER_SEED,
    "MASTER_SEED + 10000000*setting + 100000*design + 10*replication + 1",
    "MASTER_SEED + 10000000*setting + 100000*design + 10*replication + 2",
    ITERATIONS,
    BOOT,
    paste(CONFIDENCE_LEVELS, collapse = ", "),
    GRID_LEN_JEL,
    QUICK_TEST,
    "Primary CP is unconditional; invalid intervals count as noncoverage",
    "AL is conditional on successful interval construction",
    "Coverage uses the union; length is the total length of all detected components",
    "Grid-assisted inversion; repeat selected scenarios with a denser grid as a sensitivity check",
    "The chi-square cutoff presumes a valid Wilks theorem; verify the manuscript proof separately",
    "type = 1, matching the displayed order-statistic algorithm",
    "Raw manuscript-defined interval; truncated sensitivity results are also retained",
    "Each standalone method time includes the shared stratified-jackknife time",
    "Any scenario with status code 99 is flagged invalid for reporting and must be diagnosed",
    "These unequal-size designs replace, rather than directly reproduce, the manuscript's old balanced table",
    as.character(start_all),
    as.character(end_all),
    as.character(total_time),
    R.version.string,
    paste(RNGkind(), collapse = ", "),
    as.character(utils::packageVersion("openxlsx")),
    getwd(),
    OUTFILE,
    CHECKPOINT_DIR
  ),
  stringsAsFactors = FALSE
)

openxlsx::writeData(workbook, "README", readme)
openxlsx::writeData(workbook, "Table Pareto 90", pareto_table_90)
openxlsx::writeData(workbook, "Table Pareto 95", pareto_table_95)
openxlsx::writeData(
  workbook,
  "Calibration Comparison",
  calibration_comparison
)
openxlsx::writeData(workbook, "Scenario Summary", pareto_result)
openxlsx::writeData(workbook, "Parameter Pareto", pareto_settings)
openxlsx::writeData(workbook, "Sample Sizes", sample_sizes)
openxlsx::writeData(workbook, "True Values", true_value_table)
openxlsx::writeData(workbook, "Status Codes", STATUS_CODEBOOK)

if (KEEP_REPLICATION_RESULTS && !is.null(replication_result)) {
  openxlsx::writeData(
    workbook,
    "Replication Results",
    replication_result
  )
}

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (sheet in names(workbook)) {
  sheet_data <- switch(
    sheet,
    "README" = readme,
    "Table Pareto 90" = pareto_table_90,
    "Table Pareto 95" = pareto_table_95,
    "Calibration Comparison" = calibration_comparison,
    "Scenario Summary" = pareto_result,
    "Parameter Pareto" = pareto_settings,
    "Sample Sizes" = sample_sizes,
    "True Values" = true_value_table,
    "Status Codes" = STATUS_CODEBOOK,
    "Replication Results" = replication_result
  )
  
  number_columns <- ncol(sheet_data)
  
  openxlsx::addStyle(
    workbook,
    sheet,
    header_style,
    rows = 1,
    cols = seq_len(number_columns),
    gridExpand = TRUE
  )
  
  openxlsx::freezePane(workbook, sheet, firstRow = TRUE)
  openxlsx::setColWidths(
    workbook,
    sheet,
    cols = seq_len(number_columns),
    widths = "auto"
  )
}

openxlsx::saveWorkbook(
  workbook,
  OUTFILE,
  overwrite = TRUE
)

cat("\nSimulation completed successfully.\n")
cat("Excel output:", OUTFILE, "\n")
cat("Checkpoint directory:", CHECKPOINT_DIR, "\n")
cat("Total time:", as.character(total_time), "\n")
