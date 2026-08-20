## ============================================================================
## Corrected simulation for the difference of two dependent VUS values
## Marshall-Olkin bivariate exponential model (MOBVE)
##
## Methods:
##   1. Jackknife empirical likelihood (JEL)
##   2. Normal approximation with stratified jackknife variance
##   3. Kernel-smoothed bootstrap percentile confidence interval
##
## Run this file directly in RStudio. Set OUTPUT_DIR below if required.
##
## IMPORTANT: this code evaluates the stated procedures. The chi-square cutoff
## used for JEL still requires a correct Wilks theorem under the assumptions of
## the manuscript; simulation code cannot replace that proof.
## ============================================================================

## ----------------------------- user controls -------------------------------

OUTPUT_DIR <- getwd()

MASTER_SEED <- 2026L
ITERATIONS <- 1000L
BOOT <- 1000L
CONFIDENCE_LEVELS <- c(0.90, 0.95)
MC_COVERAGE_INTERVAL_LEVEL <- 0.95
GRID_LEN_JEL <- 401L

## Descriptive metadata for the single prespecified simulation suite. These
## labels are deliberately kept outside mobve_settings so that metadata-only
## edits do not invalidate otherwise compatible replication checkpoints.
SCENARIO_SUITE <- "Targeted nonzero-theta MOBVE grid"
SCENARIO_PROFILE <- paste0(
  "Four prespecified nonzero contrasts (0.95, -0.95, 0.50, -0.50) crossed ",
  "with six unequal-allocation sample-size designs"
)
SCENARIO_SOURCE <- "Author-defined manuscript simulation setting"
SCENARIO_TAG <- "theta_p095_m095_p050_m050"

## Set to TRUE for a small smoke test before the full run.
QUICK_TEST <- FALSE

## Retain replication-level results in the Excel workbook. They are always
## retained in scenario checkpoint files, irrespective of this option.
KEEP_REPLICATION_RESULTS <- TRUE

## Reuse a completed scenario only when its saved configuration signature
## exactly matches the current configuration.
RESUME_FROM_CHECKPOINTS <- TRUE

## During a smoke test, stop at the first unexpected error. During a full run,
## record the error as a failed replication and continue.
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
  "_S", SCENARIO_TAG,
  "_L", LEVEL_TAG,
  "_R", ITERATIONS,
  "_B", BOOT,
  "_G", GRID_LEN_JEL
)

OUTFILE <- file.path(
  OUTPUT_DIR,
  paste0("VUS_MOBVE_simulation_corrected_", RUN_TAG, ".xlsx")
)

CHECKPOINT_DIR <- file.path(
  OUTPUT_DIR,
  paste0("VUS_MOBVE_checkpoints_", RUN_TAG)
)

## ----------------------------- package and RNG setup -----------------------

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
  stop("CHECKPOINT_DIR does not exist or is not writable: ", CHECKPOINT_DIR)
}

RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

## ----------------------------- validation helpers --------------------------

assert_scalar_integer <- function(x, name, lower = 1L) {
  if (length(x) != 1L ||
      !is.finite(x) ||
      x != floor(x) ||
      x < lower) {
    stop(
      name,
      " must be a single integer greater than or equal to ",
      lower,
      "."
    )
  }
  invisible(TRUE)
}

assert_probability <- function(x, name, open = FALSE) {
  valid <- length(x) == 1L && is.finite(x)

  if (open) {
    valid <- valid && x > 0 && x < 1
  } else {
    valid <- valid && x >= 0 && x <= 1
  }

  if (!valid) {
    stop(name, if (open) " must lie in (0,1)." else " must lie in [0,1].")
  }
  invisible(TRUE)
}

safe_mean <- function(x) {
  keep <- is.finite(x)
  if (!any(keep)) return(NA_real_)
  mean(x[keep])
}

safe_median <- function(x) {
  keep <- is.finite(x)
  if (!any(keep)) return(NA_real_)
  stats::median(x[keep])
}

safe_iqr <- function(x) {
  keep <- is.finite(x)
  if (!any(keep)) return(NA_real_)
  stats::IQR(x[keep])
}

safe_quantile <- function(x, probability) {
  keep <- is.finite(x)
  if (!any(keep)) return(NA_real_)
  as.numeric(stats::quantile(
    x[keep],
    probs = probability,
    type = 8,
    names = FALSE
  ))
}

safe_mcse_mean <- function(x) {
  keep <- is.finite(x)
  if (sum(keep) <= 1L) return(NA_real_)
  stats::sd(x[keep]) / sqrt(sum(keep))
}

safe_binary_mcse <- function(x) {
  keep <- is.finite(x)
  if (!any(keep)) return(NA_real_)
  probability <- mean(x[keep])
  sqrt(probability * (1 - probability) / sum(keep))
}

safe_wilson_interval <- function(x, level = 0.95) {
  keep <- is.finite(x)
  values <- x[keep]

  if (length(values) == 0L || any(!values %in% c(0, 1))) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  number <- length(values)
  estimate <- mean(values)
  critical_value <- stats::qnorm(1 - (1 - level) / 2)
  denominator <- 1 + critical_value^2 / number
  center <- (estimate + critical_value^2 / (2 * number)) / denominator
  half_width <- critical_value / denominator * sqrt(
    estimate * (1 - estimate) / number +
      critical_value^2 / (4 * number^2)
  )

  c(
    lower = max(0, center - half_width),
    upper = min(1, center + half_width)
  )
}

collapse_intervals <- function(intervals, digits = 12L) {
  if (is.null(intervals) || nrow(intervals) == 0L) return("")

  pieces <- apply(
    intervals,
    1L,
    function(row) {
      paste0(
        "[",
        format(row[1], digits = digits, scientific = FALSE, trim = TRUE),
        ",",
        format(row[2], digits = digits, scientific = FALSE, trim = TRUE),
        "]"
      )
    }
  )

  paste(pieces, collapse = " U ")
}

## Separate deterministic data and bootstrap seeds ensure that changing BOOT
## or method order does not change the generated datasets.
make_replication_seed <- function(
    setting_id,
    design_id,
    replication_id,
    stream
) {
  assert_scalar_integer(setting_id, "setting_id")
  assert_scalar_integer(design_id, "design_id")
  assert_scalar_integer(replication_id, "replication_id")
  assert_scalar_integer(stream, "stream")

  if (setting_id > 99L) stop("At most 99 settings are supported.")
  if (design_id > 99L) stop("At most 99 designs are supported.")
  if (replication_id > 9999L) {
    stop("At most 9999 replications are supported by the seed scheme.")
  }
  if (stream > 9L) stop("stream must lie between 1 and 9.")

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

## ----------------------------- U-statistic helpers --------------------------

trip_count <- function(x, y, z) {
  if (length(x) == 0L || length(y) == 0L || length(z) == 0L ||
      anyNA(x) || anyNA(y) || anyNA(z) ||
      any(!is.finite(x)) || any(!is.finite(y)) || any(!is.finite(z))) {
    stop("trip_count requires nonempty, finite samples without missing values.")
  }

  sorted_x <- sort(x)
  sorted_z <- sort(z)

  ## Strict inequalities: number of x values below y and z values above y.
  number_x_below_y <- findInterval(y, sorted_x, left.open = TRUE)
  number_z_above_y <- length(sorted_z) - findInterval(y, sorted_z)

  sum(number_x_below_y * number_z_above_y)
}

u_stat <- function(x, y, z) {
  if (!is.matrix(x) || !is.matrix(y) || !is.matrix(z) ||
      ncol(x) != 2L || ncol(y) != 2L || ncol(z) != 2L) {
    stop("x, y, and z must be two-column matrices.")
  }

  n1 <- nrow(x)
  n2 <- nrow(y)
  n3 <- nrow(z)

  if (any(c(n1, n2, n3) < 1L)) {
    stop("All three samples must be nonempty.")
  }

  marker1 <- trip_count(x[, 1], y[, 1], z[, 1])
  marker2 <- trip_count(x[, 2], y[, 2], z[, 2])

  (marker1 - marker2) / (n1 * n2 * n3)
}

## ----------------------------- kernel estimator ----------------------------

bandwidth <- function(values, sample_size) {
  assert_scalar_integer(sample_size, "bandwidth sample size", lower = 2L)

  robust_scale <- min(
    stats::sd(values),
    stats::IQR(values) / 1.34,
    na.rm = TRUE
  )

  if (!is.finite(robust_scale) || robust_scale <= 0) {
    robust_scale <- stats::sd(values, na.rm = TRUE)
  }

  if (!is.finite(robust_scale) || robust_scale <= 0) {
    robust_scale <- 1e-8
  }

  0.9 * robust_scale * sample_size^(-0.2)
}

kernel_trip_sum <- function(x, y, z, scale_xy, scale_yz) {
  if (!is.finite(scale_xy) || !is.finite(scale_yz) ||
      scale_xy <= 0 || scale_yz <= 0) {
    stop("Kernel bandwidth combinations must be finite and positive.")
  }

  left_part <- rowSums(stats::pnorm(outer(y, x, "-") / scale_xy))
  right_part <- colSums(stats::pnorm(outer(z, y, "-") / scale_yz))

  sum(left_part * right_part)
}

kernel_stat <- function(x, y, z) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  n3 <- nrow(z)

  h_x1 <- bandwidth(x[, 1], n1)
  h_y1 <- bandwidth(y[, 1], n2)
  h_z1 <- bandwidth(z[, 1], n3)
  h_x2 <- bandwidth(x[, 2], n1)
  h_y2 <- bandwidth(y[, 2], n2)
  h_z2 <- bandwidth(z[, 2], n3)

  marker1 <- kernel_trip_sum(
    x[, 1],
    y[, 1],
    z[, 1],
    sqrt(h_x1^2 + h_y1^2),
    sqrt(h_y1^2 + h_z1^2)
  )

  marker2 <- kernel_trip_sum(
    x[, 2],
    y[, 2],
    z[, 2],
    sqrt(h_x2^2 + h_y2^2),
    sqrt(h_y2^2 + h_z2^2)
  )

  (marker1 - marker2) / (n1 * n2 * n3)
}

kernel_ci <- function(
    x0,
    y0,
    z0,
    boot = 1000L,
    level = 0.95,
    seed = NULL
) {
  assert_scalar_integer(boot, "boot", lower = 20L)
  assert_probability(level, "level", open = TRUE)

  n1 <- nrow(x0)
  n2 <- nrow(y0)
  n3 <- nrow(z0)

  if (!is.null(seed)) {
    assert_scalar_integer(seed, "bootstrap seed", lower = 0L)
    set.seed(seed)
  }

  bootstrap_statistics <- numeric(boot)

  for (bootstrap_index in seq_len(boot)) {
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

    bootstrap_statistics[bootstrap_index] <- kernel_stat(x, y, z)
  }

  finite_bootstrap <- is.finite(bootstrap_statistics)
  number_valid <- sum(finite_bootstrap)

  if (number_valid != boot) {
    return(list(
      ci = c(NA_real_, NA_real_),
      valid = FALSE,
      boot_requested = boot,
      boot_valid = number_valid,
      zero_width = FALSE,
      status = "nonfinite_bootstrap_statistics"
    ))
  }

  alpha <- 1 - level

  ## type = 1 implements the manuscript's order-statistic percentile rule.
  interval <- as.numeric(stats::quantile(
    bootstrap_statistics,
    probs = c(alpha / 2, 1 - alpha / 2),
    type = 1,
    names = FALSE,
    na.rm = FALSE
  ))

  valid_interval <- length(interval) == 2L &&
    all(is.finite(interval)) &&
    interval[1] <= interval[2]

  list(
    ci = interval,
    valid = valid_interval,
    boot_requested = boot,
    boot_valid = number_valid,
    zero_width = valid_interval && diff(interval) <= 1e-14,
    status = if (valid_interval) "ok" else "invalid_quantile"
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

  estimate <- u_stat(x, y, z)

  pseudo_x <- numeric(n1)
  pseudo_y <- numeric(n2)
  pseudo_z <- numeric(n3)

  for (index in seq_len(n1)) {
    leave_one_out <- u_stat(x[-index, , drop = FALSE], y, z)
    pseudo_x[index] <- n1 * estimate - (n1 - 1) * leave_one_out
  }

  for (index in seq_len(n2)) {
    leave_one_out <- u_stat(x, y[-index, , drop = FALSE], z)
    pseudo_y[index] <- n2 * estimate - (n2 - 1) * leave_one_out
  }

  for (index in seq_len(n3)) {
    leave_one_out <- u_stat(x, y, z[-index, , drop = FALSE])
    pseudo_z[index] <- n3 * estimate - (n3 - 1) * leave_one_out
  }

  tolerance <- 1e-8

  if (abs(mean(pseudo_x) - estimate) > tolerance ||
      abs(mean(pseudo_y) - estimate) > tolerance ||
      abs(mean(pseudo_z) - estimate) > tolerance) {
    stop("Stratified pseudo-values do not average to the U-statistic.")
  }

  list(
    U0 = estimate,
    vx = pseudo_x,
    vy = pseudo_y,
    vz = pseudo_z
  )
}

normal_ci_from_jackknife <- function(jackknife, level = 0.95) {
  assert_probability(level, "level", open = TRUE)

  n1 <- length(jackknife$vx)
  n2 <- length(jackknife$vy)
  n3 <- length(jackknife$vz)

  group_variances <- c(
    stats::var(jackknife$vx),
    stats::var(jackknife$vy),
    stats::var(jackknife$vz)
  )

  variance_estimate <-
    group_variances[1] / n1 +
    group_variances[2] / n2 +
    group_variances[3] / n3

  ## A zero estimated variance violates the nondegeneracy condition for the
  ## normal approximation. It is reported rather than silently treated as a
  ## valid zero-width interval.
  valid_variance <- is.finite(variance_estimate) && variance_estimate > 1e-14

  if (!valid_variance) {
    return(list(
      ci_raw = c(NA_real_, NA_real_),
      ci_truncated = c(NA_real_, NA_real_),
      valid = FALSE,
      variance = variance_estimate,
      group_variances = group_variances,
      near_degenerate = TRUE,
      status = "invalid_variance"
    ))
  }

  critical_value <- stats::qnorm(1 - (1 - level) / 2)
  standard_error <- sqrt(variance_estimate)

  interval_raw <- c(
    jackknife$U0 - critical_value * standard_error,
    jackknife$U0 + critical_value * standard_error
  )

  interval_truncated <- c(
    max(-1, interval_raw[1]),
    min(1, interval_raw[2])
  )

  list(
    ci_raw = interval_raw,
    ci_truncated = interval_truncated,
    valid = all(is.finite(interval_raw)) && interval_raw[1] <= interval_raw[2],
    variance = variance_estimate,
    group_variances = group_variances,
    near_degenerate = any(group_variances < 1e-12),
    status = "ok"
  )
}

## ----------------------------- pooled JEL ----------------------------------

pooled_jackknife_from_stratified <- function(jackknife) {
  n1 <- length(jackknife$vx)
  n2 <- length(jackknife$vy)
  n3 <- length(jackknife$vz)
  total_n <- n1 + n2 + n3

  if (total_n <= 3L) stop("The pooled sample size must exceed three.")

  ## Exact fixed-kernel pooled pseudo-values for a degree-(1,1,1)
  ## three-sample U-statistic. This transformation remains valid when the
  ## group sizes are unequal.
  pooled_x <- (total_n / (total_n - 3)) * (
    ((total_n - 1) / n1) * jackknife$vx - 2 * jackknife$U0
  )
  pooled_y <- (total_n / (total_n - 3)) * (
    ((total_n - 1) / n2) * jackknife$vy - 2 * jackknife$U0
  )
  pooled_z <- (total_n / (total_n - 3)) * (
    ((total_n - 1) / n3) * jackknife$vz - 2 * jackknife$U0
  )

  pooled_values <- c(pooled_x, pooled_y, pooled_z)

  coefficient <- function(group_size) {
    (total_n / (total_n - 3)) * (
      (total_n - 1) / group_size - 2
    )
  }

  centering_coefficients <- c(
    rep(coefficient(n1), n1),
    rep(coefficient(n2), n2),
    rep(coefficient(n3), n3)
  )

  if (abs(mean(pooled_values) - jackknife$U0) > 1e-8) {
    stop("Pooled pseudo-values do not average to the U-statistic.")
  }
  if (abs(mean(centering_coefficients) - 1) > 1e-10) {
    stop("Pooled centering coefficients do not average to one.")
  }

  list(
    U0 = jackknife$U0,
    v = pooled_values,
    a = centering_coefficients,
    group_a = c(
      X = coefficient(n1),
      Y = coefficient(n2),
      Z = coefficient(n3)
    )
  )
}

el_neg2log <- function(v, a, theta) {
  if (length(v) != length(a) ||
      any(!is.finite(v)) ||
      any(!is.finite(a)) ||
      !is.finite(theta)) {
    return(Inf)
  }

  estimating_values <- v - a * theta
  numerical_scale <- max(1, max(abs(v)), max(abs(a * theta)))

  ## Uniform empirical weights satisfy the constraint exactly in this case.
  if (max(abs(estimating_values)) <= 1e-12 * numerical_scale) {
    return(0)
  }

  ## Zero must lie in the interior of the scalar convex hull.
  if (!(min(estimating_values) < 0 && max(estimating_values) > 0)) {
    return(Inf)
  }

  lower_boundary <- max(-1 / estimating_values[estimating_values > 0])
  upper_boundary <- min(-1 / estimating_values[estimating_values < 0])
  interval_width <- upper_boundary - lower_boundary

  if (!is.finite(interval_width) || interval_width <= 0) return(Inf)

  inward_step <- max(
    interval_width * 1e-12,
    .Machine$double.eps * max(
      1,
      abs(lower_boundary),
      abs(upper_boundary)
    )
  )
  inward_step <- min(inward_step, interval_width / 4)

  lower <- lower_boundary + inward_step
  upper <- upper_boundary - inward_step

  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) return(Inf)

  score <- function(lambda) {
    denominator <- 1 + lambda * estimating_values
    if (any(denominator <= 0) || any(!is.finite(denominator))) {
      return(NA_real_)
    }
    sum(estimating_values / denominator)
  }

  score_at_zero <- sum(estimating_values)
  score_scale <- max(1, sum(abs(estimating_values)))

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
      error = function(error_condition) NA_real_
    )
  }

  if (!is.finite(lambda)) return(Inf)

  log_arguments <- lambda * estimating_values
  if (any(1 + log_arguments <= 0) || any(!is.finite(log_arguments))) {
    return(Inf)
  }

  value <- 2 * sum(log1p(log_arguments))
  if (!is.finite(value)) return(Inf)

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

    if (abs(inside - outside) <= tolerance) break
  }

  (inside + outside) / 2
}

jel_confidence_set_from_jackknife <- function(
    jackknife,
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

  pooled <- pooled_jackknife_from_stratified(jackknife)
  v <- pooled$v
  a <- pooled$a
  estimate <- pooled$U0
  critical_value <- stats::qchisq(level, df = 1)

  estimating_at_estimate <- v - a * estimate
  estimate_scale <- max(1, max(abs(v)), max(abs(a * estimate)))
  degenerate_at_estimate <-
    max(abs(estimating_at_estimate)) <= 1e-12 * estimate_scale

  ## The likelihood calculation itself defines LR(U0)=0 in the all-zero case,
  ## but the chi-square calibration requires nondegeneracy. We therefore report
  ## this replication as invalid instead of silently using a zero-width set.
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
      U0 = estimate,
      group_a = pooled$group_a,
      status = "degenerate_estimating_equations"
    ))
  }

  likelihood_ratio <- function(theta) el_neg2log(v, a, theta)

  is_inside <- function(theta) {
    value <- likelihood_ratio(theta)
    is.finite(value) && value <= critical_value + 1e-10
  }

  ## Inserting U0 prevents a false empty set when the accepted central
  ## component is narrower than the base grid spacing.
  grid <- sort(unique(c(
    seq(lower, upper, length.out = grid_len),
    min(max(estimate, lower), upper)
  )))

  likelihood_values <- vapply(grid, likelihood_ratio, numeric(1))
  inside <- is.finite(likelihood_values) &
    likelihood_values <= critical_value + 1e-10

  center_index <- which.min(abs(grid - estimate))
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
      degenerate = FALSE,
      lr_at_estimate = lr_at_estimate,
      U0 = estimate,
      group_a = pooled$group_a,
      status = "likelihood_not_zero_at_U0"
    ))
  }

  component_starts <- which(inside & !c(FALSE, head(inside, -1L)))
  component_ends <- which(inside & !c(tail(inside, -1L), FALSE))
  number_components <- length(component_starts)

  intervals <- matrix(
    NA_real_,
    nrow = number_components,
    ncol = 2L,
    dimnames = list(NULL, c("lower", "upper"))
  )

  for (component in seq_len(number_components)) {
    start_index <- component_starts[component]
    end_index <- component_ends[component]

    lower_endpoint <- if (start_index == 1L) {
      lower
    } else {
      refine_inside_outside_boundary(
        outside_point = grid[start_index - 1L],
        inside_point = grid[start_index],
        is_inside = is_inside
      )
    }

    upper_endpoint <- if (end_index == length(grid)) {
      upper
    } else {
      refine_inside_outside_boundary(
        outside_point = grid[end_index + 1L],
        inside_point = grid[end_index],
        is_inside = is_inside
      )
    }

    intervals[component, ] <- c(lower_endpoint, upper_endpoint)
  }

  central_component <- which(
    component_starts <= center_index & component_ends >= center_index
  )

  if (length(central_component) != 1L) {
    return(list(
      intervals = intervals,
      central_interval = c(NA_real_, NA_real_),
      hull = c(intervals[1, 1], intervals[nrow(intervals), 2]),
      total_length = sum(intervals[, 2] - intervals[, 1]),
      valid = FALSE,
      components = number_components,
      degenerate = FALSE,
      lr_at_estimate = lr_at_estimate,
      U0 = estimate,
      group_a = pooled$group_a,
      status = "central_component_not_unique"
    ))
  }

  central_interval <- as.numeric(intervals[central_component, ])
  valid_set <- all(is.finite(intervals)) &&
    all(intervals[, 1] <= intervals[, 2])

  list(
    intervals = intervals,
    central_interval = central_interval,
    hull = c(intervals[1, 1], intervals[nrow(intervals), 2]),
    total_length = sum(intervals[, 2] - intervals[, 1]),
    valid = valid_set,
    components = number_components,
    degenerate = FALSE,
    lr_at_estimate = lr_at_estimate,
    U0 = estimate,
    group_a = pooled$group_a,
    status = if (valid_set) "ok" else "invalid_interval_set"
  )
}

value_in_intervals <- function(value, intervals, tolerance = 1e-12) {
  if (!is.finite(value) || is.null(intervals) || nrow(intervals) == 0L) {
    return(FALSE)
  }

  any(
    value >= intervals[, 1] - tolerance &
      value <= intervals[, 2] + tolerance
  )
}

## ----------------------------- MOBVE model ---------------------------------

exponential_with_zero_rate <- function(n, rate) {
  assert_scalar_integer(n, "n")

  if (length(rate) != 1L || !is.finite(rate) || rate < 0) {
    stop("Every exponential rate must be a finite nonnegative scalar.")
  }

  if (rate == 0) return(rep(Inf, n))
  stats::rexp(n, rate = rate)
}

mobve_sample <- function(n, lambda1, lambda2, lambda3) {
  assert_scalar_integer(n, "n")

  rates <- c(lambda1, lambda2, lambda3)
  if (any(!is.finite(rates)) || any(rates < 0)) {
    stop("MOBVE rates must be finite and nonnegative.")
  }
  if (lambda1 + lambda3 <= 0 || lambda2 + lambda3 <= 0) {
    stop("Both MOBVE marginal exponential rates must be positive.")
  }
  independent_shock_1 <- exponential_with_zero_rate(n, lambda1)
  independent_shock_2 <- exponential_with_zero_rate(n, lambda2)
  common_shock <- exponential_with_zero_rate(n, lambda3)

  sample_matrix <- cbind(
    pmin(independent_shock_1, common_shock),
    pmin(independent_shock_2, common_shock)
  )

  if (any(!is.finite(sample_matrix))) {
    stop("MOBVE generation produced a nonfinite observation.")
  }

  colnames(sample_matrix) <- c("marker1", "marker2")
  sample_matrix
}

vus_exponential <- function(rate_x, rate_y, rate_z) {
  rates <- c(rate_x, rate_y, rate_z)
  if (any(!is.finite(rates)) || any(rates <= 0)) {
    stop("Marginal exponential rates must be finite and positive.")
  }

  ## P(X < Y < Z) for independent exponential variables from the three
  ## samples. Dependence between the two markers within a sample does not alter
  ## either marginal VUS; it affects their covariance and hence interval width.
  (rate_x * rate_y) /
    ((rate_z + rate_y) * (rate_z + rate_y + rate_x))
}

true_mobve_values <- function(
    lx1,
    lx2,
    lx3,
    ly1,
    ly2,
    ly3,
    lz1,
    lz2,
    lz3
) {
  vus1 <- vus_exponential(
    rate_x = lx1 + lx3,
    rate_y = ly1 + ly3,
    rate_z = lz1 + lz3
  )

  vus2 <- vus_exponential(
    rate_x = lx2 + lx3,
    rate_y = ly2 + ly3,
    rate_z = lz2 + lz3
  )

  list(
    vus1 = vus1,
    vus2 = vus2,
    theta = vus1 - vus2
  )
}

## ----------------------------- simulation settings -------------------------

## For fixed marginal rates r_Y and r_Z, solve
##   P(X < Y < Z) = r_X r_Y / {(r_Y+r_Z)(r_X+r_Y+r_Z)}
## for r_X. Finite positive exponential rates imply VUS in (0,1), so theta
## cannot equal -1 or 1 exactly. The values +/-0.95 are used as prespecified
## near-boundary contrasts, while +/-0.50 are moderate nonzero contrasts.
rate_x_for_target_vus <- function(target_vus, rate_y, rate_z) {
  denominator <- rate_y - target_vus * (rate_y + rate_z)
  if (denominator <= 0) {
    stop("The requested VUS is infeasible for the selected Y and Z rates.")
  }
  target_vus * (rate_y + rate_z)^2 / denominator
}

extreme_x_high <- rate_x_for_target_vus(0.96, 1.00, 0.01)
extreme_x_low <- rate_x_for_target_vus(0.01, 1.00, 0.01)
moderate_x_high <- rate_x_for_target_vus(0.60, 1.00, 0.10)
moderate_x_low <- rate_x_for_target_vus(0.10, 1.00, 0.10)

mobve_settings <- data.frame(
  setting_id = 1:4,
  label = c(
    "theta=0.95 (VUS1=0.96, VUS2=0.01)",
    "theta=-0.95 (VUS1=0.01, VUS2=0.96)",
    "theta=0.50 (VUS1=0.60, VUS2=0.10)",
    "theta=-0.50 (VUS1=0.10, VUS2=0.60)"
  ),
  target_vus1 = c(0.96, 0.01, 0.60, 0.10),
  target_vus2 = c(0.01, 0.96, 0.10, 0.60),
  target_theta = c(0.95, -0.95, 0.50, -0.50),
  lx1 = c(
    extreme_x_high - 0.009,
    extreme_x_low - 0.009,
    moderate_x_high - 0.100,
    moderate_x_low - 0.100
  ),
  lx2 = c(
    extreme_x_low - 0.009,
    extreme_x_high - 0.009,
    moderate_x_low - 0.100,
    moderate_x_high - 0.100
  ),
  lx3 = c(0.009, 0.009, 0.100, 0.100),
  ly1 = c(0.500, 0.500, 0.500, 0.500),
  ly2 = c(0.500, 0.500, 0.500, 0.500),
  ly3 = c(0.500, 0.500, 0.500, 0.500),
  lz1 = c(0.005, 0.005, 0.050, 0.050),
  lz2 = c(0.005, 0.005, 0.050, 0.050),
  lz3 = c(0.005, 0.005, 0.050, 0.050),
  stringsAsFactors = FALSE
)

mobve_settings$common_shock_probability_X <- with(
  mobve_settings,
  lx3 / (lx1 + lx2 + lx3)
)
mobve_settings$common_shock_probability_Y <- with(
  mobve_settings,
  ly3 / (ly1 + ly2 + ly3)
)
mobve_settings$common_shock_probability_Z <- with(
  mobve_settings,
  lz3 / (lz1 + lz2 + lz3)
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
  assert_probability(
    MC_COVERAGE_INTERVAL_LEVEL,
    "MC_COVERAGE_INTERVAL_LEVEL",
    open = TRUE
  )
  assert_scalar_integer(GRID_LEN_JEL, "GRID_LEN_JEL", lower = 101L)

  if (ITERATIONS > 9999L) {
    stop("ITERATIONS must not exceed 9999 under the seed scheme.")
  }

  if (anyDuplicated(mobve_settings$setting_id) ||
      anyDuplicated(sample_sizes$design_id)) {
    stop("setting_id and design_id values must be unique.")
  }

  target_values <- as.matrix(mobve_settings[, c(
    "target_vus1", "target_vus2", "target_theta"
  )])
  if (any(!is.finite(target_values)) ||
      any(target_values[, "target_vus1"] <= 0) ||
      any(target_values[, "target_vus1"] >= 1) ||
      any(target_values[, "target_vus2"] <= 0) ||
      any(target_values[, "target_vus2"] >= 1) ||
      any(abs(target_values[, "target_theta"]) >= 1) ||
      any(target_values[, "target_theta"] == 0) ||
      any(abs(
        target_values[, "target_vus1"] -
          target_values[, "target_vus2"] -
          target_values[, "target_theta"]
      ) > 1e-12)) {
    stop("The prespecified VUS targets or nonzero theta targets are invalid.")
  }

  invisible(lapply(
    mobve_settings$setting_id,
    assert_scalar_integer,
    name = "setting_id"
  ))
  invisible(lapply(
    sample_sizes$design_id,
    assert_scalar_integer,
    name = "design_id"
  ))

  size_matrix <- as.matrix(sample_sizes[, c("n1", "n2", "n3")])
  if (any(!is.finite(size_matrix)) ||
      any(size_matrix <= 1) ||
      any(size_matrix != floor(size_matrix))) {
    stop("Every sample size must be an integer greater than one.")
  }

  rate_columns <- c(
    "lx1", "lx2", "lx3",
    "ly1", "ly2", "ly3",
    "lz1", "lz2", "lz3"
  )
  rate_values <- unlist(mobve_settings[, rate_columns], use.names = FALSE)
  if (any(!is.finite(rate_values)) || any(rate_values < 0)) {
    stop("Every MOBVE rate must be finite and nonnegative.")
  }

  marginal_rates <- c(
    mobve_settings$lx1 + mobve_settings$lx3,
    mobve_settings$lx2 + mobve_settings$lx3,
    mobve_settings$ly1 + mobve_settings$ly3,
    mobve_settings$ly2 + mobve_settings$ly3,
    mobve_settings$lz1 + mobve_settings$lz3,
    mobve_settings$lz2 + mobve_settings$lz3
  )
  if (any(marginal_rates <= 0)) {
    stop("Every MOBVE marginal exponential rate must be positive.")
  }

  invisible(TRUE)
}

STATUS_CODEBOOK <- data.frame(
  Method = c(
    rep("JEL", 6L),
    rep("Normal", 3L),
    rep("Kernel", 4L)
  ),
  Status = c(
    "ok",
    "degenerate_estimating_equations",
    "likelihood_not_zero_at_U0",
    "central_component_not_unique",
    "invalid_interval_set",
    "unhandled_replication_error",
    "ok",
    "invalid_variance",
    "unhandled_replication_error",
    "ok",
    "nonfinite_bootstrap_statistics",
    "invalid_quantile",
    "unhandled_replication_error"
  ),
  Meaning = c(
    "Valid union confidence set",
    "All JEL estimating equations vanish at the estimate; nondegeneracy fails",
    "Numerical likelihood check failed at the U-statistic",
    "A unique accepted component containing the U-statistic was not identified",
    "Returned JEL set contains an invalid endpoint",
    "Unexpected JEL or shared-jackknife error was caught",
    "Valid normal interval",
    "Jackknife variance is nonfinite or numerically zero",
    "Unexpected Normal or shared-jackknife error was caught",
    "Valid kernel percentile interval",
    "At least one requested bootstrap statistic was nonfinite",
    "Bootstrap quantile endpoints were invalid",
    "Unexpected Kernel error was caught"
  ),
  stringsAsFactors = FALSE
)

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

  timed_method <- function(method_name, expression, fallback) {
    start_time <- proc.time()[["elapsed"]]

    value <- tryCatch(
      force(expression),
      error = function(error_condition) {
        if (FAIL_FAST) stop(error_condition)

        method_errors <<- c(
          method_errors,
          paste0(method_name, ": ", conditionMessage(error_condition))
        )
        fallback
      }
    )

    list(
      value = value,
      elapsed = proc.time()[["elapsed"]] - start_time
    )
  }

  jel_fallback <- list(
    intervals = matrix(numeric(0), ncol = 2L),
    central_interval = c(NA_real_, NA_real_),
    hull = c(NA_real_, NA_real_),
    total_length = NA_real_,
    valid = FALSE,
    components = 0L,
    degenerate = FALSE,
    group_a = c(X = NA_real_, Y = NA_real_, Z = NA_real_),
    status = "unhandled_replication_error"
  )

  normal_fallback <- list(
    ci_raw = c(NA_real_, NA_real_),
    ci_truncated = c(NA_real_, NA_real_),
    valid = FALSE,
    variance = NA_real_,
    group_variances = c(NA_real_, NA_real_, NA_real_),
    near_degenerate = FALSE,
    status = "unhandled_replication_error"
  )

  kernel_fallback <- list(
    ci = c(NA_real_, NA_real_),
    valid = FALSE,
    boot_requested = boot,
    boot_valid = 0L,
    zero_width = FALSE,
    status = "unhandled_replication_error"
  )

  ## JEL and Normal share the same stratified pseudo-values. Kernel remains
  ## independent and is still attempted if this shared calculation fails.
  jackknife_start <- proc.time()[["elapsed"]]
  jackknife <- tryCatch(
    jackknife_values_stratified(x, y, z),
    error = function(error_condition) {
      if (FAIL_FAST) stop(error_condition)

      method_errors <<- c(
        method_errors,
        paste0("Shared Jackknife: ", conditionMessage(error_condition))
      )
      NULL
    }
  )
  jackknife_time <- proc.time()[["elapsed"]] - jackknife_start

  jel_run <- if (is.null(jackknife)) {
    list(value = jel_fallback, elapsed = 0)
  } else {
    timed_method(
      "JEL",
      jel_confidence_set_from_jackknife(
        jackknife = jackknife,
        level = level,
        lower = -1,
        upper = 1,
        grid_len = grid_len_jel
      ),
      jel_fallback
    )
  }
  jel <- jel_run$value

  normal_run <- if (is.null(jackknife)) {
    list(value = normal_fallback, elapsed = 0)
  } else {
    timed_method(
      "Normal",
      normal_ci_from_jackknife(jackknife, level = level),
      normal_fallback
    )
  }
  normal <- normal_run$value

  kernel_run <- timed_method(
    "Kernel",
    kernel_ci(
      x0 = x,
      y0 = y,
      z0 = z,
      boot = boot,
      level = level,
      seed = bootstrap_seed
    ),
    kernel_fallback
  )
  kernel <- kernel_run$value

  jel_valid <- isTRUE(jel$valid)
  normal_valid <- isTRUE(normal$valid)
  kernel_valid <- isTRUE(kernel$valid)

  jel_central <- jel$central_interval
  jel_hull <- jel$hull
  normal_raw <- normal$ci_raw
  normal_truncated <- normal$ci_truncated
  kernel_interval <- kernel$ci

  data.frame(
    U0 = if (is.null(jackknife)) NA_real_ else jackknife$U0,
    Shared_jackknife_time_seconds = jackknife_time,

    a_X = as.numeric(jel$group_a["X"]),
    a_Y = as.numeric(jel$group_a["Y"]),
    a_Z = as.numeric(jel$group_a["Z"]),

    JEL_valid = jel_valid,
    JEL_cover = if (jel_valid) {
      as.numeric(value_in_intervals(theta, jel$intervals))
    } else {
      0
    },
    JEL_len = if (jel_valid) jel$total_length else NA_real_,
    JEL_central_lower = if (jel_valid) jel_central[1] else NA_real_,
    JEL_central_upper = if (jel_valid) jel_central[2] else NA_real_,
    JEL_hull_lower = if (jel_valid) jel_hull[1] else NA_real_,
    JEL_hull_upper = if (jel_valid) jel_hull[2] else NA_real_,
    JEL_components = if (jel_valid) jel$components else NA_integer_,
    JEL_multicomponent = jel_valid && jel$components > 1L,
    JEL_degenerate = isTRUE(jel$degenerate),
    JEL_status = as.character(jel$status),
    JEL_intervals = if (jel_valid) {
      collapse_intervals(jel$intervals)
    } else {
      "not_available"
    },
    JEL_standalone_time_seconds = jackknife_time + jel_run$elapsed,

    Normal_valid = normal_valid,
    Normal_cover = if (normal_valid) {
      as.numeric(theta >= normal_raw[1] && theta <= normal_raw[2])
    } else {
      0
    },
    Normal_len = if (normal_valid) diff(normal_raw) else NA_real_,
    Normal_lower = if (normal_valid) normal_raw[1] else NA_real_,
    Normal_upper = if (normal_valid) normal_raw[2] else NA_real_,
    Normal_cover_truncated = if (normal_valid) {
      as.numeric(
        theta >= normal_truncated[1] && theta <= normal_truncated[2]
      )
    } else {
      0
    },
    Normal_len_truncated = if (normal_valid) {
      diff(normal_truncated)
    } else {
      NA_real_
    },
    Normal_lower_truncated = if (normal_valid) {
      normal_truncated[1]
    } else {
      NA_real_
    },
    Normal_upper_truncated = if (normal_valid) {
      normal_truncated[2]
    } else {
      NA_real_
    },
    Normal_pseudovariance_X = normal$group_variances[1],
    Normal_pseudovariance_Y = normal$group_variances[2],
    Normal_pseudovariance_Z = normal$group_variances[3],
    Normal_variance = normal$variance,
    Normal_numerically_tiny_group_variance = isTRUE(
      normal$near_degenerate
    ),
    Normal_status = as.character(normal$status),
    Normal_standalone_time_seconds = jackknife_time + normal_run$elapsed,

    Kernel_valid = kernel_valid,
    Kernel_cover = if (kernel_valid) {
      as.numeric(
        theta >= kernel_interval[1] && theta <= kernel_interval[2]
      )
    } else {
      0
    },
    Kernel_len = if (kernel_valid) diff(kernel_interval) else NA_real_,
    Kernel_lower = if (kernel_valid) kernel_interval[1] else NA_real_,
    Kernel_upper = if (kernel_valid) kernel_interval[2] else NA_real_,
    Kernel_boot_requested = kernel$boot_requested,
    Kernel_boot_valid = kernel$boot_valid,
    Kernel_zero_width = isTRUE(kernel$zero_width),
    Kernel_status = as.character(kernel$status),
    Kernel_standalone_time_seconds = kernel_run$elapsed,

    Method_error = if (length(method_errors) == 0L) {
      "none"
    } else {
      paste(method_errors, collapse = "; ")
    },
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

failed_replication_result <- function(error_message = "") {
  data.frame(
    U0 = NA_real_,
    Shared_jackknife_time_seconds = NA_real_,
    a_X = NA_real_,
    a_Y = NA_real_,
    a_Z = NA_real_,
    JEL_valid = FALSE,
    JEL_cover = 0,
    JEL_len = NA_real_,
    JEL_central_lower = NA_real_,
    JEL_central_upper = NA_real_,
    JEL_hull_lower = NA_real_,
    JEL_hull_upper = NA_real_,
    JEL_components = NA_integer_,
    JEL_multicomponent = FALSE,
    JEL_degenerate = FALSE,
    JEL_status = "unhandled_replication_error",
    JEL_intervals = "not_available",
    JEL_standalone_time_seconds = NA_real_,
    Normal_valid = FALSE,
    Normal_cover = 0,
    Normal_len = NA_real_,
    Normal_lower = NA_real_,
    Normal_upper = NA_real_,
    Normal_cover_truncated = 0,
    Normal_len_truncated = NA_real_,
    Normal_lower_truncated = NA_real_,
    Normal_upper_truncated = NA_real_,
    Normal_pseudovariance_X = NA_real_,
    Normal_pseudovariance_Y = NA_real_,
    Normal_pseudovariance_Z = NA_real_,
    Normal_variance = NA_real_,
    Normal_numerically_tiny_group_variance = FALSE,
    Normal_status = "unhandled_replication_error",
    Normal_standalone_time_seconds = NA_real_,
    Kernel_valid = FALSE,
    Kernel_cover = 0,
    Kernel_len = NA_real_,
    Kernel_lower = NA_real_,
    Kernel_upper = NA_real_,
    Kernel_boot_requested = BOOT,
    Kernel_boot_valid = 0L,
    Kernel_zero_width = FALSE,
    Kernel_status = "unhandled_replication_error",
    Kernel_standalone_time_seconds = NA_real_,
    Method_error = as.character(error_message),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

## ----------------------------- scenario summary ----------------------------

summarize_scenario <- function(
    replication_data,
    setting,
    design,
    true_values,
    first_data_seed,
    level
) {
  jel_valid <- replication_data$JEL_valid
  normal_valid <- replication_data$Normal_valid
  kernel_valid <- replication_data$Kernel_valid

  unhandled_error <-
    replication_data$JEL_status == "unhandled_replication_error" |
    replication_data$Normal_status == "unhandled_replication_error" |
    replication_data$Kernel_status == "unhandled_replication_error"

  jel_wilson <- safe_wilson_interval(
    replication_data$JEL_cover,
    MC_COVERAGE_INTERVAL_LEVEL
  )
  normal_wilson <- safe_wilson_interval(
    replication_data$Normal_cover,
    MC_COVERAGE_INTERVAL_LEVEL
  )
  normal_truncated_wilson <- safe_wilson_interval(
    replication_data$Normal_cover_truncated,
    MC_COVERAGE_INTERVAL_LEVEL
  )
  kernel_wilson <- safe_wilson_interval(
    replication_data$Kernel_cover,
    MC_COVERAGE_INTERVAL_LEVEL
  )

  shared_variance <- replication_data$Normal_variance
  numerically_tiny_total_variance <- is.finite(shared_variance) &
    shared_variance < 1e-12

  total_n <- design$n1 + design$n2 + design$n3
  coefficient <- function(group_size) {
    (total_n / (total_n - 3)) * ((total_n - 1) / group_size - 2)
  }

  data.frame(
    Nominal_Level = level,
    Scenario_Suite = SCENARIO_SUITE,
    Setting_ID = setting$setting_id,
    Setting = setting$label,
    Target_VUS1 = setting$target_vus1,
    Target_VUS2 = setting$target_vus2,
    Target_theta = setting$target_theta,
    Design_ID = design$design_id,
    n1 = design$n1,
    n2 = design$n2,
    n3 = design$n3,
    Total_n = total_n,
    n1_fraction = design$n1 / total_n,
    n2_fraction = design$n2 / total_n,
    n3_fraction = design$n3 / total_n,
    a_X = coefficient(design$n1),
    a_Y = coefficient(design$n2),
    a_Z = coefficient(design$n3),
    True_VUS1 = true_values$vus1,
    True_VUS2 = true_values$vus2,
    True_theta = true_values$theta,
    Common_shock_probability_X = setting$common_shock_probability_X,
    Common_shock_probability_Y = setting$common_shock_probability_Y,
    Common_shock_probability_Z = setting$common_shock_probability_Z,
    R_total = nrow(replication_data),
    B_requested = BOOT,
    First_data_seed = first_data_seed,
    Unhandled_error_n = sum(unhandled_error),
    No_unhandled_errors = sum(unhandled_error) == 0L,
    Shared_variance_min = safe_quantile(shared_variance, 0),
    Shared_variance_q05 = safe_quantile(shared_variance, 0.05),
    Shared_variance_median = safe_median(shared_variance),
    Numerically_tiny_total_variance_percent = if (
      any(is.finite(shared_variance))
    ) {
      100 * mean(
        numerically_tiny_total_variance[is.finite(shared_variance)]
      )
    } else {
      NA_real_
    },

    JEL_valid_n = sum(jel_valid),
    JEL_failure_percent = 100 * mean(!jel_valid),
    JEL_CP = 100 * mean(replication_data$JEL_cover),
    JEL_CP_MCSE = 100 * safe_binary_mcse(replication_data$JEL_cover),
    JEL_CP_Wilson_lower = 100 * jel_wilson["lower"],
    JEL_CP_Wilson_upper = 100 * jel_wilson["upper"],
    JEL_CP_valid = if (any(jel_valid)) {
      100 * mean(replication_data$JEL_cover[jel_valid])
    } else {
      NA_real_
    },
    JEL_AL = safe_mean(replication_data$JEL_len),
    JEL_AL_MCSE = safe_mcse_mean(replication_data$JEL_len),
    JEL_length_median = safe_median(replication_data$JEL_len),
    JEL_length_IQR = safe_iqr(replication_data$JEL_len),
    JEL_multicomponent_percent = if (any(jel_valid)) {
      100 * mean(replication_data$JEL_multicomponent[jel_valid])
    } else {
      NA_real_
    },
    JEL_degenerate_percent = 100 * mean(replication_data$JEL_degenerate),
    JEL_time_observed_n = sum(
      is.finite(replication_data$JEL_standalone_time_seconds)
    ),
    JEL_total_standalone_time_seconds = sum(
      replication_data$JEL_standalone_time_seconds,
      na.rm = TRUE
    ),
    JEL_mean_standalone_time_seconds = safe_mean(
      replication_data$JEL_standalone_time_seconds
    ),
    JEL_median_standalone_time_seconds = safe_median(
      replication_data$JEL_standalone_time_seconds
    ),
    JEL_IQR_standalone_time_seconds = safe_iqr(
      replication_data$JEL_standalone_time_seconds
    ),

    Normal_valid_n = sum(normal_valid),
    Normal_failure_percent = 100 * mean(!normal_valid),
    Normal_CP = 100 * mean(replication_data$Normal_cover),
    Normal_CP_MCSE = 100 * safe_binary_mcse(replication_data$Normal_cover),
    Normal_CP_Wilson_lower = 100 * normal_wilson["lower"],
    Normal_CP_Wilson_upper = 100 * normal_wilson["upper"],
    Normal_CP_valid = if (any(normal_valid)) {
      100 * mean(replication_data$Normal_cover[normal_valid])
    } else {
      NA_real_
    },
    Normal_AL = safe_mean(replication_data$Normal_len),
    Normal_AL_MCSE = safe_mcse_mean(replication_data$Normal_len),
    Normal_length_median = safe_median(replication_data$Normal_len),
    Normal_length_IQR = safe_iqr(replication_data$Normal_len),
    Normal_CP_truncated = 100 * mean(
      replication_data$Normal_cover_truncated
    ),
    Normal_CP_truncated_MCSE = 100 * safe_binary_mcse(
      replication_data$Normal_cover_truncated
    ),
    Normal_CP_truncated_Wilson_lower =
      100 * normal_truncated_wilson["lower"],
    Normal_CP_truncated_Wilson_upper =
      100 * normal_truncated_wilson["upper"],
    Normal_AL_truncated = safe_mean(
      replication_data$Normal_len_truncated
    ),
    Normal_AL_truncated_MCSE = safe_mcse_mean(
      replication_data$Normal_len_truncated
    ),
    Normal_numerically_tiny_group_variance_percent = 100 * mean(
      replication_data$Normal_numerically_tiny_group_variance
    ),
    Normal_time_observed_n = sum(
      is.finite(replication_data$Normal_standalone_time_seconds)
    ),
    Normal_total_standalone_time_seconds = sum(
      replication_data$Normal_standalone_time_seconds,
      na.rm = TRUE
    ),
    Normal_mean_standalone_time_seconds = safe_mean(
      replication_data$Normal_standalone_time_seconds
    ),
    Normal_median_standalone_time_seconds = safe_median(
      replication_data$Normal_standalone_time_seconds
    ),
    Normal_IQR_standalone_time_seconds = safe_iqr(
      replication_data$Normal_standalone_time_seconds
    ),

    Kernel_valid_n = sum(kernel_valid),
    Kernel_failure_percent = 100 * mean(!kernel_valid),
    Kernel_CP = 100 * mean(replication_data$Kernel_cover),
    Kernel_CP_MCSE = 100 * safe_binary_mcse(replication_data$Kernel_cover),
    Kernel_CP_Wilson_lower = 100 * kernel_wilson["lower"],
    Kernel_CP_Wilson_upper = 100 * kernel_wilson["upper"],
    Kernel_CP_valid = if (any(kernel_valid)) {
      100 * mean(replication_data$Kernel_cover[kernel_valid])
    } else {
      NA_real_
    },
    Kernel_AL = safe_mean(replication_data$Kernel_len),
    Kernel_AL_MCSE = safe_mcse_mean(replication_data$Kernel_len),
    Kernel_length_median = safe_median(replication_data$Kernel_len),
    Kernel_length_IQR = safe_iqr(replication_data$Kernel_len),
    Kernel_boot_valid_min = min(replication_data$Kernel_boot_valid),
    Kernel_zero_width_percent = 100 * mean(
      replication_data$Kernel_zero_width
    ),
    Kernel_time_observed_n = sum(
      is.finite(replication_data$Kernel_standalone_time_seconds)
    ),
    Kernel_total_standalone_time_seconds = sum(
      replication_data$Kernel_standalone_time_seconds,
      na.rm = TRUE
    ),
    Kernel_mean_standalone_time_seconds = safe_mean(
      replication_data$Kernel_standalone_time_seconds
    ),
    Kernel_median_standalone_time_seconds = safe_median(
      replication_data$Kernel_standalone_time_seconds
    ),
    Kernel_IQR_standalone_time_seconds = safe_iqr(
      replication_data$Kernel_standalone_time_seconds
    ),
    stringsAsFactors = FALSE
  )
}

## ----------------------------- checkpoint helpers --------------------------

current_script_path <- function() {
  command_line <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_line, value = TRUE)

  if (length(file_argument) > 0L) {
    candidate <- sub("^--file=", "", file_argument[1])
    if (file.exists(candidate)) return(normalizePath(candidate))
  }

  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) as.character(frame$ofile)[1] else ""
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]

  if (length(frame_files) > 0L && file.exists(tail(frame_files, 1L))) {
    return(normalizePath(tail(frame_files, 1L)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    candidate <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(error_condition) ""
    )
    if (nzchar(candidate) && file.exists(candidate)) {
      return(normalizePath(candidate))
    }
  }

  NA_character_
}

compute_algorithm_signature <- function() {
  signature_file <- tempfile(fileext = ".rds")
  on.exit(unlink(signature_file), add = TRUE)

  saveRDS(
    list(
      trip_count = body(trip_count),
      u_stat = body(u_stat),
      bandwidth = body(bandwidth),
      kernel_trip_sum = body(kernel_trip_sum),
      kernel_stat = body(kernel_stat),
      kernel_ci = body(kernel_ci),
      jackknife = body(jackknife_values_stratified),
      normal = body(normal_ci_from_jackknife),
      pooled = body(pooled_jackknife_from_stratified),
      el = body(el_neg2log),
      boundary_refinement = body(refine_inside_outside_boundary),
      jel = body(jel_confidence_set_from_jackknife),
      zero_rate_exponential = body(exponential_with_zero_rate),
      generator = body(mobve_sample),
      truth = body(true_mobve_values),
      one_replication = body(run_one_dataset),
      settings = mobve_settings,
      sample_sizes = sample_sizes
    ),
    signature_file,
    version = 2
  )

  unname(tools::md5sum(signature_file))
}

atomic_save_rds <- function(object, path) {
  temporary_path <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary_path), add = TRUE)

  saveRDS(object, temporary_path)

  backup_path <- paste0(path, ".bak")
  if (file.exists(backup_path) && unlink(backup_path) != 0L) {
    stop("Could not remove a stale checkpoint backup: ", backup_path)
  }

  had_previous <- file.exists(path)
  if (had_previous && !file.rename(path, backup_path)) {
    stop("Could not move the previous checkpoint to backup: ", path)
  }

  if (!file.rename(temporary_path, path)) {
    if (had_previous && file.exists(backup_path)) {
      file.rename(backup_path, path)
    }
    stop("Could not commit checkpoint: ", path)
  }

  if (file.exists(backup_path)) unlink(backup_path)
  invisible(TRUE)
}

promote_file_with_backup <- function(source, destination) {
  if (!file.exists(source)) stop("Source file does not exist: ", source)

  backup <- paste0(destination, ".bak")
  if (file.exists(backup) && unlink(backup) != 0L) {
    stop("Could not remove stale output backup: ", backup)
  }

  had_previous <- file.exists(destination)
  if (had_previous && !file.rename(destination, backup)) {
    stop("Could not back up the previous output file: ", destination)
  }

  if (!file.rename(source, destination)) {
    if (had_previous && file.exists(backup)) {
      file.rename(backup, destination)
    }
    stop("Could not promote the new output file: ", destination)
  }

  if (file.exists(backup)) unlink(backup)
  invisible(TRUE)
}

read_checkpoint_safely <- function(path) {
  candidates <- c(path, paste0(path, ".bak"))
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0L) return(NULL)

  errors <- character(0)
  for (candidate in candidates) {
    value <- tryCatch(
      readRDS(candidate),
      error = function(error_condition) {
        errors <<- c(
          errors,
          paste0(candidate, ": ", conditionMessage(error_condition))
        )
        NULL
      }
    )

    if (!is.null(value)) return(value)
  }

  stop("No readable checkpoint was found. ", paste(errors, collapse = "; "))
}

scenario_configuration <- function(setting, design, level) {
  list(
    version = "MOBVE-corrected-2026-08-19-v5",
    script_md5 = SCRIPT_MD5_AT_START,
    algorithm_signature = ALGORITHM_SIGNATURE,
    r_version = R.version.string,
    openxlsx_version = as.character(utils::packageVersion("openxlsx")),
    operating_system = paste(
      unname(Sys.info()[c("sysname", "release")]),
      collapse = " "
    ),
    machine = unname(as.character(Sys.info()["machine"])),
    logical_cores = as.integer(parallel::detectCores(logical = TRUE)),
    rng_kind = RNGkind(),
    master_seed = MASTER_SEED,
    iterations = ITERATIONS,
    bootstrap_replications = BOOT,
    confidence_level = level,
    jel_grid_length = GRID_LEN_JEL,
    setting = as.list(setting),
    design = as.list(design)
  )
}

validate_checkpoint <- function(
    checkpoint,
    configuration,
    setting_id,
    design_id,
    n1,
    n2,
    n3
) {
  if (!is.list(checkpoint) ||
      !all(c("configuration", "completed", "replication_rows") %in%
           names(checkpoint))) {
    stop("Checkpoint is not a recognized MOBVE checkpoint object.")
  }

  ## Retain the script checksum for auditability, but do not reject a
  ## checkpoint solely because an output-only correction changed that
  ## checksum. Computational compatibility is enforced by the algorithm
  ## signature together with all simulation controls below.
  saved_configuration <- checkpoint$configuration
  current_configuration <- configuration
  saved_configuration$script_md5 <- NULL
  current_configuration$script_md5 <- NULL

  if (!identical(saved_configuration, current_configuration)) {
    stop(
      "Checkpoint configuration differs from the current algorithm or controls. ",
      "Use a new checkpoint directory or remove the stale scenario checkpoint."
    )
  }

  completed <- checkpoint$completed
  if (length(completed) != 1L ||
      !is.finite(completed) ||
      completed != floor(completed) ||
      completed < 0L ||
      completed > ITERATIONS) {
    stop("Checkpoint has an invalid completed-replication count.")
  }

  rows <- checkpoint$replication_rows
  if (!is.list(rows) || length(rows) != ITERATIONS) {
    stop("Checkpoint replication storage has the wrong length.")
  }

  expected_columns <- c(
    "Nominal_Level",
    "Setting_ID",
    "Setting",
    "Design_ID",
    "n1",
    "n2",
    "n3",
    "Replication",
    "Data_seed",
    "Bootstrap_seed",
    "True_theta",
    names(failed_replication_result())
  )

  if (completed > 0L) {
    for (replication in seq_len(completed)) {
      row <- rows[[replication]]

      if (!is.data.frame(row) ||
          nrow(row) != 1L ||
          !identical(names(row), expected_columns)) {
        stop("Checkpoint row schema is invalid at replication ", replication, ".")
      }

      expected_data_seed <- make_replication_seed(
        setting_id,
        design_id,
        replication,
        1L
      )
      expected_bootstrap_seed <- make_replication_seed(
        setting_id,
        design_id,
        replication,
        2L
      )

      identifiers_are_valid <-
        row$Nominal_Level == configuration$confidence_level &&
        row$Setting_ID == setting_id &&
        row$Design_ID == design_id &&
        row$n1 == n1 &&
        row$n2 == n2 &&
        row$n3 == n3 &&
        row$Replication == replication &&
        row$Data_seed == expected_data_seed &&
        row$Bootstrap_seed == expected_bootstrap_seed

      if (!isTRUE(identifiers_are_valid)) {
        stop("Checkpoint identifiers or seeds are invalid at replication ",
             replication, ".")
      }
    }
  }

  if (completed < ITERATIONS) {
    unfinished <- rows[seq.int(completed + 1L, ITERATIONS)]
    if (any(!vapply(unfinished, is.null, logical(1)))) {
      stop("Checkpoint contains rows beyond its completed prefix.")
    }
  }

  invisible(TRUE)
}

## ----------------------------- full simulation -----------------------------

run_mobve_all <- function(level) {
  validate_simulation_settings()
  assert_probability(level, "level", open = TRUE)

  number_scenarios <- nrow(mobve_settings) * nrow(sample_sizes)
  scenario_summaries <- vector("list", number_scenarios)
  scenario_status_counts <- vector("list", number_scenarios)
  retained_replications <- if (KEEP_REPLICATION_RESULTS) {
    vector("list", number_scenarios)
  } else {
    NULL
  }

  scenario_index <- 0L

  for (setting_index in seq_len(nrow(mobve_settings))) {
    setting <- mobve_settings[setting_index, , drop = FALSE]
    setting_id <- as.integer(setting$setting_id)

    true_values <- true_mobve_values(
      lx1 = setting$lx1,
      lx2 = setting$lx2,
      lx3 = setting$lx3,
      ly1 = setting$ly1,
      ly2 = setting$ly2,
      ly3 = setting$ly3,
      lz1 = setting$lz1,
      lz2 = setting$lz2,
      lz3 = setting$lz3
    )

    for (design_index in seq_len(nrow(sample_sizes))) {
      design <- sample_sizes[design_index, , drop = FALSE]
      design_id <- as.integer(design$design_id)
      n1 <- as.integer(design$n1)
      n2 <- as.integer(design$n2)
      n3 <- as.integer(design$n3)

      scenario_index <- scenario_index + 1L

      cat(sprintf(
        paste0(
          "MOBVE setting %d/%d, design %d/%d, ",
          "(n1,n2,n3) = (%d,%d,%d)\n"
        ),
        setting_index,
        nrow(mobve_settings),
        design_index,
        nrow(sample_sizes),
        n1,
        n2,
        n3
      ))

      configuration <- scenario_configuration(setting, design, level)
      checkpoint_path <- file.path(
        CHECKPOINT_DIR,
        sprintf(
          "MOBVE_L%02d_setting_%02d_design_%02d.rds",
          round(100 * level),
          setting_id,
          design_id
        )
      )

      replication_rows <- vector("list", ITERATIONS)
      completed <- 0L

      checkpoint_exists <- file.exists(checkpoint_path) ||
        file.exists(paste0(checkpoint_path, ".bak"))

      if (RESUME_FROM_CHECKPOINTS && checkpoint_exists) {
        saved_checkpoint <- read_checkpoint_safely(checkpoint_path)

        validate_checkpoint(
          checkpoint = saved_checkpoint,
          configuration = configuration,
          setting_id = setting_id,
          design_id = design_id,
          n1 = n1,
          n2 = n2,
          n3 = n3
        )

        replication_rows <- saved_checkpoint$replication_rows
        completed <- as.integer(saved_checkpoint$completed)

        cat(sprintf("  Resuming after replication %d.\n", completed))
      }

      if (completed < ITERATIONS) {
        first_new_replication <- completed + 1L

        for (replication in seq.int(first_new_replication, ITERATIONS)) {
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

          method_result <- tryCatch(
            {
              set.seed(data_seed)

              x <- mobve_sample(
                n = n1,
                lambda1 = setting$lx1,
                lambda2 = setting$lx2,
                lambda3 = setting$lx3
              )
              y <- mobve_sample(
                n = n2,
                lambda1 = setting$ly1,
                lambda2 = setting$ly2,
                lambda3 = setting$ly3
              )
              z <- mobve_sample(
                n = n3,
                lambda1 = setting$lz1,
                lambda2 = setting$lz2,
                lambda3 = setting$lz3
              )

              stopifnot(
                nrow(x) == n1,
                nrow(y) == n2,
                nrow(z) == n3,
                ncol(x) == 2L,
                ncol(y) == 2L,
                ncol(z) == 2L
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
              if (FAIL_FAST) stop(error_condition)

              message(sprintf(
                paste0(
                  "Recorded replication error: setting %d, design %d, ",
                  "replication %d: %s"
                ),
                setting_id,
                design_id,
                replication,
                conditionMessage(error_condition)
              ))

              failed_replication_result(conditionMessage(error_condition))
            }
          )

          replication_rows[[replication]] <- data.frame(
            Nominal_Level = level,
            Setting_ID = setting_id,
            Setting = as.character(setting$label),
            Design_ID = design_id,
            n1 = n1,
            n2 = n2,
            n3 = n3,
            Replication = replication,
            Data_seed = data_seed,
            Bootstrap_seed = bootstrap_seed,
            True_theta = true_values$theta,
            method_result,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )

          if (replication %% 25L == 0L || replication == ITERATIONS) {
            atomic_save_rds(
              list(
                configuration = configuration,
                completed = replication,
                replication_rows = replication_rows
              ),
              checkpoint_path
            )
          }
        }
      }

      if (any(vapply(replication_rows, is.null, logical(1)))) {
        stop("Scenario contains missing replication rows: ", checkpoint_path)
      }

      replication_data <- do.call(rbind, replication_rows)
      rownames(replication_data) <- NULL

      first_data_seed <- make_replication_seed(
        setting_id = setting_id,
        design_id = design_id,
        replication_id = 1L,
        stream = 1L
      )

      scenario_summaries[[scenario_index]] <- summarize_scenario(
        replication_data = replication_data,
        setting = setting,
        design = design,
        true_values = true_values,
        first_data_seed = first_data_seed,
        level = level
      )

      count_one_status <- function(method, values) {
        counts <- table(values, useNA = "ifany")
        data.frame(
          Nominal_Level = level,
          Setting_ID = setting_id,
          Design_ID = design_id,
          Method = method,
          Status = names(counts),
          Count = as.integer(counts),
          stringsAsFactors = FALSE
        )
      }

      scenario_status_counts[[scenario_index]] <- rbind(
        count_one_status("JEL", replication_data$JEL_status),
        count_one_status("Normal", replication_data$Normal_status),
        count_one_status("Kernel", replication_data$Kernel_status)
      )

      if (KEEP_REPLICATION_RESULTS) {
        retained_replications[[scenario_index]] <- replication_data
      }
    }
  }

  list(
    summary = do.call(rbind, scenario_summaries),
    status_counts = do.call(rbind, scenario_status_counts),
    replications = if (KEEP_REPLICATION_RESULTS) {
      do.call(rbind, retained_replications)
    } else {
      NULL
    }
  )
}

## ----------------------------- self-tests ----------------------------------

run_self_tests <- function() {
  ## Exact marker-specific VUS values and their differences.
  calculated_truth <- t(vapply(
    seq_len(nrow(mobve_settings)),
    function(index) {
      setting <- mobve_settings[index, , drop = FALSE]
      values <- true_mobve_values(
        lx1 = setting$lx1,
        lx2 = setting$lx2,
        lx3 = setting$lx3,
        ly1 = setting$ly1,
        ly2 = setting$ly2,
        ly3 = setting$ly3,
        lz1 = setting$lz1,
        lz2 = setting$lz2,
        lz3 = setting$lz3
      )
      c(values$vus1, values$vus2, values$theta)
    },
    numeric(3)
  ))

  expected_truth <- cbind(
    vus1 = mobve_settings$target_vus1,
    vus2 = mobve_settings$target_vus2,
    theta = mobve_settings$target_theta
  )
  if (max(abs(calculated_truth - expected_truth)) > 1e-12) {
    stop("MOBVE true-value self-test failed.")
  }

  ## Check the fast strict triplet counter against direct enumeration.
  x_small <- cbind(c(0.1, 0.4, 0.9), c(0.2, 0.5, 0.7))
  y_small <- cbind(c(0.3, 0.6), c(0.1, 0.6))
  z_small <- cbind(c(0.5, 0.8, 1.1, 1.3), c(0.4, 0.8, 1.0, 1.2))

  direct_count <- function(x_values, y_values, z_values) {
    total <- 0L
    for (i in seq_along(x_values)) {
      for (j in seq_along(y_values)) {
        for (k in seq_along(z_values)) {
          total <- total + as.integer(
            x_values[i] < y_values[j] && y_values[j] < z_values[k]
          )
        }
      }
    }
    total
  }

  direct_statistic <- (
    direct_count(x_small[, 1], y_small[, 1], z_small[, 1]) -
      direct_count(x_small[, 2], y_small[, 2], z_small[, 2])
  ) / (nrow(x_small) * nrow(y_small) * nrow(z_small))

  if (abs(u_stat(x_small, y_small, z_small) - direct_statistic) > 1e-12) {
    stop("Fast U-statistic self-test failed.")
  }

  ## Validate the pooled pseudo-values against the direct fixed-kernel pooled
  ## leave-one-out construction under unequal sample sizes.
  set.seed(MASTER_SEED + 31415L)
  x <- cbind(stats::runif(4), stats::runif(4))
  y <- cbind(stats::runif(5), stats::runif(5))
  z <- cbind(stats::runif(6), stats::runif(6))

  jackknife <- jackknife_values_stratified(x, y, z)
  pooled <- pooled_jackknife_from_stratified(jackknife)

  n1 <- nrow(x)
  n2 <- nrow(y)
  n3 <- nrow(z)
  total_n <- n1 + n2 + n3
  cstar <- choose(total_n, 3) / (n1 * n2 * n3)
  direct_pooled_values <- numeric(total_n)

  for (pooled_index in seq_len(total_n)) {
    if (pooled_index <= n1) {
      triplet_sum <- u_stat(x[-pooled_index, , drop = FALSE], y, z) *
        (n1 - 1) * n2 * n3
    } else if (pooled_index <= n1 + n2) {
      group_index <- pooled_index - n1
      triplet_sum <- u_stat(x, y[-group_index, , drop = FALSE], z) *
        n1 * (n2 - 1) * n3
    } else {
      group_index <- pooled_index - n1 - n2
      triplet_sum <- u_stat(x, y, z[-group_index, , drop = FALSE]) *
        n1 * n2 * (n3 - 1)
    }

    leave_one_out_pooled <-
      cstar * triplet_sum / choose(total_n - 1, 3)

    direct_pooled_values[pooled_index] <-
      total_n * jackknife$U0 -
      (total_n - 1) * leave_one_out_pooled
  }

  if (max(abs(pooled$v - direct_pooled_values)) > 1e-10 ||
      abs(mean(pooled$v) - jackknife$U0) > 1e-10 ||
      abs(mean(pooled$a) - 1) > 1e-12) {
    stop("Unequal-sample pooled-pseudo-value self-test failed.")
  }

  likelihood_at_estimate <- el_neg2log(
    pooled$v,
    pooled$a,
    jackknife$U0
  )
  if (!is.finite(likelihood_at_estimate) ||
      abs(likelihood_at_estimate) > 1e-8) {
    stop("JEL likelihood is not zero at the U-statistic in the self-test.")
  }

  if (el_neg2log(c(0, 0), c(1, 1), 0) != 0) {
    stop("All-zero JEL estimating-equation self-test failed.")
  }

  ## Exercise a strongly unbalanced design for which one centering coefficient
  ## is negative. This checks the mixed-sign path used by the confidence-set
  ## inversion, although it does not replace a full grid-sensitivity study.
  x_mixed <- cbind(
    seq(0.05, 0.25, length.out = 3),
    seq(0.15, 0.35, length.out = 3)
  )
  y_mixed <- cbind(
    seq(0.20, 0.80, length.out = 12),
    seq(0.10, 0.70, length.out = 12)
  )
  z_mixed <- cbind(
    seq(0.50, 0.95, length.out = 5),
    seq(0.40, 0.90, length.out = 5)
  )
  jackknife_mixed <- jackknife_values_stratified(
    x_mixed,
    y_mixed,
    z_mixed
  )
  pooled_mixed <- pooled_jackknife_from_stratified(jackknife_mixed)

  if (!(pooled_mixed$group_a["Y"] < 0) ||
      abs(mean(pooled_mixed$v) - jackknife_mixed$U0) > 1e-10 ||
      abs(mean(pooled_mixed$a) - 1) > 1e-12 ||
      abs(el_neg2log(
        pooled_mixed$v,
        pooled_mixed$a,
        jackknife_mixed$U0
      )) > 1e-8) {
    stop("Mixed-sign pooled-JEL self-test failed.")
  }

  confidence_set_mixed <- jel_confidence_set_from_jackknife(
    jackknife = jackknife_mixed,
    level = max(CONFIDENCE_LEVELS),
    lower = -1,
    upper = 1,
    grid_len = 201L
  )

  if (!isTRUE(confidence_set_mixed$valid) ||
      !value_in_intervals(
        jackknife_mixed$U0,
        confidence_set_mixed$intervals
      ) ||
      abs(
        confidence_set_mixed$total_length -
          sum(
            confidence_set_mixed$intervals[, 2] -
              confidence_set_mixed$intervals[, 1]
          )
      ) > 1e-10) {
    stop("Mixed-sign JEL confidence-set self-test failed.")
  }

  ## Verify uniqueness of all data and bootstrap seeds used by this run.
  seed_grid <- expand.grid(
    setting_id = mobve_settings$setting_id,
    design_id = sample_sizes$design_id,
    replication_id = seq_len(ITERATIONS),
    stream = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )

  seeds <- mapply(
    make_replication_seed,
    setting_id = seed_grid$setting_id,
    design_id = seed_grid$design_id,
    replication_id = seed_grid$replication_id,
    stream = seed_grid$stream
  )

  if (anyDuplicated(seeds)) stop("Replication seed collision detected.")

  invisible(TRUE)
}

## ----------------------------- execute simulation --------------------------

SCRIPT_PATH_AT_START <- current_script_path()
SCRIPT_MD5_AT_START <- if (
  is.character(SCRIPT_PATH_AT_START) &&
  length(SCRIPT_PATH_AT_START) == 1L &&
  !is.na(SCRIPT_PATH_AT_START) &&
  file.exists(SCRIPT_PATH_AT_START)
) {
  unname(tools::md5sum(SCRIPT_PATH_AT_START))
} else {
  "not_detected"
}
ALGORITHM_SIGNATURE <- compute_algorithm_signature()

validate_simulation_settings()
run_self_tests()

start_all <- Sys.time()
run_lock_directory <- file.path(CHECKPOINT_DIR, ".active_run_lock")
if (!dir.create(run_lock_directory, showWarnings = FALSE)) {
  stop(
    "Another run may be using this checkpoint directory, or a stale lock ",
    "remains after an interrupted process: ",
    run_lock_directory
  )
}

simulation_outputs <- lapply(CONFIDENCE_LEVELS, run_mobve_all)
end_all <- Sys.time()

mobve_result <- do.call(
  rbind,
  lapply(simulation_outputs, function(output) output$summary)
)
status_counts <- do.call(
  rbind,
  lapply(simulation_outputs, function(output) output$status_counts)
)
replication_result <- if (KEEP_REPLICATION_RESULTS) {
  do.call(
    rbind,
    lapply(simulation_outputs, function(output) output$replications)
  )
} else {
  NULL
}
total_time_seconds <- as.numeric(
  difftime(end_all, start_all, units = "secs")
)

## ----------------------------- derived output tables -----------------------

make_paper_table <- function(data) {
  reportable <- !is.na(data$No_unhandled_errors) &
    as.logical(data$No_unhandled_errors)
  mask_unreportable <- function(values) {
    values[!reportable] <- NA
    values
  }

  data.frame(
    Nominal_Level = data$Nominal_Level,
    Scenario_Suite = data$Scenario_Suite,
    Parameter_setting = data$Setting,
    Sample_size = paste0("(", data$n1, ",", data$n2, ",", data$n3, ")"),
    Common_shock_probability = round(data$Common_shock_probability_X, 3),
    R_total = data$R_total,
    Computationally_reportable = reportable,
    Unhandled_error_n = data$Unhandled_error_n,
    Numerically_tiny_total_variance_percent = round(
      data$Numerically_tiny_total_variance_percent,
      2
    ),

    JEL_valid_n = data$JEL_valid_n,
    JEL_CP_unconditional_percent = round(
      mask_unreportable(data$JEL_CP),
      2
    ),
    JEL_CP_MCSE = round(mask_unreportable(data$JEL_CP_MCSE), 2),
    JEL_CP_conditional_valid_percent = round(
      mask_unreportable(data$JEL_CP_valid),
      2
    ),
    JEL_AL_union_conditional_valid = round(
      mask_unreportable(data$JEL_AL),
      4
    ),
    JEL_AL_MCSE = round(mask_unreportable(data$JEL_AL_MCSE), 4),
    JEL_failure_percent = round(data$JEL_failure_percent, 2),
    JEL_multicomponent_percent = round(
      data$JEL_multicomponent_percent,
      2
    ),
    JEL_exact_degenerate_percent = round(data$JEL_degenerate_percent, 2),

    Normal_valid_n = data$Normal_valid_n,
    Normal_CP_conditional_valid_percent = round(
      mask_unreportable(data$Normal_CP_valid),
      2
    ),
    Normal_CP_truncated_unconditional_percent = round(
      mask_unreportable(data$Normal_CP_truncated),
      2
    ),
    Normal_CP_truncated_MCSE = round(
      mask_unreportable(data$Normal_CP_truncated_MCSE),
      2
    ),
    Normal_AL_truncated_conditional_valid = round(
      mask_unreportable(data$Normal_AL_truncated),
      4
    ),
    Normal_AL_truncated_MCSE = round(
      mask_unreportable(data$Normal_AL_truncated_MCSE),
      4
    ),
    Normal_CP_raw_unconditional_percent = round(
      mask_unreportable(data$Normal_CP),
      2
    ),
    Normal_CP_raw_MCSE = round(
      mask_unreportable(data$Normal_CP_MCSE),
      2
    ),
    Normal_AL_raw_conditional_valid = round(
      mask_unreportable(data$Normal_AL),
      4
    ),
    Normal_AL_raw_MCSE = round(
      mask_unreportable(data$Normal_AL_MCSE),
      4
    ),
    Normal_failure_percent = round(data$Normal_failure_percent, 2),

    Kernel_valid_n = data$Kernel_valid_n,
    Kernel_CP_unconditional_percent = round(
      mask_unreportable(data$Kernel_CP),
      2
    ),
    Kernel_CP_MCSE = round(mask_unreportable(data$Kernel_CP_MCSE), 2),
    Kernel_CP_conditional_valid_percent = round(
      mask_unreportable(data$Kernel_CP_valid),
      2
    ),
    Kernel_AL_conditional_valid = round(
      mask_unreportable(data$Kernel_AL),
      4
    ),
    Kernel_AL_MCSE = round(mask_unreportable(data$Kernel_AL_MCSE), 4),
    Kernel_failure_percent = round(data$Kernel_failure_percent, 2),
    stringsAsFactors = FALSE
  )
}

mobve_table_90 <- make_paper_table(
  mobve_result[mobve_result$Nominal_Level == 0.90, , drop = FALSE]
)
mobve_table_95 <- make_paper_table(
  mobve_result[mobve_result$Nominal_Level == 0.95, , drop = FALSE]
)

make_calibration_comparison <- function(data) {
  nominal_percent <- 100 * data$Nominal_Level
  errors <- cbind(
    JEL = abs(data$JEL_CP - nominal_percent),
    Normal = abs(data$Normal_CP_truncated - nominal_percent),
    Kernel = abs(data$Kernel_CP - nominal_percent)
  )
  errors[!data$No_unhandled_errors, ] <- NA_real_
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
    Scenario_Suite = data$Scenario_Suite,
    Setting_ID = data$Setting_ID,
    Design_ID = data$Design_ID,
    Sample_Size = paste0("(", data$n1, ",", data$n2, ",", data$n3, ")"),
    True_theta = data$True_theta,
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

calibration_comparison <- make_calibration_comparison(mobve_result)

true_value_table <- do.call(
  rbind,
  lapply(seq_len(nrow(mobve_settings)), function(index) {
    setting <- mobve_settings[index, , drop = FALSE]
    values <- true_mobve_values(
      lx1 = setting$lx1,
      lx2 = setting$lx2,
      lx3 = setting$lx3,
      ly1 = setting$ly1,
      ly2 = setting$ly2,
      ly3 = setting$ly3,
      lz1 = setting$lz1,
      lz2 = setting$lz2,
      lz3 = setting$lz3
    )

    data.frame(
      Setting_ID = setting$setting_id,
      Setting = setting$label,
      Scenario_Suite = SCENARIO_SUITE,
      Source = SCENARIO_SOURCE,
      Target_VUS1 = setting$target_vus1,
      Target_VUS2 = setting$target_vus2,
      Target_theta = setting$target_theta,
      Common_shock_probability_X =
        setting$common_shock_probability_X,
      Common_shock_probability_Y =
        setting$common_shock_probability_Y,
      Common_shock_probability_Z =
        setting$common_shock_probability_Z,
      True_VUS1 = values$vus1,
      True_VUS2 = values$vus2,
      True_theta = values$theta,
      stringsAsFactors = FALSE
    )
  })
)

## ----------------------------- provenance helpers --------------------------

script_path <- SCRIPT_PATH_AT_START
script_md5 <- SCRIPT_MD5_AT_START

rng_description <- paste(RNGkind(), collapse = "; ")

readme <- data.frame(
  Item = c(
    "Model",
    "Script version",
    "Run mode",
    "Scenario profile",
    "Master seed",
    "Iterations per scenario",
    "Bootstrap replications",
    "Confidence levels",
    "Monte Carlo Wilson level",
    "Sample-size design interpretation",
    "JEL grid length",
    "JEL confidence-set convention",
    "JEL length convention",
    "JEL grid limitation",
    "Coverage convention",
    "Average-length convention",
    "Coverage uncertainty",
    "Computational integrity flag",
    "Normal interval convention",
    "Normal invalid-variance threshold",
    "Normal tiny group-variance diagnostic",
    "Numerically tiny total-variance diagnostic",
    "JEL exact-degeneracy threshold",
    "Kernel estimator definition",
    "Bootstrap quantile convention",
    "Data seed formula",
    "Bootstrap seed formula",
    "RNG kind",
    "R version",
    "openxlsx version",
    "Script path",
    "Script MD5",
    "Algorithm signature",
    "Checkpoint directory",
    "Resume enabled",
    "Checkpoint provenance convention",
    "Replication sheet retained",
    "Operating system",
    "Machine architecture",
    "Detected logical cores",
    "Started",
    "Finished",
    "Current execution simulation time (seconds)",
    "Resume timing interpretation",
    "JEL/Normal timing convention",
    "Kernel timing convention",
    "Inference warning"
  ),
  Value = c(
    "Marshall-Olkin bivariate exponential",
    "MOBVE-corrected-2026-08-19-v5",
    RUN_TAG,
    SCENARIO_PROFILE,
    as.character(MASTER_SEED),
    as.character(ITERATIONS),
    as.character(BOOT),
    paste(CONFIDENCE_LEVELS, collapse = ", "),
    as.character(MC_COVERAGE_INTERVAL_LEVEL),
    paste0(
      "Every one of the four prespecified parameter settings is crossed with ",
      "the same six unequal-allocation sample-size designs."
    ),
    as.character(GRID_LEN_JEL),
    "Union of all grid-detected accepted components; coverage is membership in the union",
    "Sum of component lengths, not the hull length",
    "Components or gaps narrower than the base grid spacing can be missed; inspect multicomponent diagnostics and perform grid sensitivity checks",
    "Unconditional CP counts invalid intervals as noncoverage; CP_valid is conditional on valid construction",
    "AL and its MCSE are conditional on valid construction",
    "Scenario Summary includes plug-in MCSE and a fixed 95% Wilson interval; use the Wilson interval when observed CP is 0% or 100%",
    "FALSE only when an unhandled computational error occurred; ordinary recorded construction failures remain analyzable",
    "Scenario Summary reports both raw and intersection-with-[-1,1] Normal intervals",
    "Normal is invalid when the estimated total jackknife variance is nonfinite or <= 1e-14",
    "TRUE when any group-specific pseudo-value sample variance is below 1e-12; this is a numerical flag, not an asymptotic degeneracy test",
    "Percent of finite total jackknife variance estimates below 1e-12; this is shared descriptive evidence, not a proof of JEL nondegeneracy",
    "All |v_i-a_i*U0| <= 1e-12*max(1,max|v|,max|a_i*U0|)",
    "Product of two smoothed pairwise comparisons sharing the observed Y row; verify that this matches the estimator stated in the manuscript",
    "R quantile type=1; empirical order-statistic percentile endpoints",
    "MASTER_SEED + 10000000*setting_id + 100000*design_id + 10*replication_id + 1",
    "MASTER_SEED + 10000000*setting_id + 100000*design_id + 10*replication_id + 2",
    rng_description,
    R.version.string,
    as.character(utils::packageVersion("openxlsx")),
    ifelse(is.na(script_path), "not detected", script_path),
    ifelse(is.na(script_md5), "not detected", script_md5),
    ALGORITHM_SIGNATURE,
    normalizePath(CHECKPOINT_DIR),
    as.character(RESUME_FROM_CHECKPOINTS),
    paste0(
      "A checkpoint is reused only when controls, setting/design, algorithm ",
      "signature, R/openxlsx versions, RNG kind, and hardware signature match; ",
      "the script MD5 is retained as audit metadata"
    ),
    as.character(KEEP_REPLICATION_RESULTS),
    paste(Sys.info()[c("sysname", "release")], collapse = " "),
    as.character(Sys.info()["machine"]),
    as.character(parallel::detectCores(logical = TRUE)),
    format(start_all, "%Y-%m-%dT%H:%M:%S%z"),
    format(end_all, "%Y-%m-%dT%H:%M:%S%z"),
    format(total_time_seconds, digits = 12),
    "Current execution time excludes work loaded from checkpoints; scenario method-time fields include the stored per-replication attempts",
    "Standalone time for each includes the same shared stratified-jackknife time",
    "Standalone kernel bootstrap time",
    "JEL chi-square calibration remains conditional on a valid Wilks theorem and nondegeneracy"
  ),
  stringsAsFactors = FALSE
)

table_notes <- data.frame(
  Item = c(
    "Primary coverage",
    "Common-shock probability",
    "Conditional coverage",
    "Average length",
    "JEL set",
    "Normal primary comparator",
    "Computationally_reportable",
    "Monte Carlo uncertainty"
  ),
  Explanation = c(
    "CP_unconditional uses all R replications and counts invalid constructions as noncoverage",
    paste0(
      "lambda3/(lambda1+lambda2+lambda3), reported separately for X, Y, ",
      "and Z in the scenario summary and true-value sheets."
    ),
    "CP_conditional_valid uses only replications with a valid constructed interval/set",
    "Every AL and AL_MCSE is conditional on valid construction",
    "Grid-detected union; AL is total component length and multicomponent frequency is reported",
    "The truncated Normal interval is intersected with [-1,1] and is the like-for-like bounded comparison; raw Normal is a sensitivity result",
    "FALSE means an unhandled error occurred and inferential entries in Table 2 are masked",
    "MCSE is the plug-in Monte Carlo standard error; Scenario Summary also gives Wilson coverage limits"
  ),
  stringsAsFactors = FALSE
)

## ----------------------------- Excel output --------------------------------

excel_column_label <- function(index) {
  assert_scalar_integer(index, "Excel column index")
  label <- ""

  while (index > 0L) {
    remainder <- (index - 1L) %% 26L
    label <- paste0(LETTERS[remainder + 1L], label)
    index <- (index - 1L) %/% 26L
  }

  label
}

sheet_data <- list(
  README = readme,
  `Table 2 Notes` = table_notes,
  `Table MOBVE 90` = mobve_table_90,
  `Table MOBVE 95` = mobve_table_95,
  `Calibration Comparison` = calibration_comparison,
  `Scenario Summary` = mobve_result,
  `Parameter MOBVE` = mobve_settings,
  `Sample Sizes` = sample_sizes,
  `True Values` = true_value_table,
  `Status Codes` = STATUS_CODEBOOK,
  `Status Counts` = status_counts
)

if (KEEP_REPLICATION_RESULTS) {
  if (nrow(replication_result) > 1048575L) {
    stop("Replication Results exceeds the Excel worksheet row limit.")
  }
  sheet_data[["Replication Results"]] <- replication_result
}

workbook <- openxlsx::createWorkbook(
  creator = "MOBVE simulation script",
  title = "MOBVE VUS simulation",
  subject = "Unequal-sample JEL, Normal, and Kernel confidence intervals",
  category = "Simulation results"
)

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  halign = "center",
  valign = "center"
)

for (sheet_index in seq_along(sheet_data)) {
  sheet_name <- names(sheet_data)[sheet_index]
  data_to_write <- sheet_data[[sheet_name]]
  openxlsx::addWorksheet(workbook, sheet_name)
  openxlsx::writeData(workbook, sheet_name, data_to_write)

  used_columns <- seq_len(ncol(data_to_write))
  openxlsx::addStyle(
    workbook,
    sheet_name,
    header_style,
    rows = 1L,
    cols = used_columns,
    gridExpand = TRUE
  )
  openxlsx::freezePane(workbook, sheet_name, firstRow = TRUE)
  openxlsx::addFilter(
    workbook,
    sheet_name,
    rows = 1L,
    cols = used_columns
  )

  if (sheet_name == "Replication Results") {
    openxlsx::setColWidths(
      workbook,
      sheet_name,
      cols = used_columns,
      widths = pmin(18, pmax(10, nchar(names(data_to_write)) + 2))
    )
  } else {
    openxlsx::setColWidths(
      workbook,
      sheet_name,
      cols = used_columns,
      widths = "auto"
    )
  }

  ## openxlsx can otherwise retain the initializer <dimension ref="A1"/> when
  ## a worksheet is populated in a single write operation. Set the true used
  ## range explicitly so independent read-only parsers recover every cell.
  final_cell <- paste0(
    excel_column_label(ncol(data_to_write)),
    nrow(data_to_write) + 1L
  )
  workbook$worksheets[[sheet_index]]$dimension <- paste0(
    "<dimension ref=\"A1:",
    final_cell,
    "\"/>"
  )
}

temporary_workbook <- tempfile(
  pattern = "MOBVE_workbook_",
  tmpdir = OUTPUT_DIR,
  fileext = ".xlsx"
)

openxlsx::saveWorkbook(
  workbook,
  temporary_workbook,
  overwrite = TRUE
)

assert_workbook_roundtrip <- function(original, recovered, key_columns, label) {
  if (nrow(recovered) != nrow(original) ||
      !identical(names(recovered), names(original))) {
    stop(label, " failed the workbook dimension/name round-trip check.")
  }

  if (anyDuplicated(original[, key_columns, drop = FALSE]) ||
      anyDuplicated(recovered[, key_columns, drop = FALSE])) {
    stop(label, " contains duplicate keys.")
  }

  for (column_name in names(original)) {
    expected <- original[[column_name]]
    observed <- recovered[[column_name]]

    if (is.numeric(expected) || is.integer(expected)) {
      observed_numeric <- suppressWarnings(as.numeric(observed))

      if (!identical(is.na(expected), is.na(observed_numeric))) {
        stop(label, " changed missing values in column ", column_name, ".")
      }

      finite <- is.finite(expected) & is.finite(observed_numeric)
      if (any(finite)) {
        tolerance <- 1e-10 * max(1, max(abs(expected[finite])))
        if (max(abs(expected[finite] - observed_numeric[finite])) > tolerance) {
          stop(label, " changed numeric values in column ", column_name, ".")
        }
      }
    } else if (is.logical(expected)) {
      observed_logical <- as.logical(observed)
      if (!identical(expected, observed_logical)) {
        stop(label, " changed logical values in column ", column_name, ".")
      }
    } else {
      expected_character <- as.character(expected)
      observed_character <- as.character(observed)
      if (!identical(expected_character, observed_character)) {
        stop(label, " changed text values in column ", column_name, ".")
      }
    }
  }

  invisible(TRUE)
}

assert_raw_ooxml_dimensions <- function(workbook_path, written_sheets) {
  archive_listing <- utils::unzip(workbook_path, list = TRUE)$Name
  worksheet_files <- grep(
    "^xl/worksheets/sheet[0-9]+[.]xml$",
    archive_listing,
    value = TRUE
  )

  sheet_numbers <- as.integer(sub(
    "^xl/worksheets/sheet([0-9]+)[.]xml$",
    "\\1",
    worksheet_files
  ))
  worksheet_files <- worksheet_files[order(sheet_numbers)]

  if (length(worksheet_files) != length(written_sheets)) {
    stop("The OOXML archive contains an unexpected number of worksheets.")
  }

  extraction_directory <- tempfile(pattern = "MOBVE_ooxml_check_")
  dir.create(extraction_directory)
  on.exit(unlink(extraction_directory, recursive = TRUE), add = TRUE)

  utils::unzip(
    workbook_path,
    files = worksheet_files,
    exdir = extraction_directory
  )

  for (sheet_index in seq_along(worksheet_files)) {
    sheet_xml_path <- file.path(
      extraction_directory,
      worksheet_files[sheet_index]
    )
    connection <- file(sheet_xml_path, open = "rb")
    xml_text <- readChar(connection, nchars = 8192L, useBytes = TRUE)
    close(connection)
    dimension_match <- regexec(
      "<dimension[^>]*ref=\"([^\"]+)\"",
      xml_text,
      perl = TRUE
    )
    captured <- regmatches(xml_text, dimension_match)[[1]]

    if (length(captured) != 2L) {
      stop("Worksheet OOXML is missing its dimension reference: ",
           worksheet_files[sheet_index])
    }

    sheet_data_frame <- written_sheets[[sheet_index]]
    expected_reference <- paste0(
      "A1:",
      excel_column_label(ncol(sheet_data_frame)),
      nrow(sheet_data_frame) + 1L
    )

    if (!identical(captured[2], expected_reference)) {
      stop(
        "Worksheet OOXML dimension mismatch for ",
        names(written_sheets)[sheet_index],
        ": expected ",
        expected_reference,
        ", found ",
        captured[2]
      )
    }
  }

  invisible(TRUE)
}

normalize_archive_path <- function(path) {
  path <- sub("^/+", "", path)
  pieces <- strsplit(path, "/", fixed = TRUE)[[1]]
  output <- character(0)

  for (piece in pieces) {
    if (!nzchar(piece) || piece == ".") next
    if (piece == "..") {
      if (length(output) == 0L) return(NA_character_)
      output <- head(output, -1L)
    } else {
      output <- c(output, piece)
    }
  }

  paste(output, collapse = "/")
}

assert_ooxml_relationship_targets <- function(workbook_path) {
  archive_listing <- utils::unzip(workbook_path, list = TRUE)$Name
  relationship_files <- grep(
    "(^_rels/[.]rels$|(^|/)_rels/[^/]+[.]rels$)",
    archive_listing,
    value = TRUE
  )

  extraction_directory <- tempfile(pattern = "MOBVE_ooxml_rel_check_")
  dir.create(extraction_directory)
  on.exit(unlink(extraction_directory, recursive = TRUE), add = TRUE)

  files_to_extract <- c(relationship_files, "[Content_Types].xml")
  utils::unzip(
    workbook_path,
    files = files_to_extract,
    exdir = extraction_directory
  )

  for (relationship_file in relationship_files) {
    xml <- paste(
      readLines(
        file.path(extraction_directory, relationship_file),
        warn = FALSE
      ),
      collapse = ""
    )
    tags <- regmatches(
      xml,
      gregexpr("<Relationship\\b[^>]*/>", xml, perl = TRUE)
    )[[1]]

    if (length(tags) == 0L) next

    relationship_directory <- dirname(relationship_file)
    source_directory <- dirname(relationship_directory)
    if (identical(source_directory, ".")) source_directory <- ""

    for (tag in tags) {
      if (grepl("TargetMode=\"External\"", tag, fixed = TRUE)) next

      target_match <- regexec(
        "Target=\"([^\"]+)\"",
        tag,
        perl = TRUE
      )
      target_capture <- regmatches(tag, target_match)[[1]]
      if (length(target_capture) != 2L) {
        stop("OOXML relationship has no readable target in ",
             relationship_file)
      }

      target <- target_capture[2]
      combined <- if (startsWith(target, "/")) {
        target
      } else if (nzchar(source_directory)) {
        paste(source_directory, target, sep = "/")
      } else {
        target
      }
      resolved_target <- normalize_archive_path(combined)

      if (is.na(resolved_target) ||
          !resolved_target %in% archive_listing) {
        stop(
          "Dangling OOXML relationship in ",
          relationship_file,
          ": ",
          target
        )
      }
    }
  }

  content_type_path <- file.path(
    extraction_directory,
    "[Content_Types].xml"
  )
  content_type_xml <- paste(
    readLines(content_type_path, warn = FALSE),
    collapse = ""
  )
  override_tags <- regmatches(
    content_type_xml,
    gregexpr("<Override\\b[^>]*/>", content_type_xml, perl = TRUE)
  )[[1]]

  for (tag in override_tags) {
    part_match <- regexec("PartName=\"([^\"]+)\"", tag, perl = TRUE)
    part_capture <- regmatches(tag, part_match)[[1]]
    if (length(part_capture) == 2L) {
      part_name <- normalize_archive_path(part_capture[2])
      if (is.na(part_name) || !part_name %in% archive_listing) {
        stop("OOXML content-type override targets a missing part: ",
             part_capture[2])
      }
    }
  }

  invisible(TRUE)
}

## Some openxlsx releases can emit an unused worksheet-to-drawing
## relationship without writing the corresponding drawing part, even though
## this workbook contains no charts or images. Do not promote such a package.
## Rebuild the same sheets with writexl and then subject the rebuilt file to all
## numerical and raw-OOXML checks below.
relationship_preflight <- tryCatch(
  {
    assert_ooxml_relationship_targets(temporary_workbook)
    NULL
  },
  error = function(error_condition) error_condition
)

if (inherits(relationship_preflight, "error")) {
  message(
    "openxlsx produced an invalid internal relationship; ",
    "rewriting the workbook with writexl."
  )
  need_pkg("writexl")

  if (file.exists(temporary_workbook) &&
      unlink(temporary_workbook, force = TRUE) != 0L) {
    stop("Could not remove the invalid temporary workbook.")
  }

  writexl::write_xlsx(
    x = sheet_data,
    path = temporary_workbook,
    col_names = TRUE,
    format_headers = TRUE
  )
}

## Round-trip checks catch incomplete output before promotion.
summary_check <- openxlsx::read.xlsx(
  temporary_workbook,
  sheet = "Scenario Summary",
  check.names = FALSE
)

assert_workbook_roundtrip(
  original = mobve_result,
  recovered = summary_check,
  key_columns = c("Nominal_Level", "Setting_ID", "Design_ID"),
  label = "Scenario Summary"
)

if (KEEP_REPLICATION_RESULTS) {
  replication_check <- openxlsx::read.xlsx(
    temporary_workbook,
    sheet = "Replication Results",
    check.names = FALSE
  )

  assert_workbook_roundtrip(
    original = replication_result,
    recovered = replication_check,
    key_columns = c(
      "Nominal_Level",
      "Setting_ID",
      "Design_ID",
      "Replication"
    ),
    label = "Replication Results"
  )
}

## This independent raw-XML check detects the erroneous A1 dimensions found in
## the earlier workbook even if the same package can read its own output.
assert_raw_ooxml_dimensions(temporary_workbook, sheet_data)
assert_ooxml_relationship_targets(temporary_workbook)

promote_file_with_backup(temporary_workbook, OUTFILE)

if (unlink(run_lock_directory, recursive = TRUE, force = TRUE) != 0L) {
  warning("The completed-run lock could not be removed: ", run_lock_directory)
}

cat("\nCompleted successfully.\n")
cat("Excel workbook:", normalizePath(OUTFILE), "\n")
cat("Checkpoint directory:", normalizePath(CHECKPOINT_DIR), "\n")
cat("Current execution simulation time (seconds):", total_time_seconds, "\n")
