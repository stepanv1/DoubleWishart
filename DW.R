# ============================================================================
# DW.R
# Publication figure for:
# "Largest Finite Roots of Identity-Scale Doubly Singular Beta Type II Ensembles"
#
# The script
# simulates the ORIGINAL p-dimensional doubly singular problem, evaluates the
# exact reduced double-Wishart CDF with rootWishartHD, evaluates Johnstone's
# Tracy-Widom approximation.
# ============================================================================

# --------------------------- editable settings -------------------------------
REPO_DIR <- normalizePath(".", winslash = "/", mustWork = TRUE)
FIG_DIR <- file.path(REPO_DIR, "figures")
CACHE_DIR <- file.path(REPO_DIR, "figure_cache")

# ORIGINAL doubly singular matrices:
#   A = X X^T, X is p x M_DF
#   B = Y Y^T, Y is p x Q_DF
# The manuscript theorem assumes p > M_DF >= Q_DF and identity scale.
#
# Main sweep. Edit all three dimension parameters here.
M_DF <- 100L
Q_DF <- 6L
P_GRID <- c(101L, 110L, 120L, 130L, 200L, 500L)
P_SHOW <- c(101L, 200L, 500L)

# Separate near-saturated case. Edit all three dimensions here.
STRESS_P <- 100L
STRESS_M_DF <- 99L
STRESS_Q_DF <- 98L

# Monte Carlo and numerical grids.
N_SIM <- 10000L 
N_CURVE_GRID <- 801L
N_STRESS_GRID <- 801L
SEED <- 190501774L
SEED_STRESS <- 190601775L

# CDF plotting controls.
# X_SCALE choices:
#   "real"  : plot lambda_max on its original scale;
#   "log10" : plot log10(lambda_max), which is strongly recommended when
#             p is close to m because the beta-II right tail can be enormous.
MAIN_X_SCALE <- "log10"
STRESS_X_SCALE <- "log10"

# CDF_RANGE choices:
#   "full"     : include the complete observed Monte Carlo range;
#   "quantile" : restrict to CDF_QUANTILES for a main figure.
CDF_RANGE <- "full"
CDF_QUANTILES <- c(0.0005, 0.9995)
CDF_X_PADDING <- 0.03

# Backend controls.
# "auto" uses robust arbitrary precision whenever a Jacobi parameter is at or
# near the boundary (m_C=0 or n_C=0), including p=101,m=100,q=6.
EXACT_BACKEND <- "auto"
FORCE_RECOMPUTE <- FALSE
FORCE_RECOMPUTE_EXACT <- FALSE

# Fail rather than silently publishing a numerically broken "exact" curve.
CHECK_EXACT_AGREEMENT <- TRUE
EXACT_DISCREPANCY_MAX <- 0.05

.detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(.detected_cores) || .detected_cores < 1L) .detected_cores <- 1L
N_CORES <- min(8L, max(1L, .detected_cores - 1L))

FIG_WIDTH_IN <- 7.2
FIG_HEIGHT_IN <- 5.45
FIG_DPI <- 600L
# -----------------------------------------------------------------------------

required_packages <- c("rootWishartHD", "RMTstat", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the required packages before running DW.R: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (!all(P_SHOW %in% P_GRID)) {
  stop("Every value in P_SHOW must also occur in P_GRID.", call. = FALSE)
}
if (!(all(P_GRID > M_DF) && M_DF >= Q_DF)) {
  stop("The script requires p > m >= q for every p in P_GRID.", call. = FALSE)
}
if (!(STRESS_P > STRESS_M_DF && STRESS_M_DF >= STRESS_Q_DF)) {
  stop("The stress case must satisfy p > m >= q.", call. = FALSE)
}
if (N_SIM < 100L || N_SIM_STRESS < 100L) {
  stop("Monte Carlo sample sizes are too small for a publication figure.",
       call. = FALSE)
}

MAIN_X_SCALE <- match.arg(MAIN_X_SCALE, c("real", "log10"))
STRESS_X_SCALE <- match.arg(STRESS_X_SCALE, c("real", "log10"))
CDF_RANGE <- match.arg(CDF_RANGE, c("full", "quantile"))
if (length(CDF_QUANTILES) != 2L ||
    any(!is.finite(CDF_QUANTILES)) ||
    CDF_QUANTILES[1] < 0 ||
    CDF_QUANTILES[2] > 1 ||
    CDF_QUANTILES[1] >= CDF_QUANTILES[2]) {
  stop("CDF_QUANTILES must be two increasing probabilities in [0,1].",
       call. = FALSE)
}
if (!is.finite(CDF_X_PADDING) || CDF_X_PADDING < 0) {
  stop("CDF_X_PADDING must be a finite nonnegative number.", call. = FALSE)
}

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# Chiani/rootWishartHD parameters for
# W_q(m_df, I_q) W_q(p - m_df + q_df, I_q)^(-1).
dsb_params <- function(p, m_df, q_df) {
  if (!(p > m_df && m_df >= q_df && q_df >= 1L)) {
    stop("Expected p > m_df >= q_df >= 1.", call. = FALSE)
  }
  df2 <- p - m_df + q_df
  out <- list(
    s = as.integer(q_df),
    m = 0.5 * (m_df - q_df - 1),
    n = 0.5 * (df2 - q_df - 1),
    df2 = as.integer(df2)
  )
  if (out$m <= -1 || out$n <= -1) {
    stop("Invalid Jacobi parameters generated from p, m_df, q_df.", call. = FALSE)
  }
  out
}

# One draw from the original p-dimensional doubly singular construction.
# If X is p x m and Y is p x q, then A = XX' and B = YY'.  A thin QR
# decomposition X = QR gives A^+ = Q(RR')^(-1)Q', so the positive finite roots
# are the eigenvalues of U'U, U = R^(-1)Q'Y.  No p x p matrix is formed.
simulate_one_lmax <- function(seed, p, m_df, q_df) {
  set.seed(seed)
  X <- matrix(stats::rnorm(p * m_df), nrow = p, ncol = m_df)
  Y <- matrix(stats::rnorm(p * q_df), nrow = p, ncol = q_df)
  
  qr_x <- qr(X, LAPACK = TRUE)
  Q <- qr.Q(qr_x, complete = FALSE)
  R <- qr.R(qr_x, complete = FALSE)
  
  qty <- crossprod(Q, Y)
  U <- backsolve(R, qty)
  roots <- eigen(crossprod(U), symmetric = TRUE, only.values = TRUE)$values
  max(0, roots[1L])
}

simulate_chunk <- function(indices, seeds, p, m_df, q_df) {
  vapply(
    indices,
    function(i) simulate_one_lmax(seeds[i], p, m_df, q_df),
    numeric(1)
  )
}

simulate_for_p <- function(p, m_df, q_df, n_sim, n_cores, seed,
                           cache_dir, force = FALSE) {
  cache_file <- file.path(
    cache_dir,
    sprintf(
      "lmax_original_p%d_m%d_q%d_N%d_seed%d.rds",
      p, m_df, q_df, n_sim, seed
    )
  )
  
  if (file.exists(cache_file) && !isTRUE(force)) {
    cached <- readRDS(cache_file)
    message("Using cache: ", basename(cache_file))
    return(cached$lambda_max)
  }
  
  message(
    sprintf(
      "Simulating p=%d, m=%d, q=%d, N=%d with %d worker(s)...",
      p, m_df, q_df, n_sim, n_cores
    )
  )
  
  set.seed(seed + 100003L * p)
  seeds <- sample.int(.Machine$integer.max, n_sim, replace = FALSE)
  n_workers <- min(as.integer(n_cores), as.integer(n_sim))
  chunks <- if (n_workers <= 1L) {
    list(seq_len(n_sim))
  } else {
    split(
      seq_len(n_sim),
      cut(seq_len(n_sim), breaks = n_workers, labels = FALSE)
    )
  }
  
  if (n_workers <= 1L) {
    values <- lapply(
      chunks, simulate_chunk,
      seeds = seeds, p = p, m_df = m_df, q_df = q_df
    )
  } else if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      c("simulate_one_lmax", "simulate_chunk"),
      envir = .GlobalEnv
    )
    values <- parallel::parLapply(
      cl, chunks, simulate_chunk,
      seeds = seeds, p = p, m_df = m_df, q_df = q_df
    )
  } else {
    values <- parallel::mclapply(
      chunks, simulate_chunk,
      seeds = seeds, p = p, m_df = m_df, q_df = q_df,
      mc.cores = n_workers,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  }
  
  lambda_max <- as.numeric(unlist(values, use.names = FALSE))
  if (length(lambda_max) != n_sim || any(!is.finite(lambda_max))) {
    stop("Simulation returned an invalid number of finite roots.", call. = FALSE)
  }
  
  saveRDS(
    list(
      lambda_max = lambda_max,
      settings = list(
        p = p, m_df = m_df, q_df = q_df,
        n_sim = n_sim, seed = seed
      )
    ),
    cache_file
  )
  lambda_max
}

choose_exact_backend <- function(prm, requested = "auto") {
  requested <- match.arg(
    requested,
    c("auto", "fast_double", "scaled_double", "robust_arbitrary")
  )
  if (requested != "auto") return(requested)
  
  # Jacobi boundary parameters (m_C=0 or n_C=0) can make both the
  # historical and scaled double-precision paths collapse numerically.
  # Use the robust arbitrary-precision Pfaffian path even when s is small.
  if (min(prm$m, prm$n) <= 0.25) return("robust_arbitrary")
  if (prm$s >= 50L) return("robust_arbitrary")
  "fast_double"
}

exact_cdf_lambda <- function(lambda, p, m_df, q_df,
                             backend = "auto") {
  lambda <- as.numeric(lambda)
  if (any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("lambda must be finite and nonnegative.", call. = FALSE)
  }
  
  prm <- dsb_params(p, m_df, q_df)
  theta <- lambda / (1 + lambda)
  selected <- choose_exact_backend(prm, backend)
  
  if (selected == "fast_double") {
    ans <- rootWishartHD::doubleWishart(
      theta,
      s = prm$s,
      m = prm$m,
      n = prm$n,
      type = "double",
      verbose = FALSE
    )
  } else if (selected == "scaled_double") {
    ans <- rootWishartHD::doubleWishart_scaled(
      theta,
      s = prm$s,
      m = prm$m,
      n = prm$n,
      type = "double",
      verbose = FALSE,
      direct_b = TRUE,
      scale_iter = 8L,
      pf_method = "auto",
      adaptive = FALSE
    )
  } else {
    # High-s and near-saturated Jacobi parameters require the robust scaled
    # Pfaffian calculation. Work on the log-CDF scale and exponentiate only
    # after the stable evaluation.
    log_ans <- rootWishartHD::doubleWishart_robustPfaffians(
      theta,
      s = prm$s,
      m = prm$m,
      n = prm$n,
      type = "arbitrary",
      tail = "lower",
      verbose = FALSE,
      direct_b = TRUE,
      scale_iter = 8L,
      adaptive = FALSE,
      start_digits10 = 300L,
      max_digits10 = 20000L,
      tol = 1e-20
    )
    ans <- exp(log_ans)
    
    # Fallback for isolated backend failures.
    bad <- !is.finite(ans)
    if (any(bad)) {
      ans[bad] <- rootWishartHD::doubleWishart_scaled(
        theta[bad],
        s = prm$s,
        m = prm$m,
        n = prm$n,
        type = "arbitrary",
        verbose = FALSE,
        direct_b = TRUE,
        scale_iter = 8L,
        pf_method = "auto",
        adaptive = FALSE,
        start_digits10 = 300L,
        max_digits10 = 20000L,
        tol = 1e-20
      )
    }
  }
  
  ans <- as.numeric(ans)
  if (length(ans) != length(lambda) || any(!is.finite(ans))) {
    stop(
      sprintf(
        "Exact CDF evaluation failed for (p,m,q)=(%d,%d,%d) using %s.",
        p, m_df, q_df, selected
      ),
      call. = FALSE
    )
  }
  pmin(1, pmax(0, ans))
}

exact_cdf_lambda_cached <- function(lambda, p, m_df, q_df,
                                    backend, cache_dir,
                                    force = FALSE) {
  key <- sprintf(
    "exact_curve_p%d_m%d_q%d_G%d_%s.rds",
    p, m_df, q_df, length(lambda), backend
  )
  path <- file.path(cache_dir, key)
  
  if (file.exists(path) && !isTRUE(force)) {
    obj <- readRDS(path)
    if (isTRUE(all.equal(obj$lambda, lambda, tolerance = 0))) {
      message("Using exact-CDF cache: ", basename(path))
      return(obj$cdf)
    }
  }
  
  cdf <- exact_cdf_lambda(
    lambda, p = p, m_df = m_df, q_df = q_df, backend = backend
  )
  saveRDS(
    list(
      lambda = lambda,
      cdf = cdf,
      settings = list(
        p = p, m_df = m_df, q_df = q_df, backend = backend
      )
    ),
    path
  )
  cdf
}

# Johnstone (2008) Jacobi largest-root approximation on the logit(theta)=log(lambda)
# scale, followed by the real Tracy-Widom CDF F_1.
tw_cdf_lambda <- function(lambda, p, m_df, q_df) {
  lambda <- as.numeric(lambda)
  if (any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("lambda must be finite and nonnegative.", call. = FALSE)
  }
  
  df2 <- p - m_df + q_df
  denom <- df2 + m_df - 1
  gamma <- 2 * asin(sqrt((min(q_df, m_df) - 0.5) / denom))
  phi <- 2 * asin(sqrt((max(q_df, m_df) - 0.5) / denom))
  mu <- 2 * log(tan(0.5 * (phi + gamma)))
  sigma <- (16 / denom^2)^(1 / 3) *
    (sin(phi + gamma)^2 * sin(phi) * sin(gamma))^(-1 / 3)
  
  ans <- numeric(length(lambda))
  positive <- lambda > 0
  if (any(positive)) {
    z <- (log(lambda[positive]) - mu) / sigma
    ans[positive] <- RMTstat::ptw(z, beta = 1)
  }
  pmin(1, pmax(0, ans))
}

ks_discrepancy <- function(sample_values, cdf_function) {
  x <- sort(as.numeric(sample_values))
  n <- length(x)
  F_x <- as.numeric(cdf_function(x))
  if (length(F_x) != n || any(!is.finite(F_x))) {
    stop("CDF evaluation failed during discrepancy calculation.", call. = FALSE)
  }
  upper_jump <- seq_len(n) / n
  lower_jump <- (seq_len(n) - 1) / n
  max(abs(F_x - upper_jump), abs(F_x - lower_jump))
}

make_lambda_grid <- function(lambda_max, n_grid, x_scale,
                             cdf_range = "full",
                             cdf_quantiles = c(0.0005, 0.9995),
                             padding = 0.03) {
  x_scale <- match.arg(x_scale, c("real", "log10"))
  cdf_range <- match.arg(cdf_range, c("full", "quantile"))
  
  values <- as.numeric(lambda_max)
  values <- values[is.finite(values) & values >= 0]
  if (!length(values)) stop("No finite nonnegative roots were simulated.")
  
  if (cdf_range == "full") {
    limits <- range(values)
  } else {
    limits <- as.numeric(stats::quantile(
      values, probs = cdf_quantiles, names = FALSE, type = 8
    ))
  }
  
  if (x_scale == "log10") {
    positive <- values[values > 0]
    if (!length(positive)) {
      stop("A logarithmic x-axis requires at least one positive root.")
    }
    lower <- max(limits[1], min(positive))
    upper <- max(limits[2], lower * (1 + 1e-8))
    lower <- lower / (1 + padding)
    upper <- upper * (1 + padding)
    grid <- 10^seq(log10(lower), log10(upper), length.out = n_grid)
  } else {
    span <- diff(limits)
    if (!is.finite(span) || span <= 0) {
      span <- max(limits[2], 1) * 0.05
    }
    lower <- max(0, limits[1] - padding * span)
    upper <- limits[2] + padding * span
    grid <- seq(lower, upper, length.out = n_grid)
  }
  
  unique(as.numeric(grid))
}

transform_lambda_for_plot <- function(lambda, x_scale) {
  x_scale <- match.arg(x_scale, c("real", "log10"))
  if (x_scale == "real") return(as.numeric(lambda))
  log10(as.numeric(lambda))
}

lambda_axis_label <- function(x_scale) {
  x_scale <- match.arg(x_scale, c("real", "log10"))
  if (x_scale == "real") expression(lambda[max])
  else expression(log[10](lambda[max]))
}

make_curve_data <- function(lambda_max, p, m_df, q_df, n_grid,
                            exact_backend = "auto",
                            cache_dir = NULL,
                            cache_exact = FALSE,
                            force_exact = FALSE,
                            x_scale = "real",
                            cdf_range = "full",
                            cdf_quantiles = c(0.0005, 0.9995),
                            padding = 0.03) {
  grid <- make_lambda_grid(
    lambda_max = lambda_max,
    n_grid = n_grid,
    x_scale = x_scale,
    cdf_range = cdf_range,
    cdf_quantiles = cdf_quantiles,
    padding = padding
  )
  
  empirical <- stats::ecdf(lambda_max)(grid)
  exact <- if (isTRUE(cache_exact)) {
    if (is.null(cache_dir)) stop("cache_dir is required when cache_exact=TRUE.")
    exact_cdf_lambda_cached(
      grid, p, m_df, q_df,
      backend = exact_backend,
      cache_dir = cache_dir,
      force = force_exact
    )
  } else {
    exact_cdf_lambda(
      grid, p, m_df, q_df, backend = exact_backend
    )
  }
  tw <- tw_cdf_lambda(grid, p, m_df, q_df)
  
  out <- rbind(
    data.frame(p = p, m_df = m_df, q_df = q_df,
               lambda = grid, cdf = empirical,
               method = "Empirical CDF"),
    data.frame(p = p, m_df = m_df, q_df = q_df,
               lambda = grid, cdf = exact,
               method = "Exact double-Wishart"),
    data.frame(p = p, m_df = m_df, q_df = q_df,
               lambda = grid, cdf = tw,
               method = "Tracy-Widom")
  )
  out$x_plot <- transform_lambda_for_plot(out$lambda, x_scale)
  out$x_scale <- x_scale
  out$cdf_range <- cdf_range
  out
}

# ------------------------------ run analysis ---------------------------------
lmax_by_p <- setNames(vector("list", length(P_GRID)), as.character(P_GRID))
for (p in P_GRID) {
  lmax_by_p[[as.character(p)]] <- simulate_for_p(
    p = p,
    m_df = M_DF,
    q_df = Q_DF,
    n_sim = N_SIM,
    n_cores = N_CORES,
    seed = SEED,
    cache_dir = CACHE_DIR,
    force = FORCE_RECOMPUTE
  )
}

curve_data <- do.call(
  rbind,
  lapply(P_SHOW, function(p) {
    make_curve_data(
      lmax_by_p[[as.character(p)]],
      p = p,
      m_df = M_DF,
      q_df = Q_DF,
      n_grid = N_CURVE_GRID,
      exact_backend = EXACT_BACKEND,
      cache_dir = CACHE_DIR,
      cache_exact = FALSE,
      force_exact = FORCE_RECOMPUTE_EXACT,
      x_scale = MAIN_X_SCALE,
      cdf_range = CDF_RANGE,
      cdf_quantiles = CDF_QUANTILES,
      padding = CDF_X_PADDING
    )
  })
)
curve_data$p_label <- factor(
  paste0("p = ", curve_data$p),
  levels = paste0("p = ", P_SHOW)
)

# Quantitative CDF discrepancies on the full tested p-grid.
discrepancy_rows <- lapply(P_GRID, function(p) {
  sample_values <- lmax_by_p[[as.character(p)]]
  d_exact <- ks_discrepancy(
    sample_values,
    function(x) exact_cdf_lambda(
      x, p, M_DF, Q_DF, backend = EXACT_BACKEND
    )
  )
  d_tw <- ks_discrepancy(
    sample_values,
    function(x) tw_cdf_lambda(x, p, M_DF, Q_DF)
  )
  data.frame(
    p = p,
    method = c("Exact double-Wishart", "Tracy-Widom"),
    discrepancy = c(d_exact, d_tw)
  )
})
discrepancy_data <- do.call(rbind, discrepancy_rows)
mc_reference <- 1.36 / sqrt(N_SIM)

if (isTRUE(CHECK_EXACT_AGREEMENT)) {
  exact_check <- subset(
    discrepancy_data, method == "Exact double-Wishart"
  )
  bad <- exact_check$discrepancy > EXACT_DISCREPANCY_MAX
  if (any(bad)) {
    failed_p <- paste(exact_check$p[bad], collapse = ", ")
    stop(
      paste0(
        "The exact CDF failed its simulation agreement check for p = ",
        failed_p, ". This indicates a numerical backend failure, not failure ",
        "of the distributional identity. Boundary cases such as p=m+1 have ",
        "n_C=0 and must use the robust arbitrary-precision backend. ",
        "Try EXACT_BACKEND='robust_arbitrary', update rootWishartHD, and set ",
        "FORCE_RECOMPUTE_EXACT=TRUE."
      ),
      call. = FALSE
    )
  }
}

# Separate near-saturated case: p, m, and q differ by only one.
stress_lambda <- simulate_for_p(
  p = STRESS_P,
  m_df = STRESS_M_DF,
  q_df = STRESS_Q_DF,
  n_sim = N_SIM_STRESS,
  n_cores = N_CORES,
  seed = SEED_STRESS,
  cache_dir = CACHE_DIR,
  force = FORCE_RECOMPUTE
)
stress_curve_data <- make_curve_data(
  stress_lambda,
  p = STRESS_P,
  m_df = STRESS_M_DF,
  q_df = STRESS_Q_DF,
  n_grid = N_STRESS_GRID,
  exact_backend = "robust_arbitrary",
  cache_dir = CACHE_DIR,
  cache_exact = TRUE,
  force_exact = FORCE_RECOMPUTE_EXACT,
  x_scale = STRESS_X_SCALE,
  cdf_range = CDF_RANGE,
  cdf_quantiles = CDF_QUANTILES,
  padding = CDF_X_PADDING
)

# Grid-based discrepancies 
stress_emp <- subset(stress_curve_data, method == "Empirical CDF")
stress_exact <- subset(stress_curve_data, method == "Exact double-Wishart")
stress_tw <- subset(stress_curve_data, method == "Tracy-Widom")
if (!isTRUE(all.equal(stress_emp$lambda, stress_exact$lambda, tolerance = 0)) ||
    !isTRUE(all.equal(stress_emp$lambda, stress_tw$lambda, tolerance = 0))) {
  stop("Stress-case CDF grids are not aligned.", call. = FALSE)
}
stress_d_exact_grid <- max(
  abs(stress_emp$cdf - stress_exact$cdf), na.rm = TRUE
)
stress_d_tw_grid <- max(
  abs(stress_emp$cdf - stress_tw$cdf), na.rm = TRUE
)
stress_mc_reference <- 1.36 / sqrt(N_SIM_STRESS)

method_levels <- c("Empirical CDF", "Exact double-Wishart", "Tracy-Widom")
curve_data$method <- factor(curve_data$method, levels = method_levels)
stress_curve_data$method <- factor(
  stress_curve_data$method, levels = method_levels
)
discrepancy_data$method <- factor(
  discrepancy_data$method,
  levels = method_levels
)
discrepancy_data$p_label <- factor(
  as.character(discrepancy_data$p),
  levels = as.character(P_GRID)
)

method_colors <- c(
  "Empirical CDF" = "#000000",
  "Exact double-Wishart" = "#0072B2",
  "Tracy-Widom" = "#D55E00"
)
method_linetypes <- c(
  # The empirical curve is dotted.
  "Empirical CDF" = "dotted",
  "Exact double-Wishart" = "solid",
  "Tracy-Widom" = "longdash"
)

base_theme <- ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 9.5),
    axis.text = ggplot2::element_text(size = 8),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 9, face = "bold"),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 8),
    legend.key.width = grid::unit(1.5, "lines"),
    plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 7.5, hjust = 0.5),
    plot.margin = ggplot2::margin(4, 5, 4, 5)
  )

panel_a <- ggplot2::ggplot(
  curve_data,
  ggplot2::aes(x = x_plot, y = cdf, color = method, linetype = method)
) +
  # Draw the two analytic curves.
  ggplot2::geom_line(
    data = subset(curve_data, method != "Empirical CDF"),
    linewidth = 0.78
  ) +
  # Draw the empirical CDF.  
  ggplot2::geom_step(
    data = subset(curve_data, method == "Empirical CDF"),
    linewidth = 0.62,
    direction = "hv"
  ) +
  ggplot2::facet_wrap(~p_label, nrow = 1, scales = "free_x") +
  ggplot2::scale_color_manual(values = method_colors, breaks = method_levels) +
  ggplot2::scale_linetype_manual(values = method_linetypes, breaks = method_levels) +
  ggplot2::coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  ggplot2::labs(
    x = lambda_axis_label(MAIN_X_SCALE),
    y = "CDF"
  ) +
  base_theme +
  ggplot2::theme(legend.position = "top")

panel_b <- ggplot2::ggplot(
  discrepancy_data,
  ggplot2::aes(
    x = p_label, y = discrepancy,
    color = method, linetype = method, group = method
  )
) +
  ggplot2::geom_hline(
    yintercept = mc_reference,
    color = "grey35",
    linetype = "dotted",
    linewidth = 0.55
  ) +
  ggplot2::geom_line(linewidth = 0.72) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::annotate(
    "text",
    x = length(P_GRID),
    y = mc_reference,
    label = expression(1.36 / sqrt(N)),
    hjust = 1,
    vjust = -0.55,
    size = 2.7,
    color = "grey25"
  ) +
  ggplot2::scale_color_manual(values = method_colors, breaks = method_levels) +
  ggplot2::scale_linetype_manual(values = method_linetypes, breaks = method_levels) +
  ggplot2::labs(
    x = expression(p),
    y = expression(sup[x] * " | " * hat(F)[N](x) - F(x) * " |")
  ) +
  base_theme +
  ggplot2::theme(legend.position = "none")

stress_subtitle <- sprintf(
  "(p,m,q)=(%d,%d,%d); grid Dexact=%.3f, DTW=%.3f",
  STRESS_P, STRESS_M_DF, STRESS_Q_DF,
  stress_d_exact_grid, stress_d_tw_grid
)

panel_c <- ggplot2::ggplot(
  stress_curve_data,
  ggplot2::aes(x = x_plot, y = cdf, color = method, linetype = method)
) +
  # Draw analytic curves first, then the empirical CDF on top.
  ggplot2::geom_line(
    data = subset(stress_curve_data, method != "Empirical CDF"),
    linewidth = 0.78
  ) +
  ggplot2::geom_step(
    data = subset(stress_curve_data, method == "Empirical CDF"),
    linewidth = 0.62,
    direction = "hv"
  ) +
  ggplot2::scale_color_manual(values = method_colors, breaks = method_levels) +
  ggplot2::scale_linetype_manual(values = method_linetypes, breaks = method_levels) +
  ggplot2::coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  ggplot2::labs(
    title = "Near-saturated rank regime",
    subtitle = stress_subtitle,
    x = lambda_axis_label(STRESS_X_SCALE),
    y = "CDF"
  ) +
  base_theme +
  ggplot2::theme(legend.position = "none")

bottom_row <- patchwork::wrap_plots(
  panel_b, panel_c,
  nrow = 1,
  widths = c(1, 1.08)
)

figure_1 <- patchwork::wrap_plots(
  panel_a, bottom_row,
  ncol = 1,
  heights = c(1.55, 1.12)
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 11)
    )
  )

# ------------------------------ write outputs --------------------------------
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

ggplot2::ggsave(
  filename = file.path(FIG_DIR, "Figure1.pdf"),
  plot = figure_1,
  device = pdf_device,
  width = FIG_WIDTH_IN,
  height = FIG_HEIGHT_IN,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(FIG_DIR, "Figure1.tiff"),
  plot = figure_1,
  device = "tiff",
  compression = "lzw",
  dpi = FIG_DPI,
  width = FIG_WIDTH_IN,
  height = FIG_HEIGHT_IN,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(FIG_DIR, "Figure1.png"),
  plot = figure_1,
  device = "png",
  dpi = FIG_DPI,
  width = FIG_WIDTH_IN,
  height = FIG_HEIGHT_IN,
  units = "in"
)

utils::write.csv(
  curve_data,
  file.path(FIG_DIR, "Figure1_curve_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  discrepancy_data,
  file.path(FIG_DIR, "Figure1_discrepancy.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    settings = list(
      m_df = M_DF,
      q_df = Q_DF,
      p_grid = P_GRID,
      p_show = P_SHOW,
      n_sim = N_SIM,
      seed = SEED,
      exact_backend = EXACT_BACKEND,
      main_x_scale = MAIN_X_SCALE,
      stress_x_scale = STRESS_X_SCALE,
      cdf_range = CDF_RANGE,
      cdf_quantiles = CDF_QUANTILES,
      cdf_x_padding = CDF_X_PADDING,
      stress = list(
        p = STRESS_P,
        m_df = STRESS_M_DF,
        q_df = STRESS_Q_DF,
        n_sim = N_SIM_STRESS,
        seed = SEED_STRESS
      )
    ),
    lambda_max = lmax_by_p,
    stress_lambda_max = stress_lambda,
    curve_data = curve_data,
    stress_curve_data = stress_curve_data,
    discrepancy_data = discrepancy_data,
    mc_reference = mc_reference,
    stress_mc_reference = stress_mc_reference,
    stress_d_exact_grid = stress_d_exact_grid,
    stress_d_tw_grid = stress_d_tw_grid,
    figure = figure_1
  ),
  file.path(FIG_DIR, "Figure1_results.rds")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(FIG_DIR, "Figure1_sessionInfo.txt")
)
writeLines(
  c(
    sprintf("M_DF=%d", M_DF),
    sprintf("Q_DF=%d", Q_DF),
    sprintf("P_GRID=%s", paste(P_GRID, collapse = ",")),
    sprintf("P_SHOW=%s", paste(P_SHOW, collapse = ",")),
    sprintf("N_SIM=%d", N_SIM),
    sprintf("SEED=%d", SEED),
    sprintf("STRESS_P=%d", STRESS_P),
    sprintf("STRESS_M_DF=%d", STRESS_M_DF),
    sprintf("STRESS_Q_DF=%d", STRESS_Q_DF),
    sprintf("N_SIM_STRESS=%d", N_SIM_STRESS),
    sprintf("SEED_STRESS=%d", SEED_STRESS),
    sprintf("N_CORES=%d", N_CORES),
    sprintf("EXACT_BACKEND=%s", EXACT_BACKEND),
    sprintf("MAIN_X_SCALE=%s", MAIN_X_SCALE),
    sprintf("STRESS_X_SCALE=%s", STRESS_X_SCALE),
    sprintf("CDF_RANGE=%s", CDF_RANGE),
    sprintf("CDF_QUANTILES=%s", paste(CDF_QUANTILES, collapse = ",")),
    sprintf("CDF_X_PADDING=%.8g", CDF_X_PADDING),
    sprintf("CHECK_EXACT_AGREEMENT=%s", CHECK_EXACT_AGREEMENT),
    sprintf("EXACT_DISCREPANCY_MAX=%.8g", EXACT_DISCREPANCY_MAX),
    sprintf("MC_REFERENCE=%.8g", mc_reference),
    sprintf("STRESS_MC_REFERENCE=%.8g", stress_mc_reference),
    sprintf("STRESS_D_EXACT_GRID=%.8g", stress_d_exact_grid),
    sprintf("STRESS_D_TW_GRID=%.8g", stress_d_tw_grid)
  ),
  file.path(FIG_DIR, "Figure1_settings.txt")
)

# Console summary.
exact_values <- subset(discrepancy_data, method == "Exact double-Wishart")$discrepancy
tw_values <- subset(discrepancy_data, method == "Tracy-Widom")$discrepancy

cat("\nFigure 1 written to:\n")
cat("  ", file.path(FIG_DIR, "Figure1.pdf"), "\n", sep = "")
cat("  ", file.path(FIG_DIR, "Figure1.tiff"), "\n", sep = "")
cat("  ", file.path(FIG_DIR, "Figure1.png"), "\n", sep = "")
cat("\nValues for the manuscript TODOs:\n")
cat(sprintf("  N = %d replicates per p\n", N_SIM))
cat(sprintf("  representative p values = %s\n", paste(P_SHOW, collapse = ", ")))
cat(sprintf(
  "  D_exact range = %.5f to %.5f\n",
  min(exact_values), max(exact_values)
))
cat(sprintf(
  "  D_TW range = %.5f to %.5f\n",
  min(tw_values), max(tw_values)
))
cat(sprintf("  1.36/sqrt(N) = %.5f\n", mc_reference))
cat(sprintf(
  "  stress case = (p,m,q)=(%d,%d,%d), N_stress=%d\n",
  STRESS_P, STRESS_M_DF, STRESS_Q_DF, N_SIM_STRESS
))
cat(sprintf(
  "  stress grid D_exact = %.5f; stress grid D_TW = %.5f\n",
  stress_d_exact_grid, stress_d_tw_grid
))
cat(sprintf(
  "  stress 1.36/sqrt(N_stress) = %.5f\n\n",
  stress_mc_reference
))

print(discrepancy_data)
print(figure_1)

