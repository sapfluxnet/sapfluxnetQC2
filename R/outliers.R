################################################################################
#' Outliers substitution
#'
#' Outlier substitution by different methods
#'
#' \code{k} is the semi-value of the desired window for substitute value
#' calculation. The window is formed by the range \code{y[i] - k : y[i] + k}.
#'
#' @section out_tukeyline:
#'   \code{out_tukeyline} performs a robust fitting with the \code{line} function
#'
#' @section out_medianreg:
#'   \code{out_medianreg} performs a regression of medians with the quantreg package
#'
#' @section out_median:
#'   \code{out_median} performs a classic median outlier detection
#'
#' @family outliers
#'
#' @param y vector of values for outlier substitution
#'
#' @param k window semi-value (integer) for substitution value calculation. See
#'   details.
#'
#' @return a vector of the same lengh of y with the outlier values substituted
#'   by the calculated value
#'
#' @name outlier_subs
NULL

################################################################################
#' @describeIn outlier_subs
#'
#' @export

# START
# Function declaration
out_tukeyline <- function(y, k = 5L, parent_logger = 'test') {

  # Using calling handlers to manage errors
  withCallingHandlers({
    # STEP 0
    # Argument checks
    # y numeric
    if (!is.numeric(y)) {
      stop('Vector of values provided is not numeric')
    }
    # k integer
    if (!is.integer(k)) {
      stop('Window value is not an integer')
    }

    # STEP 1
    # Initiate needed vectors and values
    m <- length(y)
    x <- -k:k
    r <- vector('numeric', m)

    # STEP 2
    # Iteration loop to calculate outliers substitutions
    for (i in (k + 1):(m - k)) {
      z <- y[(i - k):(i + k)]
      r[i] <- ifelse(sum(!is.na(z)) > 2,
                     line(x = x,y = z)$fitted.values[k + 1],
                     NA)
    }

    # STEP 3
    # Return the res vector
    return(r)

    # END FUNCTION
  },

  # handlers
  warning = function(w){logging::logwarn(w$message,
                                         logger = paste(parent_logger,
                                                        'out_tukeyline', sep = '.'))},
  error = function(e){logging::logerror(e$message,
                                        logger = paste(parent_logger,
                                                       'out_tukeyline', sep = '.'))},
  message = function(m){logging::loginfo(m$message,
                                         logger = paste(parent_logger,
                                                        'out_tukeyline', sep = '.'))})
}

################################################################################
#' @describeIn outlier_subs
#'
#' @export

# START
# Function declaration
out_medianreg <- function(y, k = 5L, parent_logger = 'test') {

  # Using calling handlers to manage errors
  withCallingHandlers({
    # STEP 0
    # Argument checks
    # y numeric
    if (!is.numeric(y)) {
      stop('Vector of values provided is not numeric')
    }
    # k integer
    if (!is.integer(k)) {
      stop('Window value is not an integer')
    }

    # STEP 1
    # Initiate needed vectors and values
    m <- length(y)
    x <- -k:k
    r <- vector('numeric', m)

    # STEP 2
    # Iteration loop to calculate outliers substitutions
    for (i in (k+1):(m-k)) {
      z <- y[(i - k):(i + k)]
      r[i] <- ifelse(sum(!is.na(z)) > 5,
                     quantreg::rq(z ~ x, .5)$fitted.values[k + 1],
                     NA)
    }

    # STEP 3
    # Return the res vector
    return(r)

    # END FUNCTION
  },

  # handlers
  warning = function(w){logging::logwarn(w$message,
                                         logger = paste(parent_logger,
                                                        'out_medianreg', sep = '.'))},
  error = function(e){logging::logerror(e$message,
                                        logger = paste(parent_logger,
                                                       'out_medianreg', sep = '.'))},
  message = function(m){logging::loginfo(m$message,
                                         logger = paste(parent_logger,
                                                        'out_medianreg', sep = '.'))})
}

################################################################################
#' @describeIn outlier_subs
#'
#' @export

# START
# Function declaration
out_median <- function(y, k = 5L, parent_logger = 'test') {

  # Using calling handlers to manage errors
  withCallingHandlers({
    # STEP 0
    # Argument checks
    # y numeric
    if (!is.numeric(y)) {
      stop('Vector of values provided is not numeric')
    }
    # k integer
    if (!is.integer(k)) {
      stop('Window value is not an integer')
    }

    # STEP 1
    # Initiate needed objects
    n <- length(y)

    # STEP 2
    # results
    return(c(
      y[1:k],
      sapply((k + 1):(n - k), function(j) median(y[(j - k):(j + k)],
                                                 na.rm = TRUE)),
      y[(n-k+1):n]
    ))

    # END FUNCTION
  },

  # handlers
  warning = function(w){logging::logwarn(w$message,
                                         logger = paste(parent_logger,
                                                        'out_median', sep = '.'))},
  error = function(e){logging::logerror(e$message,
                                        logger = paste(parent_logger,
                                                       'out_median', sep = '.'))},
  message = function(m){logging::loginfo(m$message,
                                         logger = paste(parent_logger,
                                                        'out_median', sep = '.'))})
}

################################################################################
#' Hampel filter
#'
#' Hampel filter for detecting and substituting outliers in the environmental
#' and sap flow data.
#'
#' This is a modified version of the Hampel filter to avoid breaking when NAs
#' are present in the data. Options include the posibility of using the
#' median (classic Hampel, default), tukeyline or quantile regression estimations.
#'
#' @family outliers
#'
#' @param y vector of values for outlier substitution
#'
#' @param k window semi-value (integer) for substitution value calculation. See
#'   details.
#'
#' @param t0
#'
#' @param method Character vector indicating the method to use in the estimation:
#'   \code{hampel} (default), \code{tukey} and \code{quantile}.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item{res: A vector of the same length as y with the ouliers value
#'     substituted by the estimation}
#'     \item{index: A vector indicating which TIMESTAMPS have been modified}
#'   }
#'
#' @export

# START
# Function declaration
out_hampel_filter <- function(y, k = 5L, t0 = 3L,
                              method = 'hampel', parent_logger = 'test') {

  # Using calling handlers to manage errors
  withCallingHandlers({

    # STEP 0
    # Argument checks
    # y numeric
    if (!is.numeric(y)) {
      stop('Vector of values provided is not numeric')
    }
    # k integer
    if (!is.integer(k)) {
      stop('Window value is not an integer')
    }
    # t0 integer
    if (!is.integer(t0)) {
      stop('T0 value provided is not an integer')
    }

    # STEP 1
    # Initialising needed objects
    n <- length(y)
    z <- y
    y0 <- switch(method,
                 hampel = out_median(y, k, parent_logger),
                 tukey = out_tukeyline(y, k, parent_logger),
                 quantile = out_medianreg(y, k, parent_logger))
    y0.na <- !is.na(y0)
    y.na <- !is.na(y)
    ind <- NULL
    L <- 1.4826

    # STEP 2
    # Iteration loop
    for (i in (k + 1):(n - k)) {
      S0 <- L * median(abs(y[(i - k):(i + k)] - y0[i]), na.rm = TRUE)
      if (y0.na[i] & y.na[i] & !is.na(S0) & abs(y[i] - y0[i]) > t0 * S0) {
        z[i] <- y0[i]
        ind <- c(ind, i)
      }
    }

    # STEP 3
    # Returning the res list
    return(list(res = z, index = ind))

    # END FUNCTION

  },

  # handlers
  warning = function(w){logging::logwarn(w$message,
                                         logger = paste(parent_logger,
                                                        'out_hampel_filter',
                                                        sep = '.'))},
  error = function(e){logging::logerror(e$message,
                                        logger = paste(parent_logger,
                                                       'out_hampel_filter',
                                                       sep = '.'))},
  message = function(m){logging::loginfo(m$message,
                                         logger = paste(parent_logger,
                                                        'out_hampel_filter',
                                                        sep = '.'))})
}

################################################################################
#' Outliers QC
#'
#' Outliers detection, substitution and annotation for SfnData objects
#'
#' Outliers for sap flow data and environmental data are located and substituted
#' by the selected method, using the \code{\link{out_hampel_filter}} function.
#'
#' @family outliers
#'
#' @param sfn_data SfnData object with the site data and metadata
#'
#' @param k window semi-value (integer) for substitution value calculation. See
#'   \code{\link{out_hampel_filter}} for details.
#'
#' @param t0
#'
#' @param method Character vector indicating the method to use in the outlier
#'   estimation: \code{hampel} (default), \code{tukey} and \code{quantile}.
#'
#' @return a SfnData object as the one provided with the oulier values
#'   substituted and the flags slot updated.
#'
#' @export

# START
# Function declaration
out_remove <- function(sfn_data, k = 5L, t0 = 3L,
                       method = 'hampel', parent_logger = 'test') {

  # Using calling handlers to manage errors
  withCallingHandlers({

    # STEP 0
    # Argument checks
    # SfnData
    if (!is(sfn_data, 'SfnData')) {
      stop('sfn_data object provided is not an SfnData object')
    }
    # k integer
    if (!is.integer(k)) {
      stop('k is not an integer')
    }
    # t0 integer
    if (!is.integer(t0)) {
      stop('t0 is not an integer')
    }
    # method character
    if (!is.character(method)) {
      stop('method is not a character')
    }

    # STEP 1
    # get needed data (without timestamp)
    sapf_data <- get_sapf(sfn_data)[,-1]
    env_data <- get_env(sfn_data)[,-1]
    sapf_flags <- get_sapf_flags(sfn_data)[,-1]
    env_flags <- get_env_flags(sfn_data)[,-1]

    # STEP 2
    # Apply selected outlier filter
    sapf_out <- lapply(sapf_data, out_hampel_filter,
                       k = k, t0 = t0, method = method)
    env_out <- lapply(env_data, out_hampel_filter,
                      k = k, t0 = t0, method = method)

    # STEP 3
    # Iteration for substituting each column of data with the filtered data

    # 3.1 sapf
    for (i in 1:ncol(sapf_data)) {
      # 3.1.1 substitute values
      sapf_data[, i] <- sapf_out[[i]][[1]]
      # 3.1.2 flags update
      old_flag <- sapf_flags[sapf_out[[i]][[2]], i]

      new_flag <- vapply(old_flag, function(x) {
        if (x == '') {
          "OUT_REPLACED"
        } else { paste0(x, "; OUT_REPLACED") }
      }, character(1), USE.NAMES = FALSE)

      sapf_flags[sapf_out[[i]][[2]], i] <- new_flag
    }

    # 3.2 env
    for (i in 1:ncol(env_data)) {
      # 3.2.1 substitute values
      env_data[, i] <- env_out[[i]][[1]]
      # 3.2.2 flags update
      old_flag <- env_flags[env_out[[i]][[2]], i]

      new_flag <- vapply(old_flag, function(x) {
        if (x == '') {
          "OUT_REPLACED"
        } else { paste0(x, "; OUT_REPLACED") }
      }, character(1), USE.NAMES = FALSE)

      env_flags[env_out[[i]][[2]], i] <- new_flag
    }

    # STEP 4
    # Update sfn_data and return it!!
    get_sapf(sfn_data) <- sapf_data
    get_sapf_flags(sfn_data) <- sapf_flags
    get_env(sfn_data) <- env_data
    get_env_flags(sfn_data) <- env_flags

    return(sfn_data)

    # END FUNCTION

  },

  # handlers
  warning = function(w){logging::logwarn(w$message,
                                         logger = paste(parent_logger,
                                                        'out_remove',
                                                        sep = '.'))},
  error = function(e){logging::logerror(e$message,
                                        logger = paste(parent_logger,
                                                       'out_remove',
                                                       sep = '.'))},
  message = function(m){logging::loginfo(m$message,
                                         logger = paste(parent_logger,
                                                        'out_remove',
                                                        sep = '.'))})
}
