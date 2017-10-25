
check_input <- function(n){
  if(is.null(n)){
    stop("data is NULL")
  }
  if(!is.numeric(n)){
    stop("expecting numeric values for n")
  }
  if(any(is.na(n))){
    stop("data contains NA")
  }
  if(any(n < 0)){
    stop("negative abundances not permitted")
  }
}

#' Fits the double geometric distribution
#'
#' @param n Vector of species abundances
#' @param trials Number of interactions when fitting function.
#' @param min Minimum bound on total species richness
#' @param max Maximum bound on total species richness
#'
#' @return A list with the following components:
#' \describe{
#' \item{distribution}{The XXX}
#' \item{richness}{The XXX}
#' \item{k}{The parameter XXX}
#' }
#' @export
#' @references
#' Alroy, J. 2015. The shape of terrestrial abundance distributions. Science Advances 1(8): e1500082
#'
#' @examples
#' \dontrun{
#' fitDoubleGeometric(n)
#' }
fitDoubleGeometric <- function(n,
                               trials = 1000,
                               min = NA,
                               max = NA)	{
  n <- sort(n, decreasing = TRUE)
  spp <- length(n)
  if (is.na(min))
    min <- spp
  if (is.na(max))
    max <- 10 * spp
  p <- n / sum(n)
  best <- 999999
  k <- 0.5
  bestk <- 1
  r <- min
  bestr <- min - 1
  while (r <= max && r == bestr + 1)	{
    for (z in seq_len(trials))	{
      lastk <- k
      k <- k + rnorm(1, sd = 1 / z)
      if (k <= 0)
        k <- lastk
      dg <- array()
      dg[1] <- k
      for (i in 2:r)	{
        if (i > r - i + 1)	{
          dg[i] <- dg[i - 1] * k ** (1 / sqrt((r - i + 1) / r))
        } else	{
          dg[i] <- dg[i - 1] * k ** (1 / sqrt(i / r))
        }
      }
      raw <- dg
      dg <- dg / sum(dg)
      kl <- sum(p * log(p / dg[seq_len(spp)]))
      if (kl < best)	{
        best <- kl
        bestk <- k
        bestr <- r
        bestdg <- dg
      }
      if (bestk != k)	{
        k <- lastk
      }
    }
    r <- r + 1
  }
  return(list(
    distribution = bestdg,
    richness = bestr,
    k = bestk,
    fit = best
  ))
}

#' Fits the log normal distribution
#'
#' @inheritParams fitDoubleGeometric
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#'
#' @return A list with the following components:
#' \describe{
#' \item{distribution}{The XXX}
#' \item{richness}{The XXX}
#' \item{sd}{The parameter XXX}
#' }#' @export
#'
#' @examples
#' \dontrun{
#' fitLogNormal(n)
#' }
fitLogNormal <- function(n,
                         trials = 1000,
                         min = NA,
                         max = NA)	{
  n <- sort(n, decreasing = TRUE)
  spp <- length(n)
  if (is.na(min))
    min <- spp
  if (is.na(max))
    max <- 10 * spp
  p <- n / sum(n)
  best <- 999999
  s <- 2
  r <- min
  bestr <- min - 1
  while (r <= max && r == bestr + 1)	{
    for (z in 1:trials)	{
      lasts <- s
      s <- s + rnorm(1, sd = 1 / z)
      if (s <= 0)	{
        s <- lasts
      }
      np <- exp(sort(
                  qnorm(seq(0.5 / r, 1 - 0.5 / r, 1 / r), sd = s),
                  decreasing = TRUE))
      np <- np / sum(np)
      kl <- sum(p * log(p / np[1:spp]))
      if (kl < best)	{
        best <- kl
        bests <- s
        bestr <- r
        bestln <- np
      }
      if (bests != s)	{
        s <- lasts
      }
    }
    r <- r + 1
  }
  return(list(
    distribution = bestln,
    richness = bestr,
    sd = bests,
    fit = best
  ))
}

#' Fits the geometric series distribution
#'
#' @param n Vector of species abundances
#'
#' @return A list with the following components:
#' \describe{
#' \item{distribution}{The XXX}
#' \item{k}{The parameter XXX}
#' \item{fit}{The XXX}
#' }
#' @export
#' @importFrom stats lm
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' fitGeometric(n)
#' }
fitGeometric <- function(n)	{
  n <- sort(n, decreasing = TRUE)
  spp <- length(n)
  p <- n / sum(n)
  k <- exp(lm(log(n) ~ c(seq_len(spp)), weights = sqrt(n))$coef[2])
  best <- 999999
  for (i in seq_len(100))	{
    lastk <- k
    k <- k + rnorm(1, sd = 0.01)
    q <- array()
    for (i in seq_len(10 * spp))	{
      if (i > 1)
        q[i] <- q[i - 1] * k
      else
        q[i] <- k
    }
    q <- q / sum(q)
    kl <- sum(p * log(p / q[seq_len(spp)]))
    if (kl < best)	{
      best <- kl
      bestk <- k
      bestgs <- q[seq_len(spp)]
    }
    if (bestk != k)	{
      k <- lastk
    }
  }
  return(list(
    distribution = bestgs,
    k = as.numeric(bestk),
    fit = best
  ))
}

#' Fits the log series distribution
#'
#' @inheritParams fitGeometric
#'
#' @return A list with the following components:
#' \describe{
#' \item{distribution}{The XXX}
#' \item{alpha}{The XXX}
#' \item{fit}{The parameter XXX}
#' }
#' @export
#' @importFrom stats rnorm
#'
#' @examples
#' \dontrun{
#' fitLogSeries(n)
#' }
fitLogSeries <- function(n)	{
  n <- sort(n, decreasing = TRUE)
  spp <- length(n)
  N <- sum(n)
  p <- n / N
  if (N > 10000)
    return(list(
      distribution = NA,
      alpha = NA,
      fit = NA
    ))
  z <- logSeriesParams(n)
  alpha <- z[1]
  lsx <- z[2]
  best <- 999999
  besta <- alpha
  for (i in seq_len(1000))	{
    lasta <- alpha
    alpha <- alpha + rnorm(1, sd = 1)
    if (alpha <= 0)
      alpha <- lasta
    lsx <- N / (N + alpha)
    s <- alpha * lsx ^ (seq_len(3 * N)) / (seq_len(3 * N))
    cd <- round(cumsum(s))
    z <- 0
    q <- array()
    for (j in seq_len(cd[1]))	{
      z <- z + 1
      q[z] <- 1
    }
    for (breakpt in which(diff(cd) > 0))
      for (j in cd[breakpt]:(cd[breakpt + 1] - 1))	{
        z <- z + 1
        q[z] <- breakpt + 1
      }
    if (length(q) < spp)	{
      alpha <- lasta
      next
    }
    q <- sort(q, decreasing = TRUE)
    q <- q / sum(q)
    kl <- sum(p * log(p / q[seq_len(spp)]))
    if (kl < best)	{
      best <- kl
      besta <- alpha
      bestls <- q[seq_len(spp)]
    }
    if (besta != alpha)	{
      alpha <- lasta
    }
  }
  return(list(
    distribution = bestls,
    alpha = besta,
    fit = best
  ))
}

#' Helper function for fitLogSeries
#'
#' @inheritParams fitGeometric
#'
#' @return A vector
logSeriesParams <- function(n)	{
  a <- 10
  olda <- 0
  N <- sum(n)
  z <- 0
  while (abs(a - olda) > 0.0000001 && z < 1000)	{
    z <- z + 1
    olda <- a
    a <- length(n) / log(1 + N / a)
  }
  x <- N / (N + a)
  return(c(a, x))
}


#' Fits the broken stick distribution
#'
#' @inheritParams fitGeometric
#'
#' @return A list with the following components:
#' \describe{
#' \item{distribution}{The XXX}
#' \item{richness}{The XXX}
#' \item{fit}{The parameter XXX}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' fitBS(n)
#' }
fitBS <- function(n)	{
  n <- sort(n, decreasing = TRUE)
  best <- 999999
  spp <- length(n)
  bestr <- spp
  p <- n / sum(n)
  for (i in spp:(spp * 3))	{
    q <- rep(0, i)
    for (j in seq_len(i))	{
      for (k in 0:(i - j))	{
        q[j] <- q[j] + 1 / (i - k)
      }
      q[j] <- q[j] / i
    }
    q <- q / sum(q, na.rm = TRUE)
    kl <- sum(p * log(p / q[seq_len(spp)]))
    if (kl < best)	{
      best <- kl
      bestr <- i
      bestbs <- q
    }
    if (bestr != i)	{
      break
    }
  }
  return(list(
    distribution = bestbs,
    richness = bestr,
    fit = best
  ))
}

#' Fits the Zipf distribution
#'
#' @inheritParams fitGeometric
#'
#' @return A list with the following components:
#' \describe{
#' \item{distribution}{The XXX}
#' \item{exponent}{The XXX}
#' \item{fit}{The XXX}
#' }
#' @export
#' @importFrom stats lm
#' @importFrom stats rnorm
#'
#' @examples
#' \dontrun{
#' fitZipf(n)
#' }
fitZipf <- function(n)	{
  check_input(n)
  n <- sort(n, decreasing = TRUE)
  spp <- length(n)
  p <- n / sum(n)
  co <- lm(log(n) ~ log(seq_len(spp)), weights = sqrt(n))$coefficients
  z <- co[2]
  best <- 999999
  for (i in seq_len(100))	{
    lastexp <- z
    z <- z + rnorm(1, sd = 0.01)
    q <- exp(log(seq_len(1000)) * z)
    q <- q / sum(q)
    kl <- sum(p * log(p / q[seq_len(spp)]))
    if (kl < best)	{
      best <- kl
      bestexp <- z
      bestzipf <- q[seq_len(spp)]
    }
    if (bestexp != z)	{
      z <- lastexp
    }
  }
  return(list(
    distribution = bestzipf,
    exponent = as.numeric(bestexp),
    fit = best
  ))
}
