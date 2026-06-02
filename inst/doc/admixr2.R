## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 5
)

## ----data---------------------------------------------------------------------
library(admixr2)
library(rxode2)
library(nlmixr2)

data("examplomycin")
head(examplomycin[examplomycin$EVID == 0, c("ID", "TIME", "DV")], 9)

## ----aggregate----------------------------------------------------------------
obs   <- examplomycin[examplomycin$EVID == 0, ]
obs   <- obs[order(obs$ID, obs$TIME), ]
times <- sort(unique(obs$TIME))
ids   <- unique(obs$ID)
n     <- length(ids)                     # 500

dv_mat <- matrix(NA_real_, nrow = n, ncol = length(times))
for (i in seq_along(ids)) {
  sub         <- obs[obs$ID == ids[i], ]
  dv_mat[i, ] <- sub$DV[order(sub$TIME)]
}

E <- colMeans(dv_mat)
V <- cov.wt(dv_mat, method = "ML")$cov

round(E, 2)

## ----model--------------------------------------------------------------------
pk_model <- function() {
  ini({
    tcl     <- log(5)  ; label("Log clearance (L/hr)")
    tv1     <- log(10) ; label("Log central volume (L)")
    tv2     <- log(30) ; label("Log peripheral volume (L)")
    tq      <- log(10) ; label("Log inter-compartmental CL (L/hr)")
    tka     <- log(1)  ; label("Log absorption rate constant (1/hr)")
    prop.sd <- c(0, 0.2); label("Proportional residual error SD")
    eta.cl ~ 0.09
    eta.v1 ~ 0.09
    eta.v2 ~ 0.09
    eta.q  ~ 0.09
    eta.ka ~ 0.09
  })
  model({
    cl <- exp(tcl + eta.cl)
    v1 <- exp(tv1 + eta.v1)
    v2 <- exp(tv2 + eta.v2)
    q  <- exp(tq  + eta.q)
    ka <- exp(tka + eta.ka)
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <- ka * depot - (cl/v1 + q/v1) * central + (q/v2) * peripheral
    d/dt(peripheral) <- (q/v1) * central - (q/v2) * peripheral
    cp <- central / v1
    cp ~ prop(prop.sd)
  })
}

## ----study--------------------------------------------------------------------
study <- list(
  E     = E,
  V     = V,                       # full 9x9 covariance matrix
  n     = n,
  times = times,
  ev    = rxode2::et(amt = 100)    # single 100 mg oral dose
)

## ----fit----------------------------------------------------------------------
fit <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies   = list(examplomycin = study),
    n_sim     = 5000L,
    cov_n_sim = 10000L,
    maxeval   = 300L,
    seed      = 1L
  )
)

## ----print--------------------------------------------------------------------
print(fit)

## ----internals----------------------------------------------------------------
fit$objective                    # -2 log-likelihood
fit$env$admExtra$struct          # structural parameters (log scale)
fit$env$admExtra$omega           # estimated Omega matrix
fit$env$admExtra$sigma_var       # residual variance(s)

logLik(fit)
AIC(fit)

## ----plot, fig.cap="Left: observed vs predicted mean with residuals. Right: NLL convergence trace."----
plots <- plot(fit, which = c("mean", "nll"))

