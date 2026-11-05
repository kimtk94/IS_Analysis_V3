#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: 43_stage3d3_mr.R <harmonized.tsv[.gz]> <outdir>")
}

input <- args[1]
outdir <- args[2]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

dat <- read.delim(input, stringsAsFactors=FALSE, check.names=FALSE)

required <- c(
  "gene","mode","trait","domain","variant_id",
  "beta_exposure","se_exposure","beta_outcome","se_outcome"
)
miss <- setdiff(required, names(dat))
if (length(miss) > 0) stop(paste("Missing:", paste(miss, collapse=", ")))

dat <- dat[
  is.finite(dat$beta_exposure) &
  is.finite(dat$se_exposure) &
  is.finite(dat$beta_outcome) &
  is.finite(dat$se_outcome) &
  dat$se_exposure > 0 &
  dat$se_outcome > 0 &
  dat$beta_exposure != 0,
]

norm_p <- function(z) {
  2 * pnorm(abs(z), lower.tail=FALSE)
}

weighted_median <- function(theta, w) {
  ord <- order(theta)
  theta <- theta[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  theta[which(cw >= 0.5)[1]]
}

wm_boot <- function(bx, by, sey, B=2000, seed=20260926) {
  set.seed(seed)
  theta <- by / bx
  w <- 1 / ((sey / abs(bx))^2)
  est <- weighted_median(theta, w)

  boots <- numeric(B)
  for (b in seq_len(B)) {
    byb <- rnorm(length(by), mean=by, sd=sey)
    thb <- byb / bx
    boots[b] <- weighted_median(thb, w)
  }

  se <- sd(boots, na.rm=TRUE)
  z <- est / se
  c(beta=est, se=se, p=norm_p(z))
}

ivw_calc <- function(bx, by, sey) {
  w <- 1/(sey^2)
  den <- sum(w*bx^2)
  beta <- sum(w*bx*by)/den
  se_fixed <- sqrt(1/den)

  resid <- by - beta*bx
  Q <- sum(w*resid^2)
  df <- length(bx)-1
  Qp <- ifelse(df > 0, pchisq(Q, df=df, lower.tail=FALSE), NA_real_)

  phi <- ifelse(df > 0, max(1, Q/df), 1)
  se_random <- se_fixed * sqrt(phi)

  list(
    beta=beta,
    se_fixed=se_fixed,
    p_fixed=norm_p(beta/se_fixed),
    se_random=se_random,
    p_random=norm_p(beta/se_random),
    Q=Q,
    Q_df=df,
    Q_p=Qp,
    phi=phi
  )
}

egger_calc <- function(bx, by, sey) {
  fit <- lm(by ~ bx, weights=1/(sey^2))
  co <- summary(fit)$coefficients
  list(
    intercept=co["(Intercept)","Estimate"],
    intercept_se=co["(Intercept)","Std. Error"],
    intercept_p=co["(Intercept)","Pr(>|t|)"],
    slope=co["bx","Estimate"],
    slope_se=co["bx","Std. Error"],
    slope_p=co["bx","Pr(>|t|)"]
  )
}

keys <- unique(dat[c("gene","mode","trait","domain")])

mr <- list()
loo <- list()
ix <- 0
jx <- 0

for (i in seq_len(nrow(keys))) {
  g <- keys$gene[i]
  m <- keys$mode[i]
  t <- keys$trait[i]
  d <- keys$domain[i]

  x <- dat[
    dat$gene == g &
    dat$mode == m &
    dat$trait == t &
    dat$domain == d,
  ]

  # Deduplicate variants.
  x <- x[!duplicated(x$variant_id),]
  k <- nrow(x)
  if (k == 0) next

  bx <- x$beta_exposure
  by <- x$beta_outcome
  sey <- x$se_outcome

  base <- data.frame(
    gene=g, mode=m, trait=t, domain=d, k=k,
    method="", beta=NA_real_, se=NA_real_, p=NA_real_,
    Q=NA_real_, Q_df=NA_real_, Q_p=NA_real_,
    egger_intercept=NA_real_, egger_intercept_se=NA_real_,
    egger_intercept_p=NA_real_, notes="",
    stringsAsFactors=FALSE
  )

  if (k == 1) {
    beta <- by[1]/bx[1]
    se <- sey[1]/abs(bx[1])
    r <- base
    r$method <- "WALD"
    r$beta <- beta
    r$se <- se
    r$p <- norm_p(beta/se)
    r$notes <- "single independent cis-pQTL"
    ix <- ix + 1
    mr[[ix]] <- r
  }

  if (k >= 2) {
    v <- ivw_calc(bx,by,sey)

    r <- base
    r$method <- "IVW_FIXED"
    r$beta <- v$beta
    r$se <- v$se_fixed
    r$p <- v$p_fixed
    r$Q <- v$Q
    r$Q_df <- v$Q_df
    r$Q_p <- v$Q_p
    r$notes <- "fixed-effect IVW"
    ix <- ix + 1
    mr[[ix]] <- r

    r <- base
    r$method <- "IVW_MRE"
    r$beta <- v$beta
    r$se <- v$se_random
    r$p <- v$p_random
    r$Q <- v$Q
    r$Q_df <- v$Q_df
    r$Q_p <- v$Q_p
    r$notes <- paste0("multiplicative random-effects IVW; phi=", signif(v$phi,5))
    ix <- ix + 1
    mr[[ix]] <- r
  }

  if (k >= 3) {
    wm <- wm_boot(bx,by,sey)
    r <- base
    r$method <- "WEIGHTED_MEDIAN"
    r$beta <- unname(wm["beta"])
    r$se <- unname(wm["se"])
    r$p <- unname(wm["p"])
    r$notes <- "parametric bootstrap SE under NOME approximation"
    ix <- ix + 1
    mr[[ix]] <- r
  }

  if (k >= 10) {
    e <- egger_calc(bx,by,sey)
    r <- base
    r$method <- "MR_EGGER"
    r$beta <- e$slope
    r$se <- e$slope_se
    r$p <- e$slope_p
    r$egger_intercept <- e$intercept
    r$egger_intercept_se <- e$intercept_se
    r$egger_intercept_p <- e$intercept_p
    r$notes <- "Egger enabled only for k>=10"
    ix <- ix + 1
    mr[[ix]] <- r
  }

  if (k >= 3) {
    for (q in seq_len(k)) {
      xx <- x[-q,]
      v <- ivw_calc(xx$beta_exposure, xx$beta_outcome, xx$se_outcome)
      jx <- jx + 1
      loo[[jx]] <- data.frame(
        gene=g, mode=m, trait=t, domain=d,
        omitted_variant=x$variant_id[q],
        k_remaining=nrow(xx),
        beta_ivw=v$beta,
        se_ivw_mre=v$se_random,
        p_ivw_mre=v$p_random,
        stringsAsFactors=FALSE
      )
    }
  }
}

mr_df <- if (length(mr)) do.call(rbind,mr) else data.frame()
loo_df <- if (length(loo)) do.call(rbind,loo) else data.frame()

write.table(
  mr_df,
  file=file.path(outdir,"STAGE3D3_MR_RESULTS.tsv"),
  sep="\t", row.names=FALSE, quote=FALSE, na="NA"
)

write.table(
  loo_df,
  file=file.path(outdir,"STAGE3D3_LEAVE_ONE_OUT.tsv"),
  sep="\t", row.names=FALSE, quote=FALSE, na="NA"
)

cat("MR rows =", nrow(mr_df), "\n")
cat("LOO rows =", nrow(loo_df), "\n")
