# HC2-studentized permutation contrasts of the fixed model.
# Freedman-Lane residual transformations; subject blocks preserve repeated measures.
.t_permuco_followups <- function(data, formula, id, iv, within, covariate, np, seed) {
  if (!requireNamespace('emmeans', quietly = TRUE) || !requireNamespace('sandwich', quietly = TRUE))
    stop('置換下位検定にはemmeansとsandwichが必要です。', call. = FALSE)
  fixed_formula <- formula
  # Remove Error() through terms, without changing fixed-effect interactions.
  tt <- terms(formula, specials = 'Error')
  special <- attr(tt, 'specials')$Error
  if (length(special)) fixed_formula <- update(formula,
    paste('. ~ . -', deparse(attr(tt, 'variables')[[special + 1L]], width.cutoff = 500)))
  fit <- stats::lm(fixed_formula, data = data, contrasts = setNames(
    lapply(iv, function(x) contr.sum(nlevels(data[[x]]))), iv))
  X <- model.matrix(fit); y <- model.response(model.frame(fit))
  if (fit$rank != ncol(X)) stop('置換下位検定のモデルが推定不能です。', call. = FALSE)
  groups <- split(seq_len(nrow(data)), factor(data[[id]]), drop = TRUE)
  if (length(unique(lengths(groups))) != 1L)
    stop('参加者ごとに同じ反復測定条件が必要です。', call. = FALSE)
  # Data from t_permuco is ordered by ID and within factors before this call.
  A <- solve(crossprod(X)); Q <- qr.Q(qr(X)); beta <- coef(fit)
  cov <- sandwich::vcovCL(fit, cluster = data[[id]], type = 'HC2', cadjust = TRUE)
  adjusters <- lapply(groups, function(ix) {
    eg <- eigen(diag(length(ix)) - X[ix,,drop=FALSE] %*% A %*% t(X[ix,,drop=FALSE]), symmetric=TRUE)
    if (min(eg$values) <= 1e-10) stop('参加者の影響が大きすぎてHC2を計算できません。', call.=FALSE)
    eg$vectors %*% ((1/sqrt(eg$values)) * t(eg$vectors))
  })
  # This part of each participant's HC2 score does not depend on the contrast.
  score_design <- lapply(seq_along(groups), function(k)
    adjusters[[k]] %*% X[groups[[k]],,drop=FALSE] %*% A)
  rm(adjusters)
  at <- if (length(covariate)) setNames(as.list(rep(0,length(covariate))),covariate) else NULL
  compute <- function(effects, forced, main_effects) {
    if (!is.null(seed)) {
      had <- exists('.Random.seed', .GlobalEnv, inherits=FALSE)
      if (had) old <- get('.Random.seed', .GlobalEnv)
      on.exit(if(had) assign('.Random.seed',old,.GlobalEnv) else if(exists('.Random.seed',.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv))
      set.seed(seed)
    }
    result <- list()
    for (effect in effects) {
      factors <- strsplit(gsub('`','',effect,fixed=TRUE),':',fixed=TRUE)[[1]]
      if (!all(factors %in% iv) || length(factors)>2L) next
      directions <- if(length(factors)==1L) factors else if(effect %in% forced) intersect(factors,main_effects) else factors
      for (compared in directions) {
        fixed <- setdiff(factors,compared)
        spec <- as.formula(paste('~',paste0('`',compared,'`'),if(length(fixed))paste('|',paste0('`',fixed,'`'))else ''))
        grid <- emmeans::contrast(emmeans::emmeans(fit,spec,at=at,vcov.=cov), 'pairwise',adjust='none')
        tab <- as.data.frame(summary(grid,infer=c(FALSE,FALSE)))
        Ls <- grid@linfct
        tab$comparison <- as.character(tab$contrast); tab$contrast <- NULL
        tab$df <- NA_real_; tab$p.value <- NA_real_; tab$t.ratio <- NA_real_
        tab$permutation_scheme <- if(compared %in% within) 'subject_signflip' else 'subject_permutation'
        tab$np <- np
        tab$standard_error <- 'participant_HC2'
        for (j in seq_len(nrow(tab))) {
          L <- as.numeric(Ls[j,]); denom <- as.numeric(crossprod(L,A%*%L))
          if (!is.finite(denom) || denom<=0) stop('推定不能な対比です。',call.=FALSE)
          b0 <- beta - as.vector(A%*%L)*sum(L*beta)/denom
          f0 <- as.vector(X%*%b0); residual <- y-f0
          weights <- lapply(score_design,function(z) as.vector(z %*% L))
          statistic <- function(Y) {
            b <- A %*% crossprod(X,Y); E <- Y-Q%*%crossprod(Q,Y)
            scores <- lapply(seq_along(groups),function(k) as.vector(crossprod(weights[[k]],E[groups[[k]],,drop=FALSE])))
            se <- sqrt(Reduce('+',lapply(scores,function(z)z*z)))
            if(any(!is.finite(se)|se<=0)) stop('置換後のHC2標準誤差を計算できません。',call.=FALSE)
            as.vector(crossprod(L,b))/se
          }
          observed <- statistic(matrix(y,ncol=1))
          if (!isTRUE(all.equal(as.numeric(tab$estimate[j]/tab$SE[j]),observed,tolerance=1e-7)))
            stop('HC2の観測統計量がemmeansと一致しません。',call.=FALSE)
          exceed <- 0L
          # B=np-1 random transformations plus the observed statistic.
          for (start in seq.int(1L,np-1L,by=128L)) {
            count <- min(128L,np-start)
            Y <- replicate(count, {
              z <- residual
              if(compared %in% within) {
                signs <- sample(c(-1,1),length(groups),replace=TRUE)
                for(k in seq_along(groups)) z[groups[[k]]] <- residual[groups[[k]]]*signs[k]
              } else {
                order <- sample.int(length(groups))
                for(k in seq_along(groups)) z[groups[[k]]] <- residual[groups[[order[k]]]]
              }
              f0+z
            })
            exceed <- exceed + sum(abs(statistic(matrix(Y,nrow=length(y))))>=abs(observed)-1e-10)
          }
          tab$t.ratio[j] <- observed; tab$p.value[j] <- (1+exceed)/np
        }
        # Preserve each fixed level as its own multiplicity family.
        strata <- if(length(fixed)) split(seq_len(nrow(tab)),tab[[fixed]],drop=TRUE) else list(all=seq_len(nrow(tab)))
        for (level in names(strata)) {
          ix <- strata[[level]]; name <- paste(c(effect,if(length(factors)>1)compared,if(length(fixed))level),collapse='/')
          result[[name]] <- list(kind='permuco',comparison_grid=grid[ix],unadjusted=tab[ix,,drop=FALSE],
            formatted=tab[ix,,drop=FALSE],compared_factor=compared,source_effect=effect,
            source_slot='HC2',error=NULL,np=np,standard_error='participant_HC2',ci_adjustment='not_computed')
        }
      }
    }
    result
  }
  compute
}
