# Run from repository root: source('script/test/test_permuco_followups.R', encoding='UTF-8')
source('script/test/test_helper.R',encoding='UTF-8')
set.seed(391); n <- 24L
data <- expand.grid(ID=seq_len(n),Time=c('t1','t2','t3'),Condition=c('c1','c2'))
data$Group <- rep(rep(c('g1','g2'),c(10,14)),6)
data$Cov <- rnorm(n)[data$ID]; data$Cov2 <- rnorm(n)[data$ID]
data$Y <- rnorm(nrow(data)) + rnorm(n)[data$ID] + .4*data$Cov
make <- function(d,within,between) list(data=d,id='ID',dv='Y',iv=c(between,within),
  within=within,between=between,covariate=c('Cov','Cov2'),design='test')
experiment <- make(data,c('Time','Condition'),'Group')
run <- function(e,method='permuco',inter=TRUE,posthoc=TRUE,adjust='holm')
  analyze(e,method=method,select='arguments',show=FALSE,np=99,seed=10,
    covariate_interactions=inter,followup_scope='all',posthoc=posthoc,adjust=adjust)
check <- function(r) {
  c <- r$Y$contrasts; stopifnot(length(c)>0,all(vapply(c,function(z)is.null(z$error),logical(1))))
  for(z in c) {
    t <- z$formatted
    stopifnot(all(is.finite(t$SE)),all(t$SE>0),all(t$p.value>=t$p.uncorrected),
      all(t$p.uncorrected>=1/99),all(t$p.uncorrected<=1),
      max(abs(t$t.ratio-t$estimate/t$SE))<1e-7,
      isTRUE(all.equal(t$p.value,p.adjust(t$p.uncorrected,'holm'))))
  }
}
check(run(experiment))
for(inter in c(FALSE,TRUE)) check(run(experiment,'permuco_ancova',inter))
pure <- make(subset(data,Group=='g2'),c('Time','Condition'),character())
check(run(pure)); check(run(pure,'permuco_ancova',FALSE));check(run(pure,'permuco_ancova',TRUE))
between <- make(subset(data,Time=='t1' & Condition=='c1'),character(),'Group')
check(run(between));check(run(between,'permuco_ancova',FALSE));check(run(between,'permuco_ancova',TRUE))
a <- run(experiment); b <- run(experiment)
stopifnot(identical(lapply(a$Y$contrasts,`[[`,'formatted'),lapply(b$Y$contrasts,`[[`,'formatted')))
set.seed(123); before <- .Random.seed; invisible(run(experiment));stopifnot(identical(before,.Random.seed))
stopifnot(length(run(experiment,posthoc=FALSE)$Y$contrasts)==0L)
err <- try(run(experiment,adjust='tukey'),silent=TRUE);stopifnot(inherits(err,'try-error'))
cat('HC2 permutation follow-ups: all integration tests passed.\n')
