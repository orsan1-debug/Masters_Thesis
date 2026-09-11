## CV-PATCH TEST, session-only: no package or pipeline files touched ----------
source("R/dgp.R"); source("R/estimators.R")
source("R/simulate.R"); source("R/registry.R")


patch_cv_balnet <- function() {
  target <- function(e, head)
    is.call(e) && identical(e[[1]], as.name(head)) && length(e) == 2L &&
    is.call(e[[2]]) && identical(e[[2]][[1]], quote(matrix)) &&
    length(e[[2]]) == 4L &&
    identical(e[[2]][[if (head == "colMeans") 3L else 4L]], quote(nfolds))
  walk <- function(e, fun) {                  # guards inline: never store a child
    if (is.call(e)) {
      fun(e)
      for (i in seq_along(e)) {
        if (is.null(e[[i]])) next
        if (is.name(e[[i]]) && !nzchar(as.character(e[[i]]))) next
        walk(e[[i]], fun)
      }
    }
  }
  scan <- function(e, head) { n <- 0L; walk(e, \(x) if (target(x, head)) n <<- n + 1L); n }
  arity_ok <- function(e) { ok <- TRUE
  walk(e, \(x) if (is.name(x[[1]]) &&
                   as.character(x[[1]]) %in% c("<-", "=", "<<-") &&
                   length(x) != 3L) ok <<- FALSE); ok }
  rewrite <- function(e) {
    if (is.call(e)) {
      if (target(e, "colMeans")) {
        e[[1]] <- quote(rowMeans)
        m <- e[[2]]; tmp <- m[[3]]; m[[3]] <- m[[4]]; m[[4]] <- tmp
        e[[2]] <- m
        return(e)
      }
      for (i in seq_along(e)) {
        if (is.null(e[[i]])) next
        if (is.name(e[[i]]) && !nzchar(as.character(e[[i]]))) next
        e[[i]] <- rewrite(e[[i]])
      }
    }
    e
  }
  
  f <- getFromNamespace("cv.balnet", "balnet")
  b <- body(f)
  if (!arity_ok(b))
    stop("cv.balnet body is corrupted from an earlier patch attempt: restart R, then rerun.")
  if (scan(b, "rowMeans") == 2L) return(invisible(TRUE))     # already patched
  if (scan(b, "colMeans") != 2L)
    stop("cv.balnet body unexpected (version change?): nothing modified.")
  
  body(f) <- rewrite(b)
  stopifnot(scan(body(f), "colMeans") == 0L,
            scan(body(f), "rowMeans") == 2L,
            arity_ok(body(f)))                               # self-verify before install
  assignInNamespace("cv.balnet", f, ns = "balnet")
  if ("package:balnet" %in% search()) {
    pe <- as.environment("package:balnet")
    unlockBinding("cv.balnet", pe); assign("cv.balnet", f, envir = pe)
    lockBinding("cv.balnet", pe)
  }
  invisible(TRUE)
}

## 1. Patch main session, gate with the reprex
patch_cv_balnet()


## Gate: patched session must store the fold mean and select from it
set.seed(1)
n <- 400; p <- 10; nfolds <- 5
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, plogis(X[, 1]))
foldid <- sample(rep(seq_len(nfolds), length.out = n))

cvfit <- cv.balnet(X, W, nfolds = nfolds, foldid = foldid, max.imbalance = 1e-2)
lam   <- cvfit$`_lambda`
M <- sapply(seq_len(nfolds), \(k) {
  tr <- foldid != k
  f  <- balnet(X[tr, ], W[tr], standardize = ".inplace", max.imbalance = 1e-2)
  balnet:::get_balance_loss(f, X[!tr, , drop = FALSE], W[!tr], rep(1, sum(!tr)), lam)$treated
})

stopifnot(isTRUE(all.equal(cvfit$`_cv.info`$cv.mean$treated, rowMeans(M))),
          identical(cvfit$lambda.min$treated, lam$treated[which.min(rowMeans(M))]))










