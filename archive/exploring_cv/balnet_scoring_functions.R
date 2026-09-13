get_balance_loss <-
function (object, X.test, W.test, sample.weights, lambda, ...) 
{
    .balance_loss <- function(W, eta) {
        colSums(sample.weights * (W * exp(-eta) + (1 - W) * eta))/sum(sample.weights)
    }
    lambda <- validate_lambda(lambda)
    eta <- predict(object, X.test, lambda = lambda, type = "link", 
        .simplify = FALSE)
    loss0 <- loss1 <- NULL
    if (!is.null(object[["_fit"]]$control)) {
        loss0 <- .balance_loss(1 - W.test, eta$control)
    }
    if (!is.null(object[["_fit"]]$treated)) {
        loss1 <- .balance_loss(W.test, eta$treated)
    }
    out <- list(control = loss0, treated = loss1)
    out[!vapply(out, is.null, logical(1))]
}
get_imbalance <-
function (object, X.test, W.test, sample.weights, lambda, X.stats, 
    type.measure = "imbalance.mean", W.hat = NULL, ...) 
{
    .imbalance <- function(W, W.hat) {
        W.hat[W == 1, ] <- pmax(W.hat[W == 1, ], 0.001)
        ipw <- matrix(0, nrow = nrow(W.hat), ncol = ncol(W.hat))
        ipw[W == 1, ] <- 1/W.hat[W == 1, ]
        smd <- col_stats(X.test, ipw * sample.weights, n_threads = object[["num.threads"]])$center
        smd <- sweep(smd, 2L, X.stats$center, `-`, check.margin = FALSE)
        smd <- sweep(smd, 2L, X.stats$scale, `/`, check.margin = FALSE)
        if (type.measure == "imbalance.mean") {
            return(rowMeans(abs(smd)))
        }
        else if (type.measure == "imbalance.inf") {
            return(apply(abs(smd), 1, max))
        }
        else {
            stop("Invalid type.measure norm")
        }
    }
    if (is.null(W.hat)) {
        lambda <- validate_lambda(lambda)
        W.hat <- predict(object, X.test, lambda = lambda, type = "response", 
            .simplify = FALSE)
    }
    loss0 <- loss1 <- NULL
    if (!is.null(object[["_fit"]]$control)) {
        loss0 <- .imbalance(1 - W.test, 1 - W.hat$control)
    }
    if (!is.null(object[["_fit"]]$treated)) {
        loss1 <- .imbalance(W.test, W.hat$treated)
    }
    out <- list(control = loss0, treated = loss1)
    out[!vapply(out, is.null, logical(1))]
}
