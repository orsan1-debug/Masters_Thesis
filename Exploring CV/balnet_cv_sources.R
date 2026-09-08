balweights.cv.balnet <-
function (object, lambda = "lambda.min", ...) 
{
    if (identical(lambda, "lambda.min")) {
        lambda <- object[["_cv.info"]]$lambda.min
    }
    balweights.balnet(object, lambda = lambda)
}
coef.cv.balnet <-
function (object, lambda = "lambda.min", ...) 
{
    if (identical(lambda, "lambda.min")) {
        lambda <- object[["_cv.info"]]$lambda.min
    }
    coef.balnet(object, lambda = lambda)
}
cv.balnet <-
function (X, W, type.measure = c("balance.loss", "imbalance.mean", 
    "imbalance.inf"), nfolds = 10, foldid = NULL, ...) 
{
    dot.args <- list(...)
    type.measure <- match.arg(type.measure)
    X.stats <- NULL
    if (type.measure == "balance.loss") {
        get_loss <- get_balance_loss
    }
    else if (type.measure %in% c("imbalance.mean", "imbalance.inf")) {
        get_loss <- get_imbalance
        X.stats <- col_stats(X, dot.args[["sample.weights"]], 
            compute_sd = TRUE)
        X.stats$scale[X.stats$scale <= 0] <- 1
    }
    if (is.null(foldid)) {
        nfolds <- max(nfolds, 3)
        foldid <- sample(rep(seq(nfolds), length.out = nrow(X)))
    }
    else {
        if (length(foldid) != length(W)) {
            stop("Invalid `foldid`.")
        }
        nfolds <- max(foldid)
    }
    if (!is.null(dot.args[["verbose"]]) && dot.args[["verbose"]]) 
        cat("Fitting full model\n")
    fit.full <- balnet(X, W, ...)
    lambda.full <- fit.full[["_lambda"]]
    sample.weights <- fit.full[["sample.weights"]]
    cv.list <- list()
    for (k in 1:nfolds) {
        if (!is.null(dot.args[["verbose"]]) && dot.args[["verbose"]]) 
            cat(sprintf("\nFold: %d/%d\n", k, nfolds))
        test <- foldid == k
        train <- !test
        X.train <- X[train, , drop = FALSE]
        W.train <- W[train]
        dot.args[["sample.weights"]] <- sample.weights[train]
        fit.train <- do.call(balnet, c(list(X = X.train, W = W.train, 
            standardize = ".inplace"), dot.args))
        X.test <- X[test, , drop = FALSE]
        W.test <- W[test]
        sample.weights.test <- sample.weights[test]
        loss <- do.call(get_loss, list(fit.train, X.test, W.test, 
            sample.weights.test, lambda.full, X.stats = X.stats, 
            type.measure = type.measure))
        cv.list[[k]] <- loss
    }
    cv.mean0 <- cv.mean1 <- NULL
    idx.min0 <- idx.min1 <- NULL
    lambda.min0 <- lambda.min1 <- NULL
    if (!is.null(cv.list[[1]][["control"]])) {
        cv.mean0 <- colMeans(matrix(unlist(lapply(cv.list, `[[`, 
            "control")), nrow = length(cv.list), ncol = length(lambda.full$control), 
            byrow = TRUE))
        idx.min0 <- which.min(cv.mean0)
        lambda.min0 <- lambda.full[["control"]][idx.min0]
    }
    if (!is.null(cv.list[[1]][["treated"]])) {
        cv.mean1 <- colMeans(matrix(unlist(lapply(cv.list, `[[`, 
            "treated")), nrow = length(cv.list), ncol = length(lambda.full$treated), 
            byrow = TRUE))
        idx.min1 <- which.min(cv.mean1)
        lambda.min1 <- lambda.full[["treated"]][idx.min1]
    }
    lambda.min <- list(control = lambda.min0, treated = lambda.min1)
    lambda.min.out <- lambda.min[!vapply(lambda.min, is.null, 
        logical(1))]
    cv.info <- list(cv.mean = list(control = cv.mean0, treated = cv.mean1), 
        idx.min = list(control = idx.min0, treated = idx.min1), 
        lambda.min = lambda.min, type.measure = type.measure)
    fit.full[["lambda.min"]] <- if (length(lambda.min.out) > 
        1) 
        lambda.min.out
    else lambda.min.out[[1]]
    fit.full[["_cv.info"]] <- cv.info
    fit.full[["call"]] <- match.call()
    class(fit.full) <- c("cv.balnet", class(fit.full))
    fit.full
}
cv.boot.balnet <-
function (X, W, type.measure = c("imbalance.mean", "imbalance.inf"), 
    B = 500, ...) 
{
    dot.args <- list(...)
    type.measure <- match.arg(type.measure)
    get_loss <- get_imbalance
    X.stats <- col_stats(X, dot.args[["sample.weights"]], compute_sd = TRUE)
    X.stats$scale[X.stats$scale <= 0] <- 1
    test.list <- replicate(B, sample.int(length(W), length(W)%/%2), 
        simplify = FALSE)
    if (!is.null(dot.args[["verbose"]]) && dot.args[["verbose"]]) 
        cat("Fitting full model\n")
    fit.full <- balnet(X, W, ...)
    lambda.full <- fit.full[["_lambda"]]
    sample.weights <- fit.full[["sample.weights"]]
    W.hat.full <- predict(fit.full, X, lambda = lambda.full, 
        type = "response", .simplify = FALSE)
    cv.list <- list()
    for (k in 1:length(test.list)) {
        test <- test.list[[k]]
        X.test <- X[test, , drop = FALSE]
        W.test <- W[test]
        sample.weights.test <- sample.weights[test]
        W.hat.test <- lapply(W.hat.full, function(m) m[test, 
            , drop = FALSE])
        loss <- do.call(get_loss, list(fit.full, X.test, W.test, 
            sample.weights.test, lambda.full, X.stats = X.stats, 
            type.measure = type.measure, W.hat = W.hat.test))
        cv.list[[k]] <- loss
    }
    cv.mean0 <- cv.mean1 <- NULL
    idx.min0 <- idx.min1 <- NULL
    lambda.min0 <- lambda.min1 <- NULL
    if (!is.null(cv.list[[1]][["control"]])) {
        cv.mean0 <- colMeans(matrix(unlist(lapply(cv.list, `[[`, 
            "control")), nrow = length(cv.list), ncol = length(lambda.full$control), 
            byrow = TRUE))
        idx.min0 <- which.min(cv.mean0)
        lambda.min0 <- lambda.full[["control"]][idx.min0]
    }
    if (!is.null(cv.list[[1]][["treated"]])) {
        cv.mean1 <- colMeans(matrix(unlist(lapply(cv.list, `[[`, 
            "treated")), nrow = length(cv.list), ncol = length(lambda.full$treated), 
            byrow = TRUE))
        idx.min1 <- which.min(cv.mean1)
        lambda.min1 <- lambda.full[["treated"]][idx.min1]
    }
    lambda.min <- list(control = lambda.min0, treated = lambda.min1)
    lambda.min.out <- lambda.min[!vapply(lambda.min, is.null, 
        logical(1))]
    cv.info <- list(cv.mean = list(control = cv.mean0, treated = cv.mean1), 
        idx.min = list(control = idx.min0, treated = idx.min1), 
        lambda.min = lambda.min, type.measure = paste0(type.measure, 
            ".boot"))
    fit.full[["lambda.min"]] <- if (length(lambda.min.out) > 
        1) 
        lambda.min.out
    else lambda.min.out[[1]]
    fit.full[["_cv.info"]] <- cv.info
    fit.full[["call"]] <- match.call()
    class(fit.full) <- c("cv.boot.balnet", "cv.balnet", class(fit.full))
    fit.full
}
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
plot.cv.balnet <-
function (x, lambda = "lambda.min", ...) 
{
    if (identical(lambda, "lambda.min")) {
        lambda <- x[["_cv.info"]]$lambda.min
    }
    plot.balnet(x, lambda = lambda, ...)
}
predict.cv.balnet <-
function (object, newdata, lambda = "lambda.min", type = c("response"), 
    ...) 
{
    if (identical(lambda, "lambda.min")) {
        lambda <- object[["_cv.info"]]$lambda.min
    }
    predict.balnet(object, newdata, lambda = lambda, type = type)
}
print.cv.balnet <-
function (x, digits = max(3L, getOption("digits") - 3L), ...) 
{
    cat("Call: ", paste(deparse(x$call), collapse = "\n"), "\n\n")
    utils::capture.output(out <- print.balnet(x, digits = digits, 
        .simplify = FALSE, ...))
    df0 <- df1 <- data.frame()
    if (!is.null(x[["_fit"]]$control)) {
        idx.min0 <- x[["_cv.info"]]$idx.min[[1]]
        df0 <- cbind(Arm = "Control", out$control[idx.min0, ], 
            Index = idx.min0)
    }
    if (!is.null(x[["_fit"]]$treated)) {
        idx.min1 <- x[["_cv.info"]]$idx.min[[2]]
        df1 <- cbind(Arm = "Treated", out$treated[idx.min1, ], 
            Index = idx.min1)
    }
    type.measure <- paste0("type.measure = ", x[["_cv.info"]]$type.measure)
    cat("Cross-validated lambda minimizing ", type.measure, ":\n", 
        sep = "")
    print(rbind(df0, df1), digits = digits, row.names = FALSE, 
        right = FALSE)
}
summary.cv.balnet <-
function (object, ...) 
{
    print(object, ...)
}
