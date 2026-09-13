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
