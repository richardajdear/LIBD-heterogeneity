library("qsvaR")
library("jaffelab")
library("limma")

match_rse_tx_samples <- function(rse_tx, rse_gene) {
    # For multi-lane samples, use first entry to match
    retained_samples <- sapply(rse_norm_samples$SAMPLE_ID, `[[`, 1)
    tx_samples <- sapply(rse_tx$SAMPLE_ID, `[[`, 1)
    rse_tx_matched <- rse_tx[, tx_samples %in% retained_samples]
    return(rse_tx_matched)
}

make_qsvs <- function(rse_tx, formula = "~ Age + Sex + Race") {
    # Merge covariates
    rse_tx <- jaffelab::merge_rse_metrics(rse_tx)
    rse_tx$RIN <- sapply(rse_tx$RIN, mean)

    # Design matrix of tx data
    mod_tx <- model.matrix(as.formula(formula),
        data = colData(rse_tx)
    )
    # qSVs of matrix
    qsvs <- qSVA(rse_tx = rse_tx,
                 type = "cell_component",
                 mod = mod_tx,
                 assayname = "tpm"
    )
    return(qsvs)
}

get_residuals <- function(rse,
                          assayname = NULL,
                          formula = "~ Age + Sex + Race",
                          qsvs, n_qsvs = NULL,
                          method = 'libd', P = 0) {
    # Design matrix of data
    mod <- model.matrix(as.formula(formula),
        data = colData(rse)
    )
    # Add n qSVs to design matrix
    if (!is.null(n_qsvs)) {
        qsvs <- qsvs[, 0:n_qsvs]
    }
    mod <- cbind(mod, qsvs)

    # Set assay
    if (!is.null(assayname)) {
        y <- assays(rse, withDimnames = FALSE)[[assayname]]
    } else {
        y <- assay(rse)
    }

    # Fit model and take residuals
    # Compare a few methods for doing this
    if (method == 'limma') {
        fit <- lmFit(y, mod)
        fit <- eBayes(fit)
        res <- residuals(fit, y)
    } else if (method == 'lm') {
        fit <- lm(t(exp) ~ mod)
        res <- t(residuals(fit))
    } else if (method == 'libd') {
        # Least squares: beta = inv(X'X) X' y
        # 'solve' gives the inverse matrix in R
        Hat <- solve(t(mod) %*% mod) %*% t(mod)
        # Set NAs in the expression data to 0
        ty <- t(y)
        ty[is.na(ty)] <- 0
        # Compute beta
        beta <- Hat %*% ty
        # Get estimates of the expression
        # optionally not including the effects of all the variables
        # (P=1, including intercept or not, makes no difference to WGCNA)
        if (P > 0) {
            yhat <- t(as.matrix(mod[, -c(seq_len(P))]) %*% beta[-seq_len(P), ])
        } else {
            yhat <- t(as.matrix(mod) %*% beta)
        }
        # Get residuals of the expression
        res <- y - yhat
    }

    # Add residuals as additional assay to the rse
    assays(rse, withDimnames = FALSE)$residuals <- res
    return(rse)
}