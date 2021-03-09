#' @import nnet
#' @import plyr
#' @import data.table
#' @import qvalue
#' @import utils
#' @import EnhancedVolcano



#' @title pA_multi_logit
#' @description Function to compare the usage of alternative poly A site(s) to a reference (canonical) site.
#' @param data Dataset containing poly A (pA) site read counts. This dataset must have a long shape, meaning that there should be only one
#' column containing read counts (and it MUST be named "count"). The first four columns must be called "transcript", "pA.site", "sample" and "count".
#' Thus, each row in \code{data} contains the read count for one pA - transcript - sample combination.
#' Other sample attributes beyond sample ID may be recorded in additional variables in this dataset, or provided separately through a design matrix
#' and a key variable (e.g. sample ID) connecting the \code{data} and \code{design} matrices.
#' @param model Regression model describing the dependence of pA site usage on sample attribute(s).
#' @param design (optional) Design matrix. A matrix describing sample attributes which can be used as predictors in the regression model.
#' @param sample_ID (optional) A key variable connecting the counts dataset (\code{data}) and the design matrix.
#' @param long_output Logical variable describing output format. FALSE: Only regression coefficients and p-values are reported.
#' TRUE: Standard error of regression coefficients, and z scores are also included in the output. Default: FALSE.
#' @param adj_method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @details
#' This function uses a multinomial logistic regression algorithm from the \emph{nnet} package. For each transcript, one poly A site (pA) is set as the
#' canonical (reference) site and the usage of alternative pA(s) is compared to this reference pA. By default, the pA that comes first alphabetically
#' is used as reference. The user can specify the reference pA for each transcript by adding a prefix like \code{0_} to its name.
#' If a transcript has n pA sites, n-1 comparisons will be made. Transcripts with only one pA site should be removed from \code{data}
#' before running this function.
#' @return
#' Log ratios (multinomial logistic regression coefficients) and p-values describing the effect of predictor(s) specified in the model on the usage ratio of each alternative to reference pA site
#' per trancript. If a long output is requested, SE and z scores are also reported.
#' @examples
#' pA.toy2 <- APAlog::pA.toy2
#' pA_design <- APAlog::pA_design
#' fit1_pA <- pA_multi_logit(pA.toy2, pA.site ~ cell_line, pA_design, "sample", adj_method='fdr')
#' @export

pA_multi_logit <- function(data, model, design = NULL, sample_ID = NULL, long_output = FALSE, adj_method){
    data <- remove_0_1_pA_transcripts(data)
    xd <- merge(data, design, by = sample_ID)
    fitx <- suppressWarnings(by(xd, xd$transcript, function(y) nnet::multinom(stats::as.formula(Reduce(paste, deparse(model)), env = new.env()), data = y, weight = count)))
    sfitx <- lapply(fitx, function(x) summary(x))

    sfitxb <- lapply(sfitx, function(x) x$coefficients)
    sfitxse <- lapply(sfitx, function(x) x$standard.errors)
    sfitxz <- lapply(sfitx, function(x) x$coefficients / x$standard.errors)
    sfitxp <- lapply(sfitxz, function(x) (1 - stats::pnorm(abs(x), 0, 1))*2)

    sfitxb.df <- plyr::ldply(sfitxb, rbind, .id = "transcript")
    names(sfitxb.df)[-1] <- paste0("b_", names(sfitxb.df)[-1])
    names(sfitxb.df)[2] <- "b_intercept"

    tr_sites_list <- by(data, data$transcript, function(x) levels(droplevels(x[,strsplit(deparse(model), " ~ ")[[1]][1]])))
    levelpairs <- function(t){
        ref_site <- t[1]
        alt_sites <- t[-1]
        pairs <- expand.grid(ref_site, alt_sites)
        return(pairs)
    }
    tr_site_pairs_list <- lapply(tr_sites_list, levelpairs)
    tr_site_pairs_df <- data.table::rbindlist(tr_site_pairs_list)
    data.table::setnames(tr_site_pairs_df, c("ref_site", "alt_site"))

    sfitxse.df <- do.call(rbind.data.frame, sfitxse)
    names(sfitxse.df) <- paste0("se_", names(sfitxse.df))
    names(sfitxse.df)[1] <- "se_intercept"

    sfitxz.df <- do.call(rbind.data.frame, sfitxz)
    names(sfitxz.df) <- paste0("z_", names(sfitxz.df))
    names(sfitxz.df)[1] <- "z_intercept"

    sfitxp.df <- do.call(rbind.data.frame, sfitxp)
    names(sfitxp.df) <- paste0("p_", names(sfitxp.df))
    names(sfitxp.df)[1] <- "p_intercept"
    if (long_output == FALSE){
        sfitx.com <- data.frame(cbind(transcript = sfitxb.df[,1], tr_site_pairs_df, sfitxb.df[,-1], sfitxp.df))
    } else if (long_output == TRUE) {
        sfitx.com <- data.frame(cbind(transcript = sfitxb.df[,1], tr_site_pairs_df, sfitxb.df[,-1], sfitxse.df, sfitxz.df, sfitxp.df))
    }

    if (adj_method != 'none') {
        sfitx.com <- adj_p(sfitx.com, pcols = c(6,7), adj_method = adj_method)
    }
    return(sfitx.com)
}



#' @title pA_logit_1tr_pairs
#' @description Function to compare the usage of all pairs of pA sites of one transcript across samples
#' @param xdd input data including read counts of all pA sites of one transcript in all samples
#' @param model Regression model
#' @details This function is used inside the pA_logit_pairwise function.
#' @return
#' A dataframe containing log ratios (logistic regression coefficients), SE, z and p-values plus the names of
#' ref and alt pA sites.
#' @examples
#' outxl <- APAlog::pA_logit_1tr_pairs(y, model))
#' @export

pA_logit_1tr_pairs <- function(xdd, model){
    sites <- levels(droplevels(xdd[,strsplit(deparse(model), " ~ ")[[1]][1]]))
    levelpairs <- utils::combn(sites, 2, simplify = FALSE)
    xddp <- lapply(levelpairs, function(x) subset(xdd, xdd[,strsplit(deparse(model), " ~ ")[[1]][1]] %in% x))
    xdd_test <- lapply(xddp, function(y) stats::glm(stats::as.formula(Reduce(paste, deparse(model)), env = new.env()), data = y, family = "binomial", weight = count))
    sxdd_test <- lapply(xdd_test, function(x) summary(x)$coefficients[,c(1,4)])
    sxdd_test_flat <- lapply(sxdd_test, function(x) c(t(x)))
    sxdd_test_flat.df <- plyr::ldply(sxdd_test_flat, rbind)
    d1.names <- c("b", "p")
    d2.names <- rownames(sxdd_test[[1]])
    d2.names[1] <- "intercept"
    names(sxdd_test_flat.df) <- apply(expand.grid(d1.names, d2.names), 1, paste, collapse="_")
    ref_site <- plyr::ldply(lapply(levelpairs, function(x) x[1]), rbind)
    alt_site <- plyr::ldply(lapply(levelpairs, function(x) x[2]), rbind)
    out.df <- data.frame(ref_site, alt_site, sxdd_test_flat.df)
    names(out.df)[seq_len(2)] <- c("ref_site", "alt_site")
    return(out.df)
}


#' @title pA_logit_pairwise
#' @description Function to compare the usage of all pairs of poly sites of a transcript
#' @param data Dataset containing poly A (pA) site read counts. This dataset must have a long shape, meaning that there should be only one
#' column containing read counts (and it MUST be named "count"). The first four columns must be called "transcript", "pA.site", "sample" and "count".
#' Thus, each row in \code{data} contains the read count for one pA - transcript - sample combination.
#' Other sample attributes beyond sample ID may be recorded in additional variables in this dataset, or provided separately through a design matrix
#' and a key variable (e.g. sample ID) connecting the \code{data} and \code{design} matrices.
#' @param model Regression model describing the dependence of pA site usage on sample attribute(s).
#' @param design (optional) Design matrix. A matrix describing sample attributes which can be used as predictors in the regression model.
#' @param sample_ID (optional) A key variable connecting the counts dataset (\code{data}) and the design matrix.
#' @return
#' Log ratios (called "Estimate", logistic regression coefficients) and p-values describing the effect of predictor(s) specified in the model on the usage
#' ratio of pairs of pA sites per trancript. Standard error and z score of estimates are also given.
#' @examples
#' pA.toy2 <- APAlog::pA.toy2
#' pA_design <- APAlog::pA_design
#' fit2_pA <- pA_logit_pairwise(pA.toy2, pA.site~cell_line, pA_design, "sample")
#' @export

pA_logit_pairwise <- function(data, model, design = NULL, sample_ID = NULL){
    data <- remove_0_1_pA_transcripts(data)
    xd <- merge(data, design, by = sample_ID)
    outxl <- suppressWarnings(by(xd, xd$transcript, function(y) pA_logit_1tr_pairs(y, model)))
    outx.df <- plyr::ldply(outxl, rbind, .id = "transcript")
    return(outx.df)
}



#' @title adj_p
#' @description Function to adjust p-values for multiple testing correction.
#' @param x Dataframe containing a column or columns of p-values.
#' @param pcols A vector specifying unadjusted p-value columns in x
#' @param adj_method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @return A data frame containing all the original information with the vector(s) of adjusted p-values
#' appended to the end.
#' @details P-values e.g. coming from differential pA site usage on individual genes/transcripts
#' are correctd for multiple testing.
#' @examples
#' fit3_pA_fdr <- adj_p(fit3_pA, pcols = 2, adj_method = "fdr")
#' @export

adj_p <- function(x, pcols, adj_method){
    y <- x[, pcols, drop = FALSE]
    if (adj_method == "qvalue"){
        y <- apply(y, 2, function(t) qvalue::qvalue(t)$qvalues)
    } else {
        y <- apply(y, 2, function(t) stats::p.adjust(t, method = adj_method))
    }

    newnames <- paste0(adj_method, "_", colnames(y))
    z <- data.frame(x,y)
    colnames(z)[(NCOL(x)+1) : NCOL(z)] <- newnames
    return(z)
}



#' @title glm_deviance_test_p
#' @description Calculate the p-value from a deviance test comparing a model to its corresponding null.
#' @param x Output of a glm run
#' @return P-value calculated from the chisq test of deviance between the model and its corresponding null.
#' @examples
#' dtest_gfitx <- glm_deviance_test_p(x)
#' @export

glm_deviance_test_p <- function(x){
    p <- stats::pchisq((x$null.deviance - x$deviance), df = (x$df.null - x$df.residual), lower.tail = FALSE)
    return(p)
}



#' @title pA_logit_dev
#' @description Function to evaluate the overall effect of predictors on pA usage per transcript through a deviance test.
#' @param data Dataset containing poly A (pA) site read counts. This dataset must have a long shape, meaning that there should be only one
#' column containing read counts (and it MUST be named "count"). The first four columns must be called "transcript", "pA.site", "sample" and "count".
#' Thus, each row in \code{data} contains the read count for one pA - transcript - sample combination.
#' Other sample attributes beyond sample ID may be recorded in additional variables in this dataset, or provided separately through a design matrix
#' and a key variable (e.g. sample ID) connecting the \code{data} and \code{design} matrices.
#' @param model Regression model describing the dependence of pA site usage on sample attribute(s).
#' @param design (optional) Design matrix. A matrix describing sample attributes which can be used as predictors in the regression model.
#' @param sample_ID (optional) A key variable connecting the counts dataset (\code{data}) and the design matrix.
#' @param adj_method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @details
#' A deviance test compares the likelihood of the fitted model with its corresponding null. In other words, it tests how much the prediction of response
#' is improved by including the covariates, compared to a model with no covariates (the intercept only or null model). In the case of poly A site usage,
#' a logistic regression model is run first (but the outcome is not reported); then, instead of reporting the effects of individual predictors (covariates)
#' on the ratio of specific pairs of pA sites, a deviance test reveals the overall
#' relevance or informativeness of all the predictors in the model towards the pA site usage pattern for each transcript across samples.
#' @examples
#' pA.toy2 <- APAlog::pA.toy2
#' pA_design <- APAlog::pA_design
#' fit3_pA <- pA_logit_dev(pA.toy2, pA.site ~ cell_line, pA_design, "sample")
#' @return
#' Deviance test p-values (one per transcript).
#' @export

pA_logit_dev <- function(data, model, design = NULL, sample_ID = NULL, adj_method){
    data <- remove_0_1_pA_transcripts(data)
    xd <- merge(data, design, by = sample_ID)

    gfitx <- suppressWarnings(by(xd, xd$transcript, function(y) stats::glm(stats::as.formula(Reduce(paste, deparse(model)), env = new.env()), data = y, family = "binomial", weight = count)))

    dtest_gfitx <- lapply(gfitx, function(x) glm_deviance_test_p(x))
    dtest_gfitx.df <- plyr::ldply(dtest_gfitx, rbind, .id = "transcript")
    names(dtest_gfitx.df)[2] <- "p_devtest"

    if (adj_method != 'none'){
        dtest_gfitx.df <- adj_p(dtest_gfitx.df, pcols = 2, adj_method = adj_method)
    }

    return(dtest_gfitx.df)
}



#' @title remove_0_1_pA_transcripts
#' @description Function to remove transcripts with fewer than two active pA sites from dataset
#' @param data Dataset containing poly A (pA) site read counts. This dataset must have a long shape, meaning that there should be only one
#' column containing read counts (and it MUST be named "count"). The first four columns must be called "transcript", "pA.site", "sample" and "count".
#' Thus, each row in \code{data} contains the read count for one pA - transcript - sample combination.
#' Other sample attributes beyond sample ID may be recorded in additional variables in this dataset, or provided separately through a design matrix
#' and a key variable (e.g. sample ID) connecting the \code{data} and \code{design} matrices.
#' @details
#' This function counts the number of pA sites with non-zero read counts for each transcripts and removes transripts
#' with fewer than two active pA sites. This is essential to avoid errors when running the regression models.
#' @return
#' A subset of the input dataset where all transcripts are guaranteed to have two or more active pA sites.
#' @export

remove_0_1_pA_transcripts <- function(data){
    tt <- table(data$transcript, data$pA.site)
    ts <- rownames(tt[rowSums(tt>0) >= 2, , drop = FALSE])
    x2 <- subset(data, transcript %in% ts)
    rate <- 100 * (1 - (length(ts)/dim(tt)[1]))
    print(paste0(rate, "% of transcripts had <2 active pA sites and were removed"))
    return(x2)
}


#' @title volcano_plot
#' @description Function to make the volcano plot from the p-values matrix
#' @param fit Output dataframe of the pA_logit_dev or pA_multi_logit function
#' @param x Name of the column in the dataframe that contains the log fold change
#' @param y Name of the column in the dataframe that contains the corrected p-values
#' @return A volcano plot visualisation
#' @details A visualisation of the volcano plot resulting from the p-values and log fold change
#' @examples
#' vplot <- volcano_plot(fit1, x='logfdr', y='p_values')
#' @export

volcano_plot <- function(fit, x, xlab = "Ln fold change", y, ylab = "-Log10 FDR",
    title = "LMCN data, metastatic vs non-metastatic", titleLabSize = 12, border = "full",
    pCutoff = 0.001, FCcutoff = 1.5, xlim = c(-5, 5), ylim = c(0, 10)) {

    if (! x %in% names(fit)){
        stop(print(paste('The column', x, 'does not exist in the given dataframe.')))
    }

    if (! y %in% names(fit)){
        stop(print(paste('The column', y, 'does not exist in the given dataframe.')))
    }

    return(EnhancedVolcano::EnhancedVolcano(fit, lab = rownames(fit),x=x, xlab=xlab, y=y, ylab=ylab, title=title,
    titleLabSize=titleLabSize, border=border, pCutoff=pCutoff, FCcutoff=FCcutoff, xlim=xlim, ylim=ylim))
}



#' Toy transcript expression data
#'
#' Transcript expression data for testing the functions
#'
#' @name pA.toy2
#' @format An object of class \code{data.frame}
#' @examples
#' data <- Apalog:pA.toy2
"pA.toy2"



#' Toy experimental design
#'
#' experimental design for testing the functions
#'
#' @name pA_design
#' @format An object of class \code{data.frame}
#' @examples
#' data <- Apalog:pA_design
"pA_design"
