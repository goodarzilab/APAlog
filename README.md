![Logo-r](https://github.com/goodarzilab/APAlog/blob/master/data/APALog%20_logo3.png)

# APAlog: A tool for quantification of alternative poly A site usage


Many trancripts in human and other organisms have multiple potential poly A sites. Using different poly A sites affects the composition of 3' UTR and the the regulatory elements it contains, and may impact important aspects of mRNA life cycle including stability and translation rate. The __APAlog__ package tests the significance of differential poly A site usage for each transcript across samples. The source of input data can be a 3' sequencing protocol e.g. Tag-Seq or QuantSeq, or regular RNA-seq reads mapped to annotated poly A sites. __APAlog__ offers three tests to evaluate the differential use of poly A sites of each transcript among samples:

- Overall transcript-wise test
- Multinomial test 
- Pairwise test

The overall test evaluates the null hypothesis of no difference in poly A site usage of a transcript among samples. The multinomial test sets one of the poly A sites to canonical (reference) and compares the usage of all other poly A sites (one or more alternative sites) to this reference. The pairwise test compares the usage of all pairs of poly A sites of a transcript.  

## Installing APAlog

To install __APAlog__ directly from GitHub, you need the *devtools* package. If not already installed on your system, run
    
`install.packages("devtools")`
	
Then, load _devtools_, install and load __APAlog__ by
	
`library(devtools)`  
`install_github("Goodarzilab/APAlog", dependencies = TRUE)`  
`library(APAlog)`  

## Input data

The input to all three tests are two tables: a count table and a design table. Each row of the count table contains normalized RNA read counts pertaining to a poly A site of a transcript in one sample.  The design table describes each sample with covariates that can be used as predictors by __APAlog__. Predictors in the __APAlog__ model can be sample labels or one or more sample attributes (categorical or continuous variables) provided by the design matrix. 

Examples of a count table (only first six rows printed here):

| transcript | pA.site | sample | count |
|-|-|-|-|
| Hs.525527.1 | site.09 | MDA_sgCTRL.r1 | 1.0635544 |
| Hs.525527.1 | site.09 | MDA_sgCTRL.r2 | 0.9875983 |
| Hs.525527.1 | site.09 | MDA_sgHNRNPC.r1 | 2.6289722 |
| Hs.525527.1 | site.09 | MDA_sgHNRNPC.r2 | 4.1650605 |
| Hs.525527.1 | site.16 | MDA_sgCTRL.r1 | 38.2879574 |
| Hs.525527.1 | site.16 | MDA_sgCTRL.r2 | 48.3923165 |

This is a toy dataset with only eight transcripts and four samples (two cell lines, each replicated x2). Seven transcripts have two pA sites, one (Hs.465374.1) has four pA sites. The count table is expected to have at least four columns. The first four columns of this data frame must be named *transcript*, *pA.site*, *sample* and *count*. Additional covariates may be optionally provided. We recommend sample-specific covariates to be provided via the design table and gene- or transcript-specific covariates (optional) as additional columns in the count table. Each row in the count table contains normalized read counts at a pA site of a transcript in oe sample. A common variable named *sample* connects the design matrix and the count matrix. Example of a design matrix:

| sample | cell_line | rep |
|-|-|-|
| MDA_sgCTRL.r1 | MDA_sgCTRL | 1 |
| MDA_sgCTRL.r2 | MDA_sgCTRL | 2 |
| MDA_sgHNRNPC.r1 | MDA_sgHNRNPC | 1 |
| MDA_sgHNRNPC.r2 | MDA_sgHNRNPC | 2 |
  
## Overall transcript-wise test

The aim of this test is to indentify genes or transcripts which show differential poly A site selection among samples without specifying which pA sites are used more or less, or which covariates contribute to the difference. This is achieved through a deviance test which compares goodness-of-fit of the fitted model to the null model to answer the following question: Does the fitted model explain the data significantly better than the null model? (Null: equal usage of poly A sites across all samples.) Check the `pA_logit_dev` function documentation for description of arguments and options. 

`fit.o_HNRNPC <- APAlog::pA_logit_dev(pA.toy2, pA.site ~ cell_line, pA_design, "sample")`

| transcript | p_devtest |
|-|-|
| Hs.29665.1 | 0.6936659 |
| Hs.432760.1 | 0.0920899 |
| Hs.465374.1 | 0.6802614 |
| Hs.469154.1 | 0.5379891 |
| Hs.515329.1 | 0.7983591 |
| Hs.515688.1 | 0.7858670 |
| Hs.523054.1 | 0.3581712 |
| Hs.525527.1 | 0.0582738 |

Correct the p-values for multiple testing:

`fit.o_HNRNPC_fdr <- APAlog::adj_p(fit.o_HNRNPC, pcols = 2, adj_method = "fdr")`

| transcript | p_devtest | fdr_p_devtest |
|-|-|-|
| Hs.29665.1 | 0.6936659 | 0.7983591 |
| Hs.432760.1 | 0.0920899 | 0.3683594 |
| Hs.465374.1 | 0.6802614 | 0.7983591 |
| Hs.469154.1 | 0.5379891 | 0.7983591 |
| Hs.515329.1 | 0.7983591 | 0.7983591 |
| Hs.515688.1 | 0.7858670 | 0.7983591 |
| Hs.523054.1 | 0.3581712 | 0.7983591 |
| Hs.525527.1 | 0.0582738 | 0.3683594 |

Check the `adj_p` function documentation for description of arguments and options.  

## Pairwise test

This test compares all pairs of pA sites of a gene or transcripts and identifies those pairs whose usage ratios varies by the predictors in the model. Check the `pA_logit_pairwise` function documentation for description of arguments and options. 

`fit.p_HNRNPC <- APAlog::pA_logit_pairwise(pA.toy2, pA.site~cell_line, pA_design, "sample")`

| transcript | ref_site | alt_site | b_intercept | p_intercept | b_cell_lineMDA_sgHNRNPC | p_cell_lineMDA_sgHNRNPC |
|-|-|-|-|-|-|-|
| Hs.29665.1 | site.01 | site.04 | -0.9287158 | 0.0337176 | 0.2638369 | 0.6931826 |
| Hs.432760.1 | site.01 | site.04 | 0.1163127 | 0.7253476 | 0.8196840 | 0.0961106 |
| Hs.465374.1 | site.07 | site.14 | 1.9501648 | 0.0000384 | 0.0470609 | 0.9473727 |
| Hs.465374.1 | site.07 | site.15 | 1.1085117 | 0.0301198 | -0.3687804 | 0.6423649 |
| Hs.465374.1 | site.07 | site.16 | 2.1771230 | 0.0000032 | 0.5657150 | 0.4165262 |
| Hs.465374.1 | site.14 | site.15 | -0.8416530 | 0.0057274 | -0.4158414 | 0.4018149 |
| Hs.465374.1 | site.14 | site.16 | 0.2269582 | 0.3111764 | 0.5186541 | 0.1013934 |
| Hs.465374.1 | site.15 | site.16 | 1.0686113 | 0.0002940 | 0.9344955 | 0.0475871 |
| Hs.469154.1 | site.38 | site.39 | 0.2783812 | 0.6023843 | -0.4496357 | 0.5393186 |
| Hs.515329.1 | site.03 | site.05 | 2.6876379 | 0.0000000 | -0.1529467 | 0.7990884 |
| Hs.515688.1 | site.01 | site.05 | -0.1915948 | 0.4830983 | -0.0999251 | 0.7858553 |
| Hs.523054.1 | site.11 | site.17 | 1.0864190 | 0.0548246 | 1.0186876 | 0.3888369 |
| Hs.525527.1 | site.09 | site.16 | 3.7438244 | 0.0000001 | -1.4088383 | 0.0830058 |


Due to the mutual non-independence of p-values from testing pairs of poly A sites in transcripts with more than two sites, a major assumption of multiple testing correction procedures i.e. independence of all p-values is violated. Therefore, adjusting the p-value columns of fit.p_HNRNPC directly may be problematic. To circumvent this problem, we can merge fit.o_HNRNPC and fit.p_HNRNPC. There is only one deviance test per transcript, and the p-values from different transcripts are independent. Adjusted p-values from the deviance test can be used (with a cutoff) to determine which transcripts exhibit differential poly A site usage across the modeled conditions. For those transcripts which pass this threshold, one can select the most important shift(s) in poly A site usage based on the effect size ("b_", regression coefficient), significance (uncorrected "p_") or a combination of them from the output of the pairwise test.

`fit.op_HNRNPC <- merge(fit.o_HNRNPC_fdr, fit.p_HNRNPC, by = "transcript")`

| transcript | p_devtest | fdr_p_devtest | ref_site | alt_site | b_intercept | p_intercept | b_cell_lineMDA_sgHNRNPC | p_cell_lineMDA_sgHNRNPC |
|-|-|-|-|-|-|-|-|-|
| Hs.29665.1 | 0.6936659 | 0.7983591 | site.01 | site.04 | -0.9287158 | 0.0337176 | 0.2638369 | 0.6931826 |
| Hs.432760.1 | 0.0920899 | 0.3683594 | site.01 | site.04 | 0.1163127 | 0.7253476 | 0.8196840 | 0.0961106 |
| Hs.465374.1 | 0.6802614 | 0.7983591 | site.07 | site.14 | 1.9501648 | 0.0000384 | 0.0470609 | 0.9473727 |
| Hs.465374.1 | 0.6802614 | 0.7983591 | site.07 | site.15 | 1.1085117 | 0.0301198 | -0.3687804 | 0.6423649 |
| Hs.465374.1 | 0.6802614 | 0.7983591 | site.07 | site.16 | 2.1771230 | 0.0000032 | 0.5657150 | 0.4165262 |
| Hs.465374.1 | 0.6802614 | 0.7983591 | site.14 | site.15 | -0.8416530 | 0.0057274 | -0.4158414 | 0.4018149 |
| Hs.465374.1 | 0.6802614 | 0.7983591 | site.14 | site.16 | 0.2269582 | 0.3111764 | 0.5186541 | 0.1013934 |
| Hs.465374.1 | 0.6802614 | 0.7983591 | site.15 | site.16 | 1.0686113 | 0.0002940 | 0.9344955 | 0.0475871 |
| Hs.469154.1 | 0.5379891 | 0.7983591 | site.38 | site.39 | 0.2783812 | 0.6023843 | -0.4496357 | 0.5393186 |
| Hs.515329.1 | 0.7983591 | 0.7983591 | site.03 | site.05 | 2.6876379 | 0.0000000 | -0.1529467 | 0.7990884 |
| Hs.515688.1 | 0.7858670 | 0.7983591 | site.01 | site.05 | -0.1915948 | 0.4830983 | -0.0999251 | 0.7858553 |
| Hs.523054.1 | 0.3581712 | 0.7983591 | site.11 | site.17 | 1.0864190 | 0.0548246 | 1.0186876 | 0.3888369 |
| Hs.525527.1 | 0.0582738 | 0.3683594 | site.09 | site.16 | 3.7438244 | 0.0000001 | -1.4088383 | 0.0830058 |
  
## Multinomial test

This test is the preferred choice when one of the pA sites of each transcript e.g. the most proximal site is meant to serve as a baseline or default site vs. the alternative sites that may be activated under certain physiological conditions, as a result of 3' UTR mutations etc. Check the `pA_logit_dev` function documentation for description of arguments and options. 

Note: By deafult, the poly A site that comes first alphabetically is set as the reference level. The user can employ a naming convention using prefixes to mark the reference poly A site for each transcript. 

`fit.m_HNRNPC <- APAlog::pA_multi_logit(pA.toy2, pA.site ~ cell_line, pA_design, "sample")`

| transcript | ref_site | alt_site | b_intercept | b_cell_lineMDA_sgHNRNPC | p_intercept | p_cell_lineMDA_sgHNRNPC |
|-|-|-|-|-|-|-|
| Hs.29665.1 | site.01 | site.04 | -0.9287923 | 0.2639347 | 0.0337059 | 0.6930761 |
| Hs.432760.1 | site.01 | site.04 | 0.1163107 | 0.8196871 | 0.7253521 | 0.0961101 |
| Hs.465374.1 | site.07 | site.14 | 1.9501247 | 0.0471591 | 0.0000384 | 0.9472635 |
| Hs.465374.1 | site.07 | site.15 | 1.1084713 | -0.3686039 | 0.0301238 | 0.6425240 |
| Hs.465374.1 | site.07 | site.16 | 2.1770842 | 0.5658090 | 0.0000032 | 0.4164527 |
| Hs.469154.1 | site.38 | site.39 | 0.2783809 | -0.4496376 | 0.6023881 | 0.5393197 |
| Hs.515329.1 | site.03 | site.05 | 2.6876015 | -0.1529092 | 0.0000000 | 0.7991348 |
| Hs.515688.1 | site.01 | site.05 | -0.1915948 | -0.0999251 | 0.4830983 | 0.7858553 |
| Hs.523054.1 | site.11 | site.17 | 1.0864160 | 1.0186867 | 0.0548251 | 0.3888367 |
| Hs.525527.1 | site.09 | site.16 | 3.7438351 | -1.4088503 | 0.0000001 | 0.0830059 |

Next, adjust the p values for multiple testing:

`fit.m_HNRNPC_fdr <- APAlog::adj_p(fit.m_HNRNPC, pcols = c(6, 7), adj_method = "fdr")`

| transcript | ref_site | alt_site | b_intercept | b_cell_lineMDA_sgHNRNPC | p_intercept | p_cell_lineMDA_sgHNRNPC | fdr_p_intercept | fdr_p_cell_lineMDA_sgHNRNPC |
|-|-|-|-|-|-|-|-|-|
| Hs.29665.1 | site.01 | site.04 | -0.9287923 | 0.2639347 | 0.0337059 | 0.6930761 | 0.0561766 | 0.8879276 |
| Hs.432760.1 | site.01 | site.04 | 0.1163107 | 0.8196871 | 0.7253521 | 0.0961101 | 0.7253521 | 0.4805503 |
| Hs.465374.1 | site.07 | site.14 | 1.9501247 | 0.0471591 | 0.0000384 | 0.9472635 | 0.0000961 | 0.9472635 |
| Hs.465374.1 | site.07 | site.15 | 1.1084713 | -0.3686039 | 0.0301238 | 0.6425240 | 0.0561766 | 0.8879276 |
| Hs.465374.1 | site.07 | site.16 | 2.1770842 | 0.5658090 | 0.0000032 | 0.4164527 | 0.0000108 | 0.8879276 |
| Hs.469154.1 | site.38 | site.39 | 0.2783809 | -0.4496376 | 0.6023881 | 0.5393197 | 0.6693201 | 0.8879276 |
| Hs.515329.1 | site.03 | site.05 | 2.6876015 | -0.1529092 | 0.0000000 | 0.7991348 | 0.0000001 | 0.8879276 |
| Hs.515688.1 | site.01 | site.05 | -0.1915948 | -0.0999251 | 0.4830983 | 0.7858553 | 0.6038729 | 0.8879276 |
| Hs.523054.1 | site.11 | site.17 | 1.0864160 | 1.0186867 | 0.0548251 | 0.3888367 | 0.0783216 | 0.8879276 |
| Hs.525527.1 | site.09 | site.16 | 3.7438351 | -1.4088503 | 0.0000001 | 0.0830059 | 0.0000006 | 0.4805503 |

## Output interpretation

Regression coefficients (b's in the output column names) are equal to log of the alternative/reference poly A site (APA) ratio per covariate. Corresponding adjusted p-values mark the significance of differential poly A usage. For example, for transcript Hs.29665.1:

log APA (site.04/site.01) MDA-sgHNRNPC vs. MDA-sgCTRL = 0.264

APA (site.04/site.01) MDA-sgHNRNPC vs. MDA-sgCTRL = exp(0.264)= 1.302

But this difference is not significant (FDR=0.89). Note that in the case of categorical variables, the name of the non-reference level(s) is appended to variable name in the __APAlog__ output. In the example above, the predictor variable is "cell line", the reference level is "MDA-sgCTRL" and the non-reference level is "MDA-sgHNRNPC".

Note: The multinomial and pairwise tests use different fitting algorithms. That is why the log APA and p-values of similar comparisons in the two outputs are numerically very close but not identical.

Note: Regression intercept in the __APAlog__ output represents the alternative/reference pA site ratio in the reference sample (the sample that has the reference values for all predictors), in this case MDA-sgCTRL. Note that such a sample may not actually exist in some datasets, for example in a multivariate dataset with several categorical predictors and partial factorial design, or a dataset containing continuous variables whose ranges exclude 0. Although the reference or baseline value for continuous covariates like age is automatically set at 0, real samples almost always have non-zero values. Therefore, whether or not the intercept term has a biological meaning depends on the dataset and  sample and variable types. 

Note: __APAlog__ automatically removes transcripts that are not represented by two or more active pA sites (>=2 pA sites with non-zero counts) from the analysis because there is no comparison to make in those  cases. However, if the count of a pA site is zero in the reference sample, it can cause errors of division by zero. A simple fix that is commonly used in bioinformatics is adding a small value e.g. 0.5 to all counts before running the test to avoid this error. 
 

__APAlog__ was developed at UCSF by Hossein Asgharian under supervision of Hani Goodarzi.  

For questions and comments, email us at:  
hosseinali.asgharian@ucsf.edu  
hani.goodarzi@ucsf.edu
