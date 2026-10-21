# ================================================================
# DLBCL: COHORT- AND STROMAL-ADJUSTED CORRELATION ANALYSIS
# STEP 1: LOAD EXPRESSION DATA AND PREPARE SAMPLES
# ================================================================

library(dplyr)
library(tidyr)
library(purrr)

# ------------------------------------------------
# 1. Load the expression matrix containing lncRNAs
#    and mRNAs
# ------------------------------------------------

tpm_all <- read.csv(
  "TPM_DLBCL2_lncRNA_mRNA_boxplot.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check dimensions
cat("Expression data dimensions:", dim(tpm_all), "\n")

# Confirm identifier column
stopifnot("TranscriptSymbol" %in% colnames(tpm_all))

# ------------------------------------------------
# 2. Define sample groups
# ------------------------------------------------

normal_samples <- paste0("C", 1:8)
sra_samples    <- paste0("T", 1:10)
tcga_samples   <- paste0("T", 11:58)

tumour_samples <- c(sra_samples, tcga_samples)

# Verify sample columns
required_samples <- c(normal_samples, tumour_samples)

missing_samples <- setdiff(
  required_samples,
  colnames(tpm_all)
)

if (length(missing_samples) > 0) {
  stop(
    "Missing sample columns: ",
    paste(missing_samples, collapse = ", ")
  )
}

# ------------------------------------------------
# 3. Define the candidate transcripts
# ------------------------------------------------

lncRNAs <- c(
  "lnc-ADAMTS1-2:1",
  "lnc-FXYD2-1:5",
  "lnc-SPARC-2:3"
)

mRNAs <- c(
  "ADAMTS1",
  "ADAMTS5",
  "ATOX1",
  "BCL6",
  "FXYD2",
  "FXYD6",
  "FXYD6-FXYD2",
  "INO80",
  "SPARC"
)

# Generate 27 candidate combinations
candidate_pairs <- tidyr::expand_grid(
  lncRNA = lncRNAs,
  mRNA = mRNAs
)

cat("Number of candidate pairs:", nrow(candidate_pairs), "\n")

print(candidate_pairs)

# ------------------------------------------------
# 4. Extract expression values
# ------------------------------------------------

# Convert sample columns to numeric
for (sample in required_samples) {
  tpm_all[[sample]] <- as.numeric(
    as.character(tpm_all[[sample]])
  )
}

# Check for duplicated transcript identifiers
duplicate_transcripts <- unique(
  tpm_all$TranscriptSymbol[
    duplicated(tpm_all$TranscriptSymbol)
  ]
)

if (length(duplicate_transcripts) > 0) {
  cat("\nDuplicated identifiers in the full matrix:\n")
  print(duplicate_transcripts)
}

# Do not silently average duplicated identifiers.
# Stop if any candidate transcript has multiple rows.
candidate_ids <- unique(c(lncRNAs, mRNAs))

candidate_counts <- table(
  tpm_all$TranscriptSymbol[
    tpm_all$TranscriptSymbol %in% candidate_ids
  ]
)

if (any(candidate_counts > 1)) {
  stop(
    "At least one candidate transcript has duplicate rows. ",
    "Resolve the candidate identifiers before analysis."
  )
}

# Confirm that every candidate is present
missing_candidates <- setdiff(
  candidate_ids,
  tpm_all$TranscriptSymbol
)

if (length(missing_candidates) > 0) {
  stop(
    "These candidate transcripts were not found: ",
    paste(missing_candidates, collapse = ", ")
  )
}

# Extract candidate expression
candidate_expression <- tpm_all[
  match(candidate_ids, tpm_all$TranscriptSymbol),
  c("TranscriptSymbol", tumour_samples),
  drop = FALSE
]

# Create a named matrix
candidate_matrix <- as.matrix(
  candidate_expression[, tumour_samples, drop = FALSE]
)

storage.mode(candidate_matrix) <- "numeric"

rownames(candidate_matrix) <- candidate_expression$TranscriptSymbol

# Verify
cat("\nCandidate expression matrix dimensions:\n")
print(dim(candidate_matrix))

print(candidate_matrix[, 1:5, drop = FALSE])

# ------------------------------------------------
# 5. Prepare tumour metadata
# ------------------------------------------------
stromal_scores <- read.csv(
  "DLBCL_stromal_scores_Results.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

tumour_metadata <- stromal_scores %>%
  filter(Sample %in% tumour_samples) %>%
  select(
    Sample,
    Cohort,
    StromalScore
  )

# Arrange metadata in the same order as the expression matrix
tumour_metadata <- tumour_metadata[
  match(tumour_samples, tumour_metadata$Sample),
]

# Verify alignment
stopifnot(
  identical(tumour_metadata$Sample, tumour_samples)
)

stopifnot(
  !anyNA(tumour_metadata$Cohort),
  !anyNA(tumour_metadata$StromalScore)
)

cat("\nTumour sample counts by cohort:\n")
print(table(tumour_metadata$Cohort))

cat("\nMetadata prepared successfully.\n")

################################################################################
# =====================================================================
# DLBCL: COHORT- AND STROMAL-ADJUSTED CORRELATION ANALYSIS
# STEP 2: STATISTICAL ANALYSIS OF 27 lncRNA-mRNA PAIRS
# =====================================================================
################################################################################
library(dplyr)
library(tidyr)
library(purrr)

set.seed(2026)

# ---------------------------------------------------------------------
# 1. PREPARE EXPRESSION AND METADATA
# ---------------------------------------------------------------------

# Confirm the expected sample numbers
stopifnot(ncol(candidate_matrix) == 58)
stopifnot(nrow(tumour_metadata) == 58)

stopifnot(sum(tumour_metadata$Cohort == "SRA") == 10)
stopifnot(sum(tumour_metadata$Cohort == "TCGA") == 48)

# Ensure metadata and expression columns have identical ordering
stopifnot(
  identical(colnames(candidate_matrix), tumour_metadata$Sample)
)

# ---------------------------------------------------------------------
# 2. DEFINE SPEARMAN CORRELATION FUNCTION
# ---------------------------------------------------------------------

calculate_spearman <- function(x, y) {
  
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  n <- length(x)
  
  if (n < 4 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(
      data.frame(
        n = n,
        rho = NA_real_,
        p_value = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_
      )
    )
  }
  
  test <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    n = n,
    rho = unname(test$estimate),
    p_value = test$p.value,
    CI_lower = NA_real_,
    CI_upper = NA_real_
  )
}

# ---------------------------------------------------------------------
# 3. ANALYSIS 1: SRA TUMOUR CORRELATIONS
# ---------------------------------------------------------------------

sra_samples <- tumour_metadata$Sample[
  tumour_metadata$Cohort == "SRA"
]

sra_results <- purrr::map_dfr(
  seq_len(nrow(candidate_pairs)),
  function(i) {
    
    lnc <- candidate_pairs$lncRNA[i]
    mrna <- candidate_pairs$mRNA[i]
    
    result <- calculate_spearman(
      candidate_matrix[lnc, sra_samples],
      candidate_matrix[mrna, sra_samples]
    )
    
    data.frame(
      lncRNA = lnc,
      mRNA = mrna,
      result
    )
  }
)

sra_results$FDR_BH <- p.adjust(
  sra_results$p_value,
  method = "BH"
)

sra_results <- sra_results %>%
  arrange(FDR_BH)

cat("\nSRA tumour correlations:\n")
print(sra_results)

# ---------------------------------------------------------------------
# 4. ANALYSIS 2: TCGA TUMOUR CORRELATIONS
# ---------------------------------------------------------------------

tcga_samples <- tumour_metadata$Sample[
  tumour_metadata$Cohort == "TCGA"
]

tcga_results <- purrr::map_dfr(
  seq_len(nrow(candidate_pairs)),
  function(i) {
    
    lnc <- candidate_pairs$lncRNA[i]
    mrna <- candidate_pairs$mRNA[i]
    
    result <- calculate_spearman(
      candidate_matrix[lnc, tcga_samples],
      candidate_matrix[mrna, tcga_samples]
    )
    
    data.frame(
      lncRNA = lnc,
      mRNA = mrna,
      result
    )
  }
)

tcga_results$FDR_BH <- p.adjust(
  tcga_results$p_value,
  method = "BH"
)

tcga_results <- tcga_results %>%
  arrange(FDR_BH)

cat("\nTCGA tumour correlations:\n")
print(tcga_results)

# ---------------------------------------------------------------------
# 5. ANALYSIS 3: ADJUSTED ASSOCIATIONS
#
# Rank-based linear regression:
# rank(mRNA) ~ rank(lncRNA) + Cohort + rank(StromalScore)
#
# The rank_lncRNA coefficient estimates the association between lncRNA
# and mRNA expression after adjusting for cohort and stromal score.
# ---------------------------------------------------------------------

# Reset R's missing-value printing option
options(na.print = "NA")

# Run the adjusted regression for all candidate pairs
adjusted_results <- purrr::map_dfr(
  seq_len(nrow(candidate_pairs)),
  function(i) {
    
    # Identify the current transcript pair
    lnc <- candidate_pairs$lncRNA[i]
    mrna <- candidate_pairs$mRNA[i]
    
    # Prepare the data
    dat <- data.frame(
      lncRNA = as.numeric(
        candidate_matrix[lnc, tumour_samples]
      ),
      mRNA = as.numeric(
        candidate_matrix[mrna, tumour_samples]
      ),
      Cohort = factor(
        tumour_metadata$Cohort
      ),
      StromalScore = as.numeric(
        tumour_metadata$StromalScore
      )
    )
    
    # Remove missing and non-finite observations
    dat <- dat[
      complete.cases(dat) &
        is.finite(dat$lncRNA) &
        is.finite(dat$mRNA) &
        is.finite(dat$StromalScore),
      ,
      drop = FALSE
    ]
    
    # Return NA if there are insufficient observations or variation
    if (
      nrow(dat) < 10 ||
      length(unique(dat$lncRNA)) < 2 ||
      length(unique(dat$mRNA)) < 2 ||
      length(unique(dat$StromalScore)) < 2 ||
      length(unique(dat$Cohort)) < 2
    ) {
      
      return(
        data.frame(
          lncRNA = lnc,
          mRNA = mrna,
          n = nrow(dat),
          beta = NA_real_,
          SE = NA_real_,
          t_value = NA_real_,
          p_value = NA_real_
        )
      )
    }
    
    # Rank-transform the variables
    dat$rank_lncRNA <- rank(
      dat$lncRNA,
      ties.method = "average"
    )
    
    dat$rank_mRNA <- rank(
      dat$mRNA,
      ties.method = "average"
    )
    
    dat$rank_Stromal <- rank(
      dat$StromalScore,
      ties.method = "average"
    )
    
    # Fit the adjusted regression model
    model <- lm(
      rank_mRNA ~ rank_lncRNA + Cohort + rank_Stromal,
      data = dat
    )
    
    # Extract regression statistics
    coefficients <- summary(model)$coefficients
    
    # Check whether the lncRNA coefficient is available
    if (!"rank_lncRNA" %in% rownames(coefficients)) {
      
      return(
        data.frame(
          lncRNA = lnc,
          mRNA = mrna,
          n = nrow(dat),
          beta = NA_real_,
          SE = NA_real_,
          t_value = NA_real_,
          p_value = NA_real_
        )
      )
    }
    
    # Return the statistics for the lncRNA term
    data.frame(
      lncRNA = lnc,
      mRNA = mrna,
      n = nrow(dat),
      beta = coefficients[
        "rank_lncRNA", "Estimate"
      ],
      SE = coefficients[
        "rank_lncRNA", "Std. Error"
      ],
      t_value = coefficients[
        "rank_lncRNA", "t value"
      ],
      p_value = coefficients[
        "rank_lncRNA", "Pr(>|t|)"
      ]
    )
  }
)

# ---------------------------------------------------------------------
# 5.1. MULTIPLE-TESTING CORRECTION
# ---------------------------------------------------------------------

adjusted_results <- adjusted_results %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(FDR)

# ---------------------------------------------------------------------
# 5.2. DISPLAY RESULTS
# ---------------------------------------------------------------------

print(
  as.data.frame(adjusted_results),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 5.3. EXPORT ALL RESULTS
# ---------------------------------------------------------------------

write.csv(
  adjusted_results,
  "DLBCL_Stromal_Cohort_Adjusted_Results.csv",
  row.names = FALSE,
  na = "NA"
)

# ---------------------------------------------------------------------
# 5.4. DISPLAY SIGNIFICANT ASSOCIATIONS
# ---------------------------------------------------------------------

significant_adjusted <- adjusted_results %>%
  filter(
    !is.na(FDR),
    FDR < 0.05
  )

cat(
  "\nNumber of significant adjusted associations (FDR < 0.05):",
  nrow(significant_adjusted),
  "\n"
)

print(
  as.data.frame(significant_adjusted),
  row.names = FALSE
)

# Export significant results separately
write.csv(
  significant_adjusted,
  "DLBCL_Stromal_Cohort_Adjusted_Significant_Results.csv",
  row.names = FALSE,
  na = "NA"
)

###########################
# Check the number of tested pairs
nrow(adjusted_results)

# View all adjusted associations
print(as.data.frame(adjusted_results), row.names = FALSE)

# Count significant associations
sum(adjusted_results$FDR < 0.05, na.rm = TRUE)

# View significant associations
adjusted_results %>%
  filter(!is.na(FDR), FDR < 0.05)
# ---------------------------------------------------------------------
# END OF SECTION 5
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# 6. ANALYSIS 4: COHORT-INTERACTION ANALYSIS
#
# Model:
# rank(mRNA) ~ rank(lncRNA) * Cohort + rank(StromalScore)
#
# Tests whether the lncRNA-mRNA association differs between
# SRA and TCGA after adjusting for stromal score.
# ---------------------------------------------------------------------

library(dplyr)
library(purrr)

# Reset printing options
options(na.print = "NA")

# Set SRA as the reference cohort
tumour_metadata$Cohort <- factor(
  tumour_metadata$Cohort,
  levels = c("SRA", "TCGA")
)

# Confirm cohort sample counts
print(table(tumour_metadata$Cohort))

# Confirm sample alignment
stopifnot(
  identical(colnames(candidate_matrix), tumour_metadata$Sample)
)

# ---------------------------------------------------------------------
# 6.1. FIT INTERACTION MODELS FOR ALL 27 PAIRS
# ---------------------------------------------------------------------

interaction_results <- purrr::map_dfr(
  seq_len(nrow(candidate_pairs)),
  function(i) {
    
    lnc <- candidate_pairs$lncRNA[i]
    mrna <- candidate_pairs$mRNA[i]
    
    # Prepare pair-specific data
    dat <- data.frame(
      lncRNA = as.numeric(
        candidate_matrix[lnc, tumour_samples]
      ),
      mRNA = as.numeric(
        candidate_matrix[mrna, tumour_samples]
      ),
      Cohort = factor(
        tumour_metadata$Cohort,
        levels = c("SRA", "TCGA")
      ),
      StromalScore = as.numeric(
        tumour_metadata$StromalScore
      )
    )
    
    # Remove missing and non-finite observations
    dat <- dat[
      complete.cases(dat) &
        is.finite(dat$lncRNA) &
        is.finite(dat$mRNA) &
        is.finite(dat$StromalScore),
      ,
      drop = FALSE
    ]
    
    # Check data sufficiency
    if (
      nrow(dat) < 10 ||
      length(unique(dat$lncRNA)) < 2 ||
      length(unique(dat$mRNA)) < 2 ||
      length(unique(dat$StromalScore)) < 2 ||
      length(unique(dat$Cohort)) < 2
    ) {
      return(
        data.frame(
          lncRNA = lnc,
          mRNA = mrna,
          n = nrow(dat),
          n_SRA = sum(dat$Cohort == "SRA"),
          n_TCGA = sum(dat$Cohort == "TCGA"),
          Slope_SRA = NA_real_,
          Slope_TCGA = NA_real_,
          Interaction_Beta = NA_real_,
          Interaction_SE = NA_real_,
          Interaction_t = NA_real_,
          Interaction_P = NA_real_
        )
      )
    }
    
    # Rank-transform expression and stromal score
    dat$rank_lncRNA <- rank(
      dat$lncRNA,
      ties.method = "average"
    )
    
    dat$rank_mRNA <- rank(
      dat$mRNA,
      ties.method = "average"
    )
    
    dat$rank_Stromal <- rank(
      dat$StromalScore,
      ties.method = "average"
    )
    
    # Fit the cohort-interaction model
    model <- tryCatch(
      lm(
        rank_mRNA ~
          rank_lncRNA * Cohort +
          rank_Stromal,
        data = dat
      ),
      error = function(e) NULL
    )
    
    # Handle model-fitting errors
    if (is.null(model)) {
      return(
        data.frame(
          lncRNA = lnc,
          mRNA = mrna,
          n = nrow(dat),
          n_SRA = sum(dat$Cohort == "SRA"),
          n_TCGA = sum(dat$Cohort == "TCGA"),
          Slope_SRA = NA_real_,
          Slope_TCGA = NA_real_,
          Interaction_Beta = NA_real_,
          Interaction_SE = NA_real_,
          Interaction_t = NA_real_,
          Interaction_P = NA_real_
        )
      )
    }
    
    # Extract model coefficients
    coefficients <- summary(model)$coefficients
    
    # Identify the interaction term
    interaction_term <- grep(
      "^rank_lncRNA:CohortTCGA$|^CohortTCGA:rank_lncRNA$",
      rownames(coefficients),
      value = TRUE
    )
    
    # Check that the interaction coefficient exists
    if (length(interaction_term) != 1) {
      return(
        data.frame(
          lncRNA = lnc,
          mRNA = mrna,
          n = nrow(dat),
          n_SRA = sum(dat$Cohort == "SRA"),
          n_TCGA = sum(dat$Cohort == "TCGA"),
          Slope_SRA = NA_real_,
          Slope_TCGA = NA_real_,
          Interaction_Beta = NA_real_,
          Interaction_SE = NA_real_,
          Interaction_t = NA_real_,
          Interaction_P = NA_real_
        )
      )
    }
    
    # Extract the interaction statistics
    interaction_beta <- coefficients[
      interaction_term, "Estimate"
    ]
    
    interaction_se <- coefficients[
      interaction_term, "Std. Error"
    ]
    
    interaction_t <- coefficients[
      interaction_term, "t value"
    ]
    
    interaction_p <- coefficients[
      interaction_term, "Pr(>|t|)"
    ]
    
    # Extract cohort-specific lncRNA slopes
    slope_sra <- if ("rank_lncRNA" %in% rownames(coefficients)) {
      coefficients["rank_lncRNA", "Estimate"]
    } else {
      NA_real_
    }
    
    slope_tcga <- slope_sra + interaction_beta
    
    # Return results
    data.frame(
      lncRNA = lnc,
      mRNA = mrna,
      n = nrow(dat),
      n_SRA = sum(dat$Cohort == "SRA"),
      n_TCGA = sum(dat$Cohort == "TCGA"),
      Slope_SRA = slope_sra,
      Slope_TCGA = slope_tcga,
      Interaction_Beta = interaction_beta,
      Interaction_SE = interaction_se,
      Interaction_t = interaction_t,
      Interaction_P = interaction_p
    )
  }
)

# ---------------------------------------------------------------------
# 6.2. MULTIPLE-TESTING CORRECTION
# ---------------------------------------------------------------------

interaction_results <- interaction_results %>%
  mutate(
    Interaction_FDR = p.adjust(
      Interaction_P,
      method = "BH"
    )
  ) %>%
  arrange(Interaction_FDR)

# ---------------------------------------------------------------------
# 6.3. DISPLAY ALL RESULTS
# ---------------------------------------------------------------------

print(
  as.data.frame(interaction_results),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 6.4. EXPORT ALL RESULTS
# ---------------------------------------------------------------------

write.csv(
  interaction_results,
  "DLBCL_Cohort_Interaction_Results.csv",
  row.names = FALSE,
  na = "NA"
)

# ---------------------------------------------------------------------
# 6.5. EXTRACT SIGNIFICANT INTERACTIONS
# ---------------------------------------------------------------------

significant_interactions <- interaction_results %>%
  filter(
    !is.na(Interaction_FDR),
    Interaction_FDR < 0.05
  )

cat(
  "\nSignificant cohort interactions (FDR < 0.05):",
  nrow(significant_interactions),
  "\n"
)

print(
  as.data.frame(significant_interactions),
  row.names = FALSE
)

# Export significant interactions separately
write.csv(
  significant_interactions,
  "DLBCL_Cohort_Interaction_Significant_Results.csv",
  row.names = FALSE,
  na = "NA"
)

print(as.data.frame(interaction_results), row.names = FALSE)
# ---------------------------------------------------------------------
# END OF SECTION 6
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# 7. INTEGRATE ALL CORRELATION AND REGRESSION RESULTS
#
# Analyses:
#   1. SRA Spearman correlation
#   2. TCGA Spearman correlation
#   3. Cohort- and stromal-score-adjusted regression
#   4. SRA-versus-TCGA cohort-interaction analysis
#
# Output: 27 candidate lncRNA-mRNA pairs
# ---------------------------------------------------------------------

library(dplyr)

# ---------------------------------------------------------------------
# 7.1. VERIFY THE INPUT TABLES
# ---------------------------------------------------------------------

stopifnot(
  nrow(sra_results) == 27,
  nrow(tcga_results) == 27,
  nrow(adjusted_results) == 27,
  nrow(interaction_results) == 27
)

# Check for duplicate lncRNA-mRNA pairs
stopifnot(
  !anyDuplicated(sra_results[c("lncRNA", "mRNA")]),
  !anyDuplicated(tcga_results[c("lncRNA", "mRNA")]),
  !anyDuplicated(adjusted_results[c("lncRNA", "mRNA")]),
  !anyDuplicated(interaction_results[c("lncRNA", "mRNA")])
)

# Confirm that all four tables contain the same pairs
pair_key <- function(x) {
  paste(x$lncRNA, x$mRNA, sep = "|||")
}

reference_pairs <- sort(pair_key(sra_results))

stopifnot(
  identical(reference_pairs, sort(pair_key(tcga_results))),
  identical(reference_pairs, sort(pair_key(adjusted_results))),
  identical(reference_pairs, sort(pair_key(interaction_results)))
)

# ---------------------------------------------------------------------
# 7.2. PREPARE SRA SPEARMAN RESULTS
# ---------------------------------------------------------------------

sra_integrated <- sra_results %>%
  select(
    lncRNA,
    mRNA,
    SRA_n = n,
    SRA_rho = rho,
    SRA_p_value = p_value,
    SRA_CI_lower = CI_lower,
    SRA_CI_upper = CI_upper,
    SRA_FDR = FDR_BH
  )

# ---------------------------------------------------------------------
# 7.3. PREPARE TCGA SPEARMAN RESULTS
# ---------------------------------------------------------------------

tcga_integrated <- tcga_results %>%
  select(
    lncRNA,
    mRNA,
    TCGA_n = n,
    TCGA_rho = rho,
    TCGA_p_value = p_value,
    TCGA_CI_lower = CI_lower,
    TCGA_CI_upper = CI_upper,
    TCGA_FDR = FDR_BH
  )

# ---------------------------------------------------------------------
# 7.4. PREPARE ADJUSTED REGRESSION RESULTS
# ---------------------------------------------------------------------

adjusted_integrated <- adjusted_results %>%
  select(
    lncRNA,
    mRNA,
    Adjusted_n = n,
    Adjusted_beta = beta,
    Adjusted_SE = SE,
    Adjusted_t = t_value,
    Adjusted_p_value = p_value,
    Adjusted_FDR = FDR
  )

# ---------------------------------------------------------------------
# 7.5. PREPARE COHORT-INTERACTION RESULTS
# ---------------------------------------------------------------------

interaction_integrated <- interaction_results %>%
  select(
    lncRNA,
    mRNA,
    Interaction_n = n,
    Interaction_n_SRA = n_SRA,
    Interaction_n_TCGA = n_TCGA,
    Slope_SRA,
    Slope_TCGA,
    Interaction_beta = Interaction_Beta,
    Interaction_SE,
    Interaction_t,
    Interaction_p_value = Interaction_P,
    Interaction_FDR
  )

# ---------------------------------------------------------------------
# 7.6. MERGE THE FOUR ANALYSES
# ---------------------------------------------------------------------

integrated_results <- sra_integrated %>%
  full_join(
    tcga_integrated,
    by = c("lncRNA", "mRNA")
  ) %>%
  full_join(
    adjusted_integrated,
    by = c("lncRNA", "mRNA")
  ) %>%
  full_join(
    interaction_integrated,
    by = c("lncRNA", "mRNA")
  )

# ---------------------------------------------------------------------
# 7.7. ADD SIGNIFICANCE FLAGS
# ---------------------------------------------------------------------

integrated_results <- integrated_results %>%
  mutate(
    SRA_significant = !is.na(SRA_FDR) & SRA_FDR < 0.05,
    TCGA_significant = !is.na(TCGA_FDR) & TCGA_FDR < 0.05,
    Adjusted_significant = !is.na(Adjusted_FDR) &
      Adjusted_FDR < 0.05,
    Interaction_significant = !is.na(Interaction_FDR) &
      Interaction_FDR < 0.05
  ) %>%
  arrange(
    Interaction_FDR,
    Adjusted_FDR
  )

# ---------------------------------------------------------------------
# 7.8. VERIFY THE FINAL TABLE
# ---------------------------------------------------------------------

stopifnot(
  nrow(integrated_results) == 27,
  !anyDuplicated(integrated_results[c("lncRNA", "mRNA")])
)

cat("\nNumber of integrated candidate pairs:",
    nrow(integrated_results), "\n")

cat("\nSignificant SRA correlations (FDR < 0.05):",
    sum(integrated_results$SRA_significant), "\n")

cat("Significant TCGA correlations (FDR < 0.05):",
    sum(integrated_results$TCGA_significant), "\n")

cat("Significant adjusted associations (FDR < 0.05):",
    sum(integrated_results$Adjusted_significant), "\n")

cat("Significant cohort interactions (FDR < 0.05):",
    sum(integrated_results$Interaction_significant), "\n")

# Display the complete integrated table
options(na.print = "NA")

print(
  as.data.frame(integrated_results),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 7.9. EXPORT THE FINAL TABLE
# ---------------------------------------------------------------------

write.csv(
  integrated_results,
  "DLBCL_Integrated_Correlation_Results.csv",
  row.names = FALSE,
  na = "NA"
)

cat(
  "\nIntegrated results saved to:",
  file.path(
    getwd(),
    "DLBCL_Integrated_Correlation_Results.csv"
  ),
  "\n"
)

# ---------------------------------------------------------------------
# END OF SECTION 7
# ---------------------------------------------------------------------

write.csv(
  consistent_associations,
  "DLBCL_Consistent_TCGA_Adjusted_Associations.csv",
  row.names = FALSE,
  na = "NA"
)

# Identify the 10 prioritised pairs in the SRA results
SRA_prioritised <- sra_results %>%
  inner_join(
    consistent_associations %>%
      select(lncRNA, mRNA),
    by = c("lncRNA", "mRNA")
  ) %>%
  arrange(FDR_BH)

# Display SRA-specific results
print(
  as.data.frame(SRA_prioritised),
  row.names = FALSE
)

# Count pairs significant in SRA after FDR correction
cat(
  "\nPrioritised pairs significant in SRA (FDR < 0.05):",
  sum(SRA_prioritised$FDR_BH < 0.05, na.rm = TRUE),
  "of",
  nrow(SRA_prioritised),
  "\n"
)

# Export the results
write.csv(
  SRA_prioritised,
  "DLBCL_Prioritised_Pairs_SRA_Results.csv",
  row.names = FALSE,
  na = "NA"
)

# ---------------------------------------------------------------------
# 8. SCATTERPLOTS FOR THE THREE STRONGEST CROSS-COHORT ASSOCIATIONS
# ---------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

top_pairs <- data.frame(
  lncRNA = c(
    "lnc-ADAMTS1-2:1",
    "lnc-FXYD2-1:5",
    "lnc-SPARC-2:3"
  ),
  mRNA = c(
    "ADAMTS5",
    "FXYD6",
    "ATOX1"
  )
)

# Prepare sample metadata
plot_metadata <- tumour_metadata %>%
  select(Sample, Cohort)

# Extract the three lncRNAs and three mRNAs
plot_transcripts <- unique(
  c(top_pairs$lncRNA, top_pairs$mRNA)
)

plot_data <- as.data.frame(
  t(candidate_matrix[plot_transcripts, tumour_samples, drop = FALSE])
)

plot_data$Sample <- rownames(plot_data)

plot_data <- plot_data %>%
  left_join(
    plot_metadata,
    by = "Sample"
  )

# Create one scatterplot per pair and cohort
for (i in seq_len(nrow(top_pairs))) {
  
  lnc <- top_pairs$lncRNA[i]
  mrna <- top_pairs$mRNA[i]
  
  for (cohort_name in c("SRA", "TCGA")) {
    
    dat <- plot_data %>%
      filter(Cohort == cohort_name)
    
    p <- ggplot(
      dat,
      aes(
        x = .data[[lnc]],
        y = .data[[mrna]]
      )
    ) +
      geom_point(size = 3) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE
      ) +
      labs(
        title = paste(lnc, "vs", mrna),
        subtitle = paste(
          cohort_name,
          "tumours; n =",
          nrow(dat)
        ),
        x = paste(lnc, "expression (TPM)"),
        y = paste(mrna, "expression (TPM)")
      ) +
      theme_classic(base_size = 12)
    
    print(p)
  }
}

# =====================================================================
# STEP 9: FACETED CORRELATION PLOTS FOR ALL 27 lncRNA-mRNA PAIRS
# =====================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------
# 9.1. Define sample groups
# ---------------------------------------------------------------------

normal_samples <- paste0("C", 1:8)
sra_samples    <- paste0("T", 1:10)
tcga_samples   <- paste0("T", 11:58)

# ---------------------------------------------------------------------
# 9.2. Extract normal expression data
# ---------------------------------------------------------------------

normal_expression <- tpm_all[
  match(
    unique(c(lncRNAs, mRNAs)),
    tpm_all$TranscriptSymbol
  ),
  c("TranscriptSymbol", normal_samples),
  drop = FALSE
]

normal_matrix <- as.matrix(
  normal_expression[, normal_samples, drop = FALSE]
)

storage.mode(normal_matrix) <- "numeric"

rownames(normal_matrix) <- normal_expression$TranscriptSymbol

# ---------------------------------------------------------------------
# 9.3. Combine normal and tumour expression matrices
# ---------------------------------------------------------------------

all_samples <- c(
  normal_samples,
  sra_samples,
  tcga_samples
)

all_expression_matrix <- cbind(
  normal_matrix,
  candidate_matrix
)

# Confirm alignment
stopifnot(
  identical(colnames(all_expression_matrix), all_samples)
)

# ---------------------------------------------------------------------
# 9.4. Prepare plotting data
# ---------------------------------------------------------------------

plot_data <- as.data.frame(
  t(all_expression_matrix)
)

plot_data$Sample <- rownames(plot_data)

plot_data <- plot_data %>%
  mutate(
    SampleType = case_when(
      Sample %in% normal_samples ~ "Normal",
      Sample %in% c(sra_samples, tcga_samples) ~ "Tumour",
      TRUE ~ NA_character_
    ),
    Cohort = case_when(
      Sample %in% c(normal_samples, sra_samples) ~ "SRA",
      Sample %in% tcga_samples ~ "TCGA",
      TRUE ~ NA_character_
    ),
    Group = case_when(
      Sample %in% normal_samples ~ "Normal",
      Sample %in% sra_samples ~ "SRA Tumour",
      Sample %in% tcga_samples ~ "TCGA Tumour",
      TRUE ~ NA_character_
    )
  )

# Transform expression values for plotting
# Assumes the original data contain non-negative TPM values.

plot_data <- plot_data %>%
  mutate(
    across(
      all_of(unique(c(lncRNAs, mRNAs))),
      ~ log10(.x + 1)
    )
  )

# ---------------------------------------------------------------------
# 9.5. Calculate Spearman correlations for each pair and group
# ---------------------------------------------------------------------

calculate_plot_correlations <- function(data, group_column) {
  
  results <- purrr::map_dfr(
    seq_len(nrow(candidate_pairs)),
    function(i) {
      
      lnc <- candidate_pairs$lncRNA[i]
      mrna <- candidate_pairs$mRNA[i]
      
      purrr::map_dfr(
        unique(data[[group_column]]),
        function(group_name) {
          
          dat <- data %>%
            filter(.data[[group_column]] == group_name)
          
          x <- dat[[lnc]]
          y <- dat[[mrna]]
          
          valid <- is.finite(x) & is.finite(y)
          
          x <- x[valid]
          y <- y[valid]
          
          n <- length(x)
          
          if (
            n < 3 ||
            length(unique(x)) < 2 ||
            length(unique(y)) < 2
          ) {
            return(
              data.frame(
                lncRNA = lnc,
                mRNA = mrna,
                Group = as.character(group_name),
                n = n,
                rho = NA_real_,
                p_value = NA_real_
              )
            )
          }
          
          test <- suppressWarnings(
            cor.test(
              x,
              y,
              method = "spearman",
              exact = FALSE
            )
          )
          
          data.frame(
            lncRNA = lnc,
            mRNA = mrna,
            Group = as.character(group_name),
            n = n,
            rho = unname(test$estimate),
            p_value = test$p.value
          )
        }
      )
    }
  )
  
  results %>%
    group_by(Group) %>%
    mutate(
      FDR_BH = p.adjust(p_value, method = "BH")
    ) %>%
    ungroup()
}

# ---------------------------------------------------------------------
# 9.6. Function to generate a 27-panel correlation figure
# ---------------------------------------------------------------------

make_27_pair_plot <- function(
    data,
    groups_to_include,
    group_labels,
    plot_title,
    output_file
) {
  
  dat <- data %>%
    filter(Group %in% groups_to_include) %>%
    mutate(
      Group = factor(
        Group,
        levels = groups_to_include,
        labels = group_labels
      )
    )
  
  # Long-format data: one row per sample-pair
  pair_data <- purrr::map_dfr(
    seq_len(nrow(candidate_pairs)),
    function(i) {
      
      lnc <- candidate_pairs$lncRNA[i]
      mrna <- candidate_pairs$mRNA[i]
      
      data.frame(
        Sample = dat$Sample,
        Group = dat$Group,
        lncRNA = lnc,
        mRNA = mrna,
        Pair = paste(lnc, mrna, sep = " | "),
        x = dat[[lnc]],
        y = dat[[mrna]]
      )
    }
  )
  
  # Calculate statistics for the groups shown in this figure
  stats <- calculate_plot_correlations(
    dat,
    "Group"
  )
  
  stats <- stats %>%
    mutate(
      Pair = paste(lncRNA, mRNA, sep = " | "),
      Label = ifelse(
        is.na(rho),
        paste0(Group, "\nn = ", n, "; rho = NA"),
        paste0(
          Group,
          "\nn = ", n,
          "; rho = ", sprintf("%.2f", rho)
        )
      )
    )
  
  # Place each group's annotation near the upper-left of its panel
  positions <- pair_data %>%
    group_by(Pair) %>%
    summarise(
      x_pos = min(x, na.rm = TRUE),
      y_pos = max(y, na.rm = TRUE),
      y_range = diff(range(y, na.rm = TRUE)),
      .groups = "drop"
    )
  
  stats <- stats %>%
    left_join(positions, by = "Pair") %>%
    group_by(Pair) %>%
    mutate(
      label_y = y_pos -
        (row_number() - 1) *
        pmax(y_range, 0.1) * 0.15
    ) %>%
    ungroup()
  
  # Preserve the intended order of the 27 pairs
  pair_levels <- paste(
    candidate_pairs$lncRNA,
    candidate_pairs$mRNA,
    sep = " | "
  )
  
  pair_data$Pair <- factor(
    pair_data$Pair,
    levels = pair_levels
  )
  
  stats$Pair <- factor(
    stats$Pair,
    levels = pair_levels
  )
  
  # Create separate facet variables
  pair_data <- pair_data %>%
    mutate(
      lncRNA = factor(
        lncRNA,
        levels = unique(candidate_pairs$lncRNA)
      ),
      mRNA = factor(
        mRNA,
        levels = unique(candidate_pairs$mRNA)
      )
    )
  
  stats <- stats %>%
    mutate(
      lncRNA = factor(
        lncRNA,
        levels = unique(candidate_pairs$lncRNA)
      ),
      mRNA = factor(
        mRNA,
        levels = unique(candidate_pairs$mRNA)
      )
    )
  
  p <- ggplot(
    pair_data,
    aes(
      x = x,
      y = y,
      colour = Group,
      shape = Group
    )
  ) +
    geom_point(
      size = 1.5,
      alpha = 0.8
    ) +
    geom_smooth(
      aes(group = Group),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.45
    ) +
    geom_text(
      data = stats,
      aes(
        x = x_pos,
        y = label_y,
        label = Label,
        colour = Group
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.5,
      show.legend = FALSE
    ) +
    facet_grid(
      rows = vars(lncRNA),
      cols = vars(mRNA),
      scales = "free",
      space = "fixed",
      switch = "y"
    ) +
    labs(
      title = plot_title,
      x = expression(log[10] * "(lncRNA TPM + 1)"),
      y = expression(log[10] * "(mRNA TPM + 1)"),
      colour = "Group",
      shape = "Group"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      strip.placement = "outside",
      
      strip.text.x = element_text(
        size = 10,
        face = "bold"
      ),
      
      strip.text.y.right = element_text(
        size = 10,
        face = "bold",
        angle = 270
      ),
      
      strip.background = element_rect(
        fill = "grey90",
        colour = "black"
      ),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      legend.position = "bottom",
      panel.spacing = grid::unit(1, "mm")
    )
  
  # Convert plot to a gtable
  gt <- ggplotGrob(p)
  
  # Identify the rightmost panel column
  panel_cols <- gt$layout$l[
    grepl("^panel", gt$layout$name)
  ]
  
  rightmost_panel_col <- max(panel_cols)
  
  # Find the right-side row strips
  strip_rows <- which(
    grepl("^strip-r", gt$layout$name)
  )
  
  if (length(strip_rows) > 0) {
    
    # Add space for row labels to the right of the panels
    gt <- gtable::gtable_add_cols(
      gt,
      widths = gt$widths[rightmost_panel_col],
      pos = rightmost_panel_col
    )
    
    # Move row strips to the new right-side column
    for (i in strip_rows) {
      gt$layout$l[i] <- rightmost_panel_col + 1
      gt$layout$r[i] <- rightmost_panel_col + 1
    }
  }
  
  # Display the plot
  grid::grid.newpage()
  grid::grid.draw(gt)
  
  # Save the plot
  ggsave(
    filename = output_file,
    plot = gt,
    width = 22,
    height = 10,
    units = "in",
    dpi = 300
  )
  
  # Export correlation statistics for all 27 pairs
  write.csv(
    stats,
    sub("\\.png$", "_statistics.csv", output_file),
    row.names = FALSE
  )
  
  invisible(p)
}

# ---------------------------------------------------------------------
# 9.7. FIGURE 1: NORMAL VS SRA TUMOURS
# ---------------------------------------------------------------------

p_normal_sra <- make_27_pair_plot(
  data = plot_data,
  groups_to_include = c("Normal", "SRA Tumour"),
  group_labels = c("Normal", "SRA Tumour"),
  plot_title = "All 27 lncRNA-mRNA Pairs: Normal vs SRA Tumours",
  output_file = "DLBCL_Normal_vs_SRA_All27.png"
)

# ---------------------------------------------------------------------
# 9.8. FIGURE 2: NORMAL VS TCGA TUMOURS
# ---------------------------------------------------------------------

p_normal_tcga <- make_27_pair_plot(
  data = plot_data,
  groups_to_include = c("Normal", "TCGA Tumour"),
  group_labels = c("Normal", "TCGA Tumour"),
  plot_title = "All 27 lncRNA-mRNA Pairs: Normal vs TCGA Tumours",
  output_file = "DLBCL_Normal_vs_TCGA_All27.png"
)

# ---------------------------------------------------------------------
# 9.9. FIGURE 3: SRA VS TCGA TUMOURS
# ---------------------------------------------------------------------

p_sra_tcga <- make_27_pair_plot(
  data = plot_data,
  groups_to_include = c("SRA Tumour", "TCGA Tumour"),
  group_labels = c("SRA Tumour", "TCGA Tumour"),
  plot_title = "All 27 lncRNA-mRNA Pairs: SRA vs TCGA Tumours",
  output_file = "DLBCL_SRA_vs_TCGA_All27.png"
)

cat(
  "\nAll three 27-pair figures and their statistics have been saved to:\n",
  getwd(),
  "\n"
)
