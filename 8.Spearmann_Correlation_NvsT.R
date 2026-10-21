############################################################
# DLBCL lncRNA-mRNA correlation analysis
#
# Analyses:
# 1. Spearman correlations within groups
# 2. Normal naive/GCB subset correlations
# 3. Partial Spearman correlations adjusted for cohort
#    and stromal score in tumour samples
# 4. Differential-correlation permutation tests
# 5. Adjusted differential-association permutation tests
# 6. Analyses within annotated cell/tumour subsets
# 7. BH correction
# 8. Complete CSV exports
#
# Permutations: 9,999
############################################################


############################################################
# 1. LOAD PACKAGES
############################################################

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "purrr",
  "readr",
  "stringr",
  "tibble"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)
library(stringr)
library(tibble)


############################################################
# 2. SETTINGS
############################################################

set.seed(2026)

n_perm <- 9999

expression_file <- "TPM_DLBCL2_lncRNA_mRNA_boxplot.csv"

stromal_file <- "DLBCL_stromal_scores_Results.csv"

metadata_file <- "SRA_TCGA_DLBCL_for_PCA_meta.csv"

output_dir <- "DLBCL_Correlation_Results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 3. DEFINE CANDIDATE TRANSCRIPTS
############################################################

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

candidate_pairs <- tidyr::expand_grid(
  lncRNA = lncRNAs,
  mRNA = mRNAs
)

stopifnot(nrow(candidate_pairs) == 27)

candidate_transcripts <- unique(
  c(lncRNAs, mRNAs)
)


############################################################
# 4. READ TPM MATRIX AND EXTRACT CANDIDATE TRANSCRIPTS
############################################################

tpm_all <- read.csv(
  expression_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Identify sample columns
sample_columns <- setdiff(
  colnames(tpm_all),
  "TranscriptSymbol"
)

# Define the 12 candidate transcripts
candidate_transcripts <- unique(
  c(lncRNAs, mRNAs)
)

# Extract only the candidate transcripts
candidate_tpm <- tpm_all %>%
  dplyr::filter(
    TranscriptSymbol %in% candidate_transcripts
  )

# Create expression matrix
rownames(candidate_tpm) <- candidate_tpm$TranscriptSymbol

expression_matrix <- as.matrix(
  candidate_tpm[
    ,
    sample_columns,
    drop = FALSE
  ]
)

storage.mode(expression_matrix) <- "numeric"


############################################################
# 5. DEFINE SAMPLE GROUPS AND EXPRESSION DATA
############################################################

normal_samples <- paste0("C", 1:8)

naive_samples <- paste0("C", 1:4)

gcb_normal_samples <- paste0("C", 5:8)

sra_tumour_samples <- paste0("T", 1:10)

tcga_tumour_samples <- paste0("T", 11:58)

all_tumour_samples <- c(
  sra_tumour_samples,
  tcga_tumour_samples
)

required_samples <- c(
  normal_samples,
  all_tumour_samples
)

missing_samples <- setdiff(
  required_samples,
  colnames(expression_matrix)
)

if (length(missing_samples) > 0) {
  stop(
    "Missing samples in the expression matrix: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Convert expression matrix to sample-oriented data frame
expression_df <- as.data.frame(
  t(expression_matrix),
  check.names = FALSE
)

expression_df$Sample <- rownames(expression_df)

expression_df <- expression_df %>%
  dplyr::relocate(Sample)

############################################################
# 6. READ STROMAL SCORES AND METADATA
############################################################

stromal_scores <- read.csv(
  "DLBCL_stromal_scores_Results.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  "SRA_TCGA_DLBCL_for_PCA_meta.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Sample" %in% colnames(stromal_scores)) {
  stop(
    "The stromal score file must contain a 'Sample' column."
  )
}

if (!"StromalScore" %in% colnames(stromal_scores)) {
  stop(
    "The stromal score file must contain 'StromalScore'."
  )
}

if (anyDuplicated(stromal_scores$Sample)) {
  stop(
    "Duplicate sample IDs occur in the stromal score file. ",
    "Resolve them before continuing."
  )
}

if (!"samplename" %in% colnames(metadata)) {
  stop(
    "The metadata file must contain 'samplename'. ",
    "Inspect names(metadata) and change the code if needed."
  )
}

if (anyDuplicated(metadata$samplename)) {
  warning(
    "The metadata contains repeated sample IDs. ",
    "Only unique sample IDs will be retained for metadata joining."
  )
  
  metadata <- metadata %>%
    distinct(samplename, .keep_all = TRUE)
}

# Standardise sample identifiers.
stromal_scores$Sample <- trimws(
  as.character(stromal_scores$Sample)
)

metadata$samplename <- trimws(
  as.character(metadata$samplename)
)

# Keep the stromal score and metadata fields relevant to analysis.
stromal_keep <- intersect(
  c(
    "Sample",
    "StromalScore",
    "ImmuneScore",
    "ESTIMATEScore",
    "TumorPurity",
    "SampleType",
    "Cohort",
    "Subset"
  ),
  colnames(stromal_scores)
)

stromal_scores <- stromal_scores %>%
  select(all_of(stromal_keep))

# Avoid duplicate column names when joining.
metadata_fields <- setdiff(
  colnames(metadata),
  intersect(
    colnames(metadata),
    colnames(stromal_scores)
  )
)

metadata_for_join <- metadata %>%
  select(
    samplename,
    all_of(metadata_fields)
  ) %>%
  rename(Sample = samplename)

# The stromal file is the primary source for StromalScore,
# Cohort, SampleType and Subset when these fields exist there.
sample_metadata <- stromal_scores %>%
  left_join(
    metadata_for_join,
    by = "Sample"
  )

# Add group definitions from the sample IDs.
sample_metadata <- sample_metadata %>%
  mutate(
    Group = case_when(
      Sample %in% naive_samples ~ "Normal_Naive",
      Sample %in% gcb_normal_samples ~ "Normal_GCB",
      Sample %in% sra_tumour_samples ~ "SRA_Tumour",
      Sample %in% tcga_tumour_samples ~ "TCGA_Tumour",
      TRUE ~ NA_character_
    ),
    Tissue = case_when(
      Sample %in% normal_samples ~ "Normal",
      Sample %in% all_tumour_samples ~ "Tumour",
      TRUE ~ NA_character_
    ),
    CohortFromID = case_when(
      Sample %in% sra_tumour_samples ~ "SRA",
      Sample %in% tcga_tumour_samples ~ "TCGA",
      Sample %in% normal_samples ~ "SRA",
      TRUE ~ NA_character_
    )
  )

# For the explicitly defined groups, use CohortFromID to avoid
# accidental inconsistencies in cohort labels.
if (!"Cohort" %in% colnames(sample_metadata)) {
  sample_metadata$Cohort <- sample_metadata$CohortFromID
}

sample_metadata <- sample_metadata %>%
  mutate(
    Cohort = ifelse(
      Sample %in% required_samples,
      CohortFromID,
      as.character(Cohort)
    ),
    StromalScore = as.numeric(StromalScore)
  )

# Verify all required samples have stromal scores.
missing_stromal <- sample_metadata %>%
  filter(Sample %in% required_samples) %>%
  filter(is.na(StromalScore)) %>%
  pull(Sample)

if (length(missing_stromal) > 0) {
  stop(
    "Missing stromal scores for: ",
    paste(missing_stromal, collapse = ", "),
    ". Resolve these before running adjusted analyses."
  )
}

# Ensure each required sample has metadata.
missing_metadata <- setdiff(
  required_samples,
  sample_metadata$Sample
)

if (length(missing_metadata) > 0) {
  stop(
    "Samples missing from the metadata/stromal file: ",
    paste(missing_metadata, collapse = ", ")
  )
}

write.csv(
  sample_metadata,
  file.path(output_dir, "Analysis_Sample_Metadata.csv"),
  row.names = FALSE
)


############################################################
# 7. BUILD ANALYSIS DATA
############################################################

analysis_data <- expression_df %>%
  inner_join(
    sample_metadata,
    by = "Sample"
  ) %>%
  filter(Sample %in% required_samples)

# Export the expression matrix with aligned sample metadata.
write.csv(
  analysis_data,
  file.path(output_dir, "Expression_With_Metadata.csv"),
  row.names = FALSE
)


############################################################
# 8. SPEARMAN CORRELATION FUNCTION
############################################################

safe_spearman <- function(x, y) {
  
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
      tibble(
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
  
  tibble(
    n = n,
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}


############################################################
# 9. CORRELATIONS WITHIN EACH GROUP
############################################################

calculate_group_correlations <- function(
    data,
    samples,
    group_name
) {
  
  dat <- data %>%
    filter(Sample %in% samples)
  
  results <- purrr::map_dfr(
    seq_len(nrow(candidate_pairs)),
    function(i) {
      
      lnc <- candidate_pairs$lncRNA[i]
      mrna <- candidate_pairs$mRNA[i]
      
      stat <- safe_spearman(
        dat[[lnc]],
        dat[[mrna]]
      )
      
      tibble(
        Group = group_name,
        lncRNA = lnc,
        mRNA = mrna
      ) %>%
        bind_cols(stat)
    }
  )
  
  results %>%
    mutate(
      FDR_BH = p.adjust(
        p_value,
        method = "BH"
      )
    )
}


normal_results <- calculate_group_correlations(
  analysis_data,
  normal_samples,
  "Normal_All"
)

naive_results <- calculate_group_correlations(
  analysis_data,
  naive_samples,
  "Normal_Naive"
)

gcb_normal_results <- calculate_group_correlations(
  analysis_data,
  gcb_normal_samples,
  "Normal_GCB"
)

sra_results <- calculate_group_correlations(
  analysis_data,
  sra_tumour_samples,
  "SRA_Tumour"
)

tcga_results <- calculate_group_correlations(
  analysis_data,
  tcga_tumour_samples,
  "TCGA_Tumour"
)

write.csv(
  normal_results,
  file.path(output_dir, "Normal_All_Correlations.csv"),
  row.names = FALSE
)

write.csv(
  naive_results,
  file.path(output_dir, "Normal_Naive_Correlations.csv"),
  row.names = FALSE
)

write.csv(
  gcb_normal_results,
  file.path(output_dir, "Normal_GCB_Correlations.csv"),
  row.names = FALSE
)

write.csv(
  sra_results,
  file.path(output_dir, "SRA_Tumour_Correlations.csv"),
  row.names = FALSE
)

write.csv(
  tcga_results,
  file.path(output_dir, "TCGA_Tumour_Correlations.csv"),
  row.names = FALSE
)


############################################################
# 10. PARTIAL SPEARMAN CORRELATIONS IN TUMOURS
#
# Rank-transform lncRNA and mRNA expression.
# Residualise each rank against:
#   - Cohort
#   - StromalScore
#
# Correlate the residuals.
#
# This estimates a partial rank correlation, not the
# regression coefficient from the earlier adjusted model.
############################################################

partial_spearman <- function(
    data,
    lnc,
    mrna,
    covariate_formula
) {
  
  needed <- unique(
    c(
      lnc,
      mrna,
      all.vars(covariate_formula)
    )
  )
  
  dat <- data %>%
    select(all_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  
  n <- nrow(dat)
  
  if (n < 5) {
    return(
      tibble(
        n = n,
        partial_rho = NA_real_,
        p_value = NA_real_,
        df = NA_real_
      )
    )
  }
  
  x <- rank(dat[[lnc]], ties.method = "average")
  y <- rank(dat[[mrna]], ties.method = "average")
  
  covariate_model <- model.matrix(
    covariate_formula,
    data = dat
  )
  
  x_resid <- lm.fit(
    x = covariate_model,
    y = x
  )$residuals
  
  y_resid <- lm.fit(
    x = covariate_model,
    y = y
  )$residuals
  
  if (
    sd(x_resid) == 0 ||
    sd(y_resid) == 0
  ) {
    return(
      tibble(
        n = n,
        partial_rho = NA_real_,
        p_value = NA_real_,
        df = NA_real_
      )
    )
  }
  
  r <- cor(
    x_resid,
    y_resid,
    method = "pearson"
  )
  
  # Number of estimated nuisance parameters, excluding
  # the intercept. Cohort + StromalScore gives k = 2
  # when Cohort is a two-level factor.
  q <- qr(covariate_model)$rank - 1
  
  df <- n - q - 2
  
  if (df <= 0 || abs(r) >= 1) {
    p <- if (abs(r) >= 1 && df > 0) 0 else NA_real_
  } else {
    t_stat <- r * sqrt(df / (1 - r^2))
    p <- 2 * pt(-abs(t_stat), df = df)
  }
  
  tibble(
    n = n,
    partial_rho = r,
    p_value = p,
    df = df
  )
}


tumour_data <- analysis_data %>%
  filter(Sample %in% all_tumour_samples) %>%
  mutate(
    Cohort = factor(
      Cohort,
      levels = c("SRA", "TCGA")
    )
  )

if (nlevels(droplevels(tumour_data$Cohort)) != 2) {
  stop(
    "The tumour metadata must contain both SRA and TCGA cohorts."
  )
}

adjusted_partial_results <- purrr::map_dfr(
  seq_len(nrow(candidate_pairs)),
  function(i) {
    
    lnc <- candidate_pairs$lncRNA[i]
    mrna <- candidate_pairs$mRNA[i]
    
    stat <- partial_spearman(
      data = tumour_data,
      lnc = lnc,
      mrna = mrna,
      covariate_formula = ~ Cohort + StromalScore
    )
    
    tibble(
      lncRNA = lnc,
      mRNA = mrna
    ) %>%
      bind_cols(stat)
  }
) %>%
  mutate(
    FDR_BH = p.adjust(
      p_value,
      method = "BH"
    )
  )

write.csv(
  adjusted_partial_results,
  file.path(
    output_dir,
    "Tumour_Partial_Spearman_Adjusted_Cohort_Stromal.csv"
  ),
  row.names = FALSE
)


############################################################
# 11. RAW DIFFERENTIAL-CORRELATION PERMUTATION TEST
#
# Tests whether Spearman rho differs between two groups.
#
# Statistic:
#   delta_rho = rho(group 2) - rho(group 1)
#
# Under the null, group labels are permuted while preserving
# the original group sizes.
#
# This is an UNADJUSTED differential-Spearman test.
############################################################

permutation_delta_rho <- function(
    x,
    y,
    group,
    n_perm = 9999
) {
  
  valid <- is.finite(x) &
    is.finite(y) &
    !is.na(group)
  
  x <- x[valid]
  y <- y[valid]
  group <- droplevels(factor(group[valid]))
  
  if (nlevels(group) != 2) {
    return(
      tibble(
        n1 = NA_integer_,
        n2 = NA_integer_,
        rho1 = NA_real_,
        rho2 = NA_real_,
        delta_rho = NA_real_,
        p_permutation = NA_real_
      )
    )
  }
  
  idx1 <- which(group == levels(group)[1])
  idx2 <- which(group == levels(group)[2])
  
  n1 <- length(idx1)
  n2 <- length(idx2)
  
  if (
    n1 < 3 ||
    n2 < 3
  ) {
    return(
      tibble(
        n1 = n1,
        n2 = n2,
        rho1 = NA_real_,
        rho2 = NA_real_,
        delta_rho = NA_real_,
        p_permutation = NA_real_
      )
    )
  }
  
  rho1 <- suppressWarnings(
    cor(
      x[idx1],
      y[idx1],
      method = "spearman"
    )
  )
  
  rho2 <- suppressWarnings(
    cor(
      x[idx2],
      y[idx2],
      method = "spearman"
    )
  )
  
  observed_delta <- rho2 - rho1
  
  permuted_deltas <- replicate(
    n_perm,
    {
      perm_idx <- sample(
        seq_along(group),
        replace = FALSE
      )
      
      perm_group <- group[perm_idx]
      
      p1 <- which(perm_group == levels(group)[1])
      p2 <- which(perm_group == levels(group)[2])
      
      r1 <- suppressWarnings(
        cor(
          x[p1],
          y[p1],
          method = "spearman"
        )
      )
      
      r2 <- suppressWarnings(
        cor(
          x[p2],
          y[p2],
          method = "spearman"
        )
      )
      
      r2 - r1
    }
  )
  
  p_perm <- (
    1 + sum(
      abs(permuted_deltas) >= abs(observed_delta),
      na.rm = TRUE
    )
  ) / (
    1 + sum(is.finite(permuted_deltas))
  )
  
  tibble(
    n1 = n1,
    n2 = n2,
    rho1 = rho1,
    rho2 = rho2,
    delta_rho = observed_delta,
    p_permutation = p_perm
  )
}


############################################################
# 12. RUN RAW DIFFERENTIAL-CORRELATION TESTS
############################################################

run_delta_rho_comparison <- function(
    data,
    samples1,
    samples2,
    label1,
    label2,
    comparison_name,
    n_perm = 9999
) {
  
  dat <- data %>%
    filter(
      Sample %in% c(samples1, samples2)
    ) %>%
    mutate(
      ComparisonGroup = case_when(
        Sample %in% samples1 ~ label1,
        Sample %in% samples2 ~ label2,
        TRUE ~ NA_character_
      ),
      ComparisonGroup = factor(
        ComparisonGroup,
        levels = c(label1, label2)
      )
    )
  
  results <- purrr::map_dfr(
    seq_len(nrow(candidate_pairs)),
    function(i) {
      
      lnc <- candidate_pairs$lncRNA[i]
      mrna <- candidate_pairs$mRNA[i]
      
      stat <- permutation_delta_rho(
        x = dat[[lnc]],
        y = dat[[mrna]],
        group = dat$ComparisonGroup,
        n_perm = n_perm
      )
      
      tibble(
        Comparison = comparison_name,
        Group1 = label1,
        Group2 = label2,
        lncRNA = lnc,
        mRNA = mrna
      ) %>%
        bind_cols(stat)
    }
  )
  
  results %>%
    mutate(
      FDR_BH = p.adjust(
        p_permutation,
        method = "BH"
      )
    )
}


normal_vs_sra <- run_delta_rho_comparison(
  data = analysis_data,
  samples1 = normal_samples,
  samples2 = sra_tumour_samples,
  label1 = "Normal",
  label2 = "SRA_Tumour",
  comparison_name = "Normal_vs_SRA_Tumour",
  n_perm = n_perm
)

normal_vs_tcga <- run_delta_rho_comparison(
  data = analysis_data,
  samples1 = normal_samples,
  samples2 = tcga_tumour_samples,
  label1 = "Normal",
  label2 = "TCGA_Tumour",
  comparison_name = "Normal_vs_TCGA_Tumour",
  n_perm = n_perm
)

sra_vs_tcga <- run_delta_rho_comparison(
  data = analysis_data,
  samples1 = sra_tumour_samples,
  samples2 = tcga_tumour_samples,
  label1 = "SRA_Tumour",
  label2 = "TCGA_Tumour",
  comparison_name = "SRA_vs_TCGA_Tumour",
  n_perm = n_perm
)

write.csv(
  normal_vs_sra,
  file.path(
    output_dir,
    "Normal_vs_SRA_Differential_Correlations.csv"
  ),
  row.names = FALSE
)

write.csv(
  normal_vs_tcga,
  file.path(
    output_dir,
    "Normal_vs_TCGA_Differential_Correlations.csv"
  ),
  row.names = FALSE
)

write.csv(
  sra_vs_tcga,
  file.path(
    output_dir,
    "SRA_vs_TCGA_Differential_Correlations.csv"
  ),
  row.names = FALSE
)


############################################################
# 13. ADJUSTED DIFFERENTIAL-ASSOCIATION PERMUTATION TEST
#
# Tests whether the association between ranked lncRNA and
# ranked mRNA differs between groups after adjusting for
# specified covariates.
#
# Uses a Freedman-Lane-style residual permutation under the
# reduced model. The statistic is the t statistic for the
# lncRNA:group interaction in a rank-regression model.
#
# IMPORTANT:
# This tests a difference in adjusted rank-regression slopes,
# not a difference in Spearman rho.
############################################################

permutation_adjusted_interaction <- function(
    data,
    lnc,
    mrna,
    group_var,
    covariates = character(0),
    n_perm = 9999
) {
  
  needed <- unique(
    c(
      lnc,
      mrna,
      group_var,
      covariates
    )
  )
  
  dat <- data %>%
    select(all_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  
  dat[[group_var]] <- droplevels(
    factor(dat[[group_var]])
  )
  
  n <- nrow(dat)
  
  if (
    n < 8 ||
    nlevels(dat[[group_var]]) != 2
  ) {
    return(
      tibble(
        n = n,
        interaction_beta = NA_real_,
        interaction_t = NA_real_,
        p_permutation = NA_real_
      )
    )
  }
  
  # Rank-transform expression.
  dat$x_rank <- rank(
    dat[[lnc]],
    ties.method = "average"
  )
  
  dat$y_rank <- rank(
    dat[[mrna]],
    ties.method = "average"
  )
  
  # Construct formulas.
  reduced_rhs <- c(
    "x_rank",
    group_var,
    covariates
  )
  
  full_rhs <- c(
    paste0("x_rank * ", group_var),
    covariates
  )
  
  reduced_formula <- as.formula(
    paste(
      "y_rank ~",
      paste(reduced_rhs, collapse = " + ")
    )
  )
  
  full_formula <- as.formula(
    paste(
      "y_rank ~",
      paste(full_rhs, collapse = " + ")
    )
  )
  
  reduced_fit <- lm(
    reduced_formula,
    data = dat
  )
  
  full_fit <- lm(
    full_formula,
    data = dat
  )
  
  if (anyNA(coef(full_fit))) {
    return(
      tibble(
        n = n,
        interaction_beta = NA_real_,
        interaction_t = NA_real_,
        p_permutation = NA_real_
      )
    )
  }
  
  interaction_terms <- grep(
    ":",
    names(coef(full_fit)),
    value = TRUE
  )
  
  interaction_terms <- interaction_terms[
    grepl("x_rank", interaction_terms, fixed = TRUE)
  ]
  
  if (length(interaction_terms) != 1) {
    return(
      tibble(
        n = n,
        interaction_beta = NA_real_,
        interaction_t = NA_real_,
        p_permutation = NA_real_
      )
    )
  }
  
  interaction_term <- interaction_terms[1]
  
  observed_beta <- coef(full_fit)[interaction_term]
  
  observed_t <- summary(full_fit)$coefficients[
    interaction_term,
    "t value"
  ]
  
  fitted_reduced <- fitted(reduced_fit)
  
  residual_reduced <- residuals(reduced_fit)
  
  permuted_t <- replicate(
    n_perm,
    {
      
      # Permute residuals under the reduced model.
      y_perm <- fitted_reduced +
        sample(
          residual_reduced,
          replace = FALSE
        )
      
      dat_perm <- dat
      dat_perm$y_rank <- y_perm
      
      fit_perm <- try(
        lm(
          full_formula,
          data = dat_perm
        ),
        silent = TRUE
      )
      
      if (inherits(fit_perm, "try-error")) {
        return(NA_real_)
      }
      
      coef_perm <- coef(fit_perm)
      
      if (!interaction_term %in% names(coef_perm)) {
        return(NA_real_)
      }
      
      if (is.na(coef_perm[interaction_term])) {
        return(NA_real_)
      }
      
      summary_perm <- summary(fit_perm)$coefficients
      
      if (!interaction_term %in% rownames(summary_perm)) {
        return(NA_real_)
      }
      
      summary_perm[interaction_term, "t value"]
    }
  )
  
  valid_perm <- is.finite(permuted_t)
  
  p_perm <- (
    1 + sum(
      abs(permuted_t[valid_perm]) >= abs(observed_t)
    )
  ) / (
    1 + sum(valid_perm)
  )
  
  tibble(
    n = n,
    interaction_beta = unname(observed_beta),
    interaction_t = unname(observed_t),
    p_permutation = p_perm
  )
}


############################################################
# 14. RUN ADJUSTED INTERACTION TESTS
############################################################

# A. SRA versus TCGA tumours:
# Adjust for StromalScore.
# Cohort is the comparison group, so its main effect is
# included as part of the interaction model.

tumour_interaction_data <- analysis_data %>%
  filter(Sample %in% all_tumour_samples) %>%
  mutate(
    TumourCohort = case_when(
      Sample %in% sra_tumour_samples ~ "SRA",
      Sample %in% tcga_tumour_samples ~ "TCGA",
      TRUE ~ NA_character_
    ),
    TumourCohort = factor(
      TumourCohort,
      levels = c("SRA", "TCGA")
    )
  )

adjusted_sra_tcga <- purrr::map_dfr(
  seq_len(nrow(candidate_pairs)),
  function(i) {
    
    lnc <- candidate_pairs$lncRNA[i]
    mrna <- candidate_pairs$mRNA[i]
    
    stat <- permutation_adjusted_interaction(
      data = tumour_interaction_data,
      lnc = lnc,
      mrna = mrna,
      group_var = "TumourCohort",
      covariates = c("StromalScore"),
      n_perm = n_perm
    )
    
    tibble(
      Comparison = "SRA_vs_TCGA_Adjusted",
      lncRNA = lnc,
      mRNA = mrna
    ) %>%
      bind_cols(stat)
  }
) %>%
  mutate(
    FDR_BH = p.adjust(
      p_permutation,
      method = "BH"
    )
  )

write.csv(
  adjusted_sra_tcga,
  file.path(
    output_dir,
    "SRA_vs_TCGA_Adjusted_Interaction_Permutation.csv"
  ),
  row.names = FALSE
)


# B. Normal versus SRA tumour:
# Both groups are from SRA.
# Adjust for StromalScore.
#
# This model tests an adjusted rank-regression slope
# difference; it does not remove biological differences
# between naive and GCB normal B cells.

normal_sra_data <- analysis_data %>%
  filter(
    Sample %in% c(
      normal_samples,
      sra_tumour_samples
    )
  ) %>%
  mutate(
    ComparisonGroup = case_when(
      Sample %in% normal_samples ~ "Normal",
      Sample %in% sra_tumour_samples ~ "SRA_Tumour",
      TRUE ~ NA_character_
    ),
    ComparisonGroup = factor(
      ComparisonGroup,
      levels = c("Normal", "SRA_Tumour")
    )
  )

adjusted_normal_sra <- purrr::map_dfr(
  seq_len(nrow(candidate_pairs)),
  function(i) {
    
    lnc <- candidate_pairs$lncRNA[i]
    mrna <- candidate_pairs$mRNA[i]
    
    stat <- permutation_adjusted_interaction(
      data = normal_sra_data,
      lnc = lnc,
      mrna = mrna,
      group_var = "ComparisonGroup",
      covariates = c("StromalScore"),
      n_perm = n_perm
    )
    
    tibble(
      Comparison = "Normal_vs_SRA_Adjusted",
      lncRNA = lnc,
      mRNA = mrna
    ) %>%
      bind_cols(stat)
  }
) %>%
  mutate(
    FDR_BH = p.adjust(
      p_permutation,
      method = "BH"
    )
  )

write.csv(
  adjusted_normal_sra,
  file.path(
    output_dir,
    "Normal_vs_SRA_Adjusted_Interaction_Permutation.csv"
  ),
  row.names = FALSE
)


# C. Normal versus TCGA tumour:
#
# Cohort and tissue status are perfectly confounded:
# normals are SRA, tumours are TCGA.
# Consequently, a model cannot separately adjust for cohort
# and tissue status. This model adjusts for StromalScore only.

normal_tcga_data <- analysis_data %>%
  filter(
    Sample %in% c(
      normal_samples,
      tcga_tumour_samples
    )
  ) %>%
  mutate(
    ComparisonGroup = case_when(
      Sample %in% normal_samples ~ "Normal",
      Sample %in% tcga_tumour_samples ~ "TCGA_Tumour",
      TRUE ~ NA_character_
    ),
    ComparisonGroup = factor(
      ComparisonGroup,
      levels = c("Normal", "TCGA_Tumour")
    )
  )

adjusted_normal_tcga <- purrr::map_dfr(
  seq_len(nrow(candidate_pairs)),
  function(i) {
    
    lnc <- candidate_pairs$lncRNA[i]
    mrna <- candidate_pairs$mRNA[i]
    
    stat <- permutation_adjusted_interaction(
      data = normal_tcga_data,
      lnc = lnc,
      mrna = mrna,
      group_var = "ComparisonGroup",
      covariates = c("StromalScore"),
      n_perm = n_perm
    )
    
    tibble(
      Comparison = "Normal_vs_TCGA_Adjusted_StromalOnly",
      lncRNA = lnc,
      mRNA = mrna
    ) %>%
      bind_cols(stat)
  }
) %>%
  mutate(
    FDR_BH = p.adjust(
      p_permutation,
      method = "BH"
    )
  )

write.csv(
  adjusted_normal_tcga,
  file.path(
    output_dir,
    "Normal_vs_TCGA_Stromal_Adjusted_Interaction_Permutation.csv"
  ),
  row.names = FALSE
)


# 15. CORRELATIONS WITHIN ANNOTATED SUBSETS
#
# First analyse the two explicitly known normal subsets.
# Then, if metadata contains Subset, analyse every available
# subset represented by at least 3 samples.
#
# Subset analyses are exploratory where sample sizes are low.
############################################################

subset_results <- tibble()

# Explicit normal subsets.
subset_results <- bind_rows(
  subset_results,
  naive_results,
  gcb_normal_results
)

# Use the available Subset annotation, if present.
if ("Subset" %in% colnames(analysis_data)) {
  
  subset_data <- analysis_data %>%
    filter(
      !is.na(Subset),
      Sample %in% required_samples
    ) %>%
    mutate(
      Subset = as.character(Subset)
    )
  
  available_subsets <- subset_data %>%
    count(Subset, name = "n_samples") %>%
    filter(n_samples >= 3) %>%
    pull(Subset)
  
  for (subset_name in available_subsets) {
    
    subset_samples <- subset_data %>%
      filter(Subset == subset_name) %>%
      pull(Sample)
    
    # Avoid duplicating the explicitly defined naive/GCB
    # normal groups when the metadata labels already match.
    subset_label <- paste0(
      "Metadata_Subset_",
      make.names(subset_name)
    )
    
    subset_cor <- calculate_group_correlations(
      data = analysis_data,
      samples = subset_samples,
      group_name = subset_label
    )
    
    subset_results <- bind_rows(
      subset_results,
      subset_cor
    )
  }
}

write.csv(
  subset_results,
  file.path(
    output_dir,
    "All_Cell_and_Annotated_Subset_Correlations.csv"
  ),
  row.names = FALSE
)


############################################################
# 16. COMBINE THE MAIN GROUP CORRELATION RESULTS
############################################################

all_group_correlations <- bind_rows(
  normal_results,
  naive_results,
  gcb_normal_results,
  sra_results,
  tcga_results
)

write.csv(
  all_group_correlations,
  file.path(
    output_dir,
    "All_Main_Group_Correlations.csv"
  ),
  row.names = FALSE
)


############################################################
# 17. COMBINE DIFFERENTIAL-CORRELATION RESULTS
############################################################

all_differential_correlations <- bind_rows(
  normal_vs_sra,
  normal_vs_tcga,
  sra_vs_tcga
)

write.csv(
  all_differential_correlations,
  file.path(
    output_dir,
    "All_Differential_Correlation_Tests.csv"
  ),
  row.names = FALSE
)


############################################################
# 18. COMBINE ADJUSTED INTERACTION RESULTS
############################################################

all_adjusted_interactions <- bind_rows(
  adjusted_normal_sra,
  adjusted_normal_tcga,
  adjusted_sra_tcga
)

write.csv(
  all_adjusted_interactions,
  file.path(
    output_dir,
    "All_Adjusted_Interaction_Permutation_Tests.csv"
  ),
  row.names = FALSE
)


############################################################
# 19. COMBINED MASTER TABLE
#
# Preserve all pairs even if one of the tests is missing.
############################################################

master_table <- candidate_pairs %>%
  left_join(
    normal_results %>%
      select(
        lncRNA,
        mRNA,
        Normal_rho = rho,
        Normal_p = p_value,
        Normal_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  ) %>%
  left_join(
    sra_results %>%
      select(
        lncRNA,
        mRNA,
        SRA_rho = rho,
        SRA_p = p_value,
        SRA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  ) %>%
  left_join(
    tcga_results %>%
      select(
        lncRNA,
        mRNA,
        TCGA_rho = rho,
        TCGA_p = p_value,
        TCGA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  ) %>%
  left_join(
    adjusted_partial_results %>%
      select(
        lncRNA,
        mRNA,
        Partial_rho = partial_rho,
        Partial_p = p_value,
        Partial_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  ) %>%
  left_join(
    normal_vs_sra %>%
      select(
        lncRNA,
        mRNA,
        Normal_SRA_DeltaRho = delta_rho,
        Normal_SRA_PermutationP = p_permutation,
        Normal_SRA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  ) %>%
  left_join(
    normal_vs_tcga %>%
      select(
        lncRNA,
        mRNA,
        Normal_TCGA_DeltaRho = delta_rho,
        Normal_TCGA_PermutationP = p_permutation,
        Normal_TCGA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  ) %>%
  left_join(
    sra_vs_tcga %>%
      select(
        lncRNA,
        mRNA,
        SRA_TCGA_DeltaRho = delta_rho,
        SRA_TCGA_PermutationP = p_permutation,
        SRA_TCGA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )

write.csv(
  master_table,
  file.path(
    output_dir,
    "MASTER_DLBCL_Correlation_Results.csv"
  ),
  row.names = FALSE
)


############################################################
# 20. SUMMARY OF SIGNIFICANT RESULTS
############################################################

summary_significance <- tibble(
  Analysis = c(
    "Normal correlations",
    "SRA tumour correlations",
    "TCGA tumour correlations",
    "Adjusted partial Spearman",
    "Normal vs SRA differential rho",
    "Normal vs TCGA differential rho",
    "SRA vs TCGA differential rho",
    "Normal vs SRA adjusted interaction",
    "Normal vs TCGA adjusted interaction",
    "SRA vs TCGA adjusted interaction"
  ),
  Number_of_pairs = 27,
  Number_FDR_below_005 = c(
    sum(normal_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(sra_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(tcga_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(adjusted_partial_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(normal_vs_sra$FDR_BH < 0.05, na.rm = TRUE),
    sum(normal_vs_tcga$FDR_BH < 0.05, na.rm = TRUE),
    sum(sra_vs_tcga$FDR_BH < 0.05, na.rm = TRUE),
    sum(adjusted_normal_sra$FDR_BH < 0.05, na.rm = TRUE),
    sum(adjusted_normal_tcga$FDR_BH < 0.05, na.rm = TRUE),
    sum(adjusted_sra_tcga$FDR_BH < 0.05, na.rm = TRUE)
  )
)

write.csv(
  summary_significance,
  file.path(
    output_dir,
    "Significance_Summary.csv"
  ),
  row.names = FALSE
)

print(summary_significance)


############################################################
# 21. FINISH
############################################################

message(
  "Analysis completed. Results saved to: ",
  normalizePath(output_dir)
)

############################################################
# 22. PLOT CORRELATIONS FOR ALL 27 lncRNA-mRNA PAIRS
############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(scales)


############################################################
# 22.1 OUTPUT DIRECTORY
############################################################

plot_dir <- file.path(
  output_dir,
  "Correlation_Plots"
)

dir.create(
  plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 22.2 PREPARE EXPRESSION DATA
############################################################

plot_samples <- c(
  normal_samples,
  all_tumour_samples
)

plot_expression <- expression_df %>%
  filter(Sample %in% plot_samples) %>%
  left_join(
    sample_metadata %>%
      select(
        Sample,
        StromalScore,
        Cohort,
        Subset
      ),
    by = "Sample",
    suffix = c("", "_metadata")
  ) %>%
  mutate(
    Group = case_when(
      Sample %in% naive_samples ~ "Naive B cells",
      Sample %in% gcb_normal_samples ~ "GCB B cells",
      Sample %in% sra_tumour_samples ~ "SRA tumours",
      Sample %in% tcga_tumour_samples ~ "TCGA tumours",
      TRUE ~ NA_character_
    ),
    Tissue = ifelse(
      Sample %in% normal_samples,
      "Normal",
      "Tumour"
    )
  )

# Transform TPM values for plotting only.
# Original TPM values remain unchanged in the analysis.

for (gene in candidate_transcripts) {
  plot_expression[[paste0(gene, "_log10")]] <-
    log10(plot_expression[[gene]] + 1)
}


############################################################
# 22.3 FUNCTION TO CALCULATE CORRELATION LABELS
############################################################

calculate_plot_statistics <- function(
    dat,
    group_column = "Group"
) {
  
  group_values <- unique(
    as.character(dat[[group_column]])
  )
  
  group_values <- group_values[
    !is.na(group_values)
  ]
  
  map_dfr(
    group_values,
    function(group_value) {
      
      group_data <- dat %>%
        filter(
          as.character(.data[[group_column]]) ==
            group_value
        )
      
      map_dfr(
        seq_len(nrow(candidate_pairs)),
        function(i) {
          
          lnc <- candidate_pairs$lncRNA[i]
          mrna <- candidate_pairs$mRNA[i]
          
          x <- group_data[[lnc]]
          y <- group_data[[mrna]]
          
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
              tibble(
                Group = group_value,
                lncRNA = lnc,
                mRNA = mrna,
                n = n,
                rho = NA_real_,
                p_value = NA_real_,
                FDR_BH = NA_real_
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
          
          tibble(
            Group = group_value,
            lncRNA = lnc,
            mRNA = mrna,
            n = n,
            rho = unname(test$estimate),
            p_value = test$p.value
          )
        }
      )
    }
  ) %>%
    group_by(Group) %>%
    mutate(
      FDR_BH = p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    ungroup()
}


############################################################
# 22.4 FUNCTION TO FORMAT STATISTICAL LABELS
############################################################

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(
      p < 0.001,
      "<0.001",
      sprintf("%.3f", p)
    )
  )
}

format_rho <- function(rho) {
  ifelse(
    is.na(rho),
    "NA",
    sprintf("%.2f", rho)
  )
}

make_stat_label <- function(
    rho,
    p_value,
    FDR_BH,
    n
) {
  
  paste0(
    "\u03c1 = ", format_rho(rho),
    "\nP = ", format_p(p_value),
    "\nFDR = ", format_p(FDR_BH),
    "\nn = ", n
  )
}


############################################################
# 22.5 FUNCTION TO CREATE FACETED PLOTS
# Layout: 3 rows x 9 columns
############################################################

make_27_pair_plot <- function(
    data,
    group1_samples,
    group2_samples,
    group1_label,
    group2_label,
    filename,
    plot_title
) {
  
  # Select samples for the requested comparison.
  dat <- data %>%
    filter(
      Sample %in% c(
        group1_samples,
        group2_samples
      )
    ) %>%
    mutate(
      PlotGroup = case_when(
        Sample %in% group1_samples ~ group1_label,
        Sample %in% group2_samples ~ group2_label,
        TRUE ~ NA_character_
      ),
      PlotGroup = factor(
        PlotGroup,
        levels = c(
          group1_label,
          group2_label
        )
      )
    ) %>%
    filter(!is.na(PlotGroup))
  
  # Calculate statistics using the original TPM values.
  # calculate_plot_statistics() returns a column called Group.
  # Rename it to PlotGroup for consistency with the plotting data.
  
  stats <- calculate_plot_statistics(
    dat,
    group_column = "PlotGroup"
  ) %>%
    dplyr::rename(
      PlotGroup = Group
    ) %>%
    mutate(
      PlotGroup = factor(
        PlotGroup,
        levels = c(
          group1_label,
          group2_label
        )
      )
    )
  
  # Export the statistics shown in the figure.
  
  write.csv(
    stats,
    file.path(
      plot_dir,
      paste0(
        tools::file_path_sans_ext(filename),
        "_statistics.csv"
      )
    ),
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # Long format for all 27 pairs
  # ------------------------------------------------------------
  
  plot_long <- purrr::map_dfr(
    seq_len(nrow(candidate_pairs)),
    function(i) {
      
      lnc <- candidate_pairs$lncRNA[i]
      mrna <- candidate_pairs$mRNA[i]
      
      dat %>%
        transmute(
          Sample,
          PlotGroup,
          lncRNA = lnc,
          mRNA = mrna,
          x = log10(.data[[lnc]] + 1),
          y = log10(.data[[mrna]] + 1)
        )
    }
  ) %>%
    filter(
      is.finite(x),
      is.finite(y)
    ) %>%
    mutate(
      lncRNA = factor(
        lncRNA,
        levels = lncRNAs
      ),
      mRNA = factor(
        mRNA,
        levels = mRNAs
      )
    )
  
  # Add pair names and formatted statistics.
  
  # ------------------------------------------------------------
  # Prepare statistics and annotation positions
  # ------------------------------------------------------------
  
  stats <- stats %>%
    mutate(
      lncRNA = factor(
        lncRNA,
        levels = lncRNAs
      ),
      mRNA = factor(
        mRNA,
        levels = mRNAs
      ),
      Label = make_stat_label(
        rho,
        p_value,
        FDR_BH,
        n
      )
    )
  
  # ------------------------------------------------------------
  # Position Normal and Tumour labels with greater separation
  # ------------------------------------------------------------
  # Position labels at the top-left of every panel
  positions <- plot_long %>%
    group_by(lncRNA, mRNA) %>%
    summarise(
      x_min = min(x, na.rm = TRUE),
      x_max = max(x, na.rm = TRUE),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x_range = pmax(x_max - x_min, 0.5),
      y_range = pmax(y_max - y_min, 0.5),
      
      # Fixed position near the top-left
      x_pos = x_min + 0.03 * x_range,
      y_top = y_max + 0.85 * y_range
    )
  
  # ------------------------------------------------------------
  # Anchor labels to the top-left corner of every panel
  # ------------------------------------------------------------
  
  stats <- stats %>%
    mutate(
      label_x = -Inf,
      label_y = Inf
    )
  
  # Generate the 3 x 9 faceted plot.
  
  p <- ggplot(
    plot_long,
    aes(
      x = x,
      y = y
    )
  ) +
    
    geom_point(
      aes(
        colour = PlotGroup,
        shape = PlotGroup
      ),
      size = 1.8,
      alpha = 0.8
    ) +
    
    geom_smooth(
      aes(
        colour = PlotGroup,
        group = PlotGroup
      ),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.6
    ) +
    
    geom_text(
      data = stats %>% filter(as.character(PlotGroup) == group1_label),
      aes(
        x = -Inf,
        y = Inf,
        label = Label,
        colour = PlotGroup
      ),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.2,
      size = 3,
      lineheight = 1.0,
      show.legend = FALSE
    ) +
    
    geom_text(
      data = stats %>% filter(as.character(PlotGroup) == group2_label),
      aes(
        x = -Inf,
        y = Inf,
        label = Label,
        colour = PlotGroup
      ),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 2.5,
      size = 3,
      lineheight = 1.0,
      show.legend = FALSE
    ) +
    
    facet_grid(
      rows = vars(lncRNA),
      cols = vars(mRNA),
      scales = "free",
      space = "fixed"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.15))
    )+
    scale_x_continuous(
      breaks = c(0, 0.5, 1.0, 1.5),
      expand = expansion(mult = c(0.05, 0))
    ) +
    coord_cartesian(
      xlim = c(0, 2)
    ) +
    scale_colour_manual(
      values = setNames(
        c("#F8766D", "#00BFC4"),
        c(
          group1_label,
          group2_label
        )
      )
    ) +
    
    scale_shape_manual(
      values = setNames(
        c(16, 17),
        c(
          group1_label,
          group2_label
        )
      )
    ) +
    
    labs(
      title = plot_title,
      subtitle = paste0(
        "All 27 lncRNA-mRNA pairs | ",
        "Expression = log10(TPM + 1)"
      ),
      x = expression(
        log[10]("lncRNA TPM + 1")
      ),
      y = expression(
        log[10]("mRNA TPM + 1")
      ),
      colour = "Group",
      shape = "Group"
    ) +
    
    theme_bw(
      base_size = 9
    ) +
    
    theme(
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5
      ),
      
      strip.text = element_text(
        size = 10,
        face = "bold"
      ),
      
      axis.title = element_text(
        size = 15
      ),
      
      axis.text = element_text(
        size = 15
      ),
      
      panel.spacing = grid::unit(
        0.12,
        "cm"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 15
      ),
      
      legend.text = element_text(
        size = 15
      ),
      
      panel.grid.minor = element_blank()
    )
  
  # Save the figure at 300 dpi.
  
  ggsave(
    filename = file.path(
      plot_dir,
      filename
    ),
    plot = p,
    width = 24,
    height = 11,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
  return(p)
}


############################################################
# 22.6 PLOT 1: NORMAL VS SRA TUMOURS
############################################################

plot_normal_vs_sra <- make_27_pair_plot(
  data = plot_expression,
  group1_samples = normal_samples,
  group2_samples = sra_tumour_samples,
  group1_label = "Normal",
  group2_label = "SRA tumours",
  filename = "Normal_vs_SRA_All27.png",
  plot_title = "Normal vs SRA tumour lncRNA-mRNA correlations"
)

print(plot_normal_vs_sra)



############################################################
# 22.7 PLOT 2: NORMAL VS TCGA TUMOURS
############################################################

plot_normal_vs_tcga <- make_27_pair_plot(
  data = plot_expression,
  group1_samples = normal_samples,
  group2_samples = tcga_tumour_samples,
  group1_label = "Normal",
  group2_label = "TCGA tumours",
  filename = "Normal_vs_TCGA_All27.png",
  plot_title = "Normal vs TCGA tumour lncRNA-mRNA correlations"
)

print(plot_normal_vs_tcga)


############################################################
# 22.8 PLOT 3: SRA VS TCGA TUMOURS
############################################################

plot_sra_vs_tcga <- make_27_pair_plot(
  data = plot_expression,
  group1_samples = sra_tumour_samples,
  group2_samples = tcga_tumour_samples,
  group1_label = "SRA tumours",
  group2_label = "TCGA tumours",
  filename = "SRA_vs_TCGA_All27.png",
  plot_title = "SRA vs TCGA tumour lncRNA-mRNA correlations"
)

print(plot_sra_vs_tcga)

############################################################
# 22.9 FUNCTION FOR SINGLE-GROUP CORRELATION PLOTS
# 3 rows x 9 columns
# mRNA labels at top; lncRNA labels on far right
############################################################

make_single_group_plot <- function(
    data,
    samples,
    group_label,
    filename,
    plot_title
) {
  
  # --------------------------------------------------------
  # 22.9.1 Select samples
  # --------------------------------------------------------
  
  dat <- data %>%
    filter(Sample %in% samples) %>%
    mutate(
      PlotGroup = group_label
    )
  
  # --------------------------------------------------------
  # 22.9.2 Calculate Spearman correlation statistics
  # --------------------------------------------------------
  
  stats <- calculate_plot_statistics(
    dat,
    group_column = "PlotGroup"
  ) %>%
    dplyr::rename(
      PlotGroup = Group
    ) %>%
    mutate(
      lncRNA = factor(
        lncRNA,
        levels = lncRNAs
      ),
      mRNA = factor(
        mRNA,
        levels = mRNAs
      ),
      Label = make_stat_label(
        rho,
        p_value,
        FDR_BH,
        n
      )
    )
  
  # --------------------------------------------------------
  # 22.9.3 Export statistics
  # --------------------------------------------------------
  
  write.csv(
    stats,
    file.path(
      plot_dir,
      paste0(
        tools::file_path_sans_ext(filename),
        "_statistics.csv"
      )
    ),
    row.names = FALSE
  )
  
  # --------------------------------------------------------
  # 22.9.4 Prepare expression data for all 27 pairs
  # --------------------------------------------------------
  
  plot_long <- purrr::map_dfr(
    seq_len(nrow(candidate_pairs)),
    function(i) {
      
      lnc <- candidate_pairs$lncRNA[i]
      mrna <- candidate_pairs$mRNA[i]
      
      dat %>%
        transmute(
          Sample,
          lncRNA = lnc,
          mRNA = mrna,
          x = log10(.data[[lnc]] + 1),
          y = log10(.data[[mrna]] + 1)
        )
    }
  ) %>%
    filter(
      is.finite(x),
      is.finite(y)
    ) %>%
    mutate(
      lncRNA = factor(
        lncRNA,
        levels = lncRNAs
      ),
      mRNA = factor(
        mRNA,
        levels = mRNAs
      )
    )
  
  # --------------------------------------------------------
  # 22.9.5 Prepare statistical labels and positions
  # --------------------------------------------------------
  
  label_positions <- plot_long %>%
    group_by(lncRNA, mRNA) %>%
    summarise(
      x_min = min(x, na.rm = TRUE),
      x_max = max(x, na.rm = TRUE),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      x_range = x_max - x_min,
      y_range = y_max - y_min,
      .groups = "drop"
    ) %>%
    mutate(
      x_range = ifelse(
        x_range > 0,
        x_range,
        1
      ),
      y_range = ifelse(
        y_range > 0,
        y_range,
        1
      ),
      x_label = x_min + 0.03 * x_range,
      y_label = y_max - 0.03 * y_range
    )
  
  stats <- stats %>%
    left_join(
      label_positions %>%
        select(
          lncRNA,
          mRNA,
          x_label,
          y_label
        ),
      by = c(
        "lncRNA",
        "mRNA"
      )
    )
  
  # --------------------------------------------------------
  # 22.9.6 Generate 3 x 9 plot
  # --------------------------------------------------------
  
  p <- ggplot(
    plot_long,
    aes(
      x = x,
      y = y
    )
  ) +
    
    geom_point(
      size = 1.8,
      alpha = 0.8
    ) +
    
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.6
    ) +
    
    geom_text(
      data = stats,
      aes(
        x = x_label,
        y = y_label,
        label = Label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 1.8,
      lineheight = 0.95,
      show.legend = FALSE
    ) +
    
    facet_grid(
      rows = vars(lncRNA),
      cols = vars(mRNA),
      scales = "free",
      space = "fixed"
    ) +
    
    labs(
      title = plot_title,
      subtitle = paste0(
        "All 27 lncRNA-mRNA pairs | ",
        "Expression = log10(TPM + 1)"
      ),
      x = expression(
        log[10]("lncRNA TPM + 1")
      ),
      y = expression(
        log[10]("mRNA TPM + 1")
      )
    ) +
    
    theme_bw(
      base_size = 9
    ) +
    
    theme(
      plot.title = element_text(
        size = 13,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = 9,
        hjust = 0.5
      ),
      
      strip.text.x = element_text(
        size = 7,
        face = "bold"
      ),
      
      strip.text.y = element_text(
        size = 8,
        face = "bold",
        angle = 270,
        hjust = 0.5,
        vjust = 0.5
      ),,
      
      strip.background = element_rect(
        fill = "grey90",
        colour = "grey50"
      ),
      
      axis.title = element_text(
        size = 8
      ),
      
      axis.text = element_text(
        size = 5
      ),
      
      panel.spacing = grid::unit(
        0.12,
        "cm"
      ),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(
        5, 5, 5, 5
      )
    )
  
  # --------------------------------------------------------
  # 22.9.7 Move lncRNA row strips to the far-right side
  # --------------------------------------------------------
  
  gt <- ggplot2::ggplotGrob(p)
  
  # Identify the row-strip grobs created by facet_grid().
  strip_rows <- grep(
    "^strip-r",
    gt$layout$name
  )
  
  # Identify the panel columns.
  panel_columns <- gt$layout$l[
    grepl(
      "^panel",
      gt$layout$name
    )
  ]
  
  if (
    length(strip_rows) > 0 &&
    length(panel_columns) > 0
  ) {
    
    # Move the existing row strips after the rightmost panel.
    rightmost_panel_column <- max(panel_columns)
    
    for (i in strip_rows) {
      
      gt$layout$r[i] <- rightmost_panel_column + 1
      gt$layout$l[i] <- rightmost_panel_column + 1
      
    }
    
    # Add a column for the relocated strips.
    gt <- gtable::gtable_add_cols(
      gt,
      widths = grid::unit(
        1.5,
        "cm"
      ),
      pos = rightmost_panel_column
    )
  }
  
  # --------------------------------------------------------
  # 22.9.8 Save high-resolution plot
  # --------------------------------------------------------
  
  output_file <- file.path(
    plot_dir,
    filename
  )
  
  # Save the modified gtable rather than the original ggplot.
  grDevices::png(
    filename = output_file,
    width = 24,
    height = 11,
    units = "in",
    res = 300
  )
  
  grid::grid.newpage()
  grid::grid.draw(gt)
  
  grDevices::dev.off()
  
  return(p)
}
############################################################
# 22.10 PLOT 4: ALL NORMAL SAMPLES
############################################################

plot_normal <- make_single_group_plot(
  data = plot_expression,
  samples = normal_samples,
  group_label = "Normal",
  filename = "Normal_All27.png",
  plot_title = "Normal Sample lncRNA-mRNA Correlations"
)

print(plot_normal)


############################################################
# 22.11 PLOT 5: NAIVE B CELLS
############################################################

plot_naive <- make_single_group_plot(
  data = plot_expression,
  samples = naive_samples,
  group_label = "Naive B cells",
  filename = "Naive_B_Cells_All27.png",
  plot_title = "Naive B-Cell lncRNA-mRNA Correlations"
)

print(plot_naive)


############################################################
# 22.12 PLOT 6: GCB NORMAL B CELLS
############################################################

plot_gcb <- make_single_group_plot(
  data = plot_expression,
  samples = gcb_normal_samples,
  group_label = "GCB B cells",
  filename = "GCB_B_Cells_All27.png",
  plot_title = "GCB B-Cell lncRNA-mRNA Correlations"
)

print(plot_gcb)


############################################################
# 22.13 PLOT 7: SRA TUMOURS ONLY
############################################################

plot_sra <- make_single_group_plot(
  data = plot_expression,
  samples = sra_tumour_samples,
  group_label = "SRA tumours",
  filename = "SRA_Tumours_All27.png",
  plot_title = "SRA Tumour lncRNA-mRNA Correlations"
)

print(plot_sra)


############################################################
# 22.14 PLOT 8: TCGA TUMOURS ONLY
############################################################

plot_tcga <- make_single_group_plot(
  data = plot_expression,
  samples = tcga_tumour_samples,
  group_label = "TCGA tumours",
  filename = "TCGA_Tumours_All27.png",
  plot_title = "TCGA Tumour lncRNA-mRNA Correlations"
)

print(plot_tcga)

############################################################
# 22.15 FINISH
############################################################

message(
  "All correlation plots and corresponding statistics ",
  "have been saved in: ",
  normalizePath(plot_dir)
)

############################################################
# 23. CREATE SUPER MASTER DLBCL CORRELATION RESULTS
#
# Combines:
# 1. Normal correlations
# 2. Naive B-cell correlations
# 3. GCB normal B-cell correlations
# 4. SRA tumour correlations
# 5. TCGA tumour correlations
# 6. Adjusted partial Spearman correlations
# 7. Normal vs SRA differential correlations
# 8. Normal vs TCGA differential correlations
# 9. SRA vs TCGA differential correlations
# 10. Normal vs SRA adjusted interaction tests
# 11. Normal vs TCGA adjusted interaction tests
# 12. SRA vs TCGA adjusted interaction tests
# 13. Annotated cell/tumour subset correlations
#
# All 27 candidate pairs are retained.
############################################################


############################################################
# 23.1 OUTPUT FILE
############################################################

super_master_file <- file.path(
  output_dir,
  "SUPER_MASTER_DLBCL_Correlation_Results.csv"
)


############################################################
# 23.2 START WITH ALL 27 CANDIDATE PAIRS
############################################################

super_master <- candidate_pairs %>%
  distinct(lncRNA, mRNA)


############################################################
# 23.3 ADD MAIN GROUP CORRELATIONS
############################################################

# Normal samples
super_master <- super_master %>%
  left_join(
    normal_results %>%
      select(
        lncRNA,
        mRNA,
        Normal_n = n,
        Normal_rho = rho,
        Normal_p = p_value,
        Normal_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


# Naive B cells
super_master <- super_master %>%
  left_join(
    naive_results %>%
      select(
        lncRNA,
        mRNA,
        Naive_n = n,
        Naive_rho = rho,
        Naive_p = p_value,
        Naive_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


# GCB normal B cells
super_master <- super_master %>%
  left_join(
    gcb_normal_results %>%
      select(
        lncRNA,
        mRNA,
        GCB_Normal_n = n,
        GCB_Normal_rho = rho,
        GCB_Normal_p = p_value,
        GCB_Normal_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


# SRA tumours
super_master <- super_master %>%
  left_join(
    sra_results %>%
      select(
        lncRNA,
        mRNA,
        SRA_n = n,
        SRA_rho = rho,
        SRA_p = p_value,
        SRA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


# TCGA tumours
super_master <- super_master %>%
  left_join(
    tcga_results %>%
      select(
        lncRNA,
        mRNA,
        TCGA_n = n,
        TCGA_rho = rho,
        TCGA_p = p_value,
        TCGA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


############################################################
# 23.4 ADD ADJUSTED PARTIAL SPEARMAN CORRELATIONS
############################################################

super_master <- super_master %>%
  left_join(
    adjusted_partial_results %>%
      select(
        lncRNA,
        mRNA,
        Adjusted_Tumour_n = n,
        Adjusted_Partial_Spearman_rho = partial_rho,
        Adjusted_Partial_p = p_value,
        Adjusted_Partial_FDR = FDR_BH,
        Adjusted_Partial_df = df
      ),
    by = c("lncRNA", "mRNA")
  )


############################################################
# 23.5 ADD RAW DIFFERENTIAL-CORRELATION PERMUTATION TESTS
############################################################

# Normal vs SRA tumours

super_master <- super_master %>%
  left_join(
    normal_vs_sra %>%
      select(
        lncRNA,
        mRNA,
        Normal_vs_SRA_n1 = n1,
        Normal_vs_SRA_n2 = n2,
        Normal_vs_SRA_rho1 = rho1,
        Normal_vs_SRA_rho2 = rho2,
        Normal_vs_SRA_DeltaRho = delta_rho,
        Normal_vs_SRA_Permutation_p = p_permutation,
        Normal_vs_SRA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


# Normal vs TCGA tumours

super_master <- super_master %>%
  left_join(
    normal_vs_tcga %>%
      select(
        lncRNA,
        mRNA,
        Normal_vs_TCGA_n1 = n1,
        Normal_vs_TCGA_n2 = n2,
        Normal_vs_TCGA_rho1 = rho1,
        Normal_vs_TCGA_rho2 = rho2,
        Normal_vs_TCGA_DeltaRho = delta_rho,
        Normal_vs_TCGA_Permutation_p = p_permutation,
        Normal_vs_TCGA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


# SRA vs TCGA tumours

super_master <- super_master %>%
  left_join(
    sra_vs_tcga %>%
      select(
        lncRNA,
        mRNA,
        SRA_vs_TCGA_n1 = n1,
        SRA_vs_TCGA_n2 = n2,
        SRA_vs_TCGA_rho1 = rho1,
        SRA_vs_TCGA_rho2 = rho2,
        SRA_vs_TCGA_DeltaRho = delta_rho,
        SRA_vs_TCGA_Permutation_p = p_permutation,
        SRA_vs_TCGA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


############################################################
# 23.6 ADD ADJUSTED INTERACTION PERMUTATION TESTS
#
# These test differences in rank-regression slopes.
# They are not differential-Spearman tests.
############################################################

# Normal vs SRA: adjusted for StromalScore

super_master <- super_master %>%
  left_join(
    adjusted_normal_sra %>%
      select(
        lncRNA,
        mRNA,
        Adjusted_Normal_vs_SRA_n = n,
        Adjusted_Normal_vs_SRA_Interaction_Beta =
          interaction_beta,
        Adjusted_Normal_vs_SRA_Interaction_t =
          interaction_t,
        Adjusted_Normal_vs_SRA_Permutation_p =
          p_permutation,
        Adjusted_Normal_vs_SRA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


# Normal vs TCGA: adjusted for StromalScore only

super_master <- super_master %>%
  left_join(
    adjusted_normal_tcga %>%
      select(
        lncRNA,
        mRNA,
        Adjusted_Normal_vs_TCGA_n = n,
        Adjusted_Normal_vs_TCGA_Interaction_Beta =
          interaction_beta,
        Adjusted_Normal_vs_TCGA_Interaction_t =
          interaction_t,
        Adjusted_Normal_vs_TCGA_Permutation_p =
          p_permutation,
        Adjusted_Normal_vs_TCGA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


# SRA vs TCGA: adjusted for StromalScore

super_master <- super_master %>%
  left_join(
    adjusted_sra_tcga %>%
      select(
        lncRNA,
        mRNA,
        Adjusted_SRA_vs_TCGA_n = n,
        Adjusted_SRA_vs_TCGA_Interaction_Beta =
          interaction_beta,
        Adjusted_SRA_vs_TCGA_Interaction_t =
          interaction_t,
        Adjusted_SRA_vs_TCGA_Permutation_p =
          p_permutation,
        Adjusted_SRA_vs_TCGA_FDR = FDR_BH
      ),
    by = c("lncRNA", "mRNA")
  )


############################################################
# 23.7 ADD ANNOTATED SUBSET CORRELATIONS
#
# Keep subset-specific results in a separate long-format
# table because multiple subsets can exist for each pair.
############################################################

subset_master <- subset_results %>%
  mutate(
    Significant_FDR_0.05 = ifelse(
      !is.na(FDR_BH) & FDR_BH < 0.05,
      "Yes",
      "No"
    )
  ) %>%
  arrange(
    Group,
    lncRNA,
    mRNA
  )


write.csv(
  subset_master,
  file.path(
    output_dir,
    "SUPER_MASTER_Subset_Correlations.csv"
  ),
  row.names = FALSE,
  na = ""
)


############################################################
# 23.8 ADD SIGNIFICANCE FLAGS TO THE SUPER MASTER
############################################################

super_master <- super_master %>%
  mutate(
    Normal_Significant = ifelse(
      !is.na(Normal_FDR) & Normal_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    Naive_Significant = ifelse(
      !is.na(Naive_FDR) & Naive_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    GCB_Normal_Significant = ifelse(
      !is.na(GCB_Normal_FDR) & GCB_Normal_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    SRA_Significant = ifelse(
      !is.na(SRA_FDR) & SRA_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    TCGA_Significant = ifelse(
      !is.na(TCGA_FDR) & TCGA_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    Adjusted_Partial_Significant = ifelse(
      !is.na(Adjusted_Partial_FDR) &
        Adjusted_Partial_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    Normal_vs_SRA_Differential_Significant = ifelse(
      !is.na(Normal_vs_SRA_FDR) &
        Normal_vs_SRA_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    Normal_vs_TCGA_Differential_Significant = ifelse(
      !is.na(Normal_vs_TCGA_FDR) &
        Normal_vs_TCGA_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    SRA_vs_TCGA_Differential_Significant = ifelse(
      !is.na(SRA_vs_TCGA_FDR) &
        SRA_vs_TCGA_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    Adjusted_Normal_vs_SRA_Significant = ifelse(
      !is.na(Adjusted_Normal_vs_SRA_FDR) &
        Adjusted_Normal_vs_SRA_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    Adjusted_Normal_vs_TCGA_Significant = ifelse(
      !is.na(Adjusted_Normal_vs_TCGA_FDR) &
        Adjusted_Normal_vs_TCGA_FDR < 0.05,
      "Yes",
      "No"
    ),
    
    Adjusted_SRA_vs_TCGA_Significant = ifelse(
      !is.na(Adjusted_SRA_vs_TCGA_FDR) &
        Adjusted_SRA_vs_TCGA_FDR < 0.05,
      "Yes",
      "No"
    )
  )


############################################################
# 23.9 ADD INTERPRETATION FLAGS
############################################################

super_master <- super_master %>%
  mutate(
    Normal_and_SRA_Both_Significant =
      ifelse(
        Normal_Significant == "Yes" &
          SRA_Significant == "Yes",
        "Yes",
        "No"
      ),
    
    Normal_and_TCGA_Both_Significant =
      ifelse(
        Normal_Significant == "Yes" &
          TCGA_Significant == "Yes",
        "Yes",
        "No"
      ),
    
    SRA_and_TCGA_Both_Significant =
      ifelse(
        SRA_Significant == "Yes" &
          TCGA_Significant == "Yes",
        "Yes",
        "No"
      ),
    
    Normal_vs_SRA_Direction_Reversal =
      ifelse(
        !is.na(Normal_rho) &
          !is.na(SRA_rho) &
          Normal_rho * SRA_rho < 0,
        "Yes",
        "No"
      ),
    
    Normal_vs_TCGA_Direction_Reversal =
      ifelse(
        !is.na(Normal_rho) &
          !is.na(TCGA_rho) &
          Normal_rho * TCGA_rho < 0,
        "Yes",
        "No"
      ),
    
    SRA_vs_TCGA_Direction_Reversal =
      ifelse(
        !is.na(SRA_rho) &
          !is.na(TCGA_rho) &
          SRA_rho * TCGA_rho < 0,
        "Yes",
        "No"
      )
  )


############################################################
# 23.10 EXPORT SUPER MASTER CSV
############################################################

write.csv(
  super_master,
  super_master_file,
  row.names = FALSE,
  na = ""
)


############################################################
# 23.11 EXPORT A SUMMARY OF ALL ANALYSES
############################################################

super_master_summary <- tibble(
  Analysis = c(
    "Normal correlations",
    "Naive B-cell correlations",
    "GCB normal B-cell correlations",
    "SRA tumour correlations",
    "TCGA tumour correlations",
    "Adjusted partial Spearman correlations",
    "Normal vs SRA differential correlations",
    "Normal vs TCGA differential correlations",
    "SRA vs TCGA differential correlations",
    "Adjusted Normal vs SRA interactions",
    "Adjusted Normal vs TCGA interactions",
    "Adjusted SRA vs TCGA interactions"
  ),
  
  Number_of_pairs = 27,
  
  Number_FDR_below_005 = c(
    sum(normal_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(naive_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(gcb_normal_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(sra_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(tcga_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(adjusted_partial_results$FDR_BH < 0.05, na.rm = TRUE),
    sum(normal_vs_sra$FDR_BH < 0.05, na.rm = TRUE),
    sum(normal_vs_tcga$FDR_BH < 0.05, na.rm = TRUE),
    sum(sra_vs_tcga$FDR_BH < 0.05, na.rm = TRUE),
    sum(adjusted_normal_sra$FDR_BH < 0.05, na.rm = TRUE),
    sum(adjusted_normal_tcga$FDR_BH < 0.05, na.rm = TRUE),
    sum(adjusted_sra_tcga$FDR_BH < 0.05, na.rm = TRUE)
  )
)


write.csv(
  super_master_summary,
  file.path(
    output_dir,
    "SUPER_MASTER_Analysis_Summary.csv"
  ),
  row.names = FALSE
)


############################################################
# 23.12 VERIFY OUTPUT
############################################################

stopifnot(
  nrow(super_master) == 27
)

message(
  "\nSUPER MASTER ANALYSIS COMPLETE\n",
  "Main results: ",
  normalizePath(super_master_file),
  "\nSubset results: ",
  normalizePath(
    file.path(
      output_dir,
      "SUPER_MASTER_Subset_Correlations.csv"
    )
  ),
  "\nSummary: ",
  normalizePath(
    file.path(
      output_dir,
      "SUPER_MASTER_Analysis_Summary.csv"
    )
  ),
  "\nNumber of candidate pairs: ",
  nrow(super_master),
  "\nNumber of columns in super master: ",
  ncol(super_master)
)

print(super_master_summary)


# ============================================================
# SUPER-MASTER DLBCL CORRELATION RESULTS SUMMARY
# ============================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# ------------------------------------------------------------
# 1. INPUT AND OUTPUT
# ------------------------------------------------------------

input_file <- "DLBCL_Correlation_Results/SUPER_MASTER_DLBCL_Correlation_Results.csv"

output_dir <- "SUPER_MASTER_Results_Summary"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

master <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Rows:", nrow(master), "\n")
cat("Columns:", ncol(master), "\n")

# Inspect column names
print(names(master))

# ------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------

# Find a column using one of several possible names
find_col <- function(data, candidates, required = TRUE) {
  
  matched <- candidates[candidates %in% names(data)]
  
  if (length(matched) > 0) {
    return(matched[1])
  }
  
  if (required) {
    stop(
      "Could not find any of these columns: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  return(NA_character_)
}


# Convert a column to numeric safely
as_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}


# Extract a significant-pair table
get_significant_pairs <- function(data, fdr_col, rho_col = NULL) {
  
  if (is.na(fdr_col) || !fdr_col %in% names(data)) {
    return(data.frame())
  }
  
  result <- data %>%
    filter(
      !is.na(.data[[fdr_col]]),
      as_numeric_safe(.data[[fdr_col]]) < 0.05
    )
  
  if (!is.null(rho_col) && rho_col %in% names(result)) {
    result <- result %>%
      arrange(as_numeric_safe(.data[[fdr_col]]))
  } else {
    result <- result %>%
      arrange(as_numeric_safe(.data[[fdr_col]]))
  }
  
  result
}


# Summarise the number of significant results
summarise_fdr <- function(data, category, fdr_col) {
  
  if (is.na(fdr_col) || !fdr_col %in% names(data)) {
    
    return(data.frame(
      Category = category,
      Total_pairs = nrow(data),
      Significant_FDR_lt_0.05 = NA_integer_,
      Significant_percentage = NA_real_
    ))
  }
  
  fdr_values <- as_numeric_safe(data[[fdr_col]])
  
  data.frame(
    Category = category,
    Total_pairs = length(fdr_values),
    Significant_FDR_lt_0.05 = sum(
      !is.na(fdr_values) & fdr_values < 0.05
    ),
    Significant_percentage = round(
      100 * sum(!is.na(fdr_values) & fdr_values < 0.05) /
        max(1, sum(!is.na(fdr_values))),
      2
    )
  )
}


# ------------------------------------------------------------
# 3. IDENTIFY THE MAIN PAIR COLUMNS
# ------------------------------------------------------------

lnc_col <- find_col(
  master,
  c("lncRNA", "lncRNASymbol", "lncRNA_Symbol")
)

mrna_col <- find_col(
  master,
  c("mRNA", "Target_Gene", "mRNA_Symbol")
)

master <- master %>%
  mutate(
    lncRNA = as.character(.data[[lnc_col]]),
    mRNA = as.character(.data[[mrna_col]])
  )

# ------------------------------------------------------------
# 4. IDENTIFY CORRELATION AND FDR COLUMNS
# ------------------------------------------------------------

# Group-specific Spearman correlations
normal_rho <- find_col(
  master,
  c("Normal_rho", "Normal_Spearman_rho", "rho_Normal"),
  required = FALSE
)

normal_fdr <- find_col(
  master,
  c("Normal_FDR", "Normal_Spearman_FDR", "FDR_Normal"),
  required = FALSE
)

naive_rho <- find_col(
  master,
  c("Naive_rho", "Naive_Spearman_rho", "rho_Naive"),
  required = FALSE
)

naive_fdr <- find_col(
  master,
  c("Naive_FDR", "Naive_Spearman_FDR", "FDR_Naive"),
  required = FALSE
)

gcb_rho <- find_col(
  master,
  c("GCB_Normal_rho", "GCB_rho", "GCB_Normal_Spearman_rho"),
  required = FALSE
)

gcb_fdr <- find_col(
  master,
  c("GCB_Normal_FDR", "GCB_FDR", "GCB_Normal_Spearman_FDR"),
  required = FALSE
)

sra_rho <- find_col(
  master,
  c("SRA_Tumour_rho", "SRA_rho", "SRA_Tumour_Spearman_rho"),
  required = FALSE
)

sra_fdr <- find_col(
  master,
  c("SRA_Tumour_FDR", "SRA_FDR", "SRA_Tumour_Spearman_FDR"),
  required = FALSE
)

tcga_rho <- find_col(
  master,
  c("TCGA_Tumour_rho", "TCGA_rho", "TCGA_Tumour_Spearman_rho"),
  required = FALSE
)

tcga_fdr <- find_col(
  master,
  c("TCGA_Tumour_FDR", "TCGA_FDR", "TCGA_Tumour_Spearman_FDR"),
  required = FALSE
)

partial_rho <- find_col(
  master,
  c("Adjusted_Partial_rho", "Partial_rho", "Partial_Spearman_rho"),
  required = FALSE
)

partial_fdr <- find_col(
  master,
  c("Adjusted_Partial_FDR", "Partial_FDR", "Partial_Spearman_FDR"),
  required = FALSE
)

# ------------------------------------------------------------
# 5. GROUP-WISE CORRELATION SUMMARY
# ------------------------------------------------------------

group_summary <- bind_rows(
  
  summarise_fdr(master, "Normal", normal_fdr),
  
  summarise_fdr(master, "Naive", naive_fdr),
  
  summarise_fdr(master, "GCB normal", gcb_fdr),
  
  summarise_fdr(master, "SRA tumour", sra_fdr),
  
  summarise_fdr(master, "TCGA tumour", tcga_fdr),
  
  summarise_fdr(master, "Adjusted partial correlation", partial_fdr)
  
)

write.csv(
  group_summary,
  file.path(output_dir, "01_Groupwise_Correlation_Summary.csv"),
  row.names = FALSE
)

print(group_summary)

# ------------------------------------------------------------
# 6. EXPORT SIGNIFICANT GROUP-WISE CORRELATIONS
# ------------------------------------------------------------

export_significant <- function(
    data,
    fdr_col,
    filename,
    rho_col = NULL
) {
  
  if (is.na(fdr_col) || !fdr_col %in% names(data)) {
    warning("Skipping ", filename, ": FDR column was not found.")
    return(invisible(NULL))
  }
  
  significant <- get_significant_pairs(
    data = data,
    fdr_col = fdr_col,
    rho_col = rho_col
  )
  
  write.csv(
    significant,
    file.path(output_dir, filename),
    row.names = FALSE
  )
  
  cat(
    filename,
    ":",
    nrow(significant),
    "significant pairs\n"
  )
  
  invisible(significant)
}


export_significant(
  master, normal_fdr,
  "02_Significant_Normal_Correlations.csv",
  normal_rho
)

export_significant(
  master, naive_fdr,
  "03_Significant_Naive_Correlations.csv",
  naive_rho
)

export_significant(
  master, gcb_fdr,
  "04_Significant_GCB_Normal_Correlations.csv",
  gcb_rho
)

export_significant(
  master, sra_fdr,
  "05_Significant_SRA_Tumour_Correlations.csv",
  sra_rho
)

export_significant(
  master, tcga_fdr,
  "06_Significant_TCGA_Tumour_Correlations.csv",
  tcga_rho
)

export_significant(
  master, partial_fdr,
  "07_Significant_Adjusted_Partial_Correlations.csv",
  partial_rho
)

# ------------------------------------------------------------
# 7. DIFFERENTIAL CORRELATION ANALYSIS
# ------------------------------------------------------------

# Raw differential-correlation columns
normal_sra_fdr <- find_col(
  master,
  c("Normal_vs_SRA_FDR", "Normal_SRA_FDR",
    "FDR_Normal_vs_SRA"),
  required = FALSE
)

normal_tcga_fdr <- find_col(
  master,
  c("Normal_vs_TCGA_FDR", "Normal_TCGA_FDR",
    "FDR_Normal_vs_TCGA"),
  required = FALSE
)

sra_tcga_fdr <- find_col(
  master,
  c("SRA_vs_TCGA_FDR", "SRA_TCGA_FDR",
    "FDR_SRA_vs_TCGA"),
  required = FALSE
)

differential_summary <- bind_rows(
  
  summarise_fdr(
    master,
    "Normal vs SRA tumour",
    normal_sra_fdr
  ),
  
  summarise_fdr(
    master,
    "Normal vs TCGA tumour",
    normal_tcga_fdr
  ),
  
  summarise_fdr(
    master,
    "SRA vs TCGA tumour",
    sra_tcga_fdr
  )
  
)

write.csv(
  differential_summary,
  file.path(output_dir, "08_Differential_Correlation_Summary.csv"),
  row.names = FALSE
)

print(differential_summary)

# Export significant differential-correlation pairs
export_significant(
  master,
  normal_sra_fdr,
  "09_Significant_Differential_Normal_vs_SRA.csv"
)

export_significant(
  master,
  normal_tcga_fdr,
  "10_Significant_Differential_Normal_vs_TCGA.csv"
)

export_significant(
  master,
  sra_tcga_fdr,
  "11_Significant_Differential_SRA_vs_TCGA.csv"
)

# ------------------------------------------------------------
# 8. ADJUSTED INTERACTION-TEST SUMMARY
# ------------------------------------------------------------

# Use the exact column names present in the CSV

interaction_normal_sra_fdr <-
  "Adjusted_Normal_vs_SRA_FDR"

interaction_normal_tcga_fdr <-
  "Adjusted_Normal_vs_TCGA_FDR"

interaction_sra_tcga_fdr <-
  "Adjusted_SRA_vs_TCGA_FDR"


# ------------------------------------------------------------
# Generate the adjusted interaction-test summary
# ------------------------------------------------------------

interaction_summary <- bind_rows(
  
  summarise_fdr(
    master,
    "Adjusted interaction: Normal vs SRA",
    interaction_normal_sra_fdr
  ),
  
  summarise_fdr(
    master,
    "Adjusted interaction: Normal vs TCGA",
    interaction_normal_tcga_fdr
  ),
  
  summarise_fdr(
    master,
    "Adjusted interaction: SRA vs TCGA",
    interaction_sra_tcga_fdr
  )
  
)

# Export the summary table

write.csv(
  interaction_summary,
  file.path(
    output_dir,
    "12_Adjusted_Interaction_Summary.csv"
  ),
  row.names = FALSE
)

# Display the summary

print(interaction_summary)


# ------------------------------------------------------------
# Export significant adjusted interaction results
# ------------------------------------------------------------

export_significant(
  master,
  interaction_normal_sra_fdr,
  "13_Significant_Interaction_Normal_vs_SRA.csv"
)

export_significant(
  master,
  interaction_normal_tcga_fdr,
  "14_Significant_Interaction_Normal_vs_TCGA.csv"
)

export_significant(
  master,
  interaction_sra_tcga_fdr,
  "15_Significant_Interaction_SRA_vs_TCGA.csv"
)

# ------------------------------------------------------------
# 9. CONSOLIDATED MASTER SUMMARY
# ------------------------------------------------------------

consolidated_summary <- bind_rows(
  group_summary,
  differential_summary,
  interaction_summary
)

write.csv(
  consolidated_summary,
  file.path(output_dir, "16_Consolidated_Statistical_Summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. GENERATE A TEXT SUMMARY FOR THE RESULTS SECTION
# ------------------------------------------------------------

report_file <- file.path(
  output_dir,
  "17_Results_Section_Statistical_Summary.txt"
)

sink(report_file)

cat("DLBCL lncRNA–mRNA CORRELATION ANALYSIS\n")
cat("======================================\n\n")

cat("Total lncRNA–mRNA pairs analysed:", nrow(master), "\n\n")

cat("1. GROUP-WISE CORRELATIONS\n")
cat("--------------------------\n")

for (i in seq_len(nrow(group_summary))) {
  
  cat(
    group_summary$Category[i],
    ": ",
    group_summary$Significant_FDR_lt_0.05[i],
    " of ",
    group_summary$Total_pairs[i],
    " pairs had FDR < 0.05 (",
    group_summary$Significant_percentage[i],
    "%).\n",
    sep = ""
  )
}

cat("\n2. DIFFERENTIAL CORRELATIONS\n")
cat("---------------------------\n")

for (i in seq_len(nrow(differential_summary))) {
  
  cat(
    differential_summary$Category[i],
    ": ",
    differential_summary$Significant_FDR_lt_0.05[i],
    " of ",
    differential_summary$Total_pairs[i],
    " pairs showed a significant difference (FDR < 0.05; ",
    differential_summary$Significant_percentage[i],
    "%).\n",
    sep = ""
  )
}

cat("\n3. ADJUSTED INTERACTION TESTS\n")
cat("----------------------------\n")

for (i in seq_len(nrow(interaction_summary))) {
  
  cat(
    interaction_summary$Category[i],
    ": ",
    interaction_summary$Significant_FDR_lt_0.05[i],
    " of ",
    interaction_summary$Total_pairs[i],
    " pairs were significant (FDR < 0.05).\n",
    sep = ""
  )
}

cat("\n4. SIGNIFICANT GROUP-WISE CORRELATIONS\n")
cat("-------------------------------------\n")

for (i in seq_len(nrow(group_summary))) {
  
  category <- group_summary$Category[i]
  
  cat(
    "\n",
    category,
    ": ",
    group_summary$Significant_FDR_lt_0.05[i],
    " significant pairs.\n",
    sep = ""
  )
}

cat("\nEnd of summary.\n")

sink()

# ------------------------------------------------------------
# 11. FINAL MESSAGE
# ------------------------------------------------------------

cat("\n============================================\n")
cat("ANALYSIS SUMMARY COMPLETE\n")
cat("============================================\n")
cat("Output directory:", output_dir, "\n")
cat("Summary CSV:", "16_Consolidated_Statistical_Summary.csv", "\n")
cat("Results text:", "17_Results_Section_Statistical_Summary.txt", "\n")

