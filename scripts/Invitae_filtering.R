# VP: Load required libraries 
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(purrr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

EXCEL_PATH  <- "Invitae_Personal.xlsx"
SHEET_NAMES <- c(
  possibly_mosaic = "Possibly Mosaic",
  germline_plp   = "Germline_LP_P",
  germline_vus   = "Germline_VUS"
)

# VP: Toggle for counting unique participants (vs raw rows) in summaries
COUNT_UNIQUE_PARTICIPANTS <- TRUE

`%||%` <- function(a, b) if (!is.null(a)) a else b
trim_upper <- function(x) { if (is.null(x)) return(NA_character_); ifelse(is.na(x), NA, toupper(str_trim(as.character(x)))) }
as_bool <- function(x) { if (is.logical(x)) return(x); y <- tolower(as.character(x)); y %in% c("true","t","1","yes","y","affected","case","1.0") }

make_variant_key <- function(df) {
  has <- function(nm) nm %in% names(df)
  coalesce_chr <- function(...) { x <- list(...); out <- x[[1]]; for (i in 2:length(x)) out <- dplyr::coalesce(out, x[[i]]); out }
  coalesce_chr(
    if (has("variant_id")) df$variant_id else NA_character_,
    if (has("symbol") & has("hgvsc")) paste(df$symbol, df$hgvsc, sep="|") else NA_character_,
    if (has("symbol") & has("hgvsp")) paste(df$symbol, df$hgvsp, sep="|") else NA_character_,
    if (has("symbol") & has("REPORTED_TRANSCRIPT") & has("hgvsc"))
      paste(df$symbol, df$REPORTED_TRANSCRIPT, df$hgvsc, sep="|") else NA_character_
  )
}
manual_map <- c(
  reported_gene_name      = "REPORTED_GENE_NAME",
  variant_impact          = "VARIANT_IMPACT",
  is_individual_affected  = "IS_INDIVIDUAL_AFFECTED",
  dummy_patient_id        = "dummy_patient_id",
  dummy_report_id         = "dummy_report_id",
  sex_assigned_at_birth   = "SEX_ASSIGNED_AT_BIRTH",
  specimen_type           = "SPECIMEN_TYPE",
  reported_interpretation = "REPORTED_INTERPRETATION",
  inheritance             = "INHERITANCE",
  symbol                  = "REPORTED_GENE_NAME",
  hgvsc                   = "CNAME",
  hgvsp                   = "PNAME"
)

read_and_harmonise <- function(path, sheet) {
  raw <- readxl::read_excel(path, sheet = sheet) |> janitor::remove_empty("rows")
  if (nrow(raw) == 0) stop(paste0("Sheet '", sheet, "' is empty."))
  df <- raw

  # VP: Create columns (fill with NA if source column not present)
  for (canon in names(manual_map)) {
    src <- manual_map[[canon]]
    df[[canon]] <- if (!is.na(src) && src %in% names(df)) df[[src]] else NA
  }

  if ("REPORTED_TRANSCRIPT" %in% names(df)) df$REPORTED_TRANSCRIPT <- as.character(df$REPORTED_TRANSCRIPT)
  has_raw_eth <- "RAW_PATIENT_ETHNICITY" %in% names(df)
  has_can_eth <- "PATIENT_ETHNICITY"     %in% names(df)

  df <- df |>
    mutate(
      reported_gene_name      = trim_upper(reported_gene_name),
      symbol                  = trim_upper(symbol),
      variant_impact          = trim_upper(variant_impact),
      is_individual_affected  = as_bool(IS_INDIVIDUAL_AFFECTED %||% is_individual_affected),
      sex_assigned_at_birtH   = trim_upper(sex_assigned_at_birth %||% SEX_ASSIGNED_AT_BIRTH),
      sex_assigned_at_birth   = ifelse(!is.na(sex_assigned_at_birtH), sex_assigned_at_birtH, trim_upper(sex_assigned_at_birth)),
      specimen_type           = as.character(specimen_type %||% SPECIMEN_TYPE),
      reported_interpretation = trim_upper(reported_interpretation),
      inheritance             = trim_upper(inheritance),
      hgvsc                   = as.character(hgvsc),
      hgvsp                   = as.character(hgvsp)
    )

  # VP: Ethnicity columns copied (if present) to standardized fields used by summaries
  df$ETH_RAW   <- if (has_raw_eth) as.character(df$RAW_PATIENT_ETHNICITY) else NA_character_
  df$ETH_CANON <- if (has_can_eth) as.character(df$PATIENT_ETHNICITY)     else NA_character_

  df$variant_key <- make_variant_key(df)
  df
}

is_high_or_moderate <- function(x) {
  x <- toupper(trimws(as.character(x)))
  keep <- c(
    "HIGH","MODERATE",
    "MISSENSE",
    "FRAMESHIFT",
    "STOP GAINED",
    "SPLICE DONOR","SPLICE ACCEPTOR","SPLICE SITE",
    "IN-FRAME CODON LOSS","IN-FRAME CODON GAIN",
    "STOP LOST","INITIATOR CODON"
  )
  x %in% keep
}

apply_filters <- function(df, group_name, is_vus_group = FALSE, is_possibly_mosaic = FALSE) {
  # VP: Impact filter
  df <- df |> filter(is_high_or_moderate(variant_impact))

  # VP: Keep only affected for non-VUS groups
  if (!is_vus_group) df <- df |> filter(is_individual_affected %in% TRUE)

  # VP: Create dedup_id early; use patient ID if present, else report ID
  df <- df |> mutate(dedup_id = dplyr::coalesce(dummy_patient_id, dummy_report_id))

  # VP: Deduplicate per participant×variant (fallback key if variant_key is missing)
  df <- df |>
    mutate(
      fallback_key = ifelse(is.na(variant_key) | variant_key=="",
                            paste0(symbol %||% reported_gene_name, "|", hgvsc),
                            variant_key),
      dedup_key = paste0(dedup_id, "||", fallback_key)
    ) |>
    arrange(dedup_id, symbol, hgvsc, hgvsp) |>
    distinct(dedup_key, .keep_all = TRUE)

  # VP: X-linked filter — keep only males for variants in X-linked genes
  df <- df |>
    mutate(is_x = ifelse(!is.na(inheritance),
                         str_detect(inheritance, regex("X[- ]?LINK", ignore_case = TRUE)),
                         FALSE),
           sex_clean = trim_upper(sex_assigned_at_birth)) |>
    filter(!(is_x & !sex_clean %in% c("MALE","M"))) |>
    select(-sex_clean, -is_x)

  # VP: Possibly Mosaic 
  if (is_possibly_mosaic) {
    df <- df |> filter(reported_interpretation %in% c("PATHOGENIC","LIKELY PATHOGENIC"))
  }

  df <- df |> mutate(group = group_name)
  df
}

# VP: Read Excel, and then apply filters per group
message("\nReading Excel …")
df_poss <- read_and_harmonise(EXCEL_PATH, SHEET_NAMES["possibly_mosaic"])
df_plp  <- read_and_harmonise(EXCEL_PATH, SHEET_NAMES["germline_plp"])
df_vus  <- read_and_harmonise(EXCEL_PATH, SHEET_NAMES["germline_vus"])

message("Applying filters …")
clean_poss <- apply_filters(df_poss, "PossiblyMosaic", is_vus_group = FALSE, is_possibly_mosaic = TRUE)
clean_plp  <- apply_filters(df_plp,  "Germline_PLP",   is_vus_group = FALSE, is_possibly_mosaic = FALSE)
clean_vus  <- apply_filters(df_vus,  "Germline_VUS",   is_vus_group = TRUE,  is_possibly_mosaic = FALSE)
clean_all  <- bind_rows(clean_poss, clean_plp, clean_vus)

out_dir <- "ukb_outputs"
dir.create(out_dir, showWarnings = FALSE)
write_csv(clean_poss, file.path(out_dir, "clean_PossiblyMosaic.csv"))
write_csv(clean_plp,  file.path(out_dir, "clean_Germline_PLP.csv"))
write_csv(clean_vus,  file.path(out_dir, "clean_Germline_VUS.csv"))
write_csv(clean_all,  file.path(out_dir, "clean_ALL_groups.csv"))

specimen_breakdown <- clean_all |>
  mutate(specimen_type = ifelse(is.na(specimen_type) | specimen_type=="", "Unknown", specimen_type)) |>
  count(group, specimen_type, name = "n") |>
  group_by(group) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup() |>
  arrange(group, desc(n))
write_csv(specimen_breakdown, file.path(out_dir, "specimen_breakdown_by_group.csv"))

bucket_eth <- function(x) {
  x <- toupper(trimws(as.character(x)))
  case_when(
    str_detect(x, "WHITE|CAUCASIAN|ASHKENAZI|MEDITERRANEAN") ~ "White",
    str_detect(x, "ASIAN") ~ "Asian",
    str_detect(x, "BLACK|AFRICAN") ~ "Black or African American",
    TRUE ~ "Other"
  )
}

# VP: Build a one-column summary block for a group (participants, age stats, sex, ethnicity)
mk_summary <- function(df, label) {
  id_col <- if ("dedup_id" %in% names(df) && COUNT_UNIQUE_PARTICIPANTS) "dedup_id" else NULL
  demog <- if (!is.null(id_col)) df |> arrange(.data[[id_col]]) |> distinct(.data[[id_col]], .keep_all = TRUE) else df
  total <- nrow(demog)

  age_col <- "AGE_AT_ACCESSIONING_IN_YEARS"
  if (age_col %in% names(demog)) demog[[age_col]] <- suppressWarnings(as.numeric(demog[[age_col]]))
  age_mean   <- if (age_col %in% names(demog)) mean(demog[[age_col]], na.rm=TRUE) else NA_real_
  age_median <- if (age_col %in% names(demog)) median(demog[[age_col]], na.rm=TRUE) else NA_real_

  sex_col <- "SEX_ASSIGNED_AT_BIRTH"
  sex_counts <- if (sex_col %in% names(demog)) table(demog[[sex_col]], useNA="no") else c(Female=NA, Male=NA)

  eth_vec <- if (!all(is.na(demog$ETH_CANON))) demog$ETH_CANON else demog$ETH_RAW
  ec <- if (!all(is.na(eth_vec))) table(bucket_eth(eth_vec)) else integer()

  fmt <- function(n) if (is.na(n)) "" else sprintf("%d (%.1f%%)", n, if (total>0) 100*n/total else 0)

  tibble::tibble(
    Characteristic = c("Total participants",
                       "Mean age at recruitment (years)",
                       "Median age at recruitment (years)",
                       "Sex — Female","Sex — Male",
                       "Self-reported ethnicity",
                       "  White","  Asian","  Black or African American","  Other"),
    !!label := c(
      total,
      ifelse(is.na(age_mean),"", sprintf("%.1f", age_mean)),
      ifelse(is.na(age_median),"", sprintf("%.1f", age_median)),
      fmt(unname(sex_counts[["Female"]] %||% 0)),
      fmt(unname(sex_counts[["Male"]] %||% 0)),
      "",
      fmt(unname(ec[["White"]] %||% 0)),
      fmt(unname(ec[["Asian"]] %||% 0)),
      fmt(unname(ec[["Black or African American"]] %||% 0)),
      fmt(unname(ec[["Other"]] %||% 0))
    )
  )
}

summary_tbl <- mk_summary(clean_plp, "P/LP") |>
  left_join(mk_summary(clean_vus, "Non-P/LP (VUS)"),  by="Characteristic") |>
  left_join(mk_summary(clean_poss, "Possibly Mosaic"), by="Characteristic")

write_csv(summary_tbl, file.path(out_dir, "summary_table_ukb_style.csv"))
tg <- gridExtra::tableGrob(summary_tbl, rows = NULL)
ggsave(file.path(out_dir, "summary_table_ukb_style.png"), tg, width = 10, height = 5, dpi = 300, units = "in")

message("\nCompleted. Outputs in ", normalizePath(out_dir), "\n")
