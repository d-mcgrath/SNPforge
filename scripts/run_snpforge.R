#!/usr/bin/env Rscript

# SNPforge: configurable SNP-based genetic association testing pipeline.
# Run from the command line with:
#   Rscript run_snpforge.R /path/to/config.yaml

suppressPackageStartupMessages({
  library(broom)
  library(fs)
  library(glue)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(tidyverse)
  library(vroom)
  library(yaml)
})

# =========================================================
# Parallel settings
# =========================================================
has_future <- requireNamespace("future", quietly = TRUE)
has_furrr <- requireNamespace("furrr", quietly = TRUE)
has_vcfR <- requireNamespace("vcfR", quietly = TRUE)

DEFAULT_COVARIATES <- character()
DEFAULT_FACTOR_COLS <- character()
SUPPORTED_OUTCOME_FAMILIES <- c("gaussian", "binomial", "ordinal")

if (has_future) {
  future::plan(future::multisession, workers = future::availableCores())
  options(future.globals.maxSize = 2.5 * 1024^3)
}

# =========================================================
# Utilities
# =========================================================
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

ensure_dir <- function(path) {
  fs::dir_create(path, recurse = TRUE)
  invisible(path)
}

normalize_id_values <- function(x) {
  x |>
    as.character() |>
    trimws()
}

full_table <- function(x) {
  table(x, useNA = "always")
}

pipeline_log_env <- new.env(parent = emptyenv())
pipeline_log_env$path <- NULL

normalize_id_type <- function(x, field_name = "id_type") {
  value <- x %||% "final_id"

  allowed_values <- c("raw_id", "final_id")

  if (!is.character(value) || length(value) != 1 || !value %in% allowed_values) {
    stop(
      sprintf(
        "'%s' must be one of: raw_id, final_id.",
        field_name
      ),
      call. = FALSE
    )
  }

  value
}

normalize_flag <- function(x, field_name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", field_name), call. = FALSE)
  }
  x
}

normalize_file_format <- function(x, field_name) {
  value <- tolower(x %||% "tsv")
  aliases <- c(
    tsv = "tsv",
    txt = "tsv",
    tab = "tsv",
    csv = "csv",
    delim = "delim"
  )

  if (!value %in% names(aliases)) {
    stop(
      sprintf("'%s' must be one of: tsv, csv, delim.", field_name),
      call. = FALSE
    )
  }

  aliases[[value]]
}

normalize_sep <- function(format, sep = NULL) {
  if (!is.null(sep)) {
    return(sep)
  }

  switch(
    format,
    tsv = "\t",
    csv = ",",
    delim = stop("A separator must be provided when format = 'delim'.", call. = FALSE),
    stop(sprintf("Unsupported format '%s'.", format), call. = FALSE)
  )
}

read_delim_flex <- function(path,
                            format = "tsv",
                            sep = NULL,
                            col_names = NULL,
                            ...) {
  format <- normalize_file_format(format, "format")
  sep <- normalize_sep(format, sep)

  args <- list(
    file = path,
    delim = sep,
    show_col_types = FALSE,
    ...
  )

  if (!is.null(col_names)) {
    args$col_names <- col_names
  }

  do.call(vroom::vroom, args)
}

normalize_config <- function(config) {
  required_fields <- c(
    "variant_file",
    "variant_file_format",
    "variant_col",
    "subjects_file",
    "genotype_file_template",
    "pheno_file",
    "pheno_format",
    "pheno_id_col",
    "genotype_id_type",
    "pheno_id_type",
    "outcomes",
    "covariates",
    "output_dir"
  )

  missing_fields <- required_fields[!map_lgl(required_fields, ~ !is.null(config[[.x]]))]
  if (length(missing_fields) > 0) {
    stop(
      sprintf(
        "Config is missing required field(s): %s",
        paste(missing_fields, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  config$variant_file_format <- normalize_file_format(
    config$variant_file_format,
    "variant_file_format"
  )
  config$variant_sep <- config$variant_sep %||% NULL
  config$strip_chr_for_vcf <- config$strip_chr_for_vcf %||% FALSE
  config$strip_chr_for_file_template <- config$strip_chr_for_file_template %||% FALSE

  config$pheno_format <- normalize_file_format(
    config$pheno_format,
    "pheno_format"
  )
  config$pheno_sep <- config$pheno_sep %||% NULL

  pc_fields_present <- any(map_lgl(
    c("pc_file", "pc_format", "pc_id_col", "pc_id_type", "pc_sep"),
    ~ !is.null(config[[.x]])
  ))

  if (is.null(config$use_pc_file)) {
    if (!pc_fields_present) {
      stop(
        paste(
          "Config must either provide PC-file parameters",
          "(pc_file, pc_format, pc_id_col, pc_id_type)",
          "or explicitly set use_pc_file: false."
        ),
        call. = FALSE
      )
    }
    config$use_pc_file <- TRUE
  }
  config$use_pc_file <- normalize_flag(config$use_pc_file, "use_pc_file")

  config$genotype_id_type <- normalize_id_type(config$genotype_id_type, "genotype_id_type")
  config$pheno_id_type <- normalize_id_type(config$pheno_id_type, "pheno_id_type")

  if (isTRUE(config$use_pc_file)) {
    pc_required_fields <- c("pc_file", "pc_format", "pc_id_col", "pc_id_type")
    missing_pc_fields <- pc_required_fields[!map_lgl(pc_required_fields, ~ !is.null(config[[.x]]))]

    if (length(missing_pc_fields) > 0) {
      stop(
        sprintf(
          "Config is missing required PC field(s): %s",
          paste(missing_pc_fields, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    config$pc_format <- normalize_file_format(
      config$pc_format,
      "pc_format"
    )
    config$pc_sep <- config$pc_sep %||% NULL
    config$pc_id_type <- normalize_id_type(config$pc_id_type, "pc_id_type")
  } else {
    config$pc_file <- NULL
    config$pc_format <- NULL
    config$pc_id_col <- NULL
    config$pc_id_type <- NULL
    config$pc_sep <- NULL
  }

  uses_raw_ids <- any(c(
    config$genotype_id_type,
    config$pheno_id_type,
    if (isTRUE(config$use_pc_file)) config$pc_id_type else NULL
  ) == "raw_id")

  if (uses_raw_ids) {
    mapping_required_fields <- c(
      "mapping_file",
      "mapping_file_format",
      "mapping_id_from",
      "mapping_id_to"
    )
    missing_mapping_fields <- mapping_required_fields[!map_lgl(mapping_required_fields, ~ !is.null(config[[.x]]))]

    if (length(missing_mapping_fields) > 0) {
      stop(
        sprintf(
          "Config is missing required mapping field(s): %s",
          paste(missing_mapping_fields, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    config$mapping_file_format <- normalize_file_format(
      config$mapping_file_format,
      "mapping_file_format"
    )
    config$mapping_sep <- config$mapping_sep %||% NULL
  } else {
    config$mapping_file <- NULL
    config$mapping_file_format <- NULL
    config$mapping_id_from <- NULL
    config$mapping_id_to <- NULL
    config$mapping_sep <- NULL
  }

  config$maf_filter <- config$maf_filter %||% 0.01
  config$variant_types <- config$variant_types %||% "snps"
  config$force_samples <- config$force_samples %||% TRUE
  config$outcomes <- normalize_outcomes(config$outcomes)
  config$factor_cols <- config$factor_cols %||% DEFAULT_FACTOR_COLS
  config$seed <- config$seed %||% 111
  config$interaction_col <- normalize_interaction_col(config$interaction_col %||% NULL)

  config
}

normalize_interaction_col <- function(interaction_col) {
  if (is.null(interaction_col)) {
    return(NULL)
  }

  if (!is.character(interaction_col) || length(interaction_col) != 1 || !nzchar(interaction_col)) {
    stop(
      "'interaction_col' must be a non-empty character string when provided.",
      call. = FALSE
    )
  }

  interaction_col
}

normalize_outcomes <- function(outcomes) {
  if (is.null(outcomes) || length(outcomes) == 0) {
    stop("'outcomes' must contain at least one outcome specification.", call. = FALSE)
  }

  normalize_family <- function(x, field_name = "family") {
    value <- tolower(x %||% "")

    if (!nzchar(value) || !value %in% SUPPORTED_OUTCOME_FAMILIES) {
      stop(
        sprintf(
          "'%s' must be one of: %s.",
          field_name,
          paste(SUPPORTED_OUTCOME_FAMILIES, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    value
  }

  outcome_tbl <- purrr::map_dfr(outcomes, function(x) {
    if (is.character(x) && length(x) == 1) {
      stop(
        paste(
          "'outcomes' entries must now be objects with 'name' and 'family' fields.",
          "Example: - name: \"FEV1\" family: \"gaussian\""
        ),
        call. = FALSE
      )
    }

    if (!is.list(x) || is.null(x$name) || is.null(x$family)) {
      stop(
        paste(
          "Each entry in 'outcomes' must be a list with 'name' and 'family'.",
          "Example: - name: \"EMPH\" family: \"binomial\""
        ),
        call. = FALSE
      )
    }

    if (!is.character(x$name) || length(x$name) != 1 || !nzchar(x$name)) {
      stop("Each outcome 'name' must be a non-empty character string.", call. = FALSE)
    }

    tibble(
      name = x$name,
      family = normalize_family(x$family, sprintf("family for outcome '%s'", x$name))
    )
  })

  if (anyDuplicated(outcome_tbl$name)) {
    dupes <- unique(outcome_tbl$name[duplicated(outcome_tbl$name)])
    stop(
      sprintf("Outcome names must be unique. Duplicates: %s", paste(dupes, collapse = ", ")),
      call. = FALSE
    )
  }

  outcome_tbl
}

read_config <- function(config_file) {
  normalize_config(yaml::read_yaml(config_file))
}

timestamp_string <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

set_pipeline_log <- function(log_path) {
  pipeline_log_env$path <- log_path
  invisible(log_path)
}

append_pipeline_log <- function(level, text) {
  log_path <- pipeline_log_env$path
  if (is.null(log_path) || is.na(log_path) || !nzchar(log_path)) {
    return(invisible(NULL))
  }

  cat(
    sprintf("[%s] [%s] %s\n", timestamp_string(), level, text),
    file = log_path,
    append = TRUE
  )

  invisible(NULL)
}

initialize_pipeline_log <- function(out_dir,
                                    log_filename = "pipeline.log") {
  logs_dir <- file.path(out_dir, "logs")
  ensure_dir(logs_dir)
  log_path <- file.path(logs_dir, log_filename)
  writeLines(
    sprintf("[%s] [INFO] Pipeline log initialized", timestamp_string()),
    con = log_path
  )
  set_pipeline_log(log_path)
}

pipeline_message <- function(text) {
  message(text)
  append_pipeline_log("INFO", text)
}

pipeline_warning <- function(text) {
  warning(text, call. = FALSE)
}

with_pipeline_logging <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      append_pipeline_log("WARN", conditionMessage(w))
    }
  )
}

parallel_pmap <- function(.l, .f, seed = 111, .progress = TRUE) {
  if (has_furrr) {
    return(
      furrr::future_pmap(
        .l,
        .f,
        .options = furrr::furrr_options(seed = seed),
        .progress = .progress
      )
    )
  }

  pipeline_message("Package 'furrr' is not installed; using sequential purrr::pmap().")
  purrr::pmap(.l, .f)
}

# =========================================================
# Subject IDs
# =========================================================
read_subject_mapping <- function(mapping_file,
                                 mapping_file_format = "csv",
                                 mapping_sep = NULL,
                                 id_from,
                                 id_to,
                                 ...) {
  read_delim_flex(mapping_file, format = mapping_file_format, sep = mapping_sep, ...) |>
    transmute(
      source_id = normalize_id_values(.data[[id_from]]),
      id = normalize_id_values(.data[[id_to]])
    ) |>
    filter(!is.na(source_id), !is.na(id), source_id != "", id != "") |>
    distinct(source_id, .keep_all = TRUE)
}

harmonize_ids <- function(tbl,
                          id_col,
                          mapping_tbl = NULL,
                          input_id_type = c("final_id", "raw_id"),
                          drop_unmapped = FALSE,
                          fallback_to_raw_id = TRUE) {
  input_id_type <- match.arg(input_id_type)

  if (!id_col %in% names(tbl)) {
    stop(sprintf("ID column '%s' not found in table.", id_col), call. = FALSE)
  }

  tbl <- tbl |>
    mutate(`__raw_id__` = normalize_id_values(.data[[id_col]])) |>
    filter(!is.na(`__raw_id__`), `__raw_id__` != "")

  if (is.null(mapping_tbl)) {
    return(
      tbl |>
        select(-any_of("id")) |>
        mutate(id = `__raw_id__`) |>
        select(-`__raw_id__`)
    )
  }

  mapping_tbl <- mapping_tbl |>
    mutate(
      source_id = normalize_id_values(source_id),
      id = normalize_id_values(id)
    ) |>
    distinct(source_id, .keep_all = TRUE)

  if (identical(input_id_type, "raw_id")) {
    out <- tbl |>
      select(-any_of("id")) |>
      left_join(mapping_tbl, by = c("__raw_id__" = "source_id"))

    if (isTRUE(fallback_to_raw_id)) {
      out <- out |>
        mutate(id = coalesce(id, `__raw_id__`))
    }
  } else {
    out <- tbl |>
      select(-any_of("id")) |>
      mutate(id = `__raw_id__`)
  }

  if (drop_unmapped) {
    out <- out |>
      filter(!is.na(id), id != "")
  }

  out |>
    select(-`__raw_id__`)
}

report_id_overlap <- function(genotype_tbl, pheno_tbl, pcs_tbl = NULL) {
  genotype_ids <- unique(genotype_tbl$id)
  pheno_ids <- unique(pheno_tbl$id)

  out <- tibble(
    table1 = "genotype",
    table2 = "phenotype",
    n_table1 = length(genotype_ids),
    n_table2 = length(pheno_ids),
    n_overlap = sum(genotype_ids %in% pheno_ids)
  )

  if (!is.null(pcs_tbl)) {
    pcs_ids <- unique(pcs_tbl$id)

    out <- bind_rows(
      out,
      tibble(
        table1 = "genotype",
        table2 = "pcs",
        n_table1 = length(genotype_ids),
        n_table2 = length(pcs_ids),
        n_overlap = sum(genotype_ids %in% pcs_ids)
      ),
      tibble(
        table1 = "phenotype",
        table2 = "pcs",
        n_table1 = length(pheno_ids),
        n_table2 = length(pcs_ids),
        n_overlap = sum(pheno_ids %in% pcs_ids)
      )
    )
  }

  out
}

# =========================================================
# Variant handling
# =========================================================
normalize_variant_label <- function(x) {
  x |>
    str_replace("^chr", "") |>
    str_replace("^([^:]+):(.*)$", "chr\\1:\\2")
}

build_variant_query_tbl <- function(variant_tbl,
                                    variant_col = "variant",
                                    strip_chr_for_vcf = FALSE,
                                    strip_chr_for_file_template = FALSE) {
  names(variant_tbl) <- trimws(names(variant_tbl))
  variant_col <- trimws(variant_col)

  exact_match <- names(variant_tbl)[names(variant_tbl) == variant_col]
  ci_match <- names(variant_tbl)[tolower(names(variant_tbl)) == tolower(variant_col)]

  first_or_null <- function(x) {
    if (length(x) == 0 || all(is.na(x))) {
      return(NULL)
    }
    x[[1]]
  }

  get_matching_col <- function(candidates) {
    matches <- names(variant_tbl)[tolower(names(variant_tbl)) %in% tolower(candidates)]
    first_or_null(matches)
  }

  all_nonempty <- function(x) {
    !any(vapply(x, is.null, logical(1))) && all(nzchar(unlist(x)))
  }

  selected_variant_col <- first_or_null(exact_match) %||% first_or_null(ci_match)

  if (is.null(selected_variant_col)) {
    chr_col <- get_matching_col(c("chr", "chrom", "chromosome"))
    pos_col <- get_matching_col(c("pos", "position", "bp"))
    ref_col <- get_matching_col(c("ref", "reference", "ref_allele"))
    alt_col <- get_matching_col(c("alt", "alternate", "alt_allele", "effect_allele"))
  }

  if (!is.null(selected_variant_col)) {
    x <- variant_tbl |>
      distinct(.data[[selected_variant_col]]) |>
      rename(variant = all_of(selected_variant_col))
  } else if (all_nonempty(list(chr_col, pos_col, ref_col, alt_col))) {
    x <- variant_tbl |>
      transmute(
        chr = .data[[chr_col]],
        pos = .data[[pos_col]],
        ref = .data[[ref_col]],
        alt = .data[[alt_col]]
      ) |>
      distinct() |>
      mutate(variant = paste(chr, pos, ref, alt, sep = ":")) |>
      select(variant)
  } else {
    stop(
      paste0(
        "Could not construct variants from the provided table. ",
        "Expected a variant column or chr/pos/ref/alt columns. ",
        "Available columns: ",
        paste(names(variant_tbl), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  x |>
    mutate(
      variant = normalize_variant_label(variant),
      chr_raw = str_replace(variant, "^([^:]+):.*$", "\\1"),
      pos = str_replace(variant, "^[^:]+:([^:]+):.*$", "\\1"),
      region_chromosome = if (strip_chr_for_vcf) {
        str_replace(chr_raw, "^chr", "")
      } else {
        chr_raw
      },
      file_chromosome = if (strip_chr_for_file_template) {
        str_replace(chr_raw, "^chr", "")
      } else {
        chr_raw
      },
      region = paste0(region_chromosome, ":", pos)
    ) |>
    select(variant, region, region_chromosome, file_chromosome)
}

build_genotype_file_path <- function(file_chromosome, genotype_file_template) {
  gsub("{chromosome}", file_chromosome, genotype_file_template, fixed = TRUE)
}

# =========================================================
# bcftools extraction
# =========================================================
build_bcftools_extract_cmd <- function(region,
                                       subjects_file,
                                       input_vcf,
                                       output_vcf,
                                       maf_filter = NULL,
                                       variant_types = "snps",
                                       force_samples = TRUE) {
  force_flag <- if (isTRUE(force_samples)) "--force-samples " else ""
  maf_flag <- if (!is.null(maf_filter)) glue("-i 'MAF>={maf_filter}' ") else ""

  glue(
    "bcftools view ",
    "{force_flag}",
    "-r {shQuote(region)} ",
    "-v {variant_types} ",
    "{maf_flag}",
    "-S {shQuote(subjects_file)} ",
    "{shQuote(input_vcf)} ",
    "-O z -o {shQuote(output_vcf)}"
  )
}

run_bcftools_extract <- function(region,
                                 subjects_file,
                                 input_vcf,
                                 output_vcf,
                                 maf_filter = NULL,
                                 variant_types = "snps",
                                 force_samples = TRUE) {
  cmd <- build_bcftools_extract_cmd(
    region = region,
    subjects_file = subjects_file,
    input_vcf = input_vcf,
    output_vcf = output_vcf,
    maf_filter = maf_filter,
    variant_types = variant_types,
    force_samples = force_samples
  )

  status <- system(cmd)
  if (status != 0) {
    stop(
      glue("bcftools command failed for region {region} using input {input_vcf}"),
      call. = FALSE
    )
  }

  invisible(output_vcf)
}

count_vcf_records <- function(vcf_file) {
  cmd <- sprintf("bcftools view -H %s | wc -l", shQuote(vcf_file))
  out <- suppressWarnings(system(cmd, intern = TRUE))
  as.integer(trimws(out))
}

# =========================================================
# Genotype parsing
# =========================================================
find_variant_row <- function(vcf, requested_variant) {
  if (nrow(vcf@fix) == 0) {
    return(integer())
  }

  fix_tbl <- as_tibble(vcf@fix) |>
    transmute(
      row_index = row_number(),
      variant = normalize_variant_label(paste(CHROM, POS, REF, ALT, sep = ":"))
    )

  fix_tbl |>
    filter(variant == normalize_variant_label(requested_variant)) |>
    pull(row_index)
}

read_single_variant_vcf <- function(vcf_file, requested_variant) {
  if (!has_vcfR) {
    stop("Package 'vcfR' is required to parse extracted VCF/BCF data.", call. = FALSE)
  }

  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  if (nrow(vcf@fix) == 0) {
    warning(sprintf("No variant rows found in %s", vcf_file))
    return(NULL)
  }

  matched_rows <- find_variant_row(vcf, requested_variant)
  if (length(matched_rows) == 0) {
    warning(sprintf("Requested variant %s not found in %s", requested_variant, vcf_file))
    return(NULL)
  }

  if (length(matched_rows) > 1) {
    warning(sprintf("Multiple rows matched %s in %s; keeping the first.", requested_variant, vcf_file))
    matched_rows <- matched_rows[1]
  }

  variant_vcf <- vcf[matched_rows, ]

  genotype_df <- variant_vcf |>
    vcfR::vcfR2genlight() |>
    as.data.frame()

  if (nrow(genotype_df) == 0 || ncol(genotype_df) == 0) {
    warning(sprintf("No genotype data found in %s for %s", vcf_file, requested_variant))
    return(NULL)
  }

  original_variant_id <- names(genotype_df)[1]

  tibble(id = normalize_id_values(rownames(genotype_df))) |>
    bind_cols(as_tibble(genotype_df)) |>
    rename(!!requested_variant := all_of(original_variant_id)) |>
    mutate(across(all_of(requested_variant), as.numeric)) |>
    arrange(id)
}

extract_and_read_variant <- function(variant,
                                     region,
                                     file_chromosome,
                                     subjects_file,
                                     genotype_file_template,
                                     out_dir,
                                     maf_filter = NULL,
                                     variant_types = "snps",
                                     force_samples = TRUE) {
  ensure_dir(out_dir)

  safe_variant_name <- gsub("[/:]", "_", variant)
  output_vcf <- file.path(out_dir, paste0(safe_variant_name, ".vcf.gz"))
  input_vcf <- build_genotype_file_path(file_chromosome, genotype_file_template)
  bcftools_cmd <- build_bcftools_extract_cmd(
    region = region,
    subjects_file = subjects_file,
    input_vcf = input_vcf,
    output_vcf = output_vcf,
    maf_filter = maf_filter,
    variant_types = variant_types,
    force_samples = force_samples
  )

  run_bcftools_extract(
    region = region,
    subjects_file = subjects_file,
    input_vcf = input_vcf,
    output_vcf = output_vcf,
    maf_filter = maf_filter,
    variant_types = variant_types,
    force_samples = force_samples
  )

  n_records <- count_vcf_records(output_vcf)
  if (is.na(n_records) || n_records == 0) {
    return(
      list(
        table = NULL,
        requested_variant = variant,
        output_vcf = output_vcf,
        n_records = 0L,
        bcftools_cmd = bcftools_cmd
      )
    )
  }

  genotype_tbl <- read_single_variant_vcf(
    vcf_file = output_vcf,
    requested_variant = variant
  )

  list(
    table = genotype_tbl,
    requested_variant = variant,
    output_vcf = output_vcf,
    n_records = n_records,
    bcftools_cmd = bcftools_cmd
  )
}

safe_extract_and_read_variant <- function(...) {
  bcftools_cmd <- NULL
  warning_messages <- character()

  result <- tryCatch(
    withCallingHandlers(
      extract_and_read_variant(...),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      args <- list(...)
      safe_variant_name <- gsub("[/:]", "_", args$variant)
      output_vcf <- file.path(args$out_dir, paste0(safe_variant_name, ".vcf.gz"))
      input_vcf <- build_genotype_file_path(args$file_chromosome, args$genotype_file_template)

      bcftools_cmd <<- build_bcftools_extract_cmd(
        region = args$region,
        subjects_file = args$subjects_file,
        input_vcf = input_vcf,
        output_vcf = output_vcf,
        maf_filter = args$maf_filter %||% NULL,
        variant_types = args$variant_types %||% "snps",
        force_samples = args$force_samples %||% TRUE
      )

      structure(
        list(message = conditionMessage(e)),
        class = c("pipeline_variant_error", "error", "condition")
      )
    }
  )

  if (inherits(result, "pipeline_variant_error")) {
    return(list(result = NULL, error = result, bcftools_cmd = bcftools_cmd))
  }

  list(
    result = result,
    error = NULL,
    bcftools_cmd = result$bcftools_cmd %||% bcftools_cmd,
    warnings = warning_messages
  )
}

combine_genotypes <- function(genotype_results,
                              genotype_summary = NULL) {
  genotype_tables <- genotype_results |>
    keep(~ is.null(.x$error) && !is.null(.x$result$table)) |>
    map(~ .x$result$table)

  if (length(genotype_tables) == 0) {
    detail_msg <- NULL

    if (!is.null(genotype_summary) && nrow(genotype_summary) > 0) {
      failing_variants <- genotype_summary |>
        mutate(
          failure_reason = case_when(
            nzchar(error) ~ error,
            nzchar(warning) ~ warning,
            !is.na(n_records) & n_records == 0 ~ "bcftools extracted 0 records",
            TRUE ~ "VCF parsed but no genotype table was produced"
          )
        ) |>
        filter(!success) |>
        transmute(detail = paste0(variant, " -> ", failure_reason)) |>
        pull(detail)

      if (length(failing_variants) > 0) {
        detail_msg <- paste(head(failing_variants, 10), collapse = "; ")
      }
    }

    stop(
      paste(
        "No genotype tables were successfully extracted.",
        if (!is.null(detail_msg)) paste("Examples:", detail_msg) else NULL
      ),
      call. = FALSE
    )
  }

  reduce(genotype_tables, left_join, by = "id")
}

summarize_genotype_results <- function(genotype_results, requested_variants = NULL) {
  tibble(
    variant = requested_variants %||% rep(NA_character_, length(genotype_results)),
    success = map_lgl(
      genotype_results,
      ~ is.null(.x$error) && !is.null(.x$result$table)
    ),
    n_records = map_int(
      genotype_results,
      ~ if (is.null(.x$error) && !is.null(.x$result$n_records)) .x$result$n_records else NA_integer_
    ),
    bcftools_cmd = map_chr(
      genotype_results,
      ~ .x$bcftools_cmd %||% .x$result$bcftools_cmd %||% ""
    ),
    warning = map_chr(
      genotype_results,
      ~ paste(.x$warnings %||% character(), collapse = " | ")
    ),
    error = map_chr(
      genotype_results,
      ~ if (is.null(.x$error)) "" else conditionMessage(.x$error)
    )
  )
}

# =========================================================
# Other input readers
# =========================================================
read_pheno_table <- function(pheno_file,
                             pheno_format = "tsv",
                             pheno_sep = NULL,
                             ...) {
  read_delim_flex(pheno_file, format = pheno_format, sep = pheno_sep, ...)
}

read_pcs_table <- function(pc_file,
                           pc_format = "tsv",
                           pc_sep = NULL,
                           id_col = "id",
                           ...) {
  pc_tbl <- read_delim_flex(
    pc_file,
    format = pc_format,
    sep = pc_sep,
    ...
  )

  if (!id_col %in% names(pc_tbl)) {
    pc_tbl_no_header <- read_delim_flex(
      pc_file,
      format = pc_format,
      sep = pc_sep,
      col_names = FALSE,
      ...
    )

    if (ncol(pc_tbl_no_header) >= 2) {
      names(pc_tbl_no_header) <- c("id", paste0("PC", seq_len(ncol(pc_tbl_no_header) - 1)))
      pc_tbl <- pc_tbl_no_header
    }
  }

  if (!id_col %in% names(pc_tbl) && "id" %in% names(pc_tbl)) {
    id_col <- "id"
  }

  if (!id_col %in% names(pc_tbl)) {
    stop(
      paste0(
        "PC ID column '", id_col, "' was not found in the PC file. ",
        "Available columns: ", paste(names(pc_tbl), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  pc_tbl |>
    rename(id = all_of(id_col)) |>
    mutate(id = normalize_id_values(id))
}

# =========================================================
# Table assembly
# =========================================================
prepare_analysis_table <- function(genotype_tbl,
                                   pheno_tbl,
                                   pcs_tbl = NULL,
                                   factor_cols = DEFAULT_FACTOR_COLS) {
  out <- genotype_tbl |>
    distinct(id, .keep_all = TRUE) |>
    left_join(
      pheno_tbl |>
        distinct(id, .keep_all = TRUE),
      by = "id"
    )

  if (!is.null(pcs_tbl)) {
    out <- out |>
      left_join(
        pcs_tbl |>
          distinct(id, .keep_all = TRUE),
        by = "id"
      )
  }

  for (col in factor_cols) {
    if (col %in% names(out) && is.character(out[[col]])) {
      out[[col]] <- factor(out[[col]])
    }
  }

  out
}

# =========================================================
# Modeling
# =========================================================
prepare_model_df <- function(data,
                             vars_needed,
                             covariates,
                             interaction_term = NULL,
                             outcome = NULL,
                             family = "gaussian") {
  model_df <- data |>
    select(all_of(vars_needed)) |>
    drop_na()

  if (!is.null(outcome) &&
      outcome %in% names(model_df) &&
      identical(family, "gaussian") &&
      !is.numeric(model_df[[outcome]])) {
    model_df[[outcome]] <- suppressWarnings(as.numeric(as.character(model_df[[outcome]])))
  }

  if (!is.null(outcome) &&
      outcome %in% names(model_df) &&
      identical(family, "binomial")) {
    if (!is.numeric(model_df[[outcome]])) {
      model_df[[outcome]] <- suppressWarnings(as.numeric(as.character(model_df[[outcome]])))
    }

    model_df <- model_df |>
      filter(.data[[outcome]] %in% c(0, 1))
  }

  if (!is.null(outcome) &&
      outcome %in% names(model_df) &&
      identical(family, "ordinal")) {
    if (!is.factor(model_df[[outcome]])) {
      model_df[[outcome]] <- suppressWarnings(as.numeric(as.character(model_df[[outcome]])))
    }

    if (is.numeric(model_df[[outcome]])) {
      model_df <- model_df |>
        filter(is.finite(.data[[outcome]]))

      outcome_levels <- sort(unique(model_df[[outcome]]))
      model_df[[outcome]] <- ordered(model_df[[outcome]], levels = outcome_levels)
    } else {
      model_df[[outcome]] <- ordered(model_df[[outcome]])
    }
  }

  numeric_cols <- names(model_df)[vapply(model_df, is.numeric, logical(1))]

  if (length(numeric_cols) > 0) {
    finite_rows <- model_df |>
      transmute(across(all_of(numeric_cols), is.finite)) |>
      mutate(`__all_finite__` = if_all(everything(), identity)) |>
      pull(`__all_finite__`)

    model_df <- model_df[finite_rows, , drop = FALSE]
  }

  if (!is.null(interaction_term) && interaction_term %in% names(model_df) && is.character(model_df[[interaction_term]])) {
    model_df[[interaction_term]] <- factor(model_df[[interaction_term]])
  }

  for (cov in covariates) {
    if (cov %in% names(model_df) && is.character(model_df[[cov]])) {
      model_df[[cov]] <- factor(model_df[[cov]])
    }
  }

  droplevels(model_df)
}

assess_model_data <- function(data,
                              outcome,
                              snp,
                              covariates = DEFAULT_COVARIATES,
                              interaction_term = NULL,
                              family = "gaussian") {
  vars_needed <- unique(c(outcome, snp, covariates, interaction_term %||% character()))
  missing_vars <- setdiff(vars_needed, names(data))

  if (length(missing_vars) > 0) {
    return(tibble(
      outcome = outcome,
      snp = snp,
      status = "skip",
      reason = paste0("missing_columns: ", paste(missing_vars, collapse = ", ")),
      n_total = nrow(data),
      n_complete = NA_integer_
    ))
  }

  model_df <- prepare_model_df(
    data = data,
    vars_needed = vars_needed,
    covariates = covariates,
    interaction_term = interaction_term,
    outcome = outcome,
    family = family
  )

  if (nrow(model_df) == 0) {
    return(tibble(
      outcome = outcome,
      snp = snp,
      status = "skip",
      reason = "no_complete_finite_cases",
      n_total = nrow(data),
      n_complete = 0L
    ))
  }

  if (!is.null(interaction_term)) {
    if (is.factor(model_df[[interaction_term]]) && nlevels(model_df[[interaction_term]]) < 2) {
      return(tibble(
        outcome = outcome,
        snp = snp,
        status = "skip",
        reason = paste0("interaction_term_<2_levels: ", interaction_term),
        n_total = nrow(data),
        n_complete = nrow(model_df)
      ))
    }
  }

  bad_factors <- covariates[
    covariates %in% names(model_df) &
      map_lgl(covariates[covariates %in% names(model_df)], ~ is.factor(model_df[[.x]]) && nlevels(model_df[[.x]]) < 2)
  ]

  if (length(bad_factors) > 0) {
    return(tibble(
      outcome = outcome,
      snp = snp,
      status = "skip",
      reason = paste0("factor_<2_levels: ", paste(bad_factors, collapse = ", ")),
      n_total = nrow(data),
      n_complete = nrow(model_df)
    ))
  }

  if (identical(family, "binomial")) {
    outcome_levels <- sort(unique(model_df[[outcome]]))

    if (length(outcome_levels) < 2) {
      return(tibble(
        outcome = outcome,
        snp = snp,
        status = "skip",
        reason = paste0("binomial_outcome_<2_levels: ", outcome),
        n_total = nrow(data),
        n_complete = nrow(model_df)
      ))
    }
  }

  if (identical(family, "ordinal")) {
    if (!is.ordered(model_df[[outcome]])) {
      model_df[[outcome]] <- ordered(model_df[[outcome]])
    }

    if (nlevels(model_df[[outcome]]) < 2) {
      return(tibble(
        outcome = outcome,
        snp = snp,
        status = "skip",
        reason = paste0("ordinal_outcome_<2_levels: ", outcome),
        n_total = nrow(data),
        n_complete = nrow(model_df)
      ))
    }
  }

  tibble(
    outcome = outcome,
    snp = snp,
    status = "ok",
    reason = "ready",
    n_total = nrow(data),
    n_complete = nrow(model_df)
  )
}

fit_model <- function(data,
                      outcome,
                      snp,
                      covariates = DEFAULT_COVARIATES,
                      interaction_term = NULL,
                      family = "gaussian") {
  vars_needed <- unique(c(outcome, snp, covariates, interaction_term %||% character()))
  model_df <- prepare_model_df(
    data = data,
    vars_needed = vars_needed,
    covariates = covariates,
    interaction_term = interaction_term,
    outcome = outcome,
    family = family
  )

  if (nrow(model_df) == 0) {
    return(NULL)
  }

  snp_term <- if (is.null(interaction_term)) {
    glue("`{snp}`")
  } else {
    glue("`{snp}` * {interaction_term}")
  }

  rhs <- paste(c(snp_term, covariates), collapse = " + ")

  fml <- as.formula(glue("`{outcome}` ~ {rhs}"))

  if (identical(family, "gaussian")) {
    lm(fml, data = model_df)
  } else if (identical(family, "binomial")) {
    glm(fml, data = model_df, family = binomial())
  } else if (identical(family, "ordinal")) {
    MASS::polr(fml, data = model_df, Hess = TRUE)
  } else {
    stop("Unsupported family: ", family, call. = FALSE)
  }
}

safe_fit_model <- function(...) {
  tryCatch(
    list(model = fit_model(...), fit_error = ""),
    error = function(e) list(model = NULL, fit_error = conditionMessage(e))
  )
}

empty_model_terms_tbl <- function() {
  tibble(
    outcome = character(),
    snp = character(),
    family = character(),
    term = character(),
    estimate = numeric(),
    std.error = numeric(),
    statistic = numeric(),
    p.value = numeric(),
    conf.low = numeric(),
    conf.high = numeric(),
    n_complete = integer(),
    effect_scale = character(),
    direction = character(),
    reference_level = character(),
    odds_ratio = numeric(),
    or_conf.low = numeric(),
    or_conf.high = numeric()
  )
}

extract_reference_level <- function(term, xlevels) {
  if (is.null(xlevels) || length(xlevels) == 0 || is.na(term) || !nzchar(term)) {
    return(NA_character_)
  }

  if (term == "(Intercept)" || str_detect(term, "\\|")) {
    return(NA_character_)
  }

  term_parts <- str_split(term, ":", simplify = FALSE)[[1]] |>
    str_replace_all("`", "")

  refs <- purrr::map_chr(names(xlevels), function(var) {
    matched <- any(term_parts != var & startsWith(term_parts, var))
    if (!matched) {
      return(NA_character_)
    }

    ref_level <- xlevels[[var]][1]
    if (is.null(ref_level) || is.na(ref_level) || !nzchar(ref_level)) {
      return(NA_character_)
    }

    paste0(var, "=", ref_level)
  })

  refs <- refs[!is.na(refs) & nzchar(refs)]

  if (length(refs) == 0) {
    return(NA_character_)
  }

  paste(refs, collapse = "; ")
}

annotate_model_terms <- function(tidy_tbl,
                                 model,
                                 family,
                                 n_complete) {
  if (is.null(tidy_tbl) || nrow(tidy_tbl) == 0) {
    return(empty_model_terms_tbl())
  }

  xlevels <- model$xlevels %||% list()
  effect_scale <- case_when(
    family == "gaussian" ~ "coefficient",
    family == "binomial" ~ "log_odds",
    family == "ordinal" ~ "proportional_odds_logit",
    TRUE ~ "coefficient"
  )

  has_se <- "std.error" %in% names(tidy_tbl)

  if (!"conf.low" %in% names(tidy_tbl)) {
    tidy_tbl$conf.low <- if (has_se) tidy_tbl$estimate - 1.96 * tidy_tbl$std.error else NA_real_
  }

  if (!"conf.high" %in% names(tidy_tbl)) {
    tidy_tbl$conf.high <- if (has_se) tidy_tbl$estimate + 1.96 * tidy_tbl$std.error else NA_real_
  }

  tidy_tbl |>
    mutate(
      n_complete = as.integer(n_complete),
      effect_scale = effect_scale,
      direction = case_when(
        is.na(estimate) ~ NA_character_,
        estimate > 0 ~ "positive",
        estimate < 0 ~ "negative",
        TRUE ~ "neutral"
      ),
      reference_level = purrr::map_chr(term, ~ extract_reference_level(.x, xlevels)),
      odds_ratio = if_else(family %in% c("binomial", "ordinal") & !str_detect(term, "\\|"), exp(estimate), NA_real_),
      or_conf.low = if_else(family %in% c("binomial", "ordinal") & !str_detect(term, "\\|") & !is.na(conf.low), exp(conf.low), NA_real_),
      or_conf.high = if_else(family %in% c("binomial", "ordinal") & !str_detect(term, "\\|") & !is.na(conf.high), exp(conf.high), NA_real_)
    )
}

run_model_scan <- function(data,
                           outcomes,
                           snps,
                           covariates = DEFAULT_COVARIATES,
                           interaction_term = NULL) {
  design_tbl <- outcomes |>
    rename(outcome = name) |>
    crossing(snp = snps)

  diagnostics <- design_tbl |>
    mutate(
      diag = pmap(
        list(outcome, snp, family),
        ~ assess_model_data(
          data = data,
          outcome = ..1,
          snp = ..2,
          covariates = covariates,
          interaction_term = interaction_term,
          family = ..3
        )
      ),
      diag = map(
        diag,
        ~ select(.x, -any_of(c("outcome", "snp")))
      )
    ) |>
    unnest(diag)

  runnable <- diagnostics |>
    filter(status == "ok") |>
    select(outcome, snp, family)

  if (nrow(runnable) == 0) {
    return(list(
      diagnostics = diagnostics,
      model_terms = empty_model_terms_tbl()
    ))
  }

  model_terms <- runnable |>
    left_join(
      diagnostics |>
        select(outcome, snp, n_complete),
      by = c("outcome", "snp")
    ) |>
    mutate(
      fit = pmap(
        list(outcome, snp, family),
        ~ safe_fit_model(
          data = data,
          outcome = ..1,
          snp = ..2,
          covariates = covariates,
          interaction_term = interaction_term,
          family = ..3
        )
      ),
      model = map(fit, "model"),
      fit_error = map_chr(fit, "fit_error"),
      tidy = pmap(
        list(model, n_complete, family),
        ~ if (is.null(..1)) {
          NULL
        } else {
          tidy_tbl <- broom::tidy(..1)
          annotate_model_terms(
            tidy_tbl = tidy_tbl,
            model = ..1,
            family = ..3,
            n_complete = ..2
          )
        }
      )
    ) |>
    select(outcome, snp, family, fit_error, tidy)

  diagnostics <- diagnostics |>
    left_join(
      model_terms |>
        transmute(
          outcome,
          snp,
          fit_error,
          fit_status = if_else(fit_error == "", "fit_ok", "fit_error")
        ),
      by = c("outcome", "snp")
    ) |>
    mutate(
      fit_error = coalesce(fit_error, ""),
      fit_status = case_when(
        status != "ok" ~ "screened_out",
        fit_status == "fit_error" ~ "fit_error",
        TRUE ~ "fit_ok"
      )
    )

  model_terms <- model_terms |>
    filter(!map_lgl(tidy, is.null)) |>
    select(outcome, snp, family, tidy) |>
    unnest(tidy)

  list(
    diagnostics = diagnostics,
    model_terms = model_terms
  )
}

extract_interaction_terms <- function(model_results,
                                      interaction_term = NULL) {
  empty_terms <- empty_model_terms_tbl()

  if (is.null(interaction_term) || !"term" %in% names(model_results) || nrow(model_results) == 0) {
    return(empty_terms)
  }

  out <- model_results |>
    filter(str_detect(term, paste0(":", interaction_term)))

  if (nrow(out) == 0) {
    return(empty_terms)
  }

  out
}

# =========================================================
# Output writing
# =========================================================
save_outputs <- function(genotype_tbl,
                         analysis_tbl,
                         all_model_results,
                         interaction_results,
                         genotype_summary,
                         variant_status_tbl,
                         id_overlap_tbl,
                         model_diagnostics,
                         out_dir) {
  ensure_dir(out_dir)
  ensure_dir(file.path(out_dir, "tables"))
  ensure_dir(file.path(out_dir, "models"))
  ensure_dir(file.path(out_dir, "logs"))

  readr::write_tsv(genotype_tbl, file.path(out_dir, "tables", "genotypes.tsv.gz"))
  readr::write_tsv(analysis_tbl, file.path(out_dir, "tables", "analysis_table.tsv.gz"))
  readr::write_tsv(all_model_results, file.path(out_dir, "models", "all_model_terms.tsv.gz"))
  readr::write_tsv(interaction_results, file.path(out_dir, "models", "interaction_results.tsv.gz"))
  readr::write_tsv(genotype_summary, file.path(out_dir, "logs", "genotype_extraction_summary.tsv"))
  readr::write_tsv(variant_status_tbl, file.path(out_dir, "logs", "variant_status.tsv"))
  readr::write_tsv(id_overlap_tbl, file.path(out_dir, "logs", "id_overlap_summary.tsv"))
  readr::write_tsv(model_diagnostics, file.path(out_dir, "logs", "model_diagnostics.tsv"))

  saveRDS(genotype_tbl, file.path(out_dir, "tables", "genotypes.rds"))
  saveRDS(analysis_tbl, file.path(out_dir, "tables", "analysis_table.rds"))
  saveRDS(all_model_results, file.path(out_dir, "models", "all_model_terms.rds"))
  saveRDS(interaction_results, file.path(out_dir, "models", "interaction_results.rds"))
  saveRDS(genotype_summary, file.path(out_dir, "logs", "genotype_extraction_summary.rds"))
  saveRDS(variant_status_tbl, file.path(out_dir, "logs", "variant_status.rds"))
  saveRDS(id_overlap_tbl, file.path(out_dir, "logs", "id_overlap_summary.rds"))
  saveRDS(model_diagnostics, file.path(out_dir, "logs", "model_diagnostics.rds"))
}

# =========================================================
# Main pipeline
# =========================================================
run_pipeline <- function(config) {
  ensure_dir(config$output_dir)
  ensure_dir(file.path(config$output_dir, "vcfs"))
  if (is.null(pipeline_log_env$path)) {
    initialize_pipeline_log(config$output_dir)
  }

  with_pipeline_logging({
    pipeline_message("Reading variant table...")
    variant_tbl <- read_delim_flex(
      config$variant_file,
      format = config$variant_file_format,
      sep = config$variant_sep
    )

    variant_query_tbl <- build_variant_query_tbl(
      variant_tbl = variant_tbl,
      variant_col = config$variant_col %||% "variant",
      strip_chr_for_vcf = config$strip_chr_for_vcf %||% FALSE,
      strip_chr_for_file_template = config$strip_chr_for_file_template %||% FALSE
    )
    pipeline_message(sprintf("Prepared %s variant queries.", nrow(variant_query_tbl)))

    pipeline_message("Extracting genotype data with bcftools...")
    genotype_results <- parallel_pmap(
      list(
        variant_query_tbl$variant,
        variant_query_tbl$region,
        variant_query_tbl$file_chromosome
      ),
      function(variant, region, file_chromosome) {
        safe_extract_and_read_variant(
          variant = variant,
          region = region,
          file_chromosome = file_chromosome,
          subjects_file = config$subjects_file,
          genotype_file_template = config$genotype_file_template,
          out_dir = file.path(config$output_dir, "vcfs"),
          maf_filter = config$maf_filter %||% NULL,
          variant_types = config$variant_types %||% "snps",
          force_samples = config$force_samples %||% TRUE
        )
      },
      seed = config$seed %||% 111,
      .progress = TRUE
    )

    genotype_summary <- summarize_genotype_results(
      genotype_results = genotype_results,
      requested_variants = variant_query_tbl$variant
    )
    pipeline_message(sprintf(
      "Genotype extraction succeeded for %s of %s variants.",
      sum(genotype_summary$success, na.rm = TRUE),
      nrow(genotype_summary)
    ))
    append_pipeline_log(
      "INFO",
      paste(capture.output(print(genotype_summary, n = min(20, nrow(genotype_summary)))), collapse = "\n")
    )
    problematic_genotype_cmds <- genotype_summary |>
      filter(!success | nzchar(warning) | nzchar(error)) |>
      transmute(
        detail = paste0(
          variant,
          " | success=",
          success,
          " | n_records=",
          if_else(is.na(n_records), "NA", as.character(n_records)),
          " | warning=",
          if_else(warning == "", "none", warning),
          " | error=",
          if_else(error == "", "unknown", error),
          " | cmd=",
          bcftools_cmd
        )
      ) |>
      pull(detail)

    if (length(problematic_genotype_cmds) > 0) {
      append_pipeline_log(
        "WARN",
        paste(
          "Problematic genotype extractions:",
          paste(problematic_genotype_cmds, collapse = "\n"),
          sep = "\n"
        )
      )
    }

    genotype_tbl <- combine_genotypes(
      genotype_results,
      genotype_summary = genotype_summary
    )
    pipeline_message(sprintf("Combined genotype table has %s rows and %s SNP columns.", nrow(genotype_tbl), ncol(genotype_tbl) - 1))

    pipeline_message("Reading phenotype data...")
    pheno_tbl <- read_pheno_table(
      pheno_file = config$pheno_file,
      pheno_format = config$pheno_format,
      pheno_sep = config$pheno_sep
    )
    pipeline_message(sprintf("Phenotype table has %s rows before ID harmonization.", nrow(pheno_tbl)))

    mapping_tbl <- NULL
    if (!is.null(config$mapping_file)) {
      pipeline_message("Reading subject ID mapping...")
      mapping_tbl <- read_subject_mapping(
        mapping_file = config$mapping_file,
        mapping_file_format = config$mapping_file_format,
        mapping_sep = config$mapping_sep,
        id_from = config$mapping_id_from,
        id_to = config$mapping_id_to
      )
      pipeline_message(sprintf("Loaded %s subject ID mappings.", nrow(mapping_tbl)))
    } else {
      pipeline_message("No mapping file provided; assuming all joined tables already use final_id values.")
    }

    pipeline_message("Harmonizing genotype IDs...")
    genotype_tbl <- harmonize_ids(
      tbl = genotype_tbl,
      id_col = "id",
      mapping_tbl = mapping_tbl,
      input_id_type = config$genotype_id_type %||% "raw_id",
      drop_unmapped = FALSE
    ) |>
      distinct(id, .keep_all = TRUE)
    pipeline_message(sprintf("Genotype table has %s distinct IDs after harmonization.", nrow(genotype_tbl)))

    pipeline_message("Harmonizing phenotype IDs...")
    pheno_tbl <- harmonize_ids(
      tbl = pheno_tbl,
      id_col = config$pheno_id_col,
      mapping_tbl = mapping_tbl,
      input_id_type = config$pheno_id_type %||% "final_id",
      drop_unmapped = TRUE
    ) |>
      distinct(id, .keep_all = TRUE)
    pipeline_message(sprintf("Phenotype table has %s distinct IDs after harmonization.", nrow(pheno_tbl)))

    pcs_tbl <- NULL
    if (isTRUE(config$use_pc_file)) {
      pipeline_message("Reading PC data...")
      pcs_tbl <- read_pcs_table(
        pc_file = config$pc_file,
        pc_format = config$pc_format,
        pc_sep = config$pc_sep,
        id_col = config$pc_id_col %||% "id"
      )
      pipeline_message(sprintf("PC table has %s rows before ID harmonization.", nrow(pcs_tbl)))

      pipeline_message("Harmonizing PC IDs...")
      pcs_tbl <- harmonize_ids(
        tbl = pcs_tbl,
        id_col = "id",
        mapping_tbl = mapping_tbl,
        input_id_type = config$pc_id_type %||% "final_id",
        drop_unmapped = FALSE
      ) |>
        distinct(id, .keep_all = TRUE)
      pipeline_message(sprintf("PC table has %s distinct IDs after harmonization.", nrow(pcs_tbl)))
    } else {
      pipeline_message("No separate PC file requested; assuming any needed PC covariates are already present in the phenotype table.")
    }

    id_overlap_tbl <- report_id_overlap(
      genotype_tbl = genotype_tbl,
      pheno_tbl = pheno_tbl,
      pcs_tbl = pcs_tbl
    )
    pipeline_message("ID overlap summary:")
    append_pipeline_log("INFO", paste(capture.output(print(id_overlap_tbl)), collapse = "\n"))

    genotype_pheno_overlap <- id_overlap_tbl |>
      filter(table1 == "genotype", table2 == "phenotype") |>
      pull(n_overlap)

    if (length(genotype_pheno_overlap) == 1 && genotype_pheno_overlap == 0) {
      pipeline_warning(
        paste(
          "No genotype IDs overlapped with phenotype IDs after harmonization.",
          "Check genotype_id_type, pheno_id_type, and the subject mapping file."
        )
      )
    }

    pipeline_message("Building analysis table...")
    analysis_tbl <- prepare_analysis_table(
      genotype_tbl = genotype_tbl,
      pheno_tbl = pheno_tbl,
      pcs_tbl = pcs_tbl,
      factor_cols = config$factor_cols %||% DEFAULT_FACTOR_COLS
    )
    pipeline_message(sprintf("Analysis table has %s rows and %s columns.", nrow(analysis_tbl), ncol(analysis_tbl)))

    implicitly_factored_covariates <- setdiff(
      config$covariates[
        config$covariates %in% names(analysis_tbl) &
          map_lgl(config$covariates[config$covariates %in% names(analysis_tbl)], ~ is.character(analysis_tbl[[.x]]))
      ],
      config$factor_cols %||% character()
    )

    if (length(implicitly_factored_covariates) > 0) {
      pipeline_warning(
        paste(
          "The following character covariates were not listed in factor_cols",
          "but will be automatically converted to factors during model fitting:",
          paste(implicitly_factored_covariates, collapse = ", ")
        )
      )
    }

    interaction_term <- config$interaction_col %||% NULL
    if (!is.null(interaction_term)) {
      if (!interaction_term %in% names(analysis_tbl)) {
        stop(
          sprintf(
            "Configured interaction column '%s' was not found in the analysis table. Create it before running the pipeline.",
            interaction_term
          ),
          call. = FALSE
        )
      }

      pipeline_message(sprintf("Using pre-existing interaction column '%s'.", interaction_term))
      pipeline_message(sprintf(
        "%s distribution: %s",
        interaction_term,
        paste(capture.output(print(full_table(analysis_tbl[[interaction_term]]))), collapse = " ")
      ))
    }

    modeled_snps <- setdiff(names(genotype_tbl), "id")
    if (length(modeled_snps) == 0) {
      stop("No successfully extracted SNP genotype columns were available for modeling.", call. = FALSE)
    }

    variant_status_tbl <- variant_query_tbl |>
      mutate(extracted_successfully = variant %in% modeled_snps)

    if (is.null(interaction_term)) {
      pipeline_message("Running main-effect models...")
    } else {
      pipeline_message(sprintf("Running SNP-by-%s interaction models...", interaction_term))
    }

    model_scan <- run_model_scan(
      data = analysis_tbl,
      outcomes = config$outcomes,
      snps = modeled_snps,
      covariates = config$covariates %||% DEFAULT_COVARIATES,
      interaction_term = interaction_term
    )

    model_diagnostics <- model_scan$diagnostics
    all_model_results <- model_scan$model_terms
    interaction_results <- extract_interaction_terms(
      model_results = all_model_results,
      interaction_term = interaction_term
    )

    if ("p.value" %in% names(interaction_results)) {
      interaction_results <- interaction_results |>
        arrange(p.value)
    }

    if (all(c("fit_status", "fit_error") %in% names(model_diagnostics))) {
      fit_error_summary <- model_diagnostics |>
        filter(fit_status == "fit_error", nzchar(fit_error)) |>
        count(fit_error, sort = TRUE)

      if (nrow(fit_error_summary) > 0) {
        append_pipeline_log(
          "WARN",
          paste(
            "Model fit errors:",
            paste(capture.output(print(fit_error_summary, n = min(10, nrow(fit_error_summary)))), collapse = "\n"),
            sep = "\n"
          )
        )
      }
    }

    pipeline_message(sprintf(
      "Model scan produced %s total coefficient rows and %s interaction rows.",
      nrow(all_model_results),
      nrow(interaction_results)
    ))

    pipeline_message("Saving outputs...")
    save_outputs(
      genotype_tbl = genotype_tbl,
      analysis_tbl = analysis_tbl,
      all_model_results = all_model_results,
      interaction_results = interaction_results,
      genotype_summary = genotype_summary,
      variant_status_tbl = variant_status_tbl,
      id_overlap_tbl = id_overlap_tbl,
      model_diagnostics = model_diagnostics,
      out_dir = config$output_dir
    )
    pipeline_message("Pipeline completed successfully.")

    list(
      variant_query_tbl = variant_query_tbl,
      genotype_summary = genotype_summary,
      variant_status_tbl = variant_status_tbl,
      id_overlap_tbl = id_overlap_tbl,
      genotype_tbl = genotype_tbl,
      analysis_tbl = analysis_tbl,
      model_diagnostics = model_diagnostics,
      all_model_results = all_model_results,
      interaction_results = interaction_results,
      model_results = all_model_results
    )
  })
}

# =========================================================
# Command-line driver
# =========================================================
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) < 1 || is.na(args[1]) || !nzchar(args[1])) {
    stop("Please provide a config file path as the first command-line argument.", call. = FALSE)
  }
  config_file <- args[1]
  config <- read_config(config_file)
  initialize_pipeline_log(config$output_dir)
  results <- tryCatch(
    run_pipeline(config),
    error = function(e) {
      append_pipeline_log("ERROR", conditionMessage(e))
      stop(e)
    }
  )
  invisible(results)
}

if (sys.nframe() == 0) {
  main()
}
