library(shiny)
library(readxl)
library(dplyr)
library(glmmTMB)
library(ggplot2)

## =========================================================
## Shinylive download-button workaround
## Required for file downloads in Chromium-based browsers
## =========================================================

shinylive_downloadButton <- function(...) {
  tag <- shiny::downloadButton(...)
  tag$attribs$download <- NULL
  tag
}
## Optional base directory supplied by the launcher.
## This keeps all relative input/output paths anchored to the toolkit folder.
app_base_dir <- getOption("bioservices.base_dir", NULL)
if (!is.null(app_base_dir) && dir.exists(app_base_dir)) {
  setwd(app_base_dir)
}


## nell'app via via che aggiungo indicatori e variabili corrette log devo modificare
# log1p_indicators e biod_labels
## =========================================================
## 1. Input files
## Keep these files in the same folder as app.R
## =========================================================
indicator_file <- "multifunctionality_region_lu_indicator_sign_set.xlsx"
intercept_file <- "model_intercepts_region_lu.xlsx"
data_file      <- "df_shiny.rds"
bundle_file    <- "shiny_model_bundle.rds"

pipeline_dir            <- "Pipeline_results_MIXED"
selected_table_file      <- file.path(pipeline_dir, "selected_models_table_mixed.csv")
interaction_export_file  <- file.path(
  pipeline_dir,
  "multifunctionality_region_lu_indicator_sign_set_mixed.csv"
)
model_checks_file        <- file.path(pipeline_dir, "model_checks_mixed.csv")

if (!file.exists(data_file)) {
  stop(
    "Missing file: ", data_file,
    "\nSave df_shiny with: saveRDS(df_shiny, 'df_shiny.rds')"
  )
}

required_result_columns <- c(
  "ES", "biodiversity", "model_name", "selected_model", "slope_scope",
  "Region_LU", "trend", "SE", "p.value", "p_adjustment", "direction",
  "selection_delta_AIC", "selection_LRT_p"
)

required_intercept_columns <- c(
  "ES", "biodiversity", "model_name", "selected_model", "Region_LU",
  "biod_value_reference", "intercept_link", "intercept_SE",
  "intercept_response"
)

build_indicator_results_from_pipeline <- function() {
  required_fallback_files <- c(selected_table_file, interaction_export_file)
  missing_fallback_files <- required_fallback_files[!file.exists(required_fallback_files)]

  if (length(missing_fallback_files) > 0) {
    stop(
      "The additive-aware indicator workbook is unavailable or uses the old schema, ",
      "and the fallback pipeline files are missing: ",
      paste(missing_fallback_files, collapse = ", ")
    )
  }

  selected_table <- utils::read.csv(
    selected_table_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    dplyr::mutate(
      ES = as.character(ES),
      biodiversity = as.character(biodiversity),
      selected_model = tolower(as.character(selected_model)),
      model_name = paste(ES, biodiversity, sep = "__")
    )

  required_selected_columns <- c(
    "ES", "biodiversity", "selected_model",
    "delta_AIC_m1_vs_m0", "p_comp_m1_vs_m0_biod",
    "delta_AIC_m2_vs_m1", "p_comp_m2_vs_m1_interaction",
    "summary_est_m1_biod", "summary_SE_m1_biod",
    "summary_p_m1_biod", "additive_direction"
  )

  missing_selected_columns <- setdiff(
    required_selected_columns,
    names(selected_table)
  )

  if (length(missing_selected_columns) > 0) {
    stop(
      "Fallback selected-model table is missing columns: ",
      paste(missing_selected_columns, collapse = ", ")
    )
  }

  if (file.exists(model_checks_file)) {
    model_checks <- utils::read.csv(
      model_checks_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    ok_names <- model_checks$model_name[
      !is.na(model_checks$model_ok_numeric) &
        as.logical(model_checks$model_ok_numeric)
    ]

    selected_table <- selected_table %>%
      dplyr::filter(model_name %in% ok_names)
  }

  interaction_export <- utils::read.csv(
    interaction_export_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  interaction_lookup <- selected_table %>%
    dplyr::filter(selected_model == "interaction") %>%
    dplyr::select(
      ES, biodiversity, model_name,
      selection_delta_AIC = delta_AIC_m2_vs_m1,
      selection_LRT_p = p_comp_m2_vs_m1_interaction
    )

  interaction_rows <- interaction_export %>%
    dplyr::inner_join(
      interaction_lookup,
      by = c("ES", "biodiversity", "model_name")
    ) %>%
    dplyr::transmute(
      ES = as.character(ES),
      biodiversity = as.character(biodiversity),
      model_name = as.character(model_name),
      selected_model = "interaction",
      slope_scope = "Region_LU-specific",
      Region_LU = as.character(Region_LU),
      trend = as.numeric(trend),
      SE = as.numeric(SE),
      p.value = as.numeric(p.value),
      p_adjustment = "BH within model across Region x Land Use slopes",
      direction = tolower(as.character(direction)),
      selection_delta_AIC = as.numeric(selection_delta_AIC),
      selection_LRT_p = as.numeric(selection_LRT_p)
    )

  additive_rows <- selected_table %>%
    dplyr::filter(selected_model == "additive") %>%
    dplyr::transmute(
      ES = as.character(ES),
      biodiversity = as.character(biodiversity),
      model_name = as.character(model_name),
      selected_model = "additive",
      slope_scope = "common across Region_LU",
      Region_LU = NA_character_,
      trend = as.numeric(summary_est_m1_biod),
      SE = as.numeric(summary_SE_m1_biod),
      p.value = as.numeric(summary_p_m1_biod),
      p_adjustment = "Wald test of the common slope (unadjusted)",
      direction = tolower(as.character(additive_direction)),
      selection_delta_AIC = as.numeric(delta_AIC_m1_vs_m0),
      selection_LRT_p = as.numeric(p_comp_m1_vs_m0_biod)
    )

  dplyr::bind_rows(interaction_rows, additive_rows)
}

indicator_data_source <- "additive-aware workbook"

indicator_results <- eventReactive(input$main_tab, {
    req(input$main_tab %in% c("explorer_tab", "analysis_tab"))
    
    if (file.exists(indicator_file)) {
      readxl::read_excel(indicator_file)
    } else {
      build_indicator_results_from_pipeline()
      indicator_data_source <- "pipeline CSV fallback (old interaction-only workbook detected)"
    }
  })

if (
  is.null(indicator_results) ||
    !all(required_result_columns %in% names(indicator_results))
) {
  indicator_results <- build_indicator_results_from_pipeline()
}

missing_result_columns <- setdiff(
  required_result_columns,
  names(indicator_results)
)

if (length(missing_result_columns) > 0) {
  stop(
    "Indicator results are missing columns required by the app: ",
    paste(missing_result_columns, collapse = ", ")
  )
}

indicator_results <- indicator_results %>%
  dplyr::mutate(
    ES = as.character(ES),
    biodiversity = as.character(biodiversity),
    model_name = as.character(model_name),
    selected_model = tolower(as.character(selected_model)),
    slope_scope = as.character(slope_scope),
    Region_LU = dplyr::na_if(as.character(Region_LU), ""),
    p_adjustment = as.character(p_adjustment),
    direction = tolower(as.character(direction))
  )

if (!all(indicator_results$selected_model %in% c("interaction", "additive"))) {
  stop("Unexpected selected_model values in the indicator-result workbook.")
}

## Region_LU-specific intercepts are pre-computed by the modelling pipeline.
## The Shiny app reads these values; it does not recalculate the intercept table.
if (!file.exists(intercept_file)) {
  stop(
    "Missing file: ", intercept_file,
    "\nRun the modelling pipeline first to generate the intercept workbook."
  )
}

intercept_results <- tryCatch(
  readxl::read_excel(intercept_file),
  error = function(e) {
    stop("Could not read ", intercept_file, ": ", conditionMessage(e))
  }
)

missing_intercept_columns <- setdiff(
  required_intercept_columns,
  names(intercept_results)
)

if (length(missing_intercept_columns) > 0) {
  stop(
    "Intercept results are missing columns required by the app: ",
    paste(missing_intercept_columns, collapse = ", ")
  )
}

intercept_results <- intercept_results %>%
  dplyr::mutate(
    ES = as.character(ES),
    biodiversity = as.character(biodiversity),
    model_name = as.character(model_name),
    selected_model = tolower(as.character(selected_model)),
    Region_LU = as.character(Region_LU),
    biod_value_reference = as.numeric(biod_value_reference),
    intercept_link = as.numeric(intercept_link),
    intercept_SE = as.numeric(intercept_SE),
    intercept_response = as.numeric(intercept_response)
  )

if (any(is.na(intercept_results$Region_LU))) {
  stop("The intercept workbook contains missing Region_LU values.")
}

duplicate_intercepts <- intercept_results %>%
  dplyr::count(model_name, Region_LU, name = "n") %>%
  dplyr::filter(n > 1)

if (nrow(duplicate_intercepts) > 0) {
  stop("The intercept workbook contains duplicate model_name / Region_LU rows.")
}

#df_shiny <- readRDS(data_file)

eventReactive(input$main_tab, {
    req(input$main_tab %in% c("explorer_tab", "analysis_tab")) # Adjust to your actual tab values
    
    # Reads df_shiny.rds only when user leaves "Read me"
    readRDS(data_file)
  }, ignoreNULL = FALSE)

#### check iniziali sul file di caricamento 
#(all(unique(indicator_results$biodiversity) %in% names(df_shiny)))
#(all(unique(indicator_results$biodiversity) %in% names(biod_labels))) ## ok true ora
#(all(unique(indicator_results$ES) %in% names(df_shiny)))
#(all(log1p_indicators %in% names(df_shiny)))





## =========================================================
## 2. Friendly labels shown in the interface
## Internal variable names are not changed
## =========================================================
es_labels <- c(
  y_water_beta = "Water provision",
  y_provveg_beta = "Vegetation provision",
  Y_reg_co2_beta = "Carbon storage regulation",
  y_regtoc_beta = "TOC-based carbon regulation",
  y_nutcyc_beta = "Nutrient cycling",
  y_purif_beta = "Soil purification",
  y_reg_no_IC_beta = "Carbon storage regulation (without inorganic C)",
  y_multifunctionality_beta = "Multifunctionality",
  y_regsoc_beta = "SOC-based carbon regulation"
)
## in numero ora sono uguali a biod_vars che sono 47
biod_labels <- c(
  log_TPLFA = "Total microbial TPLFA biomass",
  log_BACT = "Bacterial TPLFA biomass",
  log_FUN = "Fungal TPLFA biomass",
  log_AMF = "Arbuscular mycorrhizal fungal TPLFA biomass",
  log_wormbio = "Earthworm biomass",
  log_macro = "Soil macrofauna biomass",
  log_mites = "Total mite abundance",
  log_coll = "Total Collembola abundance",
  bact_tax = "Bacterial taxonomic diversity (Shannon)",
  bact_tax_rich = "Bacterial observed taxonomic richness",
  fun_tax = "Fungal taxonomic diversity (Shannon)",
  fun_tax_rich = "Fungal observed taxonomic richness",
  worm_tax = "Earthworm taxonomic diversity (Shannon; biomass-weighted)",
  worm_tax_rich = "Earthworm taxonomic richness",
  macro_tax_rich = "Macrofauna order/group richness",
  macro_tax = "Macrofauna order/group diversity (Shannon)",
  mites_tax = "Mite group diversity (Shannon)",
  coll_tax = "Collembola taxonomic diversity (Shannon)",
  coll_emi = "Collembola ecomorphological index (EMI CWM)",
  nem_tax = "Nematode taxonomic diversity (Shannon)",
  nemMI = "Nematode Maturity Index",
  prot_tax = "Protist taxonomic diversity (Shannon; SRS-standardized)",
  funoverTOT = "Fungi / (fungi + bacteria) TPLFA ratio",
  funoverbact = "Fungi-to-bacteria TPLFA ratio",
  Mites_Oribatida_sapro = "Oribatid mite abundance",
  Mites_Gamasida_predators = "Gamasid mite abundance",
  Mites_Prostigmata_omnivores = "Prostigmatid mite abundance",
  Mites_Astigmata_fungivores = "Astigmatid mite abundance",
  Mites_Oribatida_sapro_prop = "Oribatid mites (relative abundance)",
  Mites_Gamasida_predators_prop = "Gamasid mites (relative abundance)",
  Mites_Prostigmata_omnivores_prop = "Prostigmatid mites (relative abundance)",
  Mites_Astigmata_fungivores_prop = "Astigmatid mites (relative abundance)",
  nematodes_Guild_Richness_srs_tax = "Nematode trophic guild richness (SRS-standardized)",
  nematodes_GuildTaxRich_srs_Animal_parasites = "Animal-parasitic nematode taxonomic richness",
  nematodes_GuildTaxRich_srs_Bacterivores = "Bacterivorous nematode taxonomic richness",
  nematodes_GuildTaxRich_srs_Fungivores = "Fungivorous nematode taxonomic richness",
  nematodes_GuildTaxRich_srs_Herbivores = "Plant-feeding nematode taxonomic richness",
  nematodes_GuildTaxRich_srs_Omnivores = "Omnivorous nematode taxonomic richness",
  nematodes_GuildTaxRich_srs_Predators = "Predatory nematode taxonomic richness",
  nematodes_prop_Animal_parasites = "Animal-parasitic nematodes (relative read abundance)",
  nematodes_prop_Bacterivores = "Bacterivorous nematodes (relative read abundance)",
  nematodes_prop_Fungivores = "Fungivorous nematodes (relative read abundance)",
  nematodes_prop_Herbivores = "Plant-feeding nematodes (relative read abundance)",
  nematodes_prop_Omnivores = "Omnivorous nematodes (relative read abundance)",
  nematodes_prop_Predators = "Predatory nematodes (relative read abundance)",
  nematodes_prop_Trophic_Shannon = "Nematode trophic diversity (Shannon)",
  nematodes_prop_Trophic_Richness = "Nematode trophic guild richness (non-SRS reads)"
)

## These variables were created with log1p() in the main pipeline.
## Users enter the original, untransformed measurement.
log1p_indicators <- c( ## 8 variabili qui anche
  "log_TPLFA",
  "log_BACT",
  "log_FUN",
  "log_AMF",
  "log_wormbio",
  "log_macro",
  "log_mites",
  "log_coll"
)

## Final consistency checks. These stop the app only when files are genuinely misaligned.
stopifnot(
  all(unique(indicator_results$biodiversity) %in% names(df_shiny)),
  all(unique(indicator_results$biodiversity) %in% names(biod_labels)),
  all(unique(indicator_results$ES) %in% names(df_shiny)),
  all(unique(indicator_results$ES) %in% names(es_labels)),
  all(log1p_indicators %in% names(df_shiny)),
  all(c("Region_LU", "Site_ID") %in% names(df_shiny)),
  all(!is.na(indicator_results$Region_LU[indicator_results$selected_model == "interaction"])),
  all(is.na(indicator_results$Region_LU[indicator_results$selected_model == "additive"])),
  all(unique(intercept_results$model_name) %in% unique(indicator_results$model_name)),
  all(unique(indicator_results$model_name) %in% unique(intercept_results$model_name)),
  all(unique(intercept_results$Region_LU) %in% as.character(unique(df_shiny$Region_LU))),
  all(!is.na(intercept_results$biod_value_reference) & intercept_results$biod_value_reference == 0)
)

friendly_es <- function(x) {
  out <- unname(es_labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

friendly_biod <- function(x) {
  out <- unname(biod_labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

friendly_region <- function(x) {
  gsub(".", " — ", x, fixed = TRUE)
}

uses_log1p <- function(indicator) {
  indicator %in% log1p_indicators
}

to_model_scale <- function(indicator, value) {
  value <- as.numeric(value)

  if (!is.finite(value)) {
    stop("Please enter a finite numeric value.")
  }

  if (value < 0) {
    stop("The measured value cannot be negative.")
  }

  if (uses_log1p(indicator)) {
    return(log1p(value))
  }

  value
}

from_model_scale <- function(indicator, value) {
  value <- as.numeric(value)

  if (uses_log1p(indicator)) {
    return(pmax(0, expm1(value)))
  }

  value
}

## =========================================================
## 3. Load the exact mixed models produced by the pipeline
##    and refit only if the bundle is unavailable
## =========================================================
## =========================================================
## 3. Load exact mixed models produced by the pipeline
##    No automatic refitting is allowed in the distributed app
## =========================================================

selected_model_names <- unique(indicator_results$model_name)
selected_model_names <- selected_model_names[!is.na(selected_model_names)]

if (!file.exists(bundle_file)) {
  stop(
    "Missing required model bundle: ", bundle_file,
    "\nPlease use the complete BIOservicES toolkit folder.",
    call. = FALSE
  )
}

#loaded_bundle <- readRDS(bundle_file)

selected_fits <- eventReactive(input$main_tab, {
    req(input$main_tab == "analysis_tab") # Adjust tab value
    
    showNotification("Loading models, please wait...", type = "message", id = "model_load_notif")
    
    loaded_bundle <- readRDS(bundle_file)
    
    removeNotification(id = "model_load_notif")
    return(loaded_bundle)
  })

missing_model_names <- setdiff(
  selected_model_names,
  names(loaded_bundle)
)

if (length(missing_model_names) > 0) {
  stop(
    "The model bundle is incomplete.\n\n",
    "Missing models:\n- ",
    paste(missing_model_names, collapse = "\n- "),
    "\n\nPlease use the complete model bundle generated by the final modelling pipeline.",
    call. = FALSE
  )
}

selected_fits <- loaded_bundle[selected_model_names]

model_loaded <- !vapply(
  selected_fits,
  is.null,
  logical(1)
)

if (!all(model_loaded)) {
  stop(
    "At least one model in the bundle is NULL or invalid.",
    call. = FALSE
  )
}

model_source <- "loaded from final pipeline bundle"


## =========================================================
## DA QUI IN POI RIPRENDE IL CODICE ORIGINALE
## NON CANCELLARE
## =========================================================

## Population-level predictions set the Site_ID random intercept to zero.
make_population_newdata <- function(fit, region_name, biod_values) {
  n_values <- length(biod_values)
  
  data.frame(
    Region_LU = factor(
      rep(region_name, n_values),
      levels = levels(fit$data$Region_LU)
    ),
    Site_ID = factor(
      rep(NA_character_, n_values),
      levels = levels(fit$data$Site_ID)
    ),
    biod_value = as.numeric(biod_values)
  )
}

predict_population <- function(fit, newdata) {
  as.numeric(
    stats::predict(
      fit$model,
      newdata = newdata,
      type = "response",
      re.form = NA,
      allow.new.levels = TRUE
    )
  )
}

es_values <- unique(indicator_results$ES)
es_choices <- stats::setNames(es_values, friendly_es(es_values))
relationship_scope_label <- function(selected_model, region_name = NULL) {
  ifelse(
    selected_model == "additive",
    "Common biodiversity slope across Region × Land Use contexts",
    ifelse(
      is.null(region_name) || is.na(region_name),
      "Region × Land Use-specific biodiversity slope",
      paste0("Slope specific to ", friendly_region(region_name))
    )
  )
}

model_regions <- function(model_name) {
  fit <- selected_fits[[model_name]]

  if (is.null(fit) || is.null(fit$data$Region_LU)) {
    return(character(0))
  }

  as.character(levels(fit$data$Region_LU))
}

model_supports_region <- function(model_name, region_name) {
  region_name %in% model_regions(model_name)
}

available_regions_for_es <- function(es_name) {
  x <- indicator_results %>%
    dplyr::filter(ES == es_name)

  regions <- x %>%
    dplyr::filter(selected_model == "interaction", !is.na(Region_LU)) %>%
    dplyr::pull(Region_LU) %>%
    as.character()

  additive_models <- x %>%
    dplyr::filter(selected_model == "additive") %>%
    dplyr::pull(model_name) %>%
    unique()

  for (model_name in additive_models) {
    regions <- union(regions, model_regions(model_name))
  }

  sort(unique(regions))
}

indicator_choice_table <- function(es_name, region_name = NULL) {
  x <- indicator_results %>%
    dplyr::filter(ES == es_name)

  if (!is.null(region_name) && nrow(x) > 0) {
    keep <- vapply(
      seq_len(nrow(x)),
      function(i) {
        if (x$selected_model[[i]] == "interaction") {
          return(!is.na(x$Region_LU[[i]]) && x$Region_LU[[i]] == region_name)
        }

        model_supports_region(x$model_name[[i]], region_name)
      },
      logical(1)
    )

    x <- x[keep, , drop = FALSE]
  }

  x %>%
    dplyr::distinct(biodiversity, selected_model) %>%
    dplyr::mutate(
      choice_label = paste0(
        friendly_biod(biodiversity),
        ifelse(
          selected_model == "additive",
          " — common slope",
          " — context-specific slopes"
        )
      )
    ) %>%
    dplyr::arrange(choice_label)
}


## =========================================================
## README content tables
## These are descriptive only and do not affect model fitting.
## =========================================================
readme_es_table <- data.frame(
  `Ecosystem service` = c(
    "Vegetation provision",
    "Water provision",
    "Carbon storage regulation",
    "Carbon storage regulation (without inorganic C)",
    "TOC-based carbon regulation",
    "SOC-based carbon regulation",
    "Nutrient cycling",
    "Soil purification",
    "Multifunctionality"
  ),
  `Construction in the toolkit` = c(
    "Mean of vegetation cover and vegetation richness SEF scores.",
    "Available water capacity (AWC) SEF score.",
    "Mean of SOC, inorganic carbon (IC), and SOC:clay ratio SEF scores.",
    "Mean of SOC and SOC:clay ratio SEF scores.",
    "TOC SEF score.",
    "SOC SEF score.",
    "Mean of the standardized scores for total N, available P, exchangeable Mg and K, CEC, bioavailable B, Cu, Fe and Mn, together with reversed scores for exchangeable Na and electrical conductivity (EC)",
    "Mean of inverse standardized contamination scores for Zn, As, Co, bioavailable Co, Cr, Pb, bioavailable Ni, bioavailable Sb, and total pesticides; higher values therefore correspond to lower measured contamination.",
    "Mean of soil purification, nutrient cycling, carbon storage regulation, water provision, and vegetation provision."
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

readme_biod_table <- data.frame(
  `Biodiversity component` = c(
    "Microbial community",
    "Earthworms",
    "Macrofauna",
    "Mites",
    "Collembola",
    "Nematodes",
    "Protists"
  ),
  `Indicators represented in the toolkit` = c(
    "Total, bacterial, fungal and AMF TPLFA biomass; bacterial and fungal Shannon diversity and observed richness; fungi:bacteria composition ratios.",
    "Total biomass, taxonomic richness and Shannon diversity.",
    "Total biomass, order/group richness and order/group Shannon diversity.",
    "Total abundance, group Shannon diversity, and absolute or relative abundance of Oribatida, Gamasida, Prostigmata and Astigmata functional groups.",
    "Total abundance, taxonomic Shannon diversity and the ecomorphological EMI community-weighted mean (CWM).",
    "Taxonomic Shannon diversity, Maturity Index, trophic-guild richness, guild-specific taxonomic richness, relative read abundance of trophic guilds, and trophic Shannon/richness metrics. SRS-standardized measures are identified in their labels.",
    "Taxonomic Shannon diversity; the toolkit label identifies the SRS-standardized measure."
  ),
  `Model input` = c(
    "TPLFA biomass variables are log1p-transformed; diversity, richness and ratio variables are used on their prepared scale.",
    "Biomass is log1p-transformed; richness/diversity measures are used on their prepared scale.",
    "Total biomass is log1p-transformed; richness/diversity measures are used on their prepared scale.",
    "Total abundance is log1p-transformed; group composition/diversity measures are used on their prepared scale.",
    "Total abundance is log1p-transformed; Shannon and EMI CWM are used on their prepared scale.",
    "Prepared taxonomic and functional metrics are used directly on their analysis scale.",
    "Prepared diversity metric is used directly on its analysis scale."
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## =========================================================
## 4. User interface
## =========================================================
ui <- shiny::fluidPage(

  shiny::tags$head(
    shiny::tags$style(
      shiny::HTML("
        .direction-positive {
          color: #18794e;
          font-weight: 600;
        }

        .direction-negative {
          color: #b42318;
          font-weight: 600;
        }

        .direction-unclear {
          color: #5f6368;
          font-weight: 600;
        }

        .summary-box {
          background: #f5f7f8;
          border-left: 5px solid #3c6e71;
          padding: 12px 16px;
          margin-bottom: 18px;
          border-radius: 4px;
        }

        .help-note {
          color: #5f6368;
          font-size: 0.95em;
        }

        .plot-note {
          margin-top: 8px;
          color: #5f6368;
          font-size: 0.9em;
        }

        .traffic-card {
          padding: 18px 20px;
          margin: 12px 0 18px 0;
          border-radius: 8px;
          border-left: 8px solid;
        }

        .traffic-red {
          background: #fde8e7;
          border-left-color: #b42318;
        }

        .traffic-yellow {
          background: #fff4bf;
          border-left-color: #d4a017;
        }

        .traffic-green {
          background: #e4f3e8;
          border-left-color: #18794e;
        }

        .traffic-dot {
          width: 22px;
          height: 22px;
          border-radius: 50%;
          display: inline-block;
          margin-right: 10px;
          vertical-align: -3px;
        }

        .dot-red {
          background: #d92d20;
        }

        .dot-yellow {
          background: #f2c94c;
          border: 1px solid #9a7b00;
        }

        .dot-green {
          background: #219653;
        }

        .prediction-score {
          font-size: 2.1em;
          font-weight: 700;
          margin: 8px 0;
        }

        .warning-box {
          background: #fff8e1;
          border-left: 5px solid #d4a017;
          padding: 10px 14px;
          margin-top: 12px;
          border-radius: 4px;
        }

        .expected-es-details {
          margin-top: 14px;
          border: 1px solid #d9e2e8;
          border-radius: 6px;
          background: #fafcfd;
        }

        .expected-es-details summary {
          cursor: pointer;
          font-weight: 600;
          padding: 10px 12px;
        }

        .expected-es-details-body {
          padding: 0 12px 10px 12px;
        }

        .readme-section {
          margin-bottom: 28px;
        }

        .readme-card {
          background: #f7f9fa;
          border: 1px solid #dfe5e8;
          border-radius: 7px;
          padding: 16px 18px;
          margin: 10px 0 18px 0;
        }

        .readme-card h4 {
          margin-top: 0;
        }

        .download-card {
          background: #f5f7f8;
          border-left: 5px solid #3c6e71;
          padding: 14px 16px;
          margin: 12px 0;
          border-radius: 4px;
        }

        .readme-small {
          color: #5f6368;
          font-size: 0.92em;
        }

        .tab-content {
          padding-top: 18px;
        }
      ")
    )
  ),

  shiny::titlePanel("BIOservicES biodiversity indicator explorer"),

  shiny::tabsetPanel(
    id = "main_tab",

    shiny::tabPanel(
      title = "Read me",
      value = "readme",

      shiny::fluidRow(
        shiny::column(
          width = 10,
          offset = 1,

          shiny::div(
            class = "readme-section",
            shiny::h2("About this toolkit"),
            shiny::p(
              "The BIOservicES biodiversity indicator explorer, links prepared soil-biodiversity indicators to ecosystem-service indicators across Region × Land Use contexts. Its purpose is to help users understand which biodiversity measurements are associated with a selected ecosystem service, where those relationships are supported, and how the fitted mixed models translate a biodiversity measurement into an expected ecosystem-service condition."
            ),
            shiny::p(
              "BIOservicES investigates relationships among soil biodiversity, ecosystem functioning, ecosystem services and land use across European biogeographical contexts. For the broader project background, objectives and consortium, see the "
            ),
            shiny::tags$a(
              href = "https://bioservices-project.eu/",
              target = "_blank",
              rel = "noopener noreferrer",
              "official BIOservicES project website"
            ),
            shiny::p(
              class = "readme-small",
              "The current toolkit reference dataset excludes the Boreal region; all scaling, model selection and displayed reference distributions therefore refer to the analysis dataset used by this application."
            )
          ),

          shiny::div(
            class = "readme-section",
            shiny::h2("How ecosystem-service indicators were constructed"),
            shiny::div(
              class = "readme-card",
              shiny::h4("Standardized ecosystem-function scores (SEF)"),
              shiny::p(
                "Individual ecosystem-function variables were transformed to a common 0–1 scale using a robust Min–Max procedure based on the 5th and 95th percentiles of the analysis dataset. Values below the 5th percentile were clamped to 0 and values above the 95th percentile to 1."
              ),
              shiny::p(
                "For variables where higher values were interpreted as greater ecosystem-service delivery, the standardized score increases with the raw measurement. For contamination-related variables and other lower-is-better measures, the score was inverted so that a higher SEF always represents the preferred direction."
              ),
              shiny::p(
                class = "readme-small",
                "These SEF scores are relative to the dataset used to calculate the quantiles. They should not be interpreted as universal ecological thresholds or fixed optimum values."
              )
            ),
            shiny::tableOutput("readme_es_table"),
            shiny::p(
              class = "readme-small",
              "For beta mixed modelling, the 0–1 ecosystem-service responses are slightly adjusted away from exact 0 and 1 so they lie inside the support of the beta distribution. This numerical adjustment does not redefine the underlying ecosystem-service indicator."
            )
          ),

          shiny::div(
            class = "readme-section",
            shiny::h2("How biodiversity indicators were constructed"),
            shiny::p(
              "The toolkit does not use one single composite biodiversity index. Instead, it tests complementary taxonomic, abundance/biomass and functional indicators prepared for the different soil-biota groups. This allows the model-selection pipeline to identify which biodiversity facet is informative for each ecosystem service."
            ),
            shiny::tableOutput("readme_biod_table"),
            shiny::div(
              class = "readme-card",
              shiny::h4("Transformation used in the mixed models"),
              shiny::p(
                "Eight abundance/biomass indicators are log1p-transformed before modelling: total microbial TPLFA biomass, bacterial TPLFA biomass, fungal TPLFA biomass, AMF TPLFA biomass, earthworm biomass, macrofauna biomass, total mite abundance and total Collembola abundance. The interface shows and accepts these measurements on their original scale and applies the transformation internally where required."
              ),
              shiny::p(
                "Other diversity, richness, compositional and functional metrics are used on the prepared analysis scale indicated by their names and definitions. DNA-derived nematode and protist metrics are identified explicitly where SRS standardization was used."
              )
            )
          ),

          shiny::div(
            class = "readme-section",
            shiny::h2("How to interpret the model output"),
            shiny::tags$ul(
              shiny::tags$li(
                shiny::strong("Additive model: "),
                "one biodiversity slope is shared across Region × Land Use contexts, while the model retains context-specific intercepts."
              ),
              shiny::tags$li(
                shiny::strong("Interaction model: "),
                "the biodiversity–ecosystem service relationship can differ among ",
                "Region × Land Use contexts. Only contexts with a statistically supported ",
                "biodiversity slope after BH correction (adjusted p < 0.01) are shown."
              ),
              
              shiny::tags$li(
                shiny::strong("Population-level predictions: "),
                "displayed predictions set the Site_ID random intercept to zero."
              ),
              shiny::tags$li(
                shiny::strong("Intercept estimates: "),
                "the toolkit reports the Region × Land Use baseline at biodiversity = 0 on the conditional logit/link scale, its standard error, and the corresponding expected ecosystem-service value on the 0–1 response scale. The Site_ID random intercept is set to zero."
              ),
              shiny::tags$li(
                shiny::strong("Association, not causation: "),
                "supported biodiversity–ecosystem-service relationships should not automatically be interpreted as causal effects."
              ),
              shiny::tags$li(
                shiny::strong("Traffic-light classes: "),
                "are based on the 33rd and 66th percentiles of the observed ecosystem-service distribution within the selected Region × Land Use context. These are relative reference categories, not validated ecological thresholds."
              )
            )
          ),

          shiny::div(
            class = "readme-section",
            shiny::h2("Download model results"),
            shiny::p(
              "The two Excel files below are the result tables distributed with the toolkit. They can be downloaded directly without modifying their contents."
            ),

            shiny::div(
              class = "download-card",
              shiny::h4("Biodiversity–ecosystem service model results"),
              shiny::p(
                "Final exported biodiversity–ES relationships, including additive common slopes and Region × Land Use-specific interaction slopes."
              ),
              shinylive_downloadButton(
                outputId = "download_model_results",
                label = "Download model results (.xlsx)",
                class = "btn-primary"
              )
            ),

            shiny::div(
              class = "download-card",
              shiny::h4("Region × Land Use expected ecosystem-service values"),
              shiny::p(
                "Model-derived Region × Land Use values evaluated at biodiversity = 0, including the link-scale estimates and their response-scale transformation where available."
              ),
              shinylive_downloadButton(
                outputId = "download_intercept_results",
                label = "Download Region × Land Use results (.xlsx)"
              )
            )
          ),

          shiny::div(
            class = "readme-section",
            
            shiny::h2("Quick references"),
            
            shiny::p(
              paste0(
                "The selection of soil biodiversity indicators and the ecosystem-service ",
                "multifunctionality framework used in this toolkit are conceptually related ",
                "to approaches developed in biodiversity–ecosystem functioning and soil ",
                "biodiversity monitoring studies, including:"
              )
            ),
            
            shiny::tags$ul(
              
              shiny::tags$li(
                shiny::strong("Byrnes et al. (2014; online 2013). "),
                "Investigating the relationship between biodiversity and ecosystem ",
                "multifunctionality: challenges and solutions. ",
                shiny::em("Methods in Ecology and Evolution"), ", 5, 111–124. ",
                shiny::tags$a(
                  href = "https://doi.org/10.1111/2041-210X.12143",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  "https://doi.org/10.1111/2041-210X.12143"
                )
              ),
              
              shiny::tags$li(
                shiny::strong("Griffiths et al. (2016). "),
                "Selecting cost effective and policy-relevant biological indicators for ",
                "European monitoring of soil biodiversity and ecosystem function. ",
                shiny::em("Ecological Indicators"), ", 69, 213–223. ",
                shiny::tags$a(
                  href = "https://doi.org/10.1016/j.ecolind.2016.04.023",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  "https://doi.org/10.1016/j.ecolind.2016.04.023"
                )
              ),
              
              shiny::tags$li(
                shiny::strong("Li et al. (2024). "),
                "Plant diversity enhances ecosystem multifunctionality via multitrophic diversity. ",
                shiny::em("Nature Ecology & Evolution"), ", 8, 2037–2047. ",
                shiny::tags$a(
                  href = "https://doi.org/10.1038/s41559-024-02517-2",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  "https://doi.org/10.1038/s41559-024-02517-2"
                )
              ),
              
              shiny::tags$li(
                shiny::strong("Manning et al. (2018). "),
                "Redefining ecosystem multifunctionality. ",
                shiny::em("Nature Ecology & Evolution"), ", 2, 427–436. ",
                shiny::tags$a(
                  href = "https://doi.org/10.1038/s41559-017-0461-7",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  "https://doi.org/10.1038/s41559-017-0461-7"
                )
              ),
              
              shiny::tags$li(
                shiny::strong("Schuldt et al. (2018). "),
                "Biodiversity across trophic levels drives multifunctionality in highly diverse forests. ",
                shiny::em("Nature Communications"), ", 9, 2989. ",
                shiny::tags$a(
                  href = "https://doi.org/10.1038/s41467-018-05421-z",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  "https://doi.org/10.1038/s41467-018-05421-z"
                )
              )
            ),
          
            shiny::p(
              class = "readme-small",
              "This Read me is a concise description of the current toolkit implementation. Detailed sampling, laboratory and bioinformatic methods remain part of the BIOservicES project documentation and the corresponding analytical workflows."
            )
          )
        )
      )
    ),

    shiny::tabPanel(
      title = "What should I measure?",
      value = "context",

      shiny::sidebarLayout(
        shiny::sidebarPanel(
          width = 4,

          shiny::selectInput(
            inputId = "context_es",
            label = "Which ecosystem service?",
            choices = es_choices
          ),

          shiny::uiOutput("context_region_selector")
        ),

        shiny::mainPanel(
          width = 8,

          shiny::h3("Candidate biodiversity indicators"),
          shiny::p(
            "Indicators supported for the selected Region × Land Use context.",
            "Additive relationships have one common slope across contexts; interaction relationships are context-specific."
          ),

          shiny::uiOutput("context_summary"),
          shiny::tableOutput("context_indicator_table")
        )
      )
    ),

    shiny::tabPanel(
      title = "Where is my indicator relevant?",
      value = "indicator",

      shiny::sidebarLayout(
        shiny::sidebarPanel(
          width = 4,

          shiny::selectInput(
            inputId = "indicator_es",
            label = "Which ecosystem service?",
            choices = es_choices
          ),

          shiny::uiOutput("indicator_selector")
        ),

        shiny::mainPanel(
          width = 8,

          shiny::h3("Where is this indicator relevant?"),
          shiny::p(
            "Scatter plot of observed data and estimated regression relationships.",
            "Common-slope models show all represented Region × Land Use contexts together."
          ),

          shiny::uiOutput("indicator_summary"),
          shiny::plotOutput("indicator_relationship_plot", height = "520px"),

          shiny::uiOutput("indicator_plot_note"),
          shiny::uiOutput("indicator_table_heading"),
          shiny::tableOutput("indicator_region_table"),
          shiny::uiOutput("indicator_additive_expected_es_ui")
        )
      )
    ),

    shiny::tabPanel(
      title = "Predict condition",
      value = "prediction",

      shiny::sidebarLayout(
        shiny::sidebarPanel(
          width = 4,

          shiny::selectInput(
            inputId = "prediction_es",
            label = "Which ecosystem service do you want to predict?",
            choices = es_choices
          ),

          shiny::uiOutput("prediction_region_selector"),
          shiny::uiOutput("prediction_indicator_selector"),
          shiny::uiOutput("prediction_value_input"),

          shiny::actionButton(
            inputId = "run_prediction",
            label = "Run prediction",
            class = "btn-primary",
            width = "100%"
          ),

          shiny::tags$br(),
          shiny::tags$br(),

          shiny::p(
            class = "help-note",
            "One biodiversity variable is used for each prediction.",
            "Additive models use a common biodiversity slope but retain context-specific intercepts.",
            "Predictions from different variables are not averaged."
          )
        ),

        shiny::mainPanel(
          width = 8,

          shiny::h3("Predicted ecosystem-service condition"),
          shiny::p(
            "The traffic-light category is relative to the observed distribution",
            "for the selected Region × Land Use context."
          ),

          shiny::uiOutput("prediction_result_card"),
          shiny::plotOutput("prediction_plot", height = "500px"),
          shiny::uiOutput("prediction_warning"),

          shiny::h4("Context-specific reference thresholds"),
          shiny::tableOutput("prediction_threshold_table"),

          shiny::p(
            class = "plot-note",
            "Red: at or below the 33rd percentile. ",
            "Yellow: between the 33rd and 66th percentiles. ",
            "Green: above the 66th percentile. ",
            "The 99th percentile is shown as an upper reference limit."
          )
        )
      )
    )

  ),

  shiny::tags$hr(),
  shiny::p(
    class = "help-note",
    "These are model-based associations and predictions. They do not necessarily imply causality.",
    "Displayed curves are population-level predictions from beta mixed models with a Site_ID random intercept.",
    "Traffic-light classes are relative reference categories rather than validated ecological thresholds."
  ),
  shiny::textOutput("app_status")
)

## =========================================================
## 5. Server logic
## =========================================================
server <- function(input, output, session) {

  ## -------------------------------------------------------
  ## Read me: descriptive tables and downloads
  ## -------------------------------------------------------
  output$readme_es_table <- shiny::renderTable({
    readme_es_table
  },
  striped = TRUE,
  bordered = TRUE,
  spacing = "s")

  output$readme_biod_table <- shiny::renderTable({
    readme_biod_table
  },
  striped = TRUE,
  bordered = TRUE,
  spacing = "s")

  output$download_model_results <- shiny::downloadHandler(
    filename = function() {
      basename(indicator_file)
    },
    content = function(file) {
      if (!file.exists(indicator_file)) {
        stop("Model-results workbook not found: ", indicator_file)
      }
      copied <- file.copy(indicator_file, file, overwrite = TRUE)
      if (!isTRUE(copied)) {
        stop("Could not copy the model-results workbook for download.")
      }
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )

  output$download_intercept_results <- shiny::downloadHandler(
    filename = function() {
      basename(intercept_file)
    },
    content = function(file) {
      if (!file.exists(intercept_file)) {
        stop("Region x Land Use workbook not found: ", intercept_file)
      }
      copied <- file.copy(intercept_file, file, overwrite = TRUE)
      if (!isTRUE(copied)) {
        stop("Could not copy the Region x Land Use workbook for download.")
      }
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )

  ## -------------------------------------------------------
  ## Context-first pathway
  ## -------------------------------------------------------
  output$context_region_selector <- shiny::renderUI({
    shiny::req(input$context_es)

    region_values <- available_regions_for_es(input$context_es)

    region_choices <- stats::setNames(
      region_values,
      friendly_region(region_values)
    )

    shiny::selectInput(
      inputId = "context_region",
      label = "Which Region × Land Use context?",
      choices = region_choices
    )
  })

  context_results <- shiny::reactive({
    shiny::req(input$context_es, input$context_region)

    x <- indicator_results %>%
      dplyr::filter(ES == input$context_es)

    keep <- vapply(
      seq_len(nrow(x)),
      function(i) {
        if (x$selected_model[[i]] == "interaction") {
          return(
            !is.na(x$Region_LU[[i]]) &&
              x$Region_LU[[i]] == input$context_region
          )
        }

        model_supports_region(
          x$model_name[[i]],
          input$context_region
        )
      },
      logical(1)
    )

    out <- x[keep, , drop = FALSE] %>%
      dplyr::mutate(
        relationship_scope = ifelse(
          selected_model == "additive",
          "Common slope across contexts",
          "Slope specific to selected context"
        ),
        intercept_region = input$context_region
      ) %>%
      dplyr::left_join(
        intercept_results %>%
          dplyr::select(
            model_name,
            intercept_region = Region_LU,
            intercept_link,
            intercept_SE,
            intercept_response
          ),
        by = c("model_name", "intercept_region")
      )

    shiny::validate(
      shiny::need(
        all(!is.na(out$intercept_response)),
        "Intercept estimates are missing for at least one displayed model/context."
      )
    )

    out %>%
      dplyr::arrange(selected_model, p.value)
  })

  output$context_summary <- shiny::renderUI({
    x <- context_results()
    n_indicators <- nrow(x)
    n_additive <- sum(x$selected_model == "additive", na.rm = TRUE)
    n_interaction <- sum(x$selected_model == "interaction", na.rm = TRUE)
    n_positive <- sum(x$direction == "positive", na.rm = TRUE)
    n_negative <- sum(x$direction == "negative", na.rm = TRUE)

    direction_text <- dplyr::case_when(
      n_positive > 0 && n_negative == 0 ~ "All identified associations are positive.",
      n_negative > 0 && n_positive == 0 ~ "All identified associations are negative.",
      n_positive > 0 && n_negative > 0 ~ "Both positive and negative associations were identified.",
      TRUE ~ ""
    )

    model_text <- paste0(
      n_interaction, " context-specific",
      if (n_additive > 0) paste0(" and ", n_additive, " common-slope") else "",
      ifelse(n_indicators == 1, " relationship", " relationships"),
      "."
    )

    shiny::div(
      class = "summary-box",
      shiny::strong(
        paste0(
          n_indicators,
          ifelse(n_indicators == 1, " candidate indicator", " candidate indicators"),
          " identified for ",
          friendly_region(input$context_region),
          "."
        )
      ),
      shiny::br(),
      model_text,
      shiny::br(),
      direction_text
    )
  })

  output$context_indicator_table <- shiny::renderTable({
    context_results() %>%
      dplyr::transmute(
        `Biodiversity indicator` = friendly_biod(biodiversity),
        `Relationship scope` = relationship_scope,
        Association = dplyr::case_when(
          direction == "positive" ~
            "<span class='direction-positive'>&uarr; Positive</span>",
          direction == "negative" ~
            "<span class='direction-negative'>&darr; Negative</span>",
          TRUE ~
            "<span class='direction-unclear'>&mdash; Unclear</span>"
        ),
        `Estimated slope (logit/link scale)` = round(trend, 3),
        `Slope standard error` = round(SE, 3),
        `Intercept estimate (logit/link scale)` = round(intercept_link, 3),
        `Intercept standard error` = round(intercept_SE, 3),
        `Expected ES at biodiversity = 0` = round(intercept_response, 3),
        `Selection Delta AIC` = round(selection_delta_AIC, 2),
        `Selection LRT p-value` = format.pval(
          selection_LRT_p,
          digits = 3,
          eps = 0.001
        ),
        `Reported slope p-value` = format.pval(
          p.value,
          digits = 3,
          eps = 0.001
        ),
        `P-value method` = p_adjustment
      )
  },
  striped = TRUE,
  bordered = TRUE,
  spacing = "s",
  sanitize.text.function = function(x) x)

  ## -------------------------------------------------------
  ## Indicator-first pathway
  ## -------------------------------------------------------
  output$indicator_selector <- shiny::renderUI({
    shiny::req(input$indicator_es)

    choice_table <- indicator_choice_table(input$indicator_es)

    shiny::selectInput(
      inputId = "selected_indicator",
      label = "Which biodiversity indicator have you measured?",
      choices = stats::setNames(
        choice_table$biodiversity,
        choice_table$choice_label
      )
    )
  })

  selected_indicator_results <- shiny::reactive({
    shiny::req(input$indicator_es, input$selected_indicator)

    indicator_results %>%
      dplyr::filter(
        ES == input$indicator_es,
        biodiversity == input$selected_indicator
      ) %>%
      dplyr::arrange(p.value)
  })

  selected_indicator_fit <- shiny::reactive({
    x <- selected_indicator_results()
    shiny::req(nrow(x) > 0)

    model_names <- unique(x$model_name)

    shiny::validate(
      shiny::need(
        length(model_names) == 1,
        "The selected ecosystem service / indicator pair maps to more than one model."
      )
    )

    model_name <- model_names[[1]]
    fit <- selected_fits[[model_name]]

    shiny::validate(
      shiny::need(
        !is.null(fit),
        paste("The selected model could not be fitted:", model_name)
      )
    )

    fit
  })

  selected_indicator_type <- shiny::reactive({
    x <- selected_indicator_results()
    shiny::req(nrow(x) > 0)
    type_values <- unique(x$selected_model)

    shiny::validate(
      shiny::need(
        length(type_values) == 1,
        "The selected indicator is associated with conflicting model types."
      )
    )

    type_values[[1]]
  })

  output$indicator_summary <- shiny::renderUI({
    x <- selected_indicator_results()
    selected_type <- selected_indicator_type()

    if (selected_type == "additive") {
      direction_text <- dplyr::case_when(
        x$direction[[1]] == "positive" ~ "The common slope is positive.",
        x$direction[[1]] == "negative" ~ "The common slope is negative.",
        TRUE ~ "The direction of the common slope is unclear."
      )

      return(
        shiny::div(
          class = "summary-box",
          shiny::strong(
            paste0(
              friendly_biod(input$selected_indicator),
              " has a supported additive association with ",
              tolower(friendly_es(input$indicator_es)),
              "."
            )
          ),
          shiny::br(),
          "The biodiversity slope is common across all Region × Land Use contexts; the model retains context-specific intercepts.",
          shiny::br(),
          paste0(
            "Model selection: Delta AIC = ",
            round(x$selection_delta_AIC[[1]], 2),
            "; LRT p = ",
            format.pval(x$selection_LRT_p[[1]], digits = 3, eps = 0.001),
            "."
          ),
          shiny::br(),
          direction_text
        )
      )
    }

    n_regions <- nrow(x)
    n_positive <- sum(x$direction == "positive", na.rm = TRUE)
    n_negative <- sum(x$direction == "negative", na.rm = TRUE)

    association_text <- paste0(
      n_positive, " positive",
      if (n_negative > 0) paste0(" and ", n_negative, " negative") else "",
      if (n_regions == 1) " association." else " associations."
    )

    shiny::div(
      class = "summary-box",
      shiny::strong(
        paste0(
          friendly_biod(input$selected_indicator),
          " is significantly associated with ",
          tolower(friendly_es(input$indicator_es)),
          " in ", n_regions,
          ifelse(n_regions == 1, " context.", " contexts.")
        )
      ),
      shiny::br(),
      association_text
    )
  })

  plot_objects <- shiny::reactive({
    x <- selected_indicator_results()
    fit <- selected_indicator_fit()
    selected_type <- selected_indicator_type()

    display_regions <- if (selected_type == "additive") {
      as.character(levels(fit$data$Region_LU))
    } else {
      sort(unique(as.character(x$Region_LU)))
    }

    display_regions <- display_regions[
      !is.na(display_regions) & nzchar(display_regions)
    ]

    observed <- fit$data %>%
      dplyr::filter(as.character(Region_LU) %in% display_regions) %>%
      dplyr::mutate(
        Region_label = friendly_region(as.character(Region_LU)),
        biod_display = from_model_scale(
          input$selected_indicator,
          biod_value
        )
      )

    prediction_list <- lapply(display_regions, function(region_name) {
      region_data <- fit$data %>%
        dplyr::filter(as.character(Region_LU) == region_name)

      raw_values <- from_model_scale(
        input$selected_indicator,
        region_data$biod_value
      )

      raw_values <- raw_values[is.finite(raw_values)]

      if (nrow(region_data) < 2 || length(raw_values) < 2) {
        return(NULL)
      }

      raw_range <- range(raw_values, na.rm = TRUE)

      if (!all(is.finite(raw_range)) || diff(raw_range) == 0) {
        return(NULL)
      }

      grid_raw <- seq(
        from = raw_range[1],
        to = raw_range[2],
        length.out = 150
      )

      model_grid <- vapply(
        grid_raw,
        function(v) to_model_scale(input$selected_indicator, v),
        numeric(1)
      )

      new_data <- make_population_newdata(
        fit = fit,
        region_name = region_name,
        biod_values = model_grid
      )
      new_data$biod_display <- grid_raw
      new_data$predicted <- predict_population(fit, new_data)
      new_data$Region_label <- friendly_region(region_name)
      new_data
    })

    list(
      observed = observed,
      predicted = dplyr::bind_rows(prediction_list),
      selected_type = selected_type,
      display_regions = display_regions
    )
  })

  output$indicator_relationship_plot <- shiny::renderPlot({
    plot_data <- plot_objects()

    shiny::validate(
      shiny::need(
        nrow(plot_data$observed) > 0,
        "No observations are available for the selected indicator."
      ),
      shiny::need(
        nrow(plot_data$predicted) > 0,
        "Predictions could not be generated for the selected indicator."
      )
    )

    subtitle_text <- if (plot_data$selected_type == "additive") {
      paste0(
        "One common biodiversity coefficient on the logit scale; ",
        "context-specific intercepts; all represented contexts shown together",
        if (uses_log1p(input$selected_indicator)) {
          "; predictor shown on the original scale (model fitted using log1p values)"
        } else {
          ""
        }
      )
    } else if (uses_log1p(input$selected_indicator)) {
      "Predictor shown on the original scale; model fitted with log1p-transformed values"
    } else {
      "Observed values and population-level mixed-model predictions for significant contexts"
    }

    legend_rows <- if (plot_data$selected_type == "additive") 4 else 2

    ggplot2::ggplot() +
      ggplot2::geom_point(
        data = plot_data$observed,
        ggplot2::aes(
          x = biod_display,
          y = ES_beta,
          color = Region_label
        ),
        alpha = 0.35,
        size = 2
      ) +
      ggplot2::geom_line(
        data = plot_data$predicted,
        ggplot2::aes(
          x = biod_display,
          y = predicted,
          color = Region_label,
          group = Region_label
        ),
        linewidth = 1.2
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste0(
          friendly_biod(input$selected_indicator),
          " and ",
          friendly_es(input$indicator_es)
        ),
        subtitle = subtitle_text,
        x = if (uses_log1p(input$selected_indicator)) {
          paste0(
            friendly_biod(input$selected_indicator),
            " — original measurement scale"
          )
        } else {
          friendly_biod(input$selected_indicator)
        },
        y = friendly_es(input$indicator_es),
        color = "Region × Land Use"
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(nrow = legend_rows, byrow = TRUE)
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold")
      )
  })

  output$indicator_plot_note <- shiny::renderUI({
    if (selected_indicator_type() == "additive") {
      shiny::p(
        class = "plot-note",
        "Points are observed values; each line is a population-level mixed-model prediction for one Region × Land Use context. ",
        "All represented contexts are shown together and each curve is restricted to the predictor range observed in that context. ",
        "The biodiversity coefficient is common on the conditional logit scale, while Region × Land Use changes the intercept. ",
        "Because predictions are displayed on the 0–1 response scale, the curves do not have to look perfectly parallel. ",
        "The Site_ID random intercept is set to zero."
      )
    } else {
      shiny::p(
        class = "plot-note",
        "Points are observed values; lines are population-level mixed-model predictions. ",
        "The Site_ID random intercept is set to zero for the displayed curves. ",
        "Only Region × Land Use slopes passing the final BH-adjusted threshold are displayed."
      )
    }
  })

  output$indicator_table_heading <- shiny::renderUI({
    if (selected_indicator_type() == "additive") {
      shiny::h4("Common slope result")
    } else {
      shiny::h4("Significant Region × Land Use contexts")
    }
  })

  output$indicator_region_table <- shiny::renderTable({
    x <- selected_indicator_results()

    if (selected_indicator_type() == "additive") {
      return(
        x %>%
          dplyr::transmute(
            Scope = "All Region × Land Use contexts",
            Association = dplyr::case_when(
              direction == "positive" ~
                "<span class='direction-positive'>&uarr; Positive</span>",
              direction == "negative" ~
                "<span class='direction-negative'>&darr; Negative</span>",
              TRUE ~
                "<span class='direction-unclear'>&mdash; Unclear</span>"
            ),
            `Estimated slope (logit/link scale)` = round(trend, 3),
            `Standard error` = round(SE, 3),
            `Selection Delta AIC` = round(selection_delta_AIC, 2),
            `Selection LRT p-value` = format.pval(
              selection_LRT_p,
              digits = 3,
              eps = 0.001
            ),
            `Reported slope p-value` = format.pval(
              p.value,
              digits = 3,
              eps = 0.001
            ),
            `P-value method` = p_adjustment
          )
      )
    }

    interaction_x <- x %>%
      dplyr::left_join(
        intercept_results %>%
          dplyr::select(
            model_name,
            Region_LU,
            intercept_link,
            intercept_SE,
            intercept_response
          ),
        by = c("model_name", "Region_LU")
      )

    shiny::validate(
      shiny::need(
        all(!is.na(interaction_x$intercept_response)),
        "Intercept estimates are missing for at least one significant context."
      )
    )

    interaction_x %>%
      dplyr::transmute(
        `Region × Land Use` = friendly_region(Region_LU),
        Association = dplyr::case_when(
          direction == "positive" ~
            "<span class='direction-positive'>&uarr; Positive</span>",
          direction == "negative" ~
            "<span class='direction-negative'>&darr; Negative</span>",
          TRUE ~
            "<span class='direction-unclear'>&mdash; Unclear</span>"
        ),
        `Estimated slope (logit/link scale)` = round(trend, 3),
        `Slope standard error` = round(SE, 3),
        `Intercept estimate (logit/link scale)` = round(intercept_link, 3),
        `Intercept standard error` = round(intercept_SE, 3),
        `Expected ES at biodiversity = 0` = round(intercept_response, 3),
        `Selection Delta AIC` = round(selection_delta_AIC, 2),
        `Selection LRT p-value` = format.pval(
          selection_LRT_p,
          digits = 3,
          eps = 0.001
        ),
        `Reported slope p-value` = format.pval(
          p.value,
          digits = 3,
          eps = 0.001
        ),
        `P-value method` = p_adjustment
      )
  },
  striped = TRUE,
  bordered = TRUE,
  spacing = "s",
  sanitize.text.function = function(x) x)

  output$indicator_additive_expected_es_ui <- shiny::renderUI({
    if (selected_indicator_type() != "additive") {
      return(NULL)
    }

    shiny::tags$details(
      class = "expected-es-details",
      shiny::tags$summary(
        "Show intercept estimates by Region × Land Use"
      ),
      shiny::tags$div(
        class = "expected-es-details-body",
        shiny::p(
          class = "help-note",
          "The additive model has one common biodiversity slope. The intercept remains context-specific and is evaluated at biodiversity = 0 with Site_ID random intercept set to zero."
        ),
        shiny::tableOutput("indicator_additive_expected_es_table")
      )
    )
  })

  output$indicator_additive_expected_es_table <- shiny::renderTable({
    shiny::req(selected_indicator_type() == "additive")

    x <- selected_indicator_results()
    model_name_now <- unique(x$model_name)

    shiny::validate(
      shiny::need(
        length(model_name_now) == 1,
        "The additive result does not map to a unique model."
      )
    )

    intercept_x <- intercept_results %>%
      dplyr::filter(.data$model_name == model_name_now[[1]]) %>%
      dplyr::filter(!is.na(.data$intercept_response)) %>%
      dplyr::arrange(.data$Region_LU)

    shiny::validate(
      shiny::need(
        nrow(intercept_x) > 0,
        "No context-specific intercept estimates are available for this additive model."
      )
    )

    intercept_x %>%
      dplyr::transmute(
        `Region × Land Use` = friendly_region(Region_LU),
        `Intercept estimate (logit/link scale)` = round(intercept_link, 3),
        `Intercept standard error` = round(intercept_SE, 3),
        `Expected ES at biodiversity = 0` = round(intercept_response, 3)
      )
  },
  striped = TRUE,
  bordered = TRUE,
  spacing = "s")

  ## -------------------------------------------------------
  ## Prediction pathway: one indicator at a time
  ## -------------------------------------------------------
  output$prediction_region_selector <- shiny::renderUI({
    shiny::req(input$prediction_es)

    region_values <- available_regions_for_es(input$prediction_es)

    shiny::selectInput(
      inputId = "prediction_region",
      label = "Which Region × Land Use context?",
      choices = stats::setNames(
        region_values,
        friendly_region(region_values)
      )
    )
  })

  output$prediction_indicator_selector <- shiny::renderUI({
    shiny::req(input$prediction_es, input$prediction_region)

    choice_table <- indicator_choice_table(
      input$prediction_es,
      input$prediction_region
    )

    shiny::selectInput(
      inputId = "prediction_indicator",
      label = "Which biodiversity variable did you measure?",
      choices = stats::setNames(
        choice_table$biodiversity,
        choice_table$choice_label
      )
    )
  })

  prediction_selected_row <- shiny::reactive({
    shiny::req(
      input$prediction_es,
      input$prediction_region,
      input$prediction_indicator
    )

    selected_row <- indicator_results %>%
      dplyr::filter(
        ES == input$prediction_es,
        biodiversity == input$prediction_indicator,
        selected_model == "additive" |
          (selected_model == "interaction" & Region_LU == input$prediction_region)
      )

    shiny::validate(
      shiny::need(
        nrow(selected_row) == 1,
        "No unique selected model is available for this combination."
      )
    )

    selected_row
  })

  prediction_fit <- shiny::reactive({
    selected_row <- prediction_selected_row()
    model_name <- selected_row$model_name[[1]]
    fit <- selected_fits[[model_name]]

    shiny::validate(
      shiny::need(
        !is.null(fit),
        paste("The selected model could not be fitted:", model_name)
      )
    )

    fit
  })

  output$prediction_value_input <- shiny::renderUI({
    shiny::req(
      input$prediction_region,
      input$prediction_indicator
    )

    fit <- prediction_fit()

    context_data <- fit$data %>%
      dplyr::filter(
        as.character(Region_LU) == input$prediction_region
      )

    raw_values <- from_model_scale(
      input$prediction_indicator,
      context_data$biod_value
    )

    raw_values <- raw_values[is.finite(raw_values)]

    default_value <- if (length(raw_values) > 0) {
      stats::median(raw_values, na.rm = TRUE)
    } else {
      0
    }

    raw_range <- range(raw_values, na.rm = TRUE)
    spread <- diff(raw_range)

    step_value <- if (all(is.finite(raw_range)) && spread > 0) {
      max(signif(spread / 100, 2), 0.001)
    } else {
      0.1
    }

    transformation_note <- if (uses_log1p(input$prediction_indicator)) {
      "Enter the original untransformed measurement. The app applies log1p() automatically."
    } else {
      "Enter the measurement in the same units used in the BIOservicES reference dataset."
    }

    shiny::tagList(
      shiny::numericInput(
        inputId = "measured_value",
        label = "Enter the single plot-measured value",
        value = signif(default_value, 5),
        min = 0,
        step = step_value
      ),
      shiny::p(
        class = "help-note",
        transformation_note
      )
    )
  })

  prediction_state <- shiny::reactiveVal(NULL)

  shiny::observeEvent(
    list(
      input$prediction_es,
      input$prediction_region,
      input$prediction_indicator,
      input$measured_value
    ),
    {
      prediction_state(NULL)
    },
    ignoreInit = TRUE
  )

  shiny::observeEvent(input$run_prediction, {
    result <- tryCatch({
      shiny::req(
        input$prediction_es,
        input$prediction_region,
        input$prediction_indicator,
        input$measured_value
      )

      fit <- prediction_fit()
      region_name <- input$prediction_region
      indicator_name <- input$prediction_indicator
      measured_raw <- as.numeric(input$measured_value)
      measured_model <- to_model_scale(indicator_name, measured_raw)

      if (!region_name %in% levels(fit$data$Region_LU)) {
        stop("The selected Region × Land Use context was not represented in the fitted model.")
      }

      context_model_data <- fit$data %>%
        dplyr::filter(as.character(Region_LU) == region_name)

      if (nrow(context_model_data) < 2) {
        stop("There are too few observations in this context to generate a prediction.")
      }

      new_data <- make_population_newdata(
        fit = fit,
        region_name = region_name,
        biod_values = measured_model
      )

      predicted_es <- predict_population(fit, new_data)

      if (length(predicted_es) != 1 || !is.finite(predicted_es)) {
        stop("The model did not return a valid prediction.")
      }

      ## Context-specific reference distribution.
      ## It uses all available ecosystem-service observations in df_shiny,
      ## not only the complete cases of the selected biodiversity variable.
      reference_values <- df_shiny[[fit$es_var]][
        as.character(df_shiny$Region_LU) == region_name
      ]
      reference_values <- as.numeric(reference_values)
      reference_values <- reference_values[is.finite(reference_values)]

      if (length(reference_values) < 3) {
        stop("There are too few ecosystem-service observations to calculate reference quantiles.")
      }

      quantiles <- stats::quantile(
        reference_values,
        probs = c(0.33, 0.66, 0.99),
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )

      traffic_class <- dplyr::case_when(
        predicted_es <= quantiles[1] ~ "Red",
        predicted_es <= quantiles[2] ~ "Yellow",
        TRUE ~ "Green"
      )

      traffic_explanation <- dplyr::case_when(
        traffic_class == "Red" ~
          "The prediction falls in the lower third of the observed distribution for this context.",
        traffic_class == "Yellow" ~
          "The prediction falls in the middle third of the observed distribution for this context.",
        TRUE ~
          "The prediction falls in the upper third of the observed distribution for this context."
      )

      observed_model_range <- range(
        context_model_data$biod_value,
        na.rm = TRUE
      )

      observed_raw_range <- from_model_scale(
        indicator_name,
        observed_model_range
      )

      outside_calibration_range <-
        measured_model < observed_model_range[1] ||
        measured_model > observed_model_range[2]

      above_99th_percentile <- predicted_es > quantiles[3]

      ## Prediction curve displayed on the original measurement scale.
      raw_context_values <- from_model_scale(
        indicator_name,
        context_model_data$biod_value
      )

      raw_x_range <- range(raw_context_values, na.rm = TRUE)

      grid_raw <- seq(
        from = raw_x_range[1],
        to = raw_x_range[2],
        length.out = 150
      )

      grid_model_values <- vapply(
        grid_raw,
        function(x) to_model_scale(indicator_name, x),
        numeric(1)
      )

      grid_newdata <- make_population_newdata(
        fit = fit,
        region_name = region_name,
        biod_values = grid_model_values
      )

      grid_prediction <- predict_population(fit, grid_newdata)

      observed_plot_data <- context_model_data %>%
        dplyr::mutate(
          measured_raw = from_model_scale(
            indicator_name,
            biod_value
          )
        )

      curve_plot_data <- data.frame(
        measured_raw = grid_raw,
        predicted = grid_prediction
      )

      list(
        error = NULL,
        fit = fit,
        region = region_name,
        indicator = indicator_name,
        selected_model = prediction_selected_row()$selected_model[[1]],
        measured_raw = measured_raw,
        measured_model = measured_model,
        predicted_es = predicted_es,
        quantiles = quantiles,
        traffic_class = traffic_class,
        traffic_explanation = traffic_explanation,
        observed_raw_range = observed_raw_range,
        outside_calibration_range = outside_calibration_range,
        above_99th_percentile = above_99th_percentile,
        observed_plot_data = observed_plot_data,
        curve_plot_data = curve_plot_data
      )
    },
    error = function(e) {
      list(
        error = conditionMessage(e)
      )
    })

    prediction_state(result)
  })

  output$prediction_result_card <- shiny::renderUI({
    result <- prediction_state()

    if (is.null(result)) {
      return(
        shiny::div(
          class = "summary-box",
          "Select a context and one measured biodiversity variable, then click Run prediction."
        )
      )
    }

    if (!is.null(result$error)) {
      return(
        shiny::div(
          class = "traffic-card traffic-red",
          shiny::strong("Prediction could not be generated."),
          shiny::br(),
          result$error
        )
      )
    }

    css_class <- paste0(
      "traffic-card traffic-",
      tolower(result$traffic_class)
    )

    dot_class <- paste0(
      "traffic-dot dot-",
      tolower(result$traffic_class)
    )

    shiny::div(
      class = css_class,
      shiny::h3(
        shiny::span(class = dot_class),
        result$traffic_class
      ),
      shiny::div(
        class = "prediction-score",
        sprintf("%.3f", result$predicted_es)
      ),
      shiny::p(
        shiny::strong("Predicted "),
        tolower(friendly_es(input$prediction_es)),
        " for ",
        friendly_region(result$region),
        "."
      ),
      shiny::p(result$traffic_explanation),
      shiny::p(
        shiny::strong("Indicator: "),
        friendly_biod(result$indicator),
        shiny::br(),
        shiny::strong("Relationship: "),
        if (result$selected_model == "additive") {
          "common biodiversity slope across contexts"
        } else {
          "Region × Land Use-specific biodiversity slope"
        },
        shiny::br(),
        shiny::strong("Entered value: "),
        format(result$measured_raw, digits = 6, trim = TRUE)
      )
    )
  })

  output$prediction_warning <- shiny::renderUI({
    result <- prediction_state()

    if (
      is.null(result) ||
      !is.null(result$error)
    ) {
      return(NULL)
    }

    warning_items <- list()

    if (isTRUE(result$outside_calibration_range)) {
      warning_items <- append(
        warning_items,
        list(
          shiny::p(
            shiny::strong("Outside the observed indicator range. "),
            "The entered value is outside the range represented for this context (",
            format(result$observed_raw_range[1], digits = 5, trim = TRUE),
            " to ",
            format(result$observed_raw_range[2], digits = 5, trim = TRUE),
            "). The result involves extrapolation and should be interpreted cautiously."
          )
        )
      )
    }

    if (isTRUE(result$above_99th_percentile)) {
      warning_items <- append(
        warning_items,
        list(
          shiny::p(
            shiny::strong("Above the 99th-percentile reference. "),
            "The prediction remains classified as Green, but it lies above the upper reference value observed for this context."
          )
        )
      )
    }

    if (length(warning_items) == 0) {
      return(NULL)
    }

    shiny::div(
      class = "warning-box",
      warning_items
    )
  })

  output$prediction_threshold_table <- shiny::renderTable({
    result <- prediction_state()

    shiny::req(
      !is.null(result),
      is.null(result$error)
    )

    data.frame(
      `Reference point` = c(
        "33rd percentile",
        "66th percentile",
        "99th percentile"
      ),
      `Ecosystem-service value` = round(result$quantiles, 3),
      `Use in the traffic light` = c(
        "Upper limit of Red",
        "Upper limit of Yellow",
        "Upper reference within Green"
      ),
      check.names = FALSE
    )
  },
  striped = TRUE,
  bordered = TRUE,
  spacing = "s")

  output$prediction_plot <- shiny::renderPlot({
    result <- prediction_state()

    shiny::req(
      !is.null(result),
      is.null(result$error)
    )

    q33 <- result$quantiles[1]
    q66 <- result$quantiles[2]
    q99 <- result$quantiles[3]

    ggplot2::ggplot() +
      ggplot2::annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = 0,
        ymax = q33,
        fill = "#fde8e7",
        alpha = 0.55
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = q33,
        ymax = q66,
        fill = "#fff4bf",
        alpha = 0.55
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = q66,
        ymax = 1,
        fill = "#e4f3e8",
        alpha = 0.55
      ) +
      ggplot2::geom_hline(
        yintercept = q33,
        linetype = "dashed",
        linewidth = 0.6
      ) +
      ggplot2::geom_hline(
        yintercept = q66,
        linetype = "dashed",
        linewidth = 0.6
      ) +
      ggplot2::geom_hline(
        yintercept = q99,
        linetype = "dotted",
        linewidth = 0.7
      ) +
      ggplot2::geom_point(
        data = result$observed_plot_data,
        ggplot2::aes(
          x = measured_raw,
          y = ES_beta
        ),
        alpha = 0.35,
        size = 2
      ) +
      ggplot2::geom_line(
        data = result$curve_plot_data,
        ggplot2::aes(
          x = measured_raw,
          y = predicted
        ),
        linewidth = 1.25
      ) +
      ggplot2::geom_point(
        data = data.frame(
          measured_raw = result$measured_raw,
          predicted = result$predicted_es
        ),
        ggplot2::aes(
          x = measured_raw,
          y = predicted
        ),
        size = 5,
        shape = 21,
        fill = "white",
        stroke = 1.4
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste0(
          friendly_biod(result$indicator),
          " in ",
          friendly_region(result$region)
        ),
        subtitle = if (result$selected_model == "additive") {
          "Common biodiversity coefficient on the logit scale; context-specific intercept and entered measurement"
        } else {
          "Context-specific fitted relationship and the entered measurement"
        },
        x = paste0(
          friendly_biod(result$indicator),
          " — original measurement scale"
        ),
        y = friendly_es(input$prediction_es)
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      )
  })

  ## -------------------------------------------------------
  ## Application status
  ## -------------------------------------------------------
  output$app_status <- shiny::renderText({
    paste0(
      "Dataset loaded: ", nrow(df_shiny), " observations. ",
      "Mixed models loaded: ", sum(model_loaded),
      " of ", length(model_loaded), ". ",
      "Interaction models in export: ",
      dplyr::n_distinct(indicator_results$model_name[indicator_results$selected_model == "interaction"]),
      "; additive models in export: ",
      dplyr::n_distinct(indicator_results$model_name[indicator_results$selected_model == "additive"]),
      ". Intercept rows loaded: ", nrow(intercept_results),
      ". Indicator-result source: ", indicator_data_source,
      ". Model source: ", model_source, "."
    )
  })
}

shiny::shinyApp(ui = ui, server = server)
