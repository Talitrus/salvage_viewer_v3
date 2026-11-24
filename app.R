# Load necessary libraries
library(shiny)
library(tidyverse)
library(lubridate)
library(plotly)
library(htmltools)
# Configure timezone database directory for Windows Conda/Pixi environments
if (.Platform$OS.type == "windows") {
  conda_tzdir <- file.path(getwd(), ".pixi/envs/default/share/zoneinfo")
  if (dir.exists(conda_tzdir)) {
    Sys.setenv(TZDIR = conda_tzdir)
  }
}

# ==============================================================================
# Security & Validation Constants
# ==============================================================================
APP_TZ <- "America/Los_Angeles"
MAX_LAD_DAYS <- 1500L
MAX_SELECT_IDS_CHARS <- 10000L
MAX_SELECT_IDS_COUNT <- 500L
ALLOWED_COLOR_COLUMNS <- c(
  "Genetic_Assignment",
  "Sherlock_Assignment",
  "Spring_tributary",
  "Facility",
  "rosa_informed",
  "PBT"
)

# ==============================================================================
# Data Discovery, Validation & Ingestion Helpers
# ==============================================================================

# HTML-escape untrusted character fields before inserting into Plotly HTML tooltips
escape_html <- function(x) {
  if (length(x) == 0L) return(character(0))
  x_chr <- as.character(x)
  out <- rep("N/A", length(x_chr))
  non_na <- !is.na(x_chr)
  if (any(non_na)) {
    out[non_na] <- htmltools::htmlEscape(x_chr[non_na], attribute = TRUE)
  }
  out
}

# Validate water_year matches ^WY\d{4}$, lies in a sensible range, and resolves safely under data_dir
validate_water_year <- function(water_year, data_dir = NULL) {
  if (!is.character(water_year) || length(water_year) != 1L || is.na(water_year) ||
      !grepl("^WY\\d{4}$", water_year)) {
    stop("Invalid water_year: must match '^WY\\d{4}$'.", call. = FALSE)
  }
  wy_num <- as.integer(sub("^WY", "", water_year))
  if (is.na(wy_num) || wy_num < 1900L || wy_num > 2100L) {
    stop("Invalid water_year: year out of range.", call. = FALSE)
  }
  if (!is.null(data_dir)) {
    if (!is.character(data_dir) || length(data_dir) != 1L || is.na(data_dir) || !dir.exists(data_dir)) {
      stop("Invalid data_dir: directory does not exist.", call. = FALSE)
    }
    norm_base <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
    wy_dir <- file.path(norm_base, water_year)
    if (!dir.exists(wy_dir)) {
      stop(sprintf("Water year directory not found: %s", water_year), call. = FALSE)
    }
    norm_wy <- normalizePath(wy_dir, winslash = "/", mustWork = TRUE)
    base_prefix <- paste0(sub("/+$", "", norm_base), "/")
    if (!startsWith(norm_wy, base_prefix) || !identical(basename(norm_wy), water_year)) {
      stop("Invalid water_year path resolution.", call. = FALSE)
    }
  }
  invisible(wy_num)
}

# Select the primary CFS_data CSV inside a water year folder, ignoring backup/temp files
select_cfs_file <- function(wy_dir, water_year) {
  if (!dir.exists(wy_dir)) return(character(0))
  files <- list.files(wy_dir, pattern = "^CFS_data.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0L) return(character(0))

  is_backup <- grepl("(^|[._-])(bak|backup|old|tmp|temp|copy)([._-]|\\.csv$)", basename(files), ignore.case = TRUE)
  files <- files[!is_backup]
  if (length(files) == 0L) return(character(0))

  canonical_wy <- tolower(paste0("CFS_data_", water_year, ".csv"))
  exact_wy <- files[tolower(basename(files)) == canonical_wy]
  if (length(exact_wy) > 0L) return(exact_wy[1])

  exact_generic <- files[tolower(basename(files)) == "cfs_data.csv"]
  if (length(exact_generic) > 0L) return(exact_generic[1])

  if (length(files) > 1L) {
    warning(sprintf("Multiple CFS_data CSV files found in %s; using %s", wy_dir, basename(sort(files)[1])), call. = FALSE)
  }
  sort(files)[1]
}

# Parse and bound user-entered fish IDs across commas, semicolons, tabs, CRLF, and spaces
parse_selected_ids <- function(ids_text, max_chars = MAX_SELECT_IDS_CHARS, max_ids = MAX_SELECT_IDS_COUNT) {
  if (is.null(ids_text) || length(ids_text) < 1L || is.na(ids_text[1])) {
    return(character(0))
  }
  raw_str <- as.character(ids_text[1])
  if (nchar(raw_str, type = "bytes") > max_chars) {
    stop(sprintf("Selected IDs input exceeds maximum length (%d characters).", max_chars), call. = FALSE)
  }
  trimmed <- trimws(raw_str)
  if (!nzchar(trimmed)) {
    return(character(0))
  }
  tokens <- stringr::str_split_1(trimmed, "[,;\\s]+")
  tokens <- unique(trimws(tokens))
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) > max_ids) {
    stop(sprintf("Too many IDs specified (%d); maximum allowed is %d.", length(tokens), max_ids), call. = FALSE)
  }
  tokens
}

# Discover available water year folders (e.g., "WY2027", "WY2026") sorted newest first
get_available_water_years <- function(data_dir = "data") {
  if (!dir.exists(data_dir)) return(character(0))
  subdirs <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
  wy_dirs <- subdirs[grepl("^WY\\d{4}$", subdirs)]
  # Keep only folders that actually contain a valid non-backup CFS_data CSV file
  wy_dirs <- wy_dirs[purrr::map_lgl(wy_dirs, function(wy) {
    length(select_cfs_file(file.path(data_dir, wy), wy)) > 0L
  })]
  sort(wy_dirs, decreasing = TRUE)
}

# Return start (Oct 1 of WY - 1) and end (Sep 30 of WY) dates for a water year string (e.g., "WY2027")
get_water_year_bounds <- function(water_year) {
  wy_num <- validate_water_year(water_year)
  list(
    wy_num = wy_num,
    start  = as.Date(sprintf("%04d-10-01", wy_num - 1L)),
    end    = as.Date(sprintf("%04d-09-30", wy_num))
  )
}

# Dynamically generate daily Length-at-Date (LAD) ribbon data for any date range
# from the year-agnostic template (`data/lad_template.tsv`)
generate_lad_data <- function(start_date, end_date, lad_template, max_days = MAX_LAD_DAYS) {
  start_dt <- suppressWarnings(as.Date(start_date))
  end_dt <- suppressWarnings(as.Date(end_date))
  if (length(start_dt) != 1L || length(end_dt) != 1L ||
      is.na(start_dt) || is.na(end_dt) ||
      !is.finite(as.numeric(start_dt)) || !is.finite(as.numeric(end_dt))) {
    stop("Invalid date range: start_date and end_date must be single valid Dates.", call. = FALSE)
  }
  if (start_dt > end_dt) {
    stop("Invalid date range: start_date must be on or before end_date.", call. = FALSE)
  }
  if (as.integer(end_dt - start_dt) > max_days) {
    stop(sprintf("Date range exceeds maximum allowed span of %d days.", max_days), call. = FALSE)
  }

  dates_df <- tibble(date = seq.Date(start_dt, end_dt, by = "day")) %>%
    mutate(
      year  = lubridate::year(date),
      month = lubridate::month(date),
      day   = lubridate::mday(date)
    )

  dates_df %>%
    inner_join(lad_template, by = c("month", "day"), relationship = "many-to-many") %>%
    mutate(
      cohort = paste0(run, sprintf("%02d", (year + brood_year_offset) %% 100))
    ) %>%
    select(month, day, year, date, run, cohort, max, min)
}

# Load and summarize PBT results for a given water year
load_pbt_summary <- function(water_year, data_dir = "data") {
  validate_water_year(water_year, data_dir = data_dir)
  pbt_dir <- file.path(data_dir, water_year, "PBT_results")
  pbt_files <- list.files(pbt_dir, pattern = "\\.tsv$", full.names = TRUE)

  if (length(pbt_files) == 0L) {
    return(tibble(D2_indiv = character(), PBT_info = character()))
  }

  expected_cols <- c(
    "D2_indiv", "D1_indiv", "logl_ratio", "num_loc",
    "total_incompat", "total_compat", "fraction_incompat"
  )

  pbt_raw <- purrr::map_dfr(unname(pbt_files), function(fpath) {
    tbl <- tryCatch(
      readr::read_tsv(fpath, col_types = cols(.default = col_character()), show_col_types = FALSE),
      error = function(e) {
        warning(sprintf("Skipping unreadable PBT file '%s': %s", basename(fpath), conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    if (is.null(tbl) || nrow(tbl) == 0L) {
      return(NULL)
    }
    if (!("D2_indiv" %in% names(tbl))) {
      warning(sprintf("Skipping PBT file '%s': missing required 'D2_indiv' column.", basename(fpath)), call. = FALSE)
      return(NULL)
    }
    for (col_nm in expected_cols) {
      if (!(col_nm %in% names(tbl))) {
        tbl[[col_nm]] <- NA_character_
      }
    }
    tbl %>%
      select(all_of(expected_cols)) %>%
      mutate(
        logl_ratio        = suppressWarnings(as.numeric(logl_ratio)),
        num_loc           = suppressWarnings(as.integer(num_loc)),
        total_incompat    = suppressWarnings(as.integer(total_incompat)),
        total_compat      = suppressWarnings(as.integer(total_compat)),
        fraction_incompat = suppressWarnings(as.numeric(fraction_incompat)),
        Source            = tools::file_path_sans_ext(basename(fpath))
      )
  })

  if (nrow(pbt_raw) == 0L) {
    return(tibble(D2_indiv = character(), PBT_info = character()))
  }

  pbt_raw %>%
    filter(!is.na(D2_indiv), nzchar(trimws(D2_indiv))) %>%
    mutate(
      parent_safe = ifelse(is.na(D1_indiv) | !nzchar(trimws(D1_indiv)), "Unknown", escape_html(D1_indiv)),
      source_safe = escape_html(Source),
      logl_ratio_num = round(logl_ratio, 2),
      logl_str = ifelse(is.na(logl_ratio_num), "N/A", as.character(logl_ratio_num)),
      loci_str = case_when(
        !is.na(total_compat) & !is.na(num_loc) ~ paste0(
          total_compat, "/", num_loc,
          ifelse(is.na(total_incompat), "", paste0(" (", total_incompat, " incompat)"))
        ),
        is.na(total_compat) & !is.na(num_loc) ~ paste0(num_loc, " loci"),
        TRUE ~ "N/A"
      ),
      match_desc = paste0(
        "• Candidate Parent: <b>", parent_safe, "</b><br>",
        "  Source: ", source_safe, "<br>",
        "  LogL Ratio: ", logl_str, "<br>",
        "  Loci (Compat/Total): ", loci_str
      )
    ) %>%
    group_by(D2_indiv) %>%
    summarise(
      PBT_info = paste(match_desc, collapse = "<br>"),
      .groups = "drop"
    )
}

# Column mapping configuration for CFS datasets.
# Maps internal standardized column names (left) to candidate raw CSV column names (right).
# If a future water year uses different column names, either:
#   1. Add the new raw column name to the candidate vector in `default_map`, OR
#   2. Add a year-specific override inside `switch(water_year, ...)` below.
get_cfs_column_map <- function(water_year) {
  default_map <- list(
    ID                  = c("CVTA_ID", "ID", "Sample_ID"),
    SampleDateTime      = c("SampleDateTime", "Sample_DateTime", "SampleDate"),
    ForkLength          = c("ForkLength", "Fork_Length", "FL"),
    LAD.Class           = c("LAD_Race", "LAD.Class", "LAD_Class"),
    Genetic_Assignment  = c("Genetic_Assignment", "GeneticAssignment"),
    Spring_tributary    = c("Spring_Pop_Assignment", "Spring_tributary"),
    Sherlock_Assignment = c("Sherlock_Assignment", "SherlockAssignment"),
    Facility            = c("Facility"),
    repunit             = c("repunit"),
    rosa_classification = c("rosa_classification")
  )

  # Optional per-year overrides if a specific water year has a unique schema
  year_overrides <- switch(water_year,
    # Example for future adaptation:
    # "WY2028" = list(ID = "New_ID_Column_Name"),
    list()
  )

  modifyList(default_map, year_overrides)
}

# Standardize raw dataframe columns using the column map and ensure expected columns exist
standardize_cfs_columns <- function(df, col_map) {
  names(df) <- make.names(names(df), unique = TRUE)

  for (target_col in names(col_map)) {
    candidates <- make.names(col_map[[target_col]])
    matched <- intersect(candidates, names(df))
    if (length(matched) > 0) {
      df[[target_col]] <- df[[matched[1]]]
    } else if (!(target_col %in% names(df))) {
      # Fill missing optional columns with NA so downstream code doesn't break
      df[[target_col]] <- NA_character_
    }
  }
  df
}

# Load and process CFS data for a given water year
load_cfs_data <- function(water_year, pbt_summary, data_dir = "data") {
  validate_water_year(water_year, data_dir = data_dir)
  wy_dir <- file.path(data_dir, water_year)
  cfs_file <- select_cfs_file(wy_dir, water_year)
  if (length(cfs_file) == 0L) {
    stop(sprintf("No valid CFS_data CSV found for %s.", water_year), call. = FALSE)
  }

  df <- read_csv(cfs_file, show_col_types = FALSE)

  # Standardize column names across water years
  col_map <- get_cfs_column_map(water_year)
  df <- standardize_cfs_columns(df, col_map)

  df %>%
    mutate(
      ID = as.character(ID),
      ForkLength = suppressWarnings(as.numeric(ForkLength)),
      SampleDate = lubridate::parse_date_time(
        SampleDateTime,
        orders = c("ymd HM", "ymd HMS", "mdy HM", "mdy HMS", "ymd", "mdy"),
        tz = APP_TZ
      ),
      Genetic_Assignment = case_when(
        Genetic_Assignment == "sacfall" ~ "Fall",
        TRUE ~ as.character(Genetic_Assignment)
      ),
      Sherlock_Assignment = case_when(
        Sherlock_Assignment %in% c(NA, "NA") ~ "Not Screened",
        TRUE ~ as.character(Sherlock_Assignment)
      ),
      rosa = case_when(
        is.na(rosa_classification) | trimws(as.character(rosa_classification)) %in% c("", "NA") ~ NA_character_,
        rosa_classification == "Unknown" ~ "Unknown",
        rosa_classification %in% c("Early", "Early-Spring", "Early-Winter", "Early-returning", "Early / Early", "Sp-Win hybrid") ~ "Early",
        rosa_classification %in% c("Late", "Late / Late") ~ "Late",
        rosa_classification %in% c("Intermediate", "Heterozygous", "Early / Late", "Sp-Fall hybrid", "Fall-Win hybrid") ~ "Intermediate",
        TRUE ~ "Unknown"
      ),
      rosa_informed = case_when(
        !is.na(repunit) & !is.na(rosa) &
          str_detect(repunit, regex("fall", ignore_case = TRUE)) & (rosa == "Early") ~ "Fall to Spring",
        !is.na(repunit) & !is.na(rosa) &
          str_detect(repunit, regex("spring", ignore_case = TRUE)) & (rosa == "Late") ~ "Spring to Fall",
        TRUE ~ "Not RoSA-informed"
      )
    ) %>%
    left_join(pbt_summary, by = c("ID" = "D2_indiv")) %>%
    mutate(
      has_pbt_match = !is.na(PBT_info),
      PBT_info = ifelse(has_pbt_match, PBT_info, "No PBT match"),
      PBT = ifelse(has_pbt_match, Genetic_Assignment, "FALSE"),
      tooltip_text = paste0(
        "<b>ID:</b> ", escape_html(ID), "<br>",
        "<b>Date:</b> ", ifelse(is.na(SampleDate), "N/A", format(SampleDate, "%Y-%m-%d %H:%M")), "<br>",
        "<b>Fork Length:</b> ", ifelse(is.na(ForkLength), "N/A", paste0(ForkLength, " mm")), "<br>",
        "<b>LAD Class:</b> ", escape_html(LAD.Class), "<br>",
        "<b>Genetic Assignment:</b> ", escape_html(Genetic_Assignment), "<br>",
        "<b>Sherlock Assignment:</b> ", escape_html(Sherlock_Assignment), "<br>",
        "<b>Facility:</b> ", escape_html(Facility), "<br>",
        "<b>PBT Match:</b> ", ifelse(has_pbt_match, "TRUE", "FALSE"),
        ifelse(has_pbt_match, paste0("<br><br><b>PBT Matches:</b><br>", PBT_info), "")
      )
    ) %>%
    select(
      ID, SampleDate, ForkLength, LAD.Class, Genetic_Assignment,
      Spring_tributary, Sherlock_Assignment, Facility, repunit,
      rosa_classification, rosa, rosa_informed,
      has_pbt_match, PBT_info, PBT, tooltip_text
    )
}

available_water_years <- get_available_water_years("data")

# Read in canonical year-agnostic LAD template once at app startup
lad_template <- read_tsv(
  file.path("data", "lad_template.tsv"),
  col_types = cols(
    month = col_integer(),
    day = col_integer(),
    run = col_character(),
    brood_year_offset = col_integer(),
    max = col_double(),
    min = col_double()
  ),
  show_col_types = FALSE
)

# Static color mappings initialized once at app startup
prerun_color_palette <- RColorBrewer::brewer.pal(9, "Set1")
base_run_color_mapping <- c(
  "Fall / Late Fall" = prerun_color_palette[5],
  "Fall or Late fall" = prerun_color_palette[5],
  "Fall" = "#c9b428",
  "Late Fall" = prerun_color_palette[7],
  "Spring" = prerun_color_palette[3],
  "Winter" = prerun_color_palette[2],
  "Non-Winter" = prerun_color_palette[9],
  "Early" = prerun_color_palette[2],
  "Early-Spring" = prerun_color_palette[3],
  "Early-Winter" = prerun_color_palette[2],
  "Sp-Win hybrid" = "#00dba4",
  "Intermediate" = "#c9b428",
  "Late" = prerun_color_palette[5],
  "Latefall" = prerun_color_palette[7],
  "Late fall" = prerun_color_palette[7],
  "female" = prerun_color_palette[1],
  "male" = prerun_color_palette[2],
  "TRUE" = prerun_color_palette[1],
  "FALSE" = prerun_color_palette[9],
  "CVP" = prerun_color_palette[2],
  "SWP" = prerun_color_palette[1],
  "Early / Early" = prerun_color_palette[2],
  "Early-returning" = prerun_color_palette[2],
  "Late / Late" = prerun_color_palette[5],
  "Early / Late" = "#c9b428",
  "Heterozygous" = "#c9b428",
  "NA" = prerun_color_palette[9],
  "Not Screened" = "#d1d1d1",
  "Other" = "#d1d1d1",
  "Unknown" = "#868686",
  "O. mykiss" = "#d1d1d1",
  "Steelhead" = "#d1d1d1",
  "Likely heterozygote" = "red",
  "Pending QA/QC" = "yellow",
  "Pending" = "yellow",
  "No Consensus" = "yellow",
  "No Tributary Assigned" = prerun_color_palette[9],
  "Butte" = prerun_color_palette[1],
  "Butte Creek Spring" = prerun_color_palette[1],
  "Deer" = prerun_color_palette[2],
  "Mill" = prerun_color_palette[3],
  "Feather River Spring" = prerun_color_palette[4],
  "Feather River" = prerun_color_palette[4],
  "Feather River Lineage" = prerun_color_palette[4],
  "Ambiguous" = "yellow",
  "Referred to GTseq: No SHERLOCK Detection" = "purple",
  "Referred to GTseq: heterozygote" = "purple",
  "Mill/Deer Creek Spring" = prerun_color_palette[3],
  "Fall to Spring" = prerun_color_palette[3],
  "Spring to Fall" = prerun_color_palette[7],
  "Not RoSA-informed" = "#868686"
)

cohort_prefix_colors <- c(
  "f"  = prerun_color_palette[6],
  "lf" = prerun_color_palette[7],
  "s"  = prerun_color_palette[3],
  "w"  = prerun_color_palette[2]
)

# Define UI for application
ui <- fluidPage(
  titlePanel("Chinook Salvage Genetic Data Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "water_year",
        "Select Water Year:",
        choices = available_water_years,
        selected = if (length(available_water_years) > 0) available_water_years[1] else NULL
      ),
      dateRangeInput("date_range", "Select date range:",
        start = NULL,
        end = NULL
      ),
      selectInput("LAD", "Select LAD size classes:", choices = NULL, multiple = TRUE),
      selectInput("color", "Select color variable:", choices = NULL),
      textAreaInput("selectIDs", "Show only selected IDs (comma-separated)")
    ),
    mainPanel(
      plotlyOutput("plot")
    ),
    fluid = TRUE
  ),
  p("Last updated: November 12, 2026. This visualization tool may not be updated in a timely manner and may not include the latest information or individuals that don't have assignments yet. Please contact Melinda Baerwald if you have questions regarding the interpretation of this data: Melinda.Baerwald@water.ca.gov")
)

# Define server logic
server <- function(input, output, session) {
  # Reactively load and join PBT + CFS dataset for the selected water year (cached across sessions)
  raw_year_data <- reactive({
    req(input$water_year)
    validate(need(
      is.character(input$water_year) && length(input$water_year) == 1L &&
        input$water_year %in% available_water_years,
      "Invalid water year selected."
    ))
    pbt <- load_pbt_summary(input$water_year)
    load_cfs_data(input$water_year, pbt)
  }) %>% bindCache(input$water_year)

  # Generate LAD ribbon data covering the entire selected water year (cached across sessions)
  lad_long <- reactive({
    req(input$water_year, raw_year_data())
    validate(need(
      is.character(input$water_year) && length(input$water_year) == 1L &&
        input$water_year %in% available_water_years,
      "Invalid water year selected."
    ))
    wy_bounds <- get_water_year_bounds(input$water_year)
    sample_dates <- as.Date(raw_year_data()$SampleDate, tz = APP_TZ)
    plausible_dates <- sample_dates[
      !is.na(sample_dates) &
        sample_dates >= (wy_bounds$start - days(366)) &
        sample_dates <= (wy_bounds$end + days(366))
    ]
    start_date <- min(c(wy_bounds$start, plausible_dates), na.rm = TRUE) - days(3)
    end_date   <- max(c(wy_bounds$end, plausible_dates), na.rm = TRUE) + days(3)
    generate_lad_data(start_date, end_date, lad_template)
  }) %>% bindCache(input$water_year)

  # Debounce ID text input so typing doesn't rebuild the plot on every keystroke
  selected_ids_debounced <- reactive({
    raw_val <- input$selectIDs
    if (is.null(raw_val) || length(raw_val) < 1L || is.na(raw_val[1])) "" else trimws(as.character(raw_val[1]))
  }) %>% debounce(400)

  parsed_ids <- reactive({
    ids_text <- selected_ids_debounced()
    res <- tryCatch(
      list(ok = TRUE, ids = parse_selected_ids(ids_text), msg = NULL),
      error = function(e) list(ok = FALSE, ids = character(0), msg = conditionMessage(e))
    )
    validate(need(res$ok, res$msg))
    res$ids
  })

  data <- reactive({
    df <- raw_year_data()
    ids <- parsed_ids()
    if (length(ids) == 0L) {
      return(df)
    }
    df %>% filter(ID %in% ids)
  })

  # Canonical cache key for ID selection (prevents unbounded shared cache pollution)
  matched_ids_cache_key <- reactive({
    ids <- parsed_ids()
    if (length(ids) == 0L) {
      return("__ALL__")
    }
    matched <- sort(unique(intersect(ids, raw_year_data()$ID)))
    if (length(matched) == 0L) "__NONE__" else matched
  })

  # Update selectInput and dateRangeInput choices only when water year data changes
  observeEvent(raw_year_data(), {
    df <- raw_year_data()
    valid_choices <- ALLOWED_COLOR_COLUMNS[ALLOWED_COLOR_COLUMNS %in% names(df)]
    current_color <- isolate(input$color)
    selected_color <- if (!is.null(current_color) && current_color %in% valid_choices) current_color else valid_choices[1]

    model_choices <- df %>%
      pull(LAD.Class) %>%
      unique() %>%
      sort()

    wy_bounds <- get_water_year_bounds(input$water_year)
    sample_dates <- as.Date(df$SampleDate, tz = APP_TZ)
    plausible_dates <- sample_dates[
      !is.na(sample_dates) &
        sample_dates >= (wy_bounds$start - days(366)) &
        sample_dates <= (wy_bounds$end + days(366))
    ]
    sample_min <- if (length(plausible_dates) > 0L) min(plausible_dates) else wy_bounds$start
    sample_max <- if (length(plausible_dates) > 0L) max(plausible_dates) else wy_bounds$end
    picker_min <- min(wy_bounds$start, sample_min)
    picker_max <- max(wy_bounds$end, sample_max)

    # Freeze date_range so p() waits for the updated dates rather than rendering twice
    freezeReactiveValue(input, "date_range")
    updateSelectInput(session, "color", choices = valid_choices, selected = selected_color)
    updateSelectInput(session, "LAD", choices = model_choices, selected = model_choices)
    updateDateRangeInput(
      session,
      "date_range",
      start = sample_min,
      end = sample_max,
      min = picker_min,
      max = picker_max
    )
  })

  p <- reactive({
    req(data(), lad_long(), input$color, length(input$date_range) >= 2L)

    color_col <- input$color
    validate(need(
      is.character(color_col) && length(color_col) == 1L &&
        color_col %in% ALLOWED_COLOR_COLUMNS && color_col %in% names(data()),
      "Invalid color variable selected."
    ))

    start_dt <- suppressWarnings(as.Date(input$date_range[1]))
    end_dt <- suppressWarnings(as.Date(input$date_range[2]))
    validate(need(
      length(start_dt) == 1L && length(end_dt) == 1L &&
        !is.na(start_dt) && !is.na(end_dt) &&
        is.finite(as.numeric(start_dt)) && is.finite(as.numeric(end_dt)),
      "Invalid date range selected."
    ))
    validate(need(start_dt <= end_dt, "Start date must be on or before end date."))
    validate(need(
      as.integer(end_dt - start_dt) <= MAX_LAD_DAYS,
      sprintf("Selected date range exceeds %d days.", MAX_LAD_DAYS)
    ))

    valid_lad <- intersect(as.character(input$LAD %||% character(0)), unique(raw_year_data()$LAD.Class))
    filtered_data <- data() %>%
      filter(
        LAD.Class %in% valid_lad,
        !is.na(SampleDate),
        as.Date(SampleDate, tz = APP_TZ) >= start_dt,
        as.Date(SampleDate, tz = APP_TZ) <= end_dt
      )
    if (color_col == "Spring_tributary") {
      filtered_data <- filtered_data %>%
        filter(!is.na(Spring_tributary))
    }

    all_cohorts <- unique(lad_long()$cohort)
    cohort_runs <- sub("\\d+$", "", all_cohorts)
    dynamic_cohort_colors <- setNames(unname(cohort_prefix_colors[cohort_runs]), all_cohorts)
    run_color_mapping <- c(base_run_color_mapping, dynamic_cohort_colors)

    lad_long_filtered <- lad_long() %>%
      filter(
        date >= start_dt - days(1),
        date <= end_dt + days(1)
      ) %>%
      mutate(
        date_dt = ISOdatetime(year, month, day, 0, 0, 0, tz = APP_TZ)
      )

    x_range_start <- ISOdatetime(
      lubridate::year(start_dt - days(3)),
      lubridate::month(start_dt - days(3)),
      lubridate::mday(start_dt - days(3)),
      0, 0, 0, tz = APP_TZ
    )
    x_range_end <- ISOdatetime(
      lubridate::year(end_dt + days(3)),
      lubridate::month(end_dt + days(3)),
      lubridate::mday(end_dt + days(3)),
      23, 59, 59, tz = APP_TZ
    )

    plt <- plot_ly(colors = run_color_mapping) %>%
      add_ribbons( # Add ribbons for min/max values overlaid on top of fork length
        data = lad_long_filtered,
        x = ~date_dt,
        ymin = ~min,
        ymax = ~max,
        color = ~cohort,
        colors = run_color_mapping,
        showlegend = TRUE,
        alpha = 0.2,
        size = 0.1,
        line = list(width = 0)
      )

    if (nrow(filtered_data) > 0L) {
      filtered_data <- filtered_data %>%
        mutate(.color_val = ifelse(is.na(.data[[color_col]]), "NA", as.character(.data[[color_col]])))

      plt <- plt %>%
        add_markers(
          data = filtered_data,
          x = ~SampleDate,
          y = ~ForkLength,
          color = ~.color_val,
          colors = run_color_mapping,
          text = ~tooltip_text,
          hoverinfo = "text",
          size = I(20),
          alpha = 0.5
        )
    } else {
      empty_msg <- if (nrow(data()) == 0L && length(parsed_ids()) > 0L) {
        "No fish records match the entered ID(s) in this water year."
      } else {
        "No fish records match the current filter criteria."
      }
      plt <- plt %>%
        layout(
          annotations = list(
            list(
              text = empty_msg,
              xref = "paper",
              yref = "paper",
              x = 0.5,
              y = 0.5,
              showarrow = FALSE,
              font = list(size = 15, color = "#555555")
            )
          )
        )
    }

    plt %>%
      layout(
        title = "Salvage Chinook Salmon",
        xaxis = list(title = "Sample Date", range = c(x_range_start, x_range_end)),
        yaxis = list(title = "Fork Length"),
        showlegend = TRUE
      )
  }) %>% bindCache(
    input$water_year,
    input$color,
    sort(unique(intersect(as.character(input$LAD %||% character(0)), unique(raw_year_data()$LAD.Class)))),
    if (length(input$date_range) >= 1L) suppressWarnings(as.Date(input$date_range[1])) else NA,
    if (length(input$date_range) >= 2L) suppressWarnings(as.Date(input$date_range[2])) else NA,
    matched_ids_cache_key(),
    cache = "session"
  )

  output$plot <- renderPlotly({
    p()
  })
}


# Run the application
shinyApp(ui = ui, server = server)
