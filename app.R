# Load necessary libraries
library(shiny)
library(tidyverse)
library(readxl)
library(lubridate)
library(plotly)
# Configure timezone database directory for Windows Conda/Pixi environments
if (.Platform$OS.type == "windows") {
  conda_tzdir <- file.path(getwd(), ".pixi/envs/default/share/zoneinfo")
  if (dir.exists(conda_tzdir)) {
    Sys.setenv(TZDIR = conda_tzdir)
  }
}
# Define UI for application
ui <- fluidPage(
  titlePanel("Chinook Salvage Genetic Data Viewer"),
  sidebarLayout(
    sidebarPanel(
      # fileInput("file", "Choose Excel File",
      #  accept = c(".xlsx")
      # ),
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
  p("Last updated: June 9, 2026. This visualization tool may not be updated in a timely manner and may not include the latest information or individuals that don't have assignments yet. Please contact Melinda Baerwald if you have questions regarding the interpretation of this data: Melinda.Baerwald@water.ca.gov")
)

# Define server logic
server <- function(input, output, session) {
  # Read in static LAD data
  tsv_files <- list.files("data", pattern = "^lad_long.*\\.txt$", full.names = TRUE)
  lad_long <- purrr::map_dfr(tsv_files, read_tsv)

  # Read in static PBT data
  pbt_files <- list.files("data/PBT_results", pattern = "\\.tsv$", full.names = TRUE)
  if (length(pbt_files) > 0) {
    pbt_raw <- purrr::map_dfr(pbt_files, read_tsv, col_types = cols(
      D2_indiv = col_character(),
      D1_indiv = col_character(),
      logl_ratio = col_double(),
      num_loc = col_integer(),
      total_incompat = col_integer(),
      total_compat = col_integer(),
      fraction_incompat = col_double()
    ))

    pbt_summary <- pbt_raw %>%
      filter(!is.na(D2_indiv), D2_indiv != "") %>%
      mutate(
        logl_ratio_num = round(as.numeric(logl_ratio), 2),
        match_desc = paste0(
          "• Candidate Parent: <b>", D1_indiv, "</b><br>",
          "  LogL Ratio: ", logl_ratio_num, "<br>",
          "  Loci (Compat/Total): ", total_compat, "/", num_loc, " (", total_incompat, " incompat)"
        )
      ) %>%
      group_by(D2_indiv) %>%
      summarise(
        PBT_info = paste(match_desc, collapse = "<br>"),
        .groups = "drop"
      )
  } else {
    pbt_summary <- tibble(D2_indiv = character(), PBT_info = character())
  }

  data <- reactive({
    df <- read_csv("data/CFS_data.csv")
    names(df) <- make.names(names(df)) # Make column names valid
    df <- df %>%
      rename(
        Spring_tributary = Spring_Pop_Assignment
      ) %>%
      mutate(
        SampleDate = lubridate::parse_date_time(SampleDateTime, orders = c("ymd HM", "mdy HM", "mdy HMS"), tz = "US/Pacific"),
        LAD.Class = LAD_Race,
        ID = CVTA_ID,
        Genetic_Assignment = case_when(
          Genetic_Assignment == "sacfall" ~ "Fall",
          TRUE ~ Genetic_Assignment
        ),
        Sherlock_Assignment = case_when(
          Sherlock_Assignment %in% c(NA, "NA") ~ "Not Screened",
          TRUE ~ Sherlock_Assignment
        ),
        rosa = case_when(
          rosa_classification == "Sp-Fall hybrid" ~ "Intermediate",
          rosa_classification == "Sp-Win hybrid" ~ "Early",
          rosa_classification == "Late" ~ "Late",
          rosa_classification == "Early-Spring" ~ "Early",
          rosa_classification == "Early-Winter" ~ "Early",
          TRUE ~ "ERROR"
        ),
        rosa_informed = case_when(
          str_detect(repunit, "fall") & (rosa == "Early") ~ "Fall to Spring",
          str_detect(repunit, "spring") & (rosa == "Late") ~ "Spring to Fall",
          TRUE ~ "Not RoSA-informed"
        )
      ) %>%
      left_join(pbt_summary, by = c("ID" = "D2_indiv")) %>%
      mutate(
        PBT = !is.na(PBT_info),
        PBT_info = ifelse(PBT, PBT_info, "No PBT match"),
        tooltip_text = paste0(
          "<b>ID:</b> ", ID, "<br>",
          "<b>Date:</b> ", ifelse(is.na(SampleDate), "N/A", format(SampleDate, "%Y-%m-%d %H:%M")), "<br>",
          "<b>Fork Length:</b> ", ifelse(is.na(ForkLength), "N/A", paste0(ForkLength, " mm")), "<br>",
          "<b>LAD Class:</b> ", ifelse(is.na(LAD.Class), "N/A", LAD.Class), "<br>",
          "<b>Genetic Assignment:</b> ", ifelse(is.na(Genetic_Assignment), "N/A", Genetic_Assignment), "<br>",
          "<b>Sherlock Assignment:</b> ", ifelse(is.na(Sherlock_Assignment), "N/A", Sherlock_Assignment), "<br>",
          "<b>Facility:</b> ", ifelse(is.na(Facility), "N/A", Facility), "<br>",
          "<b>PBT Match:</b> ", ifelse(PBT, "TRUE", "FALSE"),
          ifelse(PBT, paste0("<br><br><b>PBT Matches:</b><br>", PBT_info), "")
        )
      )
    if (length(input$selectIDs) < 1) {
      return(df)
    } else {
      parsed_IDs <- input$selectIDs %>%
        str_split_1(",\\s?")
      filtered_df <- df %>%
        filter(ID %in% parsed_IDs)
      if (nrow(filtered_df) == 0) {
        return(df)
      } else {
        return(filtered_df)
      }
    }
  })

  observe({ # Update selectInput choices when data is loaded
    req(data())
    choices <- c( # Valid and allowable choice column names
      "Genetic_Assignment",
      "Sherlock_Assignment",
      "Spring_tributary",
      "Facility",
      "rosa_informed",
      "PBT"
    )
    valid_choices <- choices[choices %in% names(data())]
    updateSelectInput(session, "color", choices = valid_choices)

    model_choices <- data() %>%
      pull(LAD.Class) %>%
      unique() %>%
      sort()
    updateSelectInput(session, "LAD", choices = model_choices, selected = model_choices)

    updateDateRangeInput(
      session,
      "date_range",
      start = as.Date(min(data()$SampleDate, na.rm = TRUE), tz = "America/Los_Angeles"),
      end = as.Date(max(data()$SampleDate, na.rm = TRUE), tz = "America/Los_Angeles"),
      min = as.Date(min(data()$SampleDate, na.rm = TRUE), tz = "America/Los_Angeles"),
      max = as.Date(max(data()$SampleDate, na.rm = TRUE), tz = "America/Los_Angeles")
    )
  })

  p <- reactive({
    req(data(), input$color)
    filtered_data <- data() %>%
      filter(
        LAD.Class %in% input$LAD
      )
    if (input$color == "Spring_tributary") {
      filtered_data <- filtered_data %>%
        filter(
          !is.na(Spring_tributary)
        )
    }
    prerun_color_palette <- RColorBrewer::brewer.pal(9, "Set1")
    run_color_mapping <- c(
      "Fall / Late Fall" = prerun_color_palette[5],
      "Fall or Late fall" = prerun_color_palette[5],
      "Fall" = "#c9b428",
      "Late Fall" = prerun_color_palette[7],
      "Spring" = prerun_color_palette[3],
      "Winter" = prerun_color_palette[2],
      "f22" = prerun_color_palette[6],
      "f23" = prerun_color_palette[6],
      "f24" = prerun_color_palette[6],
      "f25" = prerun_color_palette[6],
      "f26" = prerun_color_palette[6],
      "lf23" = prerun_color_palette[7],
      "s23" = prerun_color_palette[3],
      "s24" = prerun_color_palette[3],
      "s25" = prerun_color_palette[3],
      "s26" = prerun_color_palette[3],
      "w23" = prerun_color_palette[2],
      "w24" = prerun_color_palette[2],
      "w25" = prerun_color_palette[2],
      "w26" = prerun_color_palette[2],
      "lf24" = prerun_color_palette[7],
      "lf25" = prerun_color_palette[7],
      "lf26" = prerun_color_palette[7],
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
    lad_long_filtered <- lad_long %>%
      filter(
        date >= input$date_range[1] - days(1),
        date <= input$date_range[2] + days(1) # Filter based on date range
      )

    p <- plot_ly(colors = run_color_mapping) %>%
      add_ribbons( # Add ribbons for min/max values overlaid on top of fork length
        data = lad_long_filtered,
        x = ~date,
        ymin = ~min,
        ymax = ~max,
        color = ~cohort,
        colors = run_color_mapping,
        showlegend = TRUE,
        alpha = 0.2,
        size = 0.1,
        line = list(width = 0)
      ) %>%
      add_markers(
        data = filtered_data,
        x = ~SampleDate,
        y = ~ForkLength,
        color = ~ get(input$color),
        colors = run_color_mapping,
        text = ~tooltip_text,
        hoverinfo = "text",
        size = I(20),
        alpha = 0.5
      ) %>%
      layout(
        title = "Salvage Chinook Salmon",
        xaxis = list(title = "Sample Date", range = c(input$date_range[1] - days(3), input$date_range[2] + days(3))),
        yaxis = list(title = "Fork Length"),
        showlegend = TRUE
      )

    p
  })

  output$plot <- renderPlotly({
    p()
  })
}


# Run the application
shinyApp(ui = ui, server = server)
