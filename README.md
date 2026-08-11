# salvage_viewer_v3

Deploy with:

```r
rsconnect::deployApp(
  appName = "salvage_viewer_v3", 
  server = "shinyapps.io", 
  logLevel = "verbose",
  appFiles = c(
    "app.R",
    "data/CFS_data.csv",
    list.files("data", pattern = "^lad_long.*\\.txt$", full.names = TRUE),
    list.files("data/PBT_results", pattern = "\\.tsv$", full.names = TRUE)
  )
)
```