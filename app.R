# Full integrated Shiny app for DEG + GSEA viewing with reduced GSEA support (using merged CSVs)
# Tuned to the provided schemas.

library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(DT)
library(tibble)
library(tidyr)
library(purrr)
library(scales)
library(plotly)
library(openxlsx)
library(stringr)

`%||%` <- function(a, b) if (!is.null(a)) a else b

ui <- fluidPage(
  titlePanel("Mouse Mtb Innate Atlas DEG and GO:BP GSEA Viewer"),
  # --- Collapsible Acknowledgements / References (open by default) ---
  tags$head(tags$style(HTML("
  .ack { background:#f7f9fc; border:1px solid #e3e7ef; border-radius:8px; padding:10px 14px; margin-bottom:14px; }
  .ack > summary { font-weight:600; cursor:pointer; list-style:none; }
  .ack > summary::-webkit-details-marker { display:none; }
"))),
  tags$details(open = "open", class = "ack",
               tags$summary("Acknowledgements & References"),
               HTML(
                 "Made in collaboration between the Gerner (University of Washington, Department of Immunology) ",
                 "and Urdahl (Seattle Childrens' Research Institute, CGIDR) Laboratories.<br><br>",
                 "For use in submissions and publications, as well as further information regarding methods and analysis, please reference:<br>",
                 "<em>\"Monocytic Niches Enable Mycobacterium tuberculosis Persistence in Lymph Nodes\"</em><br>",
                 "(Manuscript in Progress)<br>",
                 "Elya A. Shamskhou, Fergal R. Duffy, Lauren M. Cross, Vitaly V. Ganusov, ",
                 "Courtney R. Plumlee, Benjamin H. Gern, Alan H. Dierks, Sara B. Cohen, ",
                 "Kevin B. Urdahl*, Michael Y. Gerner*"
               )),
  # --- Tab styling (lavender background + bold active tab) ---
  # --- Tab + Acknowledgements styling ---
  tags$head(
    tags$style(HTML("
    /* --- Tab headers --- */
    .nav-tabs > li > a {
      background-color: #f3e8ff;  /* light lavender */
      color: #333;
      border-radius: 6px 6px 0 0;
      margin-right: 4px;
    }
    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:focus,
    .nav-tabs > li.active > a:hover {
      background-color: #e9d5ff; /* darker lavender when active */
      color: #000;
      font-weight: bold;
    }
    .nav-tabs > li > a:hover {
      background-color: #ede9fe;
    }

    /* --- Acknowledgements block --- */
    details.ack {
      background-color: #f3e8ff;   /* light lavender closed */
      border: 1px solid #e3e7ef;
      border-radius: 8px;
      padding: 10px 14px;
      margin-bottom: 14px;
      transition: background-color 0.2s ease;
    }
    details.ack[open] {
      background-color: #e9d5ff;   /* darker lavender when expanded */
    }
    details.ack summary {
      font-weight: 600;
      cursor: pointer;
      list-style: none;
    }
    details.ack summary::-webkit-details-marker {
      display: none; /* hide default triangle marker */
    }
    details.ack:hover {
      background-color: #ede9fe;   /* hover shade */
    }
  "))
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput("analysis_type", "Select Analysis Type:",
                  choices = c("DEG","GSEA","GSEA_reduced")),
      
      uiOutput("comparison_info"),
      selectInput("comparison",   "Select Comparison:", choices = NULL),
      
      selectInput("group1", "Filter by Group 1:", choices = c("All")),
      selectInput("group2", "Filter by Group 2:", choices = c("All")),
      uiOutput("axis_badge"),  # <— NEW: shows which axis (time/tissue/infection) is active
      
      selectInput("celltype",     "Select Cell Type:", choices = NULL, multiple=TRUE),
      
      selectInput("plot_type",    "Select Plot Type:",
                  choices = c("Volcano Plot","Dot Plot","Bar Plot")),  # server updates per analysis_type
      
      selectInput("tissue",       "Filter by Tissue:", choices = c("All")),
      selectInput("condition",    "Filter by Infection Status:", choices = c("All")),
      selectInput("day",          "Filter by Day:", choices = c("All")),
      sliderInput("font_size",    "Font Size:", min=2, max=20, value=12),
      numericInput("pval_thresh", "P-Value Threshold:",  value=0.05, min=0, max=1, step=0.01),
      numericInput("effect_thresh","Fold Change / NES Threshold:", value=1, min=0, step=0.5),
      textInput("gene_filter",    "Filter by Gene/Pathway (comma-separated or substrings):", value=""),
      selectInput("sort_by",      "Sort by:", 
                  c("P-value"="pval","Effect"="effect","Name"="name")),
      textAreaInput("manual_order","Manual gene/pathway order (comma-separated):","", rows=3),
      downloadButton("download_plot", "Download Plot"),
      downloadButton("download_data","Download Data")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Table", DTOutput("data_table")),
        tabPanel("Plot", 
                 plotlyOutput("main_plot", height="600px"),
                 conditionalPanel("output.plotInProgress", tags$p("⏳ Generating plot…"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  last_ggplot <- reactiveVal(NULL)
  plotting <- reactiveVal(FALSE)
  
  base_folders <- list(
    DEG          = "merged_DEG_results",
    GSEA         = "merged_GSEA_results",
    GSEA_reduced = "merged_GSEA_reduced_results"
  )
  
  # Plot choices adapt to analysis type (GSEA also allows Volcano)
  observeEvent(input$analysis_type, {
    typ <- input$analysis_type %||% "DEG"
    choices <- if (typ == "DEG") {
      c("Volcano Plot","Dot Plot","Bar Plot")
    } else {
      c("Volcano Plot","Dot Plot","Bar Plot","Bubble Plot")
    }
    updateSelectInput(session, "plot_type", choices = choices, selected = choices[1])
  }, ignoreInit = FALSE)
  
  # --- populate comparison & metadata menus ---
  observe({
    folder <- base_folders[[input$analysis_type]]
    validate(need(dir.exists(folder), sprintf("Folder '%s' not found.", folder)))
    fns    <- list.files(folder, pattern="\\.csv$", full.names = FALSE)
    validate(need(length(fns) > 0, sprintf("No CSVs found in '%s'.", folder)))
    
    prefix <- paste0("^merged_", input$analysis_type, "_")
    comps  <- fns %>%
      sub(prefix, "", .) %>%
      sub("\\.csv$", "", .) %>%
      sort()
    updateSelectInput(session, "comparison", choices = comps, selected = comps[1])
    
    all_data <- purrr::map_dfr(
      file.path(folder, fns),
      ~ readr::read_csv(.x, show_col_types = FALSE) %>%
        mutate(across(any_of(
          c("celltype","tissue","Mtb","day",
            "group1","group2","comparison","comparison_family","term","parentTerm")
        ), as.character))
    )
    
    updateSelectInput(session, "celltype",
                      choices = sort(unique(all_data$celltype)))
    updateSelectInput(session, "tissue",
                      choices = c("All", sort(na.omit(unique(all_data$tissue)))),
                      selected = "All")
    updateSelectInput(session, "day",
                      choices = c("All", sort(na.omit(unique(all_data$day)))),
                      selected = "All")
    updateSelectInput(session, "condition",
                      choices = c("All", sort(na.omit(unique(all_data$Mtb)))),
                      selected = "All")
    updateSelectInput(session, "group1",
                      choices = c("All", sort(na.omit(unique(all_data$group1)))),
                      selected = "All")
    updateSelectInput(session, "group2",
                      choices = c("All", sort(na.omit(unique(all_data$group2)))),
                      selected = "All")
  })
  
  output$plotInProgress <- reactive({ plotting() })
  outputOptions(output, "plotInProgress", suspendWhenHidden = FALSE)
  
  # --- helpers to detect token types from group1/group2 values ---
  .is_numeric_level <- function(x) { all(grepl("^\\s*-?\\d+\\s*$", x)) }
  .normalize_infection <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub("mtb\\+", "pos", x)
    x <- gsub("mtb-",  "neg", x)
    x
  }
  .is_infection_levels <- function(x) {
    y <- .normalize_infection(x)
    any(y %in% c("pos","neg","infected","uninfected","positive","negative","control"))
  }
  .is_tissue_levels <- function(x) {
    any(tolower(trimws(x)) %in% c("blood","lung","mln"))
  }
  
  # Dynamic sizing helpers
  compute_plot_height <- function(n_rows, base = 380, per = 22, max_h = 2400) {
    # total height = base + per * rows, but capped so it doesn't get absurdly tall
    min(max_h, base + per * n_rows)
  }
  size_range_for <- function(n_rows) {
    # tighter sizes when there are many rows
    if (n_rows >= 60) return(c(1.5, 5))
    if (n_rows >= 40) return(c(1.8, 6))
    if (n_rows >= 25) return(c(2, 8))
    c(2, 10)
  }
  
  # --- robust classifier: use comparison_family > group1/2 tokens > filename > soft fallback ---
  classify_comparison <- function(comp_name, dat = NULL) {
    s <- tolower(comp_name %||% "")
    
    # 1) explicit comparison_family wins
    if (!is.null(dat) && "comparison_family" %in% names(dat)) {
      fam <- unique(na.omit(tolower(as.character(dat$comparison_family))))
      fam <- gsub("^mtb$", "infection", fam)
      fam <- fam[fam %in% c("time","tissue","infection")]
      if (length(fam) == 1) return(fam)
    }
    
    # 2) infer from group1/group2 token types (most reliable for DEG)
    if (!is.null(dat) && any(c("group1","group2") %in% names(dat))) {
      gvals <- unique(na.omit(as.character(c(dat$group1, dat$group2))))
      if (length(gvals)) {
        if (.is_tissue_levels(gvals)) return("tissue")
        if (.is_infection_levels(gvals)) return("infection")
        if (.is_numeric_level(gvals)) return("time")
        # mixed tokens—try to see which bucket dominates
        tis_hit <- mean(tolower(trimws(gvals)) %in% c("blood","lung","mln"))
        inf_hit <- mean(.normalize_infection(gvals) %in% c("pos","neg","infected","uninfected","positive","negative","control"))
        num_hit <- mean(grepl("^\\s*-?\\d+\\s*$", gvals))
        if (max(tis_hit, inf_hit, num_hit) > 0) {
          return(c("tissue","infection","time")[which.max(c(tis_hit, inf_hit, num_hit))])
        }
      }
    }
    
    # 3) filename-based hints
    if (grepl("^\\d+([_-]|vs)\\d+$", s) || grepl("^(\\d+[_-]\\d+)$", s)) return("time")
    if (grepl("^(blood|lung|mln)[^a-z0-9]+(blood|lung|mln)$", s)) return("tissue")
    if (grepl("gobp_time", s)) return("time")
    if (grepl("gobp_tissue", s)) return("tissue")
    if (grepl("gobp_mtb", s) || grepl("mtb", s) || grepl("neg|pos", s)) return("infection")
    
    # 4) soft fallback (last resort; much less aggressive than before)
    if (!is.null(dat)) {
      get_vals <- function(col) if (col %in% names(dat)) unique(na.omit(as.character(dat[[col]]))) else character()
      # Prefer variation in group columns over whole-table columns
      gvals <- unique(na.omit(as.character(c(dat$group1, dat$group2))))
      if (.is_tissue_levels(gvals)) return("tissue")
      if (.is_infection_levels(gvals)) return("infection")
      if (.is_numeric_level(gvals)) return("time")
    }
    "unknown"
  }
  
  
  comparison_help <- function(type) {
    switch(type,
           "time" = HTML(
             "<b>Time comparison: </b> Possible group comparisons are <b>15 v 0, 15 v 28, 28 v 0</b>. ",
             "Set <i>Filter by Day</i> to <b>All</b> since you are comparing timepoints.",
             "Tip: The first group is UP, while the second group DOWN"
           ),
           "tissue" = HTML(
             "<b>Tissue comparison:</b> for tissues affected by aerosol Mtb infection- <i>lung, mln (mediastinal lymph node), and blood</i>",
             "Possible group comparisons are <b>lung vs mln, blood v mln, blood v lung <'b>. ",
             "Set <i>Filter by Tissue</i> to <b>All</b> since you are comparing tissues",
             "Tip: The first group is UP, while the second group DOWN"
           ),
           "infection" = HTML(
             "<b>mtb (infection) comparison:</b> contrasts Mtb-infected (<b>pos</b>) vs uninfected (<b>neg</b>). ",
             "Valid for <b>lung</b> and <b>mln</b> tissues only. ",
             "Set <i>Filter by Infection Status</i> to <b>All</b> since you are comparing pos vs neg",
             "Tip: The first group is UP, while the second group DOWN"
           ),
           HTML(
             "<b>Mixed/unknown comparison:</b> If this contrasts groups on a specific axis ",
             "(time, tissue, or infection), set that axis's filter to <b>All</b>."
           )
    )
  }
  
  # More forgiving parser (supports 0_15, 0-15, 0vs15, blood-lung...)
  parse_groups_from_comparison <- function(comp, typ) {
    s <- tolower(comp %||% "")
    toks <- unlist(strsplit(s, "[^a-z0-9]+"))
    toks <- toks[nzchar(toks)]
    if (typ == "time" || grepl("gobp_time", s)) {
      dg <- toks[grepl("^-?\\d+$", toks)]
      if (length(dg) == 0) dg <- c("0","15","28")
      unique(dg)
    } else if (typ == "tissue" || grepl("gobp_tissue", s)) {
      valid <- c("blood","lung","mln")
      tt <- toks[toks %in% valid]
      if (length(tt) == 0) tt <- valid
      unique(tt)
    } else if (typ == "infection" || grepl("gobp_mtb", s)) {
      c("pos","neg")
    } else {
      character()
    }
  }
  
  # Canonicalizers for group1/group2 levels
  canon_infection <- function(v) {
    if (is.null(v)) return(character())
    v <- tolower(trimws(v))
    v <- gsub("mtb\\+", "pos", v)
    v <- gsub("mtb-",  "neg", v)
    v <- gsub("^positive$|^infected$", "pos", v)
    v <- gsub("^negative$|^uninfected$|^control$", "neg", v)
    unique(v[v %in% c("pos","neg")])
  }
  canon_tissue <- function(v) {
    if (is.null(v)) return(character())
    v <- tolower(trimws(v))
    valid <- c("blood","lung","mln")
    unique(v[v %in% valid])
  }
  canon_day <- function(v) {
    if (is.null(v)) return(character())
    v <- as.character(v)
    v <- gsub("\\s+", "", v)
    v <- v[grepl("^-?\\d+$", v)]
    sort(unique(v), na.last = NA)
  }
  
  get_raw_data <- reactive({
    req(input$comparison)
    folder <- base_folders[[input$analysis_type]]
    file   <- file.path(folder, paste0("merged_", input$analysis_type, "_", input$comparison, ".csv"))
    validate(need(file.exists(file), sprintf("No such file: %s", basename(file))))
    readr::read_csv(file, show_col_types = FALSE) %>%
      mutate(across(any_of(c("celltype","tissue","Mtb","day","group1","group2","term","parentTerm","comparison_family")), as.character))
  })
  
  valid_groups_from_dat <- function(dat, typ, comp) {
    if (typ == "time") {
      pool <- c(dat$group1, dat$group2, dat$day)
      vg <- canon_day(pool)
      if (length(vg) == 0) vg <- parse_groups_from_comparison(comp, typ)
      return(vg)
    }
    if (typ == "infection") {
      pool <- c(dat$group1, dat$group2, dat$Mtb)
      vg <- canon_infection(pool)
      if (length(vg) == 0) vg <- parse_groups_from_comparison(comp, typ)
      return(vg)
    }
    if (typ == "tissue") {
      pool <- c(dat$group1, dat$group2, dat$tissue)
      vg <- canon_tissue(pool)
      if (length(vg) == 0) vg <- parse_groups_from_comparison(comp, typ)
      return(vg)
    }
    character()
  }
  
  # --- Axis badge (small pill under Group 1/2) ---
  render_axis_badge <- function(typ) {
    col <- switch(typ, time="#6ab04c", tissue="#0984e3", infection="#e17055", "#b2bec3")
    lbl <- switch(typ, time="time", tissue="tissue", infection="infection", "unknown")
    div(
      style = sprintf("margin:-6px 0 8px 0;"),
      span(
        lbl,
        style = sprintf(
          "display:inline-block;padding:2px 8px;border-radius:999px;font-size:12px;
           color:white;background:%s;letter-spacing:.5px;", col
        )
      )
    )
  }
  output$axis_badge <- renderUI({
    req(input$comparison)
    typ <- classify_comparison(input$comparison, tryCatch(get_raw_data(), error = function(...) NULL))
    render_axis_badge(typ)
  })
  
  output$comparison_info <- renderUI({
    req(input$comparison)
    typ <- classify_comparison(input$comparison, tryCatch(get_raw_data(), error = function(...) NULL))
    div(
      style = "padding:8px; margin-bottom:6px; border-left:4px solid #3498db; background:#f5f9ff;",
      comparison_help(typ)
    )
  })
  
  observeEvent(input$comparison, {
    req(input$comparison)
    typ <- classify_comparison(input$comparison, tryCatch(get_raw_data(), error = function(...) NULL))
    if (typ == "time") {
      updateSelectInput(session, "day", selected = "All")
      showNotification("Time comparison selected — including both timepoints (Filter by Day → All).", type="message", duration=4)
    } else if (typ == "tissue") {
      updateSelectInput(session, "tissue", selected = "All")
      showNotification("Tissue comparison selected — including both tissues (Filter by Tissue → All).", type="message", duration=4)
    } else if (typ == "infection") {
      updateSelectInput(session, "condition", selected = "All")
      showNotification("Infection comparison selected — including pos & neg (Filter by Infection Status → All).", type="message", duration=4)
    }
  }, ignoreInit = FALSE)
  
  observeEvent(list(input$comparison, input$tissue, input$condition, input$day), {
    req(input$comparison)
    typ <- classify_comparison(input$comparison, tryCatch(get_raw_data(), error = function(...) NULL))
    if (typ == "time" && !identical(input$day, "All")) {
      showNotification("Reminder: for a time comparison, 'Filter by Day' should be 'All'.", type = "warning", duration = 5)
    }
    if (typ == "tissue" && !identical(input$tissue, "All")) {
      showNotification("Reminder: for a tissue comparison, 'Filter by Tissue' should be 'All'.", type = "warning", duration = 5)
    }
    if (typ == "infection") {
      if (!identical(input$condition, "All")) {
        showNotification("Reminder: for an infection comparison, 'Filter by Infection Status' should be 'All'.", type = "warning", duration = 5)
      }
      if (identical(input$tissue, "blood")) {
        showNotification("Note: mtb comparisons are only valid for lung and mln tissues.", type = "warning", duration = 6)
      }
    }
  }, ignoreInit = FALSE)
  
  # --- column detection helpers ---
  detect_col <- function(nms, candidates) {
    cand <- candidates[candidates %in% nms]
    if (length(cand)) cand[1] else NA_character_
  }
  pick_gsea_label_col <- function(nms) {
    detect_col(nms, c("parentTerm","term","Description","pathway","gs_name","ID","Name","name","title"))
  }
  
  # choose best available grouping column for dot/bar (x-axis)
  pick_group_axis <- function(nms) {
    detect_col(nms, c("celltype","group1","group2","comparison_family","tissue","day","Mtb"))
  }
  
  # --- load & filter merged data according to inputs ---
  get_data <- reactive({
    req(input$comparison)
    folder <- base_folders[[input$analysis_type]]
    file   <- file.path(folder, paste0("merged_", input$analysis_type, "_", input$comparison, ".csv"))
    validate(need(file.exists(file), sprintf("No such file: %s", basename(file))))
    
    dat <- read_csv(file, show_col_types=FALSE) %>%
      mutate(across(any_of(c("tissue","Mtb","day","comparison","celltype","group1","group2","term","parentTerm","comparison_family")), as.character))
    
    # Filters (skip if columns absent)
    if (input$group1   != "All" && "group1" %in% names(dat)) dat <- dat %>% filter(group1 == input$group1)
    if (input$group2   != "All" && "group2" %in% names(dat)) dat <- dat %>% filter(group2 == input$group2)
    if (input$tissue   != "All" && "tissue" %in% names(dat)) dat <- dat %>% filter(tissue==input$tissue)
    if (input$condition!= "All" && "Mtb"    %in% names(dat)) dat <- dat %>% filter(Mtb==input$condition)
    if (input$day      != "All" && "day"    %in% names(dat)) dat <- dat %>% filter(day==input$day)
    if (length(input$celltype)>0 && "celltype" %in% names(dat)) dat <- dat %>% filter(celltype %in% input$celltype)
    
    nms <- names(dat)
    # Column detection tuned to your schemas
    pval_col   <- if (input$analysis_type=="DEG") detect_col(nms, c("p_val_adj")) else detect_col(nms, c("padj"))
    if (is.na(pval_col)) pval_col <- detect_col(nms, c("p.adjust","qval","pval","p_value"))
    effect_col <- if (input$analysis_type=="DEG") detect_col(nms, c("avg_log2FC","log2FC","logFC")) else detect_col(nms, c("NES","enrichmentScore","ES"))
    main_col   <- if (input$analysis_type=="DEG") detect_col(nms, c("gene","Gene","symbol")) else detect_col(nms, c("pathway","term","Description","gs_name","ID"))
    parent_col <- detect_col(nms, c("parentTerm","parent","parent_term"))
    celltype_col <- detect_col(nms, c("celltype","CellType","cell_type"))
    gsea_label_col <- if (input$analysis_type=="DEG") NA_character_ else pick_gsea_label_col(nms)
    
    validate(
      need(!is.na(main_col), "Could not find a gene/pathway column."),
      need(!is.na(pval_col), "Could not find an adjusted p-value column."),
      need(!is.na(effect_col), "Could not find an effect (log2FC/NES) column.")
    )
    
    # Thresholds
    if (pval_col   %in% names(dat)) dat <- dat %>% filter(.data[[pval_col]]   <= input$pval_thresh)
    if (effect_col %in% names(dat)) dat <- dat %>% filter(abs(.data[[effect_col]]) >= input$effect_thresh)
    
    # Flexible name filter
    if (str_trim(input$gene_filter) != "") {
      tokens <- str_split(input$gene_filter, ",")[[1]] %>% str_trim() %>% discard(~ .x == "")
      rx <- paste0("(", paste(stringr::str_replace_all(tokens, "([\\W])", "\\\\\\1"), collapse="|"), ")")
      dat <- dat %>% filter(str_detect(!!sym(main_col), regex(rx, ignore_case = TRUE)))
    }
    
    # Guaranteed label for hover/axes (DEG: gene; GSEA: parentTerm -> term -> Description -> pathway)
    if (input$analysis_type == "DEG") {
      dat <- dat %>% mutate(`..label` = .data[[main_col]])
    } else {
      lab_col <- if (!is.na(parent_col) && parent_col %in% names(dat)) parent_col
      else if (!is.na(gsea_label_col) && gsea_label_col %in% names(dat)) gsea_label_col
      else main_col
      dat <- dat %>% mutate(`..label` = .data[[lab_col]])
    }
    
    # Sorting
    sort_col <- switch(input$sort_by,
                       pval   = pval_col,
                       effect = effect_col,
                       name   = if (input$analysis_type=="DEG") main_col else "..label")
    dat <- dat %>% arrange(.data[[sort_col]])
    
    # Manual order
    if (input$manual_order!="") {
      ord <- str_trim(str_split(input$manual_order,",")[[1]])
      ord_col <- if (input$analysis_type=="DEG") main_col else "..label"
      if (ord_col %in% names(dat)) {
        dat[[ord_col]] <- factor(dat[[ord_col]], levels=ord)
        dat <- dat %>% arrange(.data[[ord_col]])
      }
    }
    validate(need(nrow(dat) > 0, "No rows after filtering. Loosen filters or adjust thresholds."))
    
    attr(dat, "cols") <- list(
      pval=pval_col, effect=effect_col, main=main_col,
      parent=parent_col, celltype=celltype_col, gsea_label=gsea_label_col,
      group_axis = pick_group_axis(nms)
    )
    dat
  })
  
  # Use smarter classifier anywhere comparison "type" matters
  observeEvent(list(input$comparison, input$analysis_type), {
    req(input$comparison)
    dat <- get_raw_data()
    typ <- classify_comparison(input$comparison, dat)
    valid <- unique(valid_groups_from_dat(dat, typ, input$comparison))
    if (length(valid) == 0) {
      updateSelectInput(session, "group1", choices = c("All"), selected = "All")
      updateSelectInput(session, "group2", choices = c("All"), selected = "All")
    } else {
      updateSelectInput(session, "group1", choices = c("All", valid), selected = "All")
      updateSelectInput(session, "group2", choices = c("All", valid), selected = "All")
    }
  }, ignoreInit = FALSE)
  
  output$data_table <- renderDT({
    dat <- get_data()
    
    # 1) Drop internal label columns used only for plotting hovers
    drop_candidates <- c("..label", "label")
    drop_candidates <- drop_candidates[drop_candidates %in% names(dat)]
    if (length(drop_candidates)) {
      dat <- dat[, setdiff(names(dat), drop_candidates), drop = FALSE]
    }
    
    # 2) For GSEA / GSEA_reduced, move leading-edge style columns to the far right
    if (input$analysis_type %in% c("GSEA", "GSEA_reduced")) {
      le_candidates <- c(
        "leadingEdge", "leading_edge", "leading.edge",
        "core_enrichment", "coreEnrichment", "leadingGenes", "leading_genes"
      )
      le_present <- le_candidates[le_candidates %in% names(dat)]
      if (length(le_present)) {
        ord <- c(setdiff(names(dat), le_present), le_present)
        dat <- dat[, ord, drop = FALSE]
      }
    }
    
    datatable(
      dat,
      options = list(pageLength = 10, scrollX = TRUE),
      filter = "top"
    )
  })

  
  make_bubble <- function(dat, cols, size_range){
    pval_col   <- cols$pval
    effect_col <- cols$effect
    x_lab <- if (!is.na(effect_col) && effect_col == "NES") "Normalized Enrichment Score (NES)" else effect_col
    ggplot(dat, aes(x=.data[[effect_col]], y=`..label`)) +
      geom_point(aes(size=-log10(.data[[pval_col]]), color=.data[[effect_col]])) +
      scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
      scale_size_continuous(range=size_range, name="-log10(adj p)") +
      theme_minimal(base_size=input$font_size) +
      labs(title=unique(na.omit(dat$comparison)),
           x=x_lab, y="Term",
           color=effect_col) +
      theme(plot.title=element_text(hjust=0))
  }
  
  output$main_plot <- renderPlotly({
    plotting(TRUE)                         # turn ON when rendering starts
    on.exit(plotting(FALSE), add = TRUE)   # guarantee OFF when it finishes (success or error)
    
    validate(need(!is.null(input$plot_type), "Choose a plot type."))
    dat <- get_data(); req(nrow(dat) > 0)
    cols <- attr(dat, "cols"); req(!is.null(cols))
    pval_col     <- cols$pval
    effect_col   <- cols$effect
    main_col     <- cols$main
    celltype_col <- cols$celltype
    group_axis   <- cols$group_axis
    # figure out the Y categories for dot/bubble plots
    y_var <- if (input$analysis_type == "DEG") main_col else "..label"
    y_levels <- if (input$plot_type %in% c("Dot Plot","Bubble Plot")) {
      length(unique(na.omit(dat[[y_var]])))
    } else 0
    dyn_height <- if (y_levels > 0) compute_plot_height(y_levels) else 600
    dot_sizes  <- size_range_for(y_levels)
    
    p <- switch(input$plot_type,
                "Volcano Plot" = {
                  if (input$analysis_type=="DEG") {
                    xcol <- detect_col(names(dat), c("avg_log2FC","log2FC")) %||% "avg_log2FC"
                    ycol <- detect_col(names(dat), c("p_val_adj","padj")) %||% "p_val_adj"
                    base <- ggplot(dat, aes(x=.data[[xcol]], y=-log10(.data[[ycol]]), text=`..label`)) +
                      geom_hline(yintercept=-log10(input$pval_thresh), linetype="dashed") +
                      geom_vline(xintercept=c(-input$effect_thresh,input$effect_thresh), linetype="dashed") +
                      labs(title="Volcano Plot", x=xcol, y="-log10(adj p)") +
                      theme_minimal(base_size=input$font_size)
                    if (!is.na(celltype_col)) {
                      base + geom_point(aes(color=.data[[celltype_col]]), alpha=0.6) + labs(color="Cell Type")
                    } else {
                      base + geom_point(alpha=0.6)
                    }
                  } else {
                    ggplot(dat, aes(x=.data[[effect_col]], y=-log10(.data[[pval_col]]), text=`..label`)) +
                      geom_point(aes(color=.data[[effect_col]]), alpha=0.7) +
                      scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
                      geom_hline(yintercept=-log10(input$pval_thresh), linetype="dashed") +
                      geom_vline(xintercept=c(-input$effect_thresh,input$effect_thresh), linetype="dashed") +
                      labs(title="GSEA Volcano Plot",
                           x=if (effect_col=="NES") "Normalized Enrichment Score (NES)" else effect_col,
                           y="-log10(adj p)", color=effect_col) +
                      theme_minimal(base_size=input$font_size)
                  }
                },
                "Dot Plot" = {
                  x_var <- group_axis %||% main_col
                  ggplot(dat, aes(x=.data[[x_var]], y=if (input$analysis_type=="DEG") .data[[main_col]] else `..label`)) +
                    geom_point(aes(size=-log10(.data[[pval_col]]), color=.data[[effect_col]])) +
                    scale_color_gradient2(low="blue",mid="white",high="red",midpoint=0) +
                    scale_size_continuous(range=dot_sizes, name="-log10(adj p)") +
                    labs(x=str_to_title(x_var), y=if(input$analysis_type=="DEG") "Gene" else "Term") +
                    theme_minimal(base_size=input$font_size) +
                    theme(axis.text.x=element_text(angle=45,hjust=1))
                },
                
                "Bar Plot" = {
                  x_var <- if (input$analysis_type=="DEG") main_col else "..label"
                  grp   <- group_axis %||% x_var
                  df <- dat %>%
                    group_by(.data[[x_var]], .data[[grp]]) %>%
                    summarize(mean_eff=mean(.data[[effect_col]],na.rm=TRUE), .groups="drop")
                  ggplot(df, aes(x=.data[[x_var]], y=mean_eff, fill=.data[[grp]]))+
                    geom_col(position="dodge")+
                    theme_minimal(base_size=input$font_size)+
                    theme(axis.text.x=element_text(angle=45,hjust=1)) +
                    labs(x=if(input$analysis_type=="DEG") "Gene" else "Term",
                         y=paste("Mean", effect_col),
                         fill=str_to_title(grp))
                },
                "Bubble Plot" = {
                  validate(need(input$analysis_type!="DEG","Bubble only for GSEA/GSEA_reduced"))
                  make_bubble(dat, cols, size_range = dot_sizes)
                }
    )
    
    last_ggplot(p)
    tips <- c("text", effect_col, pval_col)
    if (!is.na(group_axis)) tips <- c(tips, group_axis)
    if (!is.na(celltype_col)) tips <- c(tips, celltype_col)
    
    ggplotly(p, tooltip = tips) %>%
      layout(height = dyn_height, margin = list(t = 60, r = 20, b = 60, l = 100))
  })
  
  output$download_plot <- downloadHandler(
    filename = function() paste0("plot_",Sys.Date(),".png"),
    content = function(file){
      ggsave(file, last_ggplot(), width=10, height=7, dpi=300)
    }
  )
  
  output$download_data <- downloadHandler(
    filename = function() paste0("data_",Sys.Date(),
                                 if(input$analysis_type=="DEG")".csv" else ".xlsx"),
    content = function(file){
      if(input$analysis_type=="DEG"){
        write_csv(get_data(), file)
      } else {
        openxlsx::write.xlsx(get_data(), file)
      }
    }
  )
}

shinyApp(ui, server)
