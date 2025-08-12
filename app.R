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
  sidebarLayout(
    sidebarPanel(
      selectInput("analysis_type", "Select Analysis Type:",
                  choices = c("DEG","GSEA","GSEA_reduced")),
      
      uiOutput("comparison_info"),
      selectInput("comparison",   "Select Comparison:", choices = NULL),
      
      selectInput("group1", "Filter by Group 1:", choices = c("All")),
      selectInput("group2", "Filter by Group 2:", choices = c("All")),
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
  
  output$plotInProgress <- reactive(TRUE)
  outputOptions(output, "plotInProgress", suspendWhenHidden=FALSE)
  
  # --- comparison helpers (classification + info + auto-set filter-to-All) ---
  classify_comparison <- function(x) {
    x <- tolower(x %||% "")
    if (grepl("^(\\d+_\\d+)$", x)) return("time")
    if (grepl("^(blood|lung|mln)_(blood|lung|mln)$", x) ||
        grepl("(blood|lung|mln)_", x) || grepl("_(blood|lung|mln)$", x)) return("tissue")
    if (grepl("mtb", x) || grepl("mtb-.*mtb\\+", x) || grepl("neg.*pos|pos.*neg", x)) return("infection")
    if (grepl("gobp_time", x)) return("time")
    if (grepl("gobp_tissue", x)) return("tissue")
    if (grepl("gobp_mtb", x)) return("infection")
    "unknown"
  }
  
  comparison_help <- function(type) {
    switch(type,
           "time" = HTML(
             "<b>Time comparison:</b> contrasts samples from two days (e.g., 0 vs 15). ",
             "Tip: set <i>Filter by Day</i> to <b>All</b> so both timepoints are included."
           ),
           "tissue" = HTML(
             "<b>Tissue comparison:</b> contrasts tissues (e.g., blood vs lung). ",
             "Tip: set <i>Filter by Tissue</i> to <b>All</b> so both tissues are included."
           ),
           "infection" = HTML(
             "<b>mtb (infection) comparison:</b> contrasts Mtb-infected (<b>pos</b>) vs uninfected (<b>neg</b>). ",
             "Valid for <b>lung</b> and <b>mln</b> tissues. ",
             "Tip: set <i>Filter by Infection Status</i> to <b>All</b> so both groups are included."
           ),
           HTML(
             "<b>Mixed/unknown comparison:</b> If this contrasts groups on a specific axis ",
             "(time, tissue, or infection), set that axis's filter to <b>All</b>."
           )
    )
  }
  
  output$comparison_info <- renderUI({
    req(input$comparison)
    typ <- classify_comparison(input$comparison)
    div(
      style = "padding:8px; margin-bottom:6px; border-left:4px solid #3498db; background:#f5f9ff;",
      comparison_help(typ)
    )
  })
  
  observeEvent(input$comparison, {
    req(input$comparison)
    typ <- classify_comparison(input$comparison)
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
    typ <- classify_comparison(input$comparison)
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
  
  parse_groups_from_comparison <- function(comp, typ) {
    comp <- tolower(comp %||% "")
    if (typ == "time" || grepl("gobp_time", comp)) {
      toks <- unlist(strsplit(comp, "_"))
      toks <- toks[grepl("^-?\\d+$", toks)]
      if (length(toks) == 0) toks <- c("0","15","28")
      unique(toks)
    } else if (typ == "tissue" || grepl("gobp_tissue", comp)) {
      valid <- c("blood","lung","mln")
      toks <- unlist(strsplit(comp, "_"))
      toks <- toks[toks %in% valid]
      if (length(toks) == 0) toks <- valid
      unique(toks)
    } else if (typ == "infection" || grepl("gobp_mtb", comp)) {
      c("pos","neg")
    } else {
      character()
    }
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
  
  observeEvent(list(input$comparison, input$analysis_type), {
    req(input$comparison)
    typ <- classify_comparison(input$comparison)
    dat <- get_raw_data()
    valid <- valid_groups_from_dat(dat, typ, input$comparison)
    if (length(valid) == 0) {
      updateSelectInput(session, "group1", choices = c("All"), selected = "All")
      updateSelectInput(session, "group2", choices = c("All"), selected = "All")
    } else {
      updateSelectInput(session, "group1", choices = c("All", valid), selected = "All")
      updateSelectInput(session, "group2", choices = c("All", valid), selected = "All")
    }
  }, ignoreInit = FALSE)
  
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
  
  output$data_table <- renderDT({
    datatable(get_data(), options=list(pageLength=10), filter='top')
  })
  
  make_bubble <- function(dat, cols){
    pval_col   <- cols$pval
    effect_col <- cols$effect
    x_lab <- if (!is.na(effect_col) && effect_col == "NES") "Normalized Enrichment Score (NES)" else effect_col
    ggplot(dat, aes(x=.data[[effect_col]], y=`..label`)) +
      geom_point(aes(size=-log10(.data[[pval_col]]), color=.data[[effect_col]])) +
      scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
      scale_size_continuous(range=c(2,10), name="-log10(adj p)") +   # size legend label
      theme_minimal(base_size=input$font_size) +
      labs(title=unique(na.omit(dat$comparison)),
           x=x_lab, y="Term",
           color=effect_col) +
      theme(plot.title=element_text(hjust=0))
  }
  
  output$main_plot <- renderPlotly({
    validate(need(!is.null(input$plot_type), "Choose a plot type."))
    dat <- get_data(); req(nrow(dat)>0)
    cols <- attr(dat, "cols"); req(!is.null(cols))
    pval_col     <- cols$pval
    effect_col   <- cols$effect
    main_col     <- cols$main
    celltype_col <- cols$celltype
    group_axis   <- cols$group_axis
    
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
                  # X axis: best available grouping (celltype->group1->group2->comparison_family->tissue/day/Mtb)
                  x_var <- group_axis %||% main_col
                  ggplot(dat, aes(x=.data[[x_var]], y=if (input$analysis_type=="DEG") .data[[main_col]] else `..label`)) +
                    geom_point(aes(size=-log10(.data[[pval_col]]), color=.data[[effect_col]])) +
                    scale_color_gradient2(low="blue",mid="white",high="red",midpoint=0) +
                    scale_size_continuous(range=c(2,10), name="-log10(adj p)") +  # size legend label
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
                  make_bubble(dat, cols)
                }
    )
    
    last_ggplot(p)
    # Use the text aesthetic (mapped to ..label) so genes/terms always show
    tips <- c("text", effect_col, pval_col)
    if (!is.na(group_axis)) tips <- c(tips, group_axis)
    if (!is.na(celltype_col)) tips <- c(tips, celltype_col)
    ggplotly(p, tooltip = tips)
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
