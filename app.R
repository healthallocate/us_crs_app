library(shiny)
library(bslib)
library(bsicons)
library(shinyWidgets)
library(tidyverse)
library(here)
library(DT)

options(shiny.maxRequestSize = 150 * 1024^2)

mapping_table <- read.csv(here::here("uscrs_mapping_table.csv"))

# ---- Helpers ----------------------------------------------------------------

# Risk tier metadata for a given integer US-CRS score (1-50)
risk_info <- function(score) {
  if (length(score) == 0 || is.na(score)) {
    return(list(
      color  = "#6c757d",
      theme  = "secondary",
      label  = "Awaiting input",
      desc   = "Enter patient characteristics to calculate score"
    ))
  }
  if (score <= 10) {
    list(color = "#2a9d8f", theme = "success",
         label = "Lower Risk",
         desc  = "Lower estimated waitlist mortality")
  } else if (score <= 25) {
    list(color = "#e9c46a", theme = "warning",
         label = "Moderate Risk",
         desc  = "Moderate estimated waitlist mortality")
  } else if (score <= 40) {
    list(color = "#f4a261", theme = "orange",
         label = "High Risk",
         desc  = "High estimated waitlist mortality")
  } else {
    list(color = "#e76f51", theme = "danger",
         label = "Very High Risk",
         desc  = "Very high estimated waitlist mortality")
  }
}

# Core score computation — identical math to the original app
compute_uscrs <- function(df) {
  df <- df %>%
    mutate(
      eGFR = case_when(
        sex == 1 ~ 142 * (pmin((creatinine / 0.7), 1)^(-0.241)) *
          (pmax((creatinine / 0.7), 1)^(-1.2)) * 0.9938^(age) * 1.012,
        sex == 0 ~ 142 * (pmin((creatinine / 0.9), 1)^(-0.302)) *
          (pmax((creatinine / 0.7), 1)^(-1.2)) * 0.9938^(age)
      ),
      bnp_reg_indicator = ifelse(BNP_NT_Pro != 1, 1, 0),
      uscrs_raw_score = 1.02 * short_MCS_ever + 0.55 * log(bilirubin + 1) -
        0.01 * eGFR +
        0.40 * log(BNP) * bnp_reg_indicator +
        0.20 * log(BNP) * BNP_NT_Pro -
        0.63 * albumin - 0.07 * sodium - 1.12 * LVAD
    )

  breaks <- sort(unique(c(
    mapping_table$min_raw_score,
    max(mapping_table$max_raw_score)
  )))

  df %>%
    mutate(
      uscrs_score = cut(
        uscrs_raw_score,
        breaks = breaks,
        labels = FALSE,
        include.lowest = TRUE
      )
    )
}

# Reference patient values for contribution analysis
# Clinically typical values for a heart failure patient on the waitlist
REF_ALBUMIN    <- 4.0    # g/dL (normal)
REF_BILIRUBIN  <- 1.0    # mg/dL (upper normal)
REF_CREATININE <- 1.0    # mg/dL (normal)
REF_SODIUM     <- 138    # mEq/L (normal)
REF_AGE        <- 50     # years
REF_BNP_REG    <- 150    # pg/mL (moderate for regular BNP)
REF_BNP_NTPRO  <- 900    # pg/mL (comparable for NT-pro BNP)

# Compute each variable's contribution to the raw score relative to
# a reference patient with typical values. Positive = risk-increasing,
# negative = protective.
compute_contributions <- function(albumin, bilirubin, sex, age, creatinine,
                                  sodium, LVAD, short_MCS_ever, BNP_NT_Pro, BNP) {
  # Patient eGFR
  patient_eGFR <- if (sex == 1) {
    142 * (pmin((creatinine / 0.7), 1)^(-0.241)) *
      (pmax((creatinine / 0.7), 1)^(-1.2)) * 0.9938^(age) * 1.012
  } else {
    142 * (pmin((creatinine / 0.9), 1)^(-0.302)) *
      (pmax((creatinine / 0.7), 1)^(-1.2)) * 0.9938^(age)
  }

  # Reference eGFR (same sex, reference age & creatinine)
  ref_eGFR <- if (sex == 1) {
    142 * (pmin((REF_CREATININE / 0.7), 1)^(-0.241)) *
      (pmax((REF_CREATININE / 0.7), 1)^(-1.2)) * 0.9938^(REF_AGE) * 1.012
  } else {
    142 * (pmin((REF_CREATININE / 0.9), 1)^(-0.302)) *
      (pmax((REF_CREATININE / 0.7), 1)^(-1.2)) * 0.9938^(REF_AGE)
  }

  # BNP contribution (same assay type for fair comparison)
  ref_BNP <- if (BNP_NT_Pro == 1) REF_BNP_NTPRO else REF_BNP_REG
  bnp_reg <- if (BNP_NT_Pro != 1) 1 else 0
  patient_bnp_term <- 0.40 * log(BNP) * bnp_reg + 0.20 * log(BNP) * BNP_NT_Pro
  ref_bnp_term     <- 0.40 * log(ref_BNP) * bnp_reg + 0.20 * log(ref_BNP) * BNP_NT_Pro

  data.frame(
    variable = c("Albumin", "Bilirubin", "Kidney Function\n(eGFR)",
                  "Sodium", "BNP", "LVAD", "Short-term MCS"),
    contribution = c(
      -0.63 * (albumin - REF_ALBUMIN),
       0.55 * (log(bilirubin + 1) - log(REF_BILIRUBIN + 1)),
      -0.01 * (patient_eGFR - ref_eGFR),
      -0.07 * (sodium - REF_SODIUM),
      patient_bnp_term - ref_bnp_term,
      -1.12 * (LVAD - 0),
       1.02 * (short_MCS_ever - 0)
    ),
    stringsAsFactors = FALSE
  )
}

# ---- Theme ------------------------------------------------------------------

app_theme <- bs_theme(
  version      = 5,
  bootswatch   = "flatly",
  primary      = "#2C6FB3",
  base_font    = font_google("Inter"),
  heading_font = font_google("Inter"),
  font_scale   = 0.95
)

app_css <- HTML("
  .navbar-brand { font-weight: 600; letter-spacing: 0.2px; }
  .bslib-value-box .value-box-title { font-size: 1rem; font-weight: 500; opacity: 0.95; }
  .bslib-value-box .value-box-value { font-size: 2.8rem; font-weight: 700; line-height: 1.1; }
  .bslib-value-box .value-box-showcase { opacity: 0.88; }
  .form-group, .shiny-input-container { margin-bottom: 1rem; }
  .irs--shiny .irs-bar { background: #2C6FB3; border-color: #2C6FB3; }
  .irs--shiny .irs-handle { border-color: #2C6FB3; }
  .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
    background: #2C6FB3;
  }
  .section-card { box-shadow: 0 1px 2px rgba(20,40,80,0.04); }
  .section-card > .card-header {
    background: #f5f8fc; font-weight: 600; color: #1f4f85;
    border-bottom: 1px solid #e3ecf6;
  }
  .card { border-radius: 10px; }
  .card-header { border-top-left-radius: 10px; border-top-right-radius: 10px; }
  .tier-chip {
    display: inline-block; padding: 2px 10px; border-radius: 999px;
    color: white; font-weight: 600; font-size: 0.8rem; margin-right: 6px;
  }
  .btn-primary { background-color: #2C6FB3; border-color: #2C6FB3; }
  .btn-primary:hover { background-color: #245a94; border-color: #245a94; }
  .navbar a.nav-link { font-weight: 500; }
  .accordion-button:not(.collapsed) {
    background-color: #eaf2fb; color: #1f4f85;
  }

  /* Mobile responsive tweaks */
  @media (max-width: 576px) {
    .navbar-brand { font-size: 1rem; }
    .bslib-value-box .value-box-value { font-size: 2.2rem; }
    .bslib-value-box .value-box-title { font-size: 0.9rem; }
    .card-header { padding: 0.6rem 0.8rem; font-size: 0.95rem; }
    .card-body { padding: 0.8rem; }
    .tier-chip { font-size: 0.75rem; padding: 2px 8px; }
    .sidebar-title { font-size: 1rem; }
  }

  /* Card header that wraps title + subtitle */
  .card-header-wrap {
    display: flex; flex-wrap: wrap; align-items: baseline;
    gap: 4px 10px;
  }
  .card-header-sub {
    font-weight: 400; font-size: 0.8rem; color: #888;
  }
")

# ---- UI ---------------------------------------------------------------------

ui <- page_navbar(
  title = tagList(
    bs_icon("heart-pulse-fill"),
    span("US-CRS Calculator", style = "margin-left: 6px;")
  ),
  theme = app_theme,
  bg = "#2C6FB3",
  inverse = TRUE,
  fillable = FALSE,
  header = tags$head(tags$style(app_css)),

  # ---- Calculator tab ----
  nav_panel(
    title = tagList(bs_icon("calculator"), " Calculator"),

    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        title = tagList(bs_icon("clipboard2-pulse"), " Patient Characteristics"),
        bg = "#ffffff",

        card(
          class = "section-card",
          card_header(tagList(bs_icon("person-fill"), " Demographics")),
          card_body(
            prettyRadioButtons(
              inputId = "sex", label = "Sex",
              choices = c("Female" = 1, "Male" = 0),
              selected = 1, inline = TRUE,
              status = "primary", icon = icon("check"),
              animation = "smooth"
            ),
            numericInput("age", "Age (years)", value = 50, min = 0, max = 120)
          )
        ),

        card(
          class = "section-card",
          card_header(tagList(bs_icon("droplet-half"), " Laboratory Values")),
          card_body(
            sliderInput("albumin", "Albumin (g/dL)",
                        min = 1, max = 5.5, step = 0.05, value = 4.5),
            sliderInput("bilirubin", "Bilirubin (mg/dL)",
                        min = 0.5, max = 10, step = 0.05, value = 0.5),
            sliderInput("creatinine", "Serum Creatinine (mg/dL)",
                        min = 0, max = 5, step = 0.05, value = 1.5),
            sliderInput("sodium", "Sodium (mEq/L)",
                        min = 110, max = 145, step = 1, value = 140)
          )
        ),

        card(
          class = "section-card",
          card_header(tagList(bs_icon("activity"), " BNP")),
          card_body(
            prettyRadioButtons(
              inputId = "BNP_type", label = "BNP Assay Type",
              choices = c("NT-pro BNP" = 1, "Regular BNP" = 0),
              selected = 0, inline = TRUE,
              status = "primary", icon = icon("check"),
              animation = "smooth"
            ),
            sliderInput("BNP_value", "Natural Log of BNP (pg/mL)",
                        min = -1, max = 10, step = 0.1, value = 1.5)
          )
        ),

        card(
          class = "section-card",
          card_header(tagList(bs_icon("heart-pulse"), " Mechanical Support")),
          card_body(
            prettyRadioButtons(
              inputId = "lvad", label = "Left ventricular assist device?",
              choices = c("Yes" = 1, "No" = 0),
              selected = 0, inline = TRUE,
              status = "primary", icon = icon("check"),
              animation = "smooth"
            ),
            prettyRadioButtons(
              inputId = "short_term_MCS", label = "Ever on short-term MCS?",
              choices = c("Yes" = 1, "No" = 0),
              selected = 0, inline = TRUE,
              status = "primary", icon = icon("check"),
              animation = "smooth"
            )
          )
        )
      ),

      # ---- Main content ----
      layout_column_wrap(
        width = 1,
        fill = FALSE,
        uiOutput("score_value_box")
      ),

      card(
        full_screen = FALSE,
        card_header(
          tagList(bs_icon("bullseye"), " Score Position on the US-CRS Scale")
        ),
        card_body(
          min_height = 210,
          plotOutput("gauge", height = "180px")
        )
      ),

      card(
        full_screen = FALSE,
        card_header(
          div(
            class = "card-header-wrap",
            span(bs_icon("bar-chart-line"), " Score Drivers"),
            span("vs. typical reference patient", class = "card-header-sub")
          )
        ),
        card_body(
          plotOutput("drivers", height = "320px"),
          p(
            tags$em("Bars show how each variable shifts the raw score relative to a ",
                    "reference patient with typical values (albumin 4.0, bilirubin 1.0, ",
                    "creatinine 1.0, sodium 138, no devices, moderate BNP)."),
            style = "font-size:0.8rem; color:#777; margin-top:8px; margin-bottom:0;"
          )
        )
      ),

      accordion(
        open = "about",
        accordion_panel(
          title = tagList(bs_icon("info-circle"), " About the US-CRS Score"),
          value = "about",
          p("The ", strong("US-CRS Score"),
            " is a continuous score that estimates the probability a heart ",
            "transplant candidate will die on the waitlist. Scores range from ",
            strong("0 - 50"), ", with higher values indicating greater urgency."),
          p(tags$em(
            "Zhang KC, Parker WF, et al. JAMA, February 2024."
          ))
        ),
        accordion_panel(
          title = tagList(bs_icon("list-check"), " Interpretation Guide"),
          value = "tiers",
          tags$div(
            tags$p(tags$span(class = "tier-chip",
                             style = "background:#2a9d8f;", "1 - 10"),
                   " Lower estimated waitlist mortality."),
            tags$p(tags$span(class = "tier-chip",
                             style = "background:#e9c46a; color:#333;", "11 - 25"),
                   " Moderate estimated waitlist mortality."),
            tags$p(tags$span(class = "tier-chip",
                             style = "background:#f4a261;", "26 - 40"),
                   " High estimated waitlist mortality."),
            tags$p(tags$span(class = "tier-chip",
                             style = "background:#e76f51;", "41 - 50"),
                   " Very high estimated waitlist mortality.")
          )
        )
      )
    )
  ),

  # ---- Bulk upload tab ----
  nav_panel(
    title = tagList(bs_icon("file-earmark-spreadsheet"), " Bulk Upload"),

    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        title = tagList(bs_icon("cloud-upload"), " Upload a Dataset"),
        bg = "#ffffff",

        card(
          class = "section-card",
          card_header(tagList(bs_icon("file-earmark-arrow-up"), " CSV File")),
          card_body(
            fileInput(
              "datafile", "Select a CSV file",
              accept = ".csv",
              buttonLabel = "Browse...",
              placeholder = "No file selected"
            ),
            downloadButton("download_scores", "Download scored CSV",
                           class = "btn-primary w-100")
          )
        ),

        card(
          class = "section-card",
          card_header(tagList(bs_icon("list-columns"), " Required Columns")),
          card_body(
            tags$ul(
              style = "padding-left: 1.1rem; margin-bottom: 0;",
              tags$li(tags$code("albumin")),
              tags$li(tags$code("bilirubin")),
              tags$li(tags$code("sex"), " (1 = Female, 0 = Male)"),
              tags$li(tags$code("age")),
              tags$li(tags$code("creatinine")),
              tags$li(tags$code("sodium")),
              tags$li(tags$code("LVAD"), " (1/0)"),
              tags$li(tags$code("short_MCS_ever"), " (1/0)"),
              tags$li(tags$code("BNP_NT_Pro"), " (1/0)"),
              tags$li(tags$code("BNP"), " (pg/mL)")
            )
          )
        )
      ),

      layout_column_wrap(
        width = 1/2,
        fill = FALSE,
        value_box(
          title = "Patients scored",
          value = textOutput("n_patients"),
          showcase = bs_icon("people-fill"),
          theme = "primary"
        ),
        value_box(
          title = "Mean US-CRS score",
          value = textOutput("mean_score"),
          showcase = bs_icon("graph-up"),
          theme = "info"
        )
      ),

      card(
        full_screen = TRUE,
        card_header(tagList(bs_icon("table"), " Scored Patients")),
        card_body(DTOutput("bulk_scores"))
      )
    )
  ),

  nav_spacer(),
  nav_item(
    tags$a(
      href = "https://jamanetwork.com/journals/jama/article-abstract/2814884",
      target = "_blank", rel = "noopener",
      style = "color: white;",
      bs_icon("journal-text"), " JAMA Paper"
    )
  )
)

# ---- Server -----------------------------------------------------------------

server <- function(input, output, session) {

  # Reactive score computation for the manual tab
  manual_score <- reactive({
    df <- data.frame(
      albumin         = input$albumin,
      bilirubin       = input$bilirubin,
      sex             = as.numeric(input$sex),
      age             = input$age,
      creatinine      = input$creatinine,
      sodium          = input$sodium,
      LVAD            = as.numeric(input$lvad),
      short_MCS_ever  = as.numeric(input$short_term_MCS),
      BNP_NT_Pro      = as.numeric(input$BNP_type),
      BNP             = exp(input$BNP_value)
    )
    compute_uscrs(df)$uscrs_score
  })

  # Dynamic value box so the theme color tracks the risk tier
  output$score_value_box <- renderUI({
    score <- tryCatch(manual_score(), error = function(e) NA_integer_)
    info  <- risk_info(score)

    display_value <- if (is.na(score)) "—" else as.character(score)

    value_box(
      title     = "US-CRS Score",
      value     = display_value,
      p(info$label, style = paste0("font-weight:600; color:", info$color, ";")),
      p(info$desc, style = "font-size:0.85rem; color:#555; margin-bottom:0;"),
      showcase  = bs_icon("heart-fill"),
      theme     = value_box_theme(bg = info$color, fg = "white")
    )
  })

  # Gradient gauge visualization
  output$gauge <- renderPlot({
    score <- tryCatch(manual_score(), error = function(e) NA_integer_)

    plot_w <- tryCatch(
      session$clientData[[paste0("output_gauge_width")]],
      error = function(e) NULL
    )
    is_narrow <- !is.null(plot_w) && plot_w < 520

    grad <- data.frame(x = seq(0, 49.8, by = 0.2))
    grad$xend <- grad$x + 0.2

    x_breaks   <- if (is_narrow) seq(0, 50, by = 10) else seq(0, 50, by = 5)
    xtext_size <- if (is_narrow) 10 else 12
    xtitle_sz  <- if (is_narrow) 11 else 13
    label_size <- if (is_narrow) 3.6 else 4.4

    # Clamp the "Score: N" label position so it stays inside the panel
    label_x <- if (!is.na(score)) max(4, min(46, score)) else NA_real_

    p <- ggplot() +
      geom_rect(
        data = grad,
        aes(xmin = x, xmax = xend, ymin = 0, ymax = 1, fill = x),
        color = NA
      ) +
      scale_fill_gradientn(
        colors = c("#2a9d8f", "#8ab17d", "#e9c46a", "#f4a261", "#e76f51"),
        limits = c(0, 50)
      ) +
      scale_x_continuous(breaks = x_breaks,
                         limits = c(-1, 51), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-0.3, 1.9), expand = c(0, 0)) +
      labs(x = "US-CRS Score", y = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = xtext_size, color = "#444"),
        axis.title.x = element_text(size = xtitle_sz, color = "#333",
                                    margin = margin(t = 10)),
        plot.margin = margin(10, 10, 5, 10)
      )

    if (!is.na(score)) {
      p <- p +
        annotate("segment", x = score, xend = score,
                 y = -0.05, yend = 1.35,
                 color = "#111", linewidth = 1.1) +
        annotate("point", x = score, y = 1.35,
                 shape = 25, size = 5,
                 fill = "#111", color = "#111") +
        annotate("label", x = label_x, y = 1.72,
                 label = paste0("Score: ", score),
                 size = label_size, fontface = "bold",
                 color = "#111", fill = "white",
                 label.size = 0.4, label.r = unit(0.25, "lines"))
    }
    p
  }, res = 96)

  # Variable contribution analysis
  manual_contributions <- reactive({
    compute_contributions(
      albumin        = input$albumin,
      bilirubin      = input$bilirubin,
      sex            = as.numeric(input$sex),
      age            = input$age,
      creatinine     = input$creatinine,
      sodium         = input$sodium,
      LVAD           = as.numeric(input$lvad),
      short_MCS_ever = as.numeric(input$short_term_MCS),
      BNP_NT_Pro     = as.numeric(input$BNP_type),
      BNP            = exp(input$BNP_value)
    )
  })

  output$drivers <- renderPlot({
    contribs <- tryCatch(manual_contributions(), error = function(e) NULL)
    if (is.null(contribs)) return(NULL)

    contribs <- contribs %>%
      filter(abs(contribution) > 0.001)

    if (nrow(contribs) == 0) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "All values near the reference \u2014 no strong drivers.",
                   size = 5, color = "#888") +
          theme_void()
      )
    }

    # Narrow-screen (mobile) adjustments
    plot_w <- tryCatch(
      session$clientData[[paste0("output_drivers_width")]],
      error = function(e) NULL
    )
    is_narrow <- !is.null(plot_w) && plot_w < 520

    # Convert to percentage of total absolute deviation
    total_abs <- sum(abs(contribs$contribution))
    contribs <- contribs %>%
      mutate(
        pct       = (contribution / total_abs) * 100,
        direction = ifelse(contribution > 0, "Risk-increasing", "Protective"),
        variable  = reorder(variable, abs(pct)),
        label     = paste0(round(abs(pct)), "%"),
        inside    = abs(pct) >= 12,
        label_y   = ifelse(inside, pct / 2, pct),
        label_hjust = case_when(
          inside       ~ 0.5,
          pct >= 0     ~ -0.2,
          TRUE         ~ 1.2
        ),
        label_color = ifelse(inside, "white", "#333")
      )

    label_size <- if (is_narrow) 3.2 else 3.8
    ytitle_size <- if (is_narrow) 10 else 11
    ytext_size  <- if (is_narrow) 9  else 10
    xtext_size  <- if (is_narrow) 10 else 11

    ggplot(contribs, aes(x = variable, y = pct, fill = direction)) +
      geom_col(width = 0.65) +
      geom_hline(yintercept = 0, linewidth = 0.5, color = "#333") +
      geom_text(
        aes(y = label_y, label = label,
            hjust = label_hjust, color = label_color),
        size = label_size, fontface = "bold", show.legend = FALSE
      ) +
      coord_flip(clip = "off") +
      scale_fill_manual(
        values = c("Risk-increasing" = "#e76f51", "Protective" = "#2a9d8f"),
        labels = c("Risk-increasing" = "Risk \u2191",
                   "Protective"      = "Protective"),
        name = NULL
      ) +
      scale_color_identity() +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      labs(x = NULL, y = "% of score deviation") +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        legend.text = element_text(size = if (is_narrow) 10 else 11),
        legend.margin = margin(b = 2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = xtext_size, color = "#333",
                                    face = "bold", lineheight = 0.85),
        axis.text.x = element_text(size = ytext_size, color = "#555"),
        axis.title.x = element_text(size = ytitle_size, color = "#555",
                                    margin = margin(t = 8)),
        plot.margin = margin(5, 15, 5, 5)
      )
  }, res = 96)

  # ---- Bulk scoring reactive ----
  bulk_data <- reactive({
    req(input$datafile)
    df <- read.csv(input$datafile$datapath)

    required_cols <- c("albumin", "bilirubin", "sex", "age", "creatinine",
                       "sodium", "LVAD", "short_MCS_ever", "BNP_NT_Pro", "BNP")
    if (!all(required_cols %in% names(df))) {
      missing <- setdiff(required_cols, names(df))
      return(list(
        error = paste0("Missing required columns: ",
                       paste(missing, collapse = ", "))
      ))
    }

    scored <- compute_uscrs(df) %>%
      select(uscrs_score, everything())
    list(data = scored)
  })

  output$n_patients <- renderText({
    res <- tryCatch(bulk_data(), error = function(e) NULL)
    if (is.null(res) || !is.null(res$error)) "—" else nrow(res$data)
  })

  output$mean_score <- renderText({
    res <- tryCatch(bulk_data(), error = function(e) NULL)
    if (is.null(res) || !is.null(res$error)) {
      "—"
    } else {
      sprintf("%.1f", mean(res$data$uscrs_score, na.rm = TRUE))
    }
  })

  output$bulk_scores <- renderDT({
    res <- bulk_data()
    if (!is.null(res$error)) {
      return(datatable(
        data.frame(Error = res$error),
        options = list(dom = "t"), rownames = FALSE
      ))
    }

    dt <- datatable(
      res$data,
      rownames = FALSE,
      class = "stripe hover compact",
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = "ftip",
        columnDefs = list(list(className = "dt-center", targets = 0))
      )
    )

    # Color-code the score column by risk tier
    dt %>%
      formatStyle(
        "uscrs_score",
        backgroundColor = styleInterval(
          c(10, 25, 40),
          c("#c8e6c9", "#fff2c2", "#ffd8a8", "#ffb4a2")
        ),
        fontWeight = "bold"
      )
  })

  output$download_scores <- downloadHandler(
    filename = function() {
      paste0("uscrs_scored_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      res <- bulk_data()
      if (is.null(res$error)) {
        write.csv(res$data, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Error = res$error), file, row.names = FALSE)
      }
    }
  )
}

shinyApp(ui = ui, server = server)
