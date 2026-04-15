library(shiny)
library(shinythemes)
library(shinyWidgets)
library(tidyverse)
library(here)
options(shiny.maxRequestSize = 150*1024^2) 

mapping_table <- read.csv(here::here("uscrs_mapping_table.csv"))

# UI ----
ui <- fluidPage(
  theme = shinytheme('cerulean'),
  titlePanel('The US-CRS Score for Heart Transplant Allocation'),
  
  tabsetPanel(
    # Manual input tab
    tabPanel("Manual Entry",
             sidebarLayout(
               sidebarPanel(
                 h2('Patient Characteristics'),
                 radioButtons('sex', 'Gender', choices = c('Female' = 1, 'Male' = 0), selected = 1),
                 numericInput('age', 'Age (years):', value = 50),
                 sliderInput('albumin', 'Albumin (g/dL):', min = 1, max = 5.5, step = 0.05, value = 4.5),
                 sliderInput('bilirubin', 'Bilirubin (mg/dL):', min = 0.5, max = 10, step = 0.05, value = 0.5),
                 sliderInput('creatinine', 'Serum Creatinine (mg/dL):', min = 0, max = 5, step = 0.05, value = 1.5),
                 sliderInput('sodium', 'Sodium (mEq/L):', min = 110, max = 145, step = 1, value = 140),
                 radioButtons('lvad', 'Left ventricular assist device?', choices = c('Yes' = 1, 'No' = 0), selected = 0),
                 radioButtons('short_term_MCS', 'Ever on short-term MCS?', choices = c('Yes' = 1, 'No' = 0), selected = 0),
                 radioButtons('BNP_type', 'Type of BNP', choices = c('NT-pro BNP' = 1, 'Regular BNP' = 0), selected = 0),
                 sliderInput('BNP_value', 'Natural Log of BNP (pg/mL)', min = -1, max = 10, step = 0.1, value = 1.5)
               ),
               mainPanel(
                 h2(textOutput('result'), style = 'text-align:center; font-size:25px;'),
                 plotOutput(outputId = 'histogram'),
                 br(),
                 h3('About the US-CRS Score'),
                 h5('The US-CRS Score is a continuous score that estimates the probability a heart transplant candidate will die on the waitlist. Scores range from 0–50; higher = greater urgency.'),
                 h5('JAMA, Feb. 2024. Zhang KC, Parker WF, et al.')
               )
             )
    ),
    
    # Upload tab
    tabPanel("Upload Dataset",
             sidebarLayout(
               sidebarPanel(
                 fileInput("datafile", "Upload CSV File", accept = ".csv"),
                 helpText("Required columns: albumin, bilirubin, sex, age, creatinine, sodium, LVAD, short_MCS_ever, BNP_NT_Pro, BNP")
               ),
               mainPanel(
                 tableOutput("bulk_scores")
               )
             )
    )
  )
)

# SERVER ----
server <- function(input, output) {
  
  # MANUAL score calculator ----
  observe({
    df <- data.frame(
      albumin = input$albumin,
      bilirubin = input$bilirubin,
      sex = input$sex,
      age = input$age,
      creatinine = input$creatinine,
      sodium = input$sodium,
      LVAD = as.numeric(input$lvad),
      short_MCS_ever = as.numeric(input$short_term_MCS),
      BNP_NT_Pro = as.numeric(input$BNP_type),
      BNP = exp(input$BNP_value)
    )
    
    df <- df %>%
      mutate(
        eGFR = case_when(
          sex == '1' ~ 142 * (pmin((creatinine / 0.7), 1)^(-0.241)) *
            (pmax((creatinine / 0.7), 1)^(-1.2)) * 0.9938^(age) * 1.012,
          sex == '0' ~ 142 * (pmin((creatinine / 0.9), 1)^(-0.302)) *
            (pmax((creatinine / 0.7), 1)^(-1.2)) * 0.9938^(age)
        ),
        bnp_reg_indicator = ifelse(BNP_NT_Pro != 1, 1, 0),
        uscrs_raw_score = 1.02*short_MCS_ever + 0.55*log(bilirubin + 1) - 0.01*eGFR +
          0.40*log(BNP)*bnp_reg_indicator + 0.20*log(BNP)*BNP_NT_Pro - 
          0.63*albumin - 0.07*sodium -1.12*LVAD
      )
    
    breaks <- sort(unique(c(mapping_table$min_raw_score, max(mapping_table$max_raw_score))))
    df <- df %>%
      mutate(uscrs_score = cut(uscrs_raw_score, breaks = breaks, labels = FALSE, include.lowest = TRUE)) 
    
    score <- df$uscrs_score
    
    output$result <- renderText({ paste("US-CRS Score: ", score) })
    
    output$histogram <- renderPlot({
      df_plot <- data.frame(
        bin = seq(2.5, 47.5, by = 5),  # midpoints
        score_range = paste(seq(0, 45, by = 5), seq(5, 50, by = 5), sep = "-")
      )
      
      # Determine which bin the score falls into
      highlight_bin <- cut(score, breaks = seq(0, 50, by = 5), labels = FALSE, include.lowest = TRUE)
      
      # Mark color = 1 for bin containing score, 0 otherwise
      df_plot$color <- ifelse(seq_len(nrow(df_plot)) == highlight_bin, 1, 0)
      
      # Plot
      ggplot(df_plot, aes(x = bin, y = 1, fill = factor(color))) +
        geom_col(width = 5, color = "black") +
        scale_fill_manual(values = c('white', '#3f88c5')) +
        scale_x_continuous(breaks = seq(0, 50, by = 5), limits = c(0, 50)) +
        labs(x = 'US-CRS Score', y = '') +
        theme_bw() +
        theme(
          legend.position = 'none',
          panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16, margin = margin(t = 15))
        )
    })
  })
  
  # BULK scoring from file ----
  output$bulk_scores <- renderTable({
    req(input$datafile)
    infile <- input$datafile$datapath
    
    df <- read.csv(infile)
    
    if (!all(c("albumin", "bilirubin", "sex", "age", "creatinine", "sodium", 
               "LVAD", "short_MCS_ever", "BNP_NT_Pro", "BNP") %in% names(df))) {
      return(data.frame(Error = "Missing required columns. Please check your input."))
    }
    
    df <- df %>%
      mutate(
        eGFR = case_when(
          sex == 1 ~ 142 * (pmin((creatinine / 0.7), 1)^(-0.241)) *
            (pmax((creatinine / 0.7), 1)^(-1.2)) * 0.9938^(age) * 1.012,
          sex == 0 ~ 142 * (pmin((creatinine / 0.9), 1)^(-0.302)) *
            (pmax((creatinine / 0.7), 1)^(-1.2)) * 0.9938^(age)
        ),
        bnp_reg_indicator = ifelse(BNP_NT_Pro != 1, 1, 0),
        uscrs_raw_score = 1.02*short_MCS_ever + 0.55*log(bilirubin + 1) - 0.01*eGFR +
          0.40*log(BNP)*bnp_reg_indicator + 0.20*log(BNP)*BNP_NT_Pro - 
          0.63*albumin - 0.07*sodium -1.12*LVAD
      )
    
    breaks <- sort(unique(c(mapping_table$min_raw_score, max(mapping_table$max_raw_score))))
    df <- df %>%
      mutate(uscrs_score = cut(uscrs_raw_score, breaks = breaks, labels = FALSE, include.lowest = TRUE))
    
    df %>% select(uscrs_score, everything())
  })
}

shinyApp(ui = ui, server = server)
