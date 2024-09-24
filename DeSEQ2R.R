#Start Deseq2 calculation and visulalization
# 필요한 라이브러리 로드
library(shiny)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(dplyr)

# 샘플 데이터를 위한 데이터 생성 함수
generate_sample_data <- function() {
  set.seed(123)
  counts <- matrix(rnbinom(100 * 6, mu = 10, size = 1), ncol = 6)
  col_data <- data.frame(
    condition = factor(rep(c("Control", "Treatment"), each = 3))
  )
  rownames(counts) <- paste0("Gene_", 1:100)
  colnames(counts) <- paste0("Sample_", 1:6)
  list(counts = counts, col_data = col_data)
}

# Shiny UI 정의
ui <- fluidPage(
  titlePanel("DESeq2 Analysis with Volcano Plot"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("count_file", "Upload Count Data (CSV)", accept = ".csv"),
      fileInput("coldata_file", "Upload Column Data (CSV)", accept = ".csv"),
      selectInput("condition_column", "Condition Column", choices = NULL),
      actionButton("run_deseq", "Run DESeq2"),
      hr(),
      sliderInput("padj_threshold", "padj Threshold", min = 0, max = 1, value = 0.05, step = 0.01),
      sliderInput("fc_threshold", "Fold Change Threshold", min = 0, max = 5, value = 1, step = 0.1),
      hr(),
      downloadButton("download_results", "Download DESeq2 Results")
    ),
    
    mainPanel(
      plotOutput("volcano_plot"),
      tableOutput("deseq_results")
    )
  )
)

# Shiny 서버 정의
server <- function(input, output, session) {
  
  # DESeq2 데이터셋을 저장할 reactive 변수
  deseq_data <- reactiveValues(dds = NULL, results = NULL)
  
  # 샘플 데이터를 로드하거나 업로드된 데이터를 읽음
  observeEvent(input$count_file, {
    counts <- read.csv(input$count_file$datapath, row.names = 1)
    coldata <- read.csv(input$coldata_file$datapath, row.names = 1)
    
    # Condition 선택을 위한 UI 업데이트
    updateSelectInput(session, "condition_column", choices = colnames(coldata))
    
    # DESeq2 객체 생성
    deseq_data$dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = coldata,
      design = as.formula(paste("~", input$condition_column))
    )
  })
  
  # DESeq2 분석 실행
  observeEvent(input$run_deseq, {
    req(deseq_data$dds)
    dds <- DESeq(deseq_data$dds)
    res <- results(dds)
    
    # DESeq2 결과를 정렬하여 저장
    deseq_data$results <- as.data.frame(res[order(res$padj),])
    
    # DESeq2 결과 출력
    output$deseq_results <- renderTable({
      deseq_data$results %>%
        mutate(Significant = ifelse(padj < input$padj_threshold, "Yes", "No"))
    }, rownames = TRUE)
  })
  
  # Volcano Plot 생성
  output$volcano_plot <- renderPlot({
    req(deseq_data$results)
    EnhancedVolcano(deseq_data$results,
                    lab = rownames(deseq_data$results),
                    x = "log2FoldChange",
                    y = "padj",
                    pCutoff = input$padj_threshold,
                    FCcutoff = input$fc_threshold,
                    pointSize = 2.0,
                    labSize = 4.0,
                    title = "Volcano Plot",
                    subtitle = "DESeq2 Analysis",
                    caption = "Threshold: padj < 0.05",
                    xlim = c(-5, 5),
                    ylim = c(0, -log10(min(deseq_data$results$padj, na.rm = TRUE))),
                    col = c("grey", "blue", "red", "green"),
                    colAlpha = 1,
                    legendPosition = "right")
  })
  
  # DESeq2 결과 다운로드 기능
  output$download_results <- downloadHandler(
    filename = function() {
      paste("deseq2_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(deseq_data$results, file)
    }
  )
}

# Shiny 앱 실행
shinyApp(ui = ui, server = server)

