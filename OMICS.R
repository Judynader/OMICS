if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

# Install shinyWidgets if not already installed
if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
  install.packages("shinyWidgets")
}

# Load necessary libraries for analysis
library(DESeq2)
library(ggplot2)
library(CCA)
library(shiny)
library(shinyWidgets)

# Load raw count data
raw_counts_data <- read.delim("/Users/judynader/Downloads/GSE199826_raw-counts-matrix.txt", header = TRUE, sep = "\t")

# Preprocess the data
rownames(raw_counts_data) <- raw_counts_data$Genes
raw_counts_data <- raw_counts_data[, -1]  # Remove the genes column
raw_counts_data[is.na(raw_counts_data)] <- 0  # Replace NA values with 0
raw_counts_data[raw_counts_data < 0] <- 0  # Ensure no negative values
raw_counts_data <- as.data.frame(lapply(raw_counts_data, as.integer))  # Convert values to integers


# Define experimental conditions
condition <- factor(rep(c("control", "treated"), each = 13))
condition <- c(condition, "treated")  # Adjust to match sample size
colData <- data.frame(condition = condition, row.names = colnames(raw_counts_data))

# Check if `dds` exists
if (exists("dds")) {
  # Extract normalized counts and apply log2 transformation
  log_normalized_data <- log2(counts(dds, normalized = TRUE) + 1)
  print(dim(log_normalized_data))  # Check dimensions of the normalized data
} else {
  stop("The `dds` object is missing. Ensure DESeq2 analysis was performed.")
}

# Separate transcriptomics and proteomics data based on column indices
transcriptomics_data <- log_normalized_data[, 1:10]  # Assuming first 10 columns are Transcriptomics
proteomics_data <- log_normalized_data[, 11:20]     # Assuming next 10 columns are Proteomics

# Check dimensions
print(dim(transcriptomics_data))
print(dim(proteomics_data))


print(dim(log_normalized_data))
print(colnames(log_normalized_data))


pca_transcriptomics <- prcomp(transcriptomics_data, scale. = TRUE)
variance_explained <- (pca_transcriptomics$sdev^2) / sum(pca_transcriptomics$sdev^2) * 100

# Create and plot Scree Plot
scree_plot_data <- data.frame(Component = seq_along(variance_explained), Variance = variance_explained)
ggplot(scree_plot_data, aes(x = Component, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(group = 1, color = "blue") +
  geom_point(size = 2, color = "red") +
  labs(title = "Scree Plot for Transcriptomics PCA", x = "Principal Component", y = "Variance Explained (%)") +
  theme_minimal()


pca_proteomics <- prcomp(proteomics_data, scale. = TRUE)
variance_explained_proteomics <- (pca_proteomics$sdev^2) / sum(pca_proteomics$sdev^2) * 100

# Create and plot Scree Plot
scree_plot_proteomics <- data.frame(Component = seq_along(variance_explained_proteomics), Variance = variance_explained_proteomics)
ggplot(scree_plot_proteomics, aes(x = Component, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(group = 1, color = "blue") +
  geom_point(size = 2, color = "red") +
  labs(title = "Scree Plot for Proteomics PCA", x = "Principal Component", y = "Variance Explained (%)") +
  theme_minimal()


# Perform CCA
cca_result <- cc(transcriptomics_data, proteomics_data)

# Create and plot Correlation Results
correlation_data <- data.frame(Canonical_Variate = seq_along(cca_result$cor), Correlation = cca_result$cor)
ggplot(correlation_data, aes(x = Canonical_Variate, y = Correlation)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  labs(title = "Canonical Correlation Analysis", x = "Canonical Variate", y = "Correlation") +
  theme_minimal()


ui <- fluidPage(
  titlePanel("Multi-Omics Data Integration Tool"),
  sidebarLayout(
    sidebarPanel(
      selectInput("analysis_type", "Choose Analysis Type:", choices = c("PCA", "CCA")),
      numericInput("num_components", "Number of PCA Components:", min = 1, max = 50, value = 2),
      sliderTextInput("color_choice", "Choose Plot Color:", choices = c("blue", "red", "green", "black", "yellow"), selected = "blue"),
      actionButton("run_analysis", "Execute Analysis")
    ),
    mainPanel(
      plotOutput("analysis_plot"),
      verbatimTextOutput("cca_summary")
    )
  )
)


server <- function(input, output) {
  observeEvent(input$run_analysis, {
    analysis_type <- input$analysis_type
    if (analysis_type == "PCA") {
      pca_result <- prcomp(log_normalized_data, scale. = TRUE)
      output$analysis_plot <- renderPlot({
        ggplot(data.frame(Component = seq_along(pca_result$sdev^2), Variance = pca_result$sdev^2), aes(x = Component, y = Variance)) +
          geom_line(color = input$color_choice) +
          geom_point(color = input$color_choice) +
          ggtitle("Explained Variance by PCA Components")
      })
    } else if (analysis_type == "CCA") {
      cca_result <- cc(transcriptomics_data, proteomics_data)
      output$cca_summary <- renderPrint({ summary(cca_result) })
    }
  })
}


shinyApp(ui = ui, server = server)
