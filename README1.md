# README

1. Project Description
I developed this project as a tool for analyzing multi-omics data using:
- Principal Component Analysis (PCA):To reduce data dimensionality and explore variance within the datasets.
- Canonical Correlation Analysis (CCA): To identify and analyze relationships between transcriptomics and proteomics datasets.

The tool is implemented in R and features an interactive user interface built with the (Shiny) package.

 2. System Requirements
To run this project, make sure please  you have the following installed:
- R version 4.0.0 or newer.
- The following R packages:
  - `DESeq2`
  - `ggplot2`
  - `CCA`
  - `shiny`
  - `shinyWidgets`

 3. Installing Required Packages
Before running the code, you’ll need to install the required R packages. 
Use the following commands in your R environment please :

```R
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")

if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
  install.packages("shinyWidgets")
}

install.packages(c("ggplot2", "CCA", "shiny"))

4. How to Run the Code

Option 1: 
Using RStudio
	1.	Open the app.R file in RStudio.
	2.	Make sure the dataset file GSE199826_raw-counts-matrix.txt is in the same directory as the script.
	3.	Click the “Run App” button in RStudio to launch the tool.

Option 2:
 Using Terminal
	1.	Open Terminal or Command Prompt.
	2.	Navigate to the project directory containing the app.R file:
cd /path/to/your/project
	3.	Run the code using the following command:
Rscript app.R

5. Project Structure

Here’s the file structure for the project:
	app.R: The main code file for running the tool.
	data/GSE199826_raw-counts-matrix.txt: The raw data file required for analysis.
	README.md: Instructions for running the project.

6. User Configuration

Analysis Options:
	Upon launching the app, you can select the type of analysis:
	PCA: Perform Principal Component Analysis.
	CCA: Perform Canonical Correlation Analysis.

Customization:
	The app allows you to specify:
	The number of PCA components to analyze.
	Plot colors using the Shiny interface.
