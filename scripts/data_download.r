library(TCGAbiolinks)
library(here)

data_dir <- here("data", "TCGA_data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# Get a list of all TCGA projects
tcga_projects <- TCGAbiolinks::getGDCprojects()$project_id

# Filter for only TCGA projects (projects starting with 'TCGA-')
tcga_projects <- tcga_projects[grep("TCGA-", tcga_projects)]

# Loop through each TCGA project
for (project in tcga_projects) {
  
  # Print the project being processed
  message("Processing: ", project)

  # Create a query for gene expression data (STAR - Counts) from the project
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  # Define the path to save the processed data
  save_path <- file.path(data_dir, paste0(project, ".rds"))

  # Check if the data has already been downloaded (avoid duplicates)
  if (!file.exists(save_path)) {
    
    # If not, download the data and save it to the specified directory
    GDCdownload(query, directory = data_dir)
    
    # Prepare the downloaded data for analysis (e.g., formatting into a data frame)
    data <- GDCprepare(query, directory = data_dir)
    
    # Save the prepared data as an RDS file
    saveRDS(data, save_path)
  } else {
    
    # If the data already exists, skip downloading and processing
    message("Skipping: ", project, " (already downloaded)")
  }
}
