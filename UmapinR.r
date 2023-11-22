# Install necessary packages if not installed
if (!req("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!req("umap", quietly = TRUE)) {
  install.packages("umap")
}

if (!req("DESeq2", quietly = TRUE)) {
  if (!req("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("DESeq2")
}

# Load libraries
library(tidyverse)
library(umap)
library(DESeq2)

# Read count matrix
counts <- read.csv('count_matrix_folder_path', row.names = 1, header = TRUE)
df <- t(counts)

# Extract the first row as column names
colnames(df) <- df[1, ]

# Drop the first row, which now contains the column names
df <- df[-1, ]

# Read sample information
metadata <- read.csv('sample_information_folder_path')

# Set the first column 'Name' as the index
metadata <- metadata %>% 
  column_to_rownames(var = "sample_ID")

# Merge data frames on the row names (sample IDs)
merged_df <- merge(df, metadata, by = 'row.names', all = TRUE)
rownames(merged_df) <- merged_df$Row.names
merged_df <- merged_df[, -1]

# Extract count matrix for DESeq2
deseq_data <- merged_df[, -c(1:6)]  # Adjust columns based on your data

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = deseq_data,
                              colData = metadata,
                              design = ~ Genotype + Sex)

# Perform DESeq analysis
dds <- DESeq(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Fit the UMAP model to the normalized count matrix
umap_model <- umap(normalized_counts, n_components = 2)
umap_embedding <- umap_transform(umap_model, normalized_counts)

# Define custom marker styles for each type
type_markers <- c('Eos' = 17, 'Mac' = 16)  # 17 for triangle, 16 for circle

# Create a data frame with the necessary information
umap_data_plot <- data.frame(
  x = umap_embedding[, 1],
  y = umap_embedding[, 2],
  Type = metadata$Type,
  Genotype = metadata$Genotype,
  Sex = metadata$Sex
)

# Create a ggplot
umap_plot <- ggplot(umap_data_plot, aes(x = x, y = y, color = Genotype, shape = Type, fill = Sex)) +
  geom_point(size = 3) +
  scale_shape_manual(values = type_markers) +
  labs(title = "UMAP Plot by Genotype and Sex",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2",
       color = "Genotype",
       shape = "Type",
       fill = "Sex") +
  theme_minimal()

# Show the plot
print(umap_plot)
