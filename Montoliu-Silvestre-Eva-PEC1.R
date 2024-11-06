# PEC 1 - ANÁLISIS DE DATOS ÓMICOS

# Instalación de BioConductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SummarizedExperiment")

# Carga de paquetes necesarios
library(SummarizedExperiment)
library(readr)  # Para leer archivos csv
library(tidyverse)  # Esto incluye dplyr y otras librerías útiles

# Leer los datos
abundance_data <- read_tsv("data/m_MTBLS10856_LC-MS_positive_hilic_metabolite_profiling_v2_maf.tsv", , show_col_types = FALSE)
sample_info <- read.table("data/s_MTBLS10856.txt", header = TRUE, sep = "\t")
additional_sample_info <- read.table("data/a_MTBLS10856_LC-MS_positive_hilic_metabolite_profiling.txt", header = TRUE, sep = "\t")
metadata <- readLines("data/i_Investigation.txt")

# Dimensiones y estructura de cada elemento 
dim(abundance_data)
str(abundance_data)

dim(sample_info)
str(sample_info)

dim(additional_sample_info)
str(additional_sample_info)

# PREPARAMOS LOS DIFERENTES ELEMENTOS DEL OBJETO 'SummarizedExperiment'

## Assay
# Guardar los nombres de las muestras
sample_names <- sample_info$Sample.Name
# Crea la matriz con todas las filas y solo las columnas correspondientes a las muestras
abundance_matrix <- as.matrix(abundance_data[, sample_names])
# Nombra las filas con el nombre del metabolito correspondiente
rownames(abundance_matrix) <- abundance_data$metabolite_identification

## rowData
# Eliminamos las columnas que no tienen información (100% de NA)
cols_borrar <- which(colMeans(is.na(abundance_data)) == 1)
rowData <- data.frame(abundance_data[, -cols_borrar])
# Eliminamos las columnas que tienen la información de abundancias (guardada en Assay)
rowData <- rowData[, 1:(ncol(rowData)-length(sample_names))]

# También podemos seleccionar las columnas que vamos a incluir manualmente
#rowData <- data.frame(mass_to_charge = abundance_data$mass_to_charge,
#                      retention_time = abundance_data$retention_time,
#                      data_base = abundance_data$database)
#rownames(rowData) <- abundance_data$metabolite_identification

## colData
# Eliminamos las columnas que no tienen información (100% de NA)
cols_borrar <- which(colMeans(is.na(sample_info)) == 1)
sample_info <- sample_info[, -cols_borrar]

cols_borrar <- which(colMeans(is.na(additional_sample_info)) == 1)
additional_sample_info <- additional_sample_info[, -cols_borrar]

## Metadata
  # no hay que hacer nada

# CREA EL OBJETO SummarizedExperiment
se <- SummarizedExperiment(assays = list(counts = abundance_matrix),
                           rowData = rowData,
                           colData = list(sample_info = sample_info, additional_sample_info = additional_sample_info),
                           metadata = list(metadata = metadata)
                           )


# EXPLORACIÓN DE LOS DATOS
# Revisamos el objeto SummarizedExperiment
se

# Y cada uno de sus elementos
dim(assay(se))
head(assay(se))
tail(assay(se))

head(rowData(se))
tail(rowData(se))

head(colData(se))
tail(colData(se))


# REPOSITORIO GitHub
# Guardar el objeto 'SummarizedExperiment' en formato .Rda
save(se, file = "data/se_dataset.Rda")

# Guardar metadata en formato markdown
writeLines(metadata, "data/metadata.md")


