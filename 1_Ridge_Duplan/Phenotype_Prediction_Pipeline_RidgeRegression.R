# Updated by: # DUPLAN Alexandre - 2025
# Contact: alexandre.duplan@inrae.fr
# INRAE, UMR BioForA, UMR 0588, F-45075 Orleans, France
#! R version 4.4.2

# Cleaning working space
rm(list = ls())
while (!is.null(dev.list())) {
  dev.off()
}

# Import packages
library(anyLib)
pkg <- c("caret","doSNOW","glmnet", "dplyr", "argparse")
suppressMessages(anyLib(pkg))
print("Library imported...")

# Define argument parser
parser <- ArgumentParser(description = "Script to perform ridge regression and phenotype prediction")

# Define arguments
parser$add_argument(
  "--indiv_names", 
  help = "Path to a file containing the individuals to exclude, one per line (e.g., TO_EXCLUDE.txt)", 
  type = "character", 
  default="/home/INRA/aduplan/Documents/5_Predictive_model/1_RidgeRegression/TO_EXCLUDE.txt"
)

parser$add_argument(
  "--working_dir", 
  help = "Directory where you work",
  default = "/home/INRA/aduplan/Documents/5_Predictive_model/"
)

parser$add_argument(
  "--output_dir", 
  help = "Directory where output will be saved",
  default = "Output_RidgeR_Prediction/"
)

parser$add_argument(
  "--nbr_core", 
  help = "Number of cores to use for parallel processing", 
  type = "integer", 
  default = 20
)
parser$add_argument(
  "--path_omics", 
  help = "Number of cores to use for parallel processing", 
  type = "character", 
  default = "/home/INRA/aduplan/Documents/5_Predictive_model/0_Data/PoplarOmicsData.Rdata" # "/home/INRA/aduplan/Documents/5_Predictive_model/0_Data/PoplarOmicsData_reduce.Rdata"
)

# Parse the arguments
args <- parser$parse_args()

# Get indiv_names from the file if provided
if (!is.null(args$indiv_names)) {
  # Read individuals to exclude from the file
  indiv_names <- readLines(args$indiv_names)
  indiv_names <- as.vector(indiv_names)
} else{
  indiv_names <- NA
}

# Define Variables
Output_dir <- args$output_dir
Working_dir <- args$working_dir
Nbr_core <- args$nbr_core
path_data_omics <- args$path_omics

#Create folder if necessary
if (!dir.exists(Output_dir)) {
  dir.create(Output_dir, recursive = TRUE)
}

# Get predictions functions
source(paste0(Working_dir, "Phenotype_Prediction_Functions.R"))
print("Functions sourced...")

# Load Omics Data ------------------------------------------------------------
PoplarData <- get(load(path_data_omics))
print("Data loaded...")

# Delete individuals that you don't want -------------------------------------
for (layer_name in names(PoplarData)) {
  # Vérifier si la couche est un data frame ou une matrice avec des rownames
  if (is.matrix(PoplarData[[layer_name]]) || is.data.frame(PoplarData[[layer_name]])) {   
    # Retirer les lignes avec les rownames spécifiés
    PoplarData[[layer_name]] <- PoplarData[[layer_name]][!rownames(PoplarData[[layer_name]]) %in% indiv_names, ]
  }
  else{
    print(paste0("Error: Omics data have a problem in the format: "), layer_name)
  }
}
print("Individuals excluded...(if necessary)")

Ridge.mod <- list() #Init result list

#Just to be sure about the class of Omics
PoplarData$Phenotype <- as.data.frame(PoplarData$Phenotype)

# Liste des colonnes à supprimer
colonnes_a_supprimer <- c("HT2009", "CIRC2009", "HT2011", "CIRC2011", "H", "G",
                          "S", "H_G", "S_G", "Klason_lignin", "Glucose", "Xyl_Glu",
			  "C5_C6", "Wet_chem_extractives", "Dia2015_sqrt",
		          "InfraDens", "Rust", "BrAnglVert",
			  "RamifSyllep", "BudFlushSlope")

colonnes_a_supprimer <- c("HT2009", "CIRC2009", "HT2011", "CIRC2011", "H", "G",
                          "S", "H_G", "S_G", "Klason_lignin", "Glucose", "Xyl_Glu",
                          "C5_C6", "Wet_chem_extractives", "Dia2015_sqrt",
                          "InfraDens")


# Supprimer les colonnes spécifiées du data frame PoplarDataORL$Phenotype
PoplarData$Phenotype <- PoplarData$Phenotype[, !(colnames(PoplarData$Phenotype) %in% colonnes_a_supprimer)]
print(colnames(PoplarData$Phenotype))
#Delete NA
PoplarData$Phenotype <- na.omit(PoplarData$Phenotype)

#Epigenetic is the new sources of information here but also the limit factor
#So let just get the individuals that we have DNA meth.
epigenetic_rownames <- rownames(PoplarData$Epigenetic_CG)
PoplarData$Phenotype <- as.data.frame(PoplarData$Phenotype[rownames(PoplarData$Phenotype) %in% epigenetic_rownames, , drop = FALSE])

# Delete omics information for individual when you don't have the phenotype associated
phenotype_rownames <- rownames(PoplarData$Phenotype)

PoplarData$Genomic <- PoplarData$Genomic[rownames(PoplarData$Genomic) %in% phenotype_rownames, , drop = FALSE]
PoplarData$Transcriptomic <- PoplarData$Transcriptomic[rownames(PoplarData$Transcriptomic) %in% phenotype_rownames, , drop = FALSE]
PoplarData$Epigenetic_CG <- PoplarData$Epigenetic_CG[rownames(PoplarData$Epigenetic_CG) %in% phenotype_rownames, , drop = FALSE]
PoplarData$Epigenetic_CHG <- PoplarData$Epigenetic_CHG[rownames(PoplarData$Epigenetic_CHG) %in% phenotype_rownames, , drop = FALSE]
PoplarData$Epigenetic_CHH <- PoplarData$Epigenetic_CHH[rownames(PoplarData$Epigenetic_CHH) %in% phenotype_rownames, , drop = FALSE]

#Combine all epigenetic in one dataframe
Epigenetic_Data <- cbind(PoplarData$Epigenetic_CG, PoplarData$Epigenetic_CHG, PoplarData$Epigenetic_CHH)

#Just to be sure about the class of Omics
PoplarData$Phenotype <- as.data.frame(PoplarData$Phenotype)
PoplarData$Genomic <- as.matrix(PoplarData$Genomic)
PoplarData$Transcriptomic <- as.matrix(PoplarData$Transcriptomic)
PoplarData$Epigenetic_CG <- as.matrix(PoplarData$Epigenetic_CG)
PoplarData$Epigenetic_CHG <- as.matrix(PoplarData$Epigenetic_CHG)
PoplarData$Epigenetic_CHH <- as.matrix(PoplarData$Epigenetic_CHH)
PoplarData$Epigenetic <- as.matrix(Epigenetic_Data)

# Nested Cross Validation Fold Sampling ---------------------------------------
PoplarData$CV.fold <- apply(PoplarData$Phenotype,2, FUN = Ncv_Fold_Sampling)

# check if data are numeric, have no NA, and share the same rownames
check_data_integrity <- function(PoplarData_temp) {
  PoplarData <- PoplarData_temp
  PoplarData$CV.fold <- NULL
  # Vérifier que toutes les colonnes ont les mêmes rownames
  rownames_list <- lapply(PoplarData, function(x) rownames(x))
  unique_rownames <- Reduce(union, rownames_list)

  # Vérifier les différences de rownames
  for (i in seq_along(rownames_list)) {
    missing_in_table <- setdiff(unique_rownames, rownames_list[[i]])
    if (length(missing_in_table) > 0) {
      cat("Les rownames suivants sont absents dans la table", names(PoplarData)[i], ":", missing_in_table, "\n")
    }
  }

  # Vérifier les valeurs manquantes
  missing_values <- sapply(PoplarData, function(x) any(is.na(x)))
  if (any(missing_values)) {
    stop("Il y a des valeurs manquantes dans les données.")
  }

  # Vérifier que toutes les colonnes dans chaque table sont numériques
  non_numeric_columns <- sapply(PoplarData, function(x) {
    if (is.data.frame(x)) {
      return(any(!sapply(x, is.numeric)))
    }
    return(!is.numeric(x))
  })

  if (any(non_numeric_columns)) {
    stop("Il y a des colonnes non numériques dans les données.")
  }

  # Afficher les données avec des rownames différents
  cat("Les données qui ont des rownames différents :\n")
  for (i in seq_along(rownames_list)) {
    extra_rows <- setdiff(rownames_list[[i]], unique_rownames)
    if (length(extra_rows) > 0) {
      cat("Table", names(PoplarData)[i], ": ", extra_rows, "\n")
    }
  }
  

  print(paste0("Phenotype available:", colnames(PoplarData$Phenotype)))
  # Message de confirmation si tout est ok
  message("Les données sont intactes : mêmes rownames, pas de valeurs manquantes et toutes les colonnes sont numériques.")
}

# Exemple d'utilisation :
check_data_integrity(PoplarData)


print("Lets start Prediction ...")

#Init result list
Ridge.mod <- list()

# Regression Genomic only -------------------------------------------------
Ridge.mod$Genomic <- sapply(colnames(PoplarData$Phenotype), 
                            function(trait){
                              Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
                                             First.Omic=PoplarData$Genomic, 
                                             CV.fold=PoplarData$CV.fold[[trait]],
                                             cores = Nbr_core
                                           )
                            }, simplify = FALSE)

save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_Genomic.RData"))

# Regression Epigenomic only -------------------------------------------------
#Ridge.mod$Epigenetic_CG <- sapply(colnames(PoplarData$Phenotype), 
#                               function(trait){
#                                 Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
#                                              First.Omic=PoplarData$Epigenetic_CG, 
#                                              CV.fold=PoplarData$CV.fold[[trait]],
#                                              cores = Nbr_core)
#                               },simplify = FALSE)
#save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_epigenomic_CG.RData"))
#
#Ridge.mod$Epigenetic_CHG <- sapply(colnames(PoplarData$Phenotype), 
#                                  function(trait){
#                                    Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
#                                                 First.Omic=PoplarData$Epigenetic_CHG, 
#                                                 CV.fold=PoplarData$CV.fold[[trait]],
#                                                 cores = Nbr_core)
#                                  },simplify = FALSE)
#
#save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_epigenomic_CHG.RData"))

#Ridge.mod$Epigenetic_CHH <- sapply(colnames(PoplarData$Phenotype), 
#                                  function(trait){
#                                    Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
#                                                 First.Omic=PoplarData$Epigenetic_CHH, 
#                                                 CV.fold=PoplarData$CV.fold[[trait]],
#                                                 cores = Nbr_core)
#                                  },simplify = FALSE)

#save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_epigenomic_CHH.RData"))

#Ridge.mod$Epigenetic <- sapply(colnames(PoplarData$Phenotype), 
#                               function(trait){
#                                 Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
#                                              First.Omic=PoplarData$Epigenetic, 
#                                              CV.fold=PoplarData$CV.fold[[trait]],
#                                              cores = Nbr_core)
#                               },simplify = FALSE)
#
#save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_epigenomic_AllContext.RData"))

# Regression Transcriptomic only ----------------------------------------------
#Ridge.mod$Expr <- sapply(colnames(PoplarData$Phenotype), 
#                         function(trait){
#                            Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
#                                         First.Omic=PoplarData$Transcriptomic, 
#                                         CV.fold=PoplarData$CV.fold[[trait]],
#                                         cores = Nbr_core)
#                            },simplify = FALSE)

#save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_Expr.RData"))

# Regression Concatenate Genomic and Transcriptomic --------------------------
#Ridge.mod$Comb_G_T <- sapply(colnames(PoplarData$Phenotype), 
#                         function(trait){
#                           Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
#                                        First.Omic=PoplarData$Genomic, 
#                                        Second.Omic=PoplarData$Transcriptomic,
#                                        CV.fold=PoplarData$CV.fold[[trait]],
#                                        cores = Nbr_core)
#                         },simplify = FALSE)

#save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_G_T.RData"))

## Regression Concatenate Epigenetic and Transcriptomic ------------------------
#Ridge.mod$Comb_E_T <- sapply(colnames(PoplarData$Phenotype), 
#                         function(trait){
#                           Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
#                                        First.Omic=PoplarData$Epigenetic, 
#                                         Second.Omic=PoplarData$Transcriptomic,
#                                         CV.fold=PoplarData$CV.fold[[trait]],
#                                        cores = Nbr_core)
#                         },simplify = FALSE)
#save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_E_T.RData"))

# Regression Concatenate Epigenetic and Genomic -------------------------------
#Ridge.mod$Comb_E_G <- sapply(colnames(PoplarData$Phenotype), 
#                         function(trait){
#                           Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
#                                        First.Omic=PoplarData$Epigenetic, 
#                                        Second.Omic=PoplarData$Genomic,
#                                        CV.fold=PoplarData$CV.fold[[trait]],
#                                        cores = Nbr_core)
#                         },simplify = FALSE)
#save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_E_G.RData"))

# Regression Concatenate Epigenetic, Transcriptomic and Genomic -------------------------------
#Ridge.mod$Comb_E_G_T <- sapply(colnames(PoplarData$Phenotype), 
#                         function(trait){
#                           Ncv_RidgeReg_2(Y=PoplarData$Phenotype[,trait],
#                                        First.Omic=PoplarData$Epigenetic, 
#                                        Second.Omic=PoplarData$Genomic,
#                                        Third.Omic=PoplarData$Trancriptomic,
#                                        CV.fold=PoplarData$CV.fold[[trait]],
#                                        cores = Nbr_core)
#                         },simplify = FALSE)
#save(Ridge.mod, file = paste0(Output_dir, "Pred_MultiOmic_Output_E_G_T.RData"))

# Save final Output -----------------------------------------------------------------
#save(Ridge.mod, file = paste0(Output_dir,"Pred_MultiOmic_Output_all.RData"))






