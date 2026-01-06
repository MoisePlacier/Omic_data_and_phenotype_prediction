## ----echo=FALSE----------------------------------------------------------------------------------
library(rMVP)
library(randomForest)
library(yaml)
library(ggplot2)
library(data.table)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript rMVP_GWAS.R config_rMVP_GWAS.yml")
}

config <- yaml::read_yaml(args[1])
message("Config chargée depuis : ", config)

## ------------------------------------------------------------------------------------------------
#config <- read_yaml("config_rMVP_GWAS.yml")
dir.create(config$settings$output_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------------------------------------------
load(config$data$phenotype)    # Phenotype
load(config$data$genomic)      # Genomic
load(config$data$localisation)
load(config$data$kinship)


## ------------------------------------------------------------------------------------------------
# QC génétique : MAF & Missing
geno <- Genomic

maf <- colMeans(geno, na.rm = TRUE) / 2
maf <- pmin(maf, 1 - maf)
missing <- colMeans(is.na(geno))

keep_snps <- which(
  maf >= config$analysis$maf_threshold &
    missing <= config$analysis$missing_threshold
)

geno_qc <- geno[, keep_snps]


qc_dir <- file.path(config$settings$output_dir, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  data.frame(MAF = maf[keep_snps], Missing = missing[keep_snps]),
  file.path(qc_dir, "maf_missing_summary.csv"),
  row.names = FALSE
)


## ------------------------------------------------------------------------------------------------
common_ids <- Reduce(intersect, list(
  rownames(Genomic),
  rownames(Phenotype),
  localisation$ID,
  rownames(Kinship_Poplar)
))

Genomic    <- Genomic[common_ids, ]
Phenotype  <- Phenotype[common_ids, , drop = FALSE]
localisation <- localisation[match(common_ids, localisation$ID), ]
Kinship_Poplar_matrix <- Kinship_Poplar[common_ids, common_ids]

message(sprintf("Alignement effectué : %d individus communs", length(common_ids)))


## ------------------------------------------------------------------------------------------------
marker_names <- colnames(Genomic)

chr <- rep(0, length(marker_names))
pos <- seq_along(marker_names)

## Cas 1 : Chromosomes principaux (Chr01_12345)
is_chr <- grepl("^Chr", marker_names)
chr[is_chr] <- as.numeric(gsub("Chr0?(\\d+)_.*", "\\1", marker_names[is_chr]))
pos[is_chr] <- as.numeric(gsub(".*_", "", marker_names[is_chr]))

## Cas 2 : Scaffolds non assignés (scaffold_104_5033)
is_scaffold <- grepl("^scaffold", marker_names)
chr[is_scaffold] <- 20
pos[is_scaffold] <- as.numeric(gsub(".*_", "", marker_names[is_scaffold]))

na_pos <- is.na(pos)
if (any(na_pos)) {
  warning(sum(na_pos), " SNP avec position non reconnue → positions séquentielles assignées")
  positions[na_pos] <- seq_len(sum(na_pos))
  chromosomes[na_pos] <- 0
}

map_info <- data.frame(
  SNP = marker_names,
  Chromosome = chr,
  Position = pos
)



## ------------------------------------------------------------------------------------------------

# rMVP exige : SNP en lignes, individus en colonnes
geno_mat <- t(as.matrix(Genomic))
storage.mode(geno_mat) <- "numeric"



## ------------------------------------------------------------------------------------------------
# Création du dossier de sortie MVP
mvp_dir <- file.path(config$settings$output_dir, "MVP_Data")
dir.create(mvp_dir, recursive = TRUE, showWarnings = FALSE)

geno_num_file <- file.path(mvp_dir, "geno_numeric.txt")
map_file      <- file.path(mvp_dir, "geno_map.txt")
phe_file      <- file.path(mvp_dir, "phenotype.txt")

write.table(
  geno_mat,
  file = geno_num_file,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

# Forcer les types AVANT l’export sinon erreur dans la matrice des SNP
map_info <- data.frame(
  SNP = as.character(marker_names),
  Chromosome = as.integer(chr),
  Position = as.numeric(pos),
  stringsAsFactors = FALSE
)

write.table(
  map_info,
  file = map_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  data.frame(ID = rownames(Phenotype), Phenotype),
  phe_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## ------------------------------------------------------------------------------------------------
MVP.Data.Numeric2MVP(
  num_file = geno_num_file,
  map_file = map_file,
  out = "mvp_poplar",
  auto_transpose = FALSE,  # SNP en lignes
  type.geno = "char",
  verbose = TRUE
)

# Import du phénotype (format MVP)
MVP.Data.Pheno(
  pheno_file = phe_file,
  out = "mvp_poplar",
  sep = "\t"
)


## ------------------------------------------------------------------------------------------------
geno_big <- attach.big.matrix("mvp_poplar.geno.desc")
map      <- read.table("mvp_poplar.geno.map", header = TRUE)
phe      <- read.table("mvp_poplar.phe", header = TRUE)

message(" Données rMVP prêtes (geno / map / phe)")


## ------------------------------------------------------------------------------------------------

geno_big <- attach.big.matrix("mvp_poplar.geno.desc")

pca <- MVP.PCA(geno_big, pcs.keep = 10, cpu = config$settings$n_cores)

# sauvegarde des resultats
pca_dir <- file.path(config$settings$output_dir, "PCA")
dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  pca,
  file = file.path(pca_dir, "pca_scores.csv"),
  row.names = FALSE
)

pca_graph <- paste0(config$settings$output_dir, "/PCA/PCA_PC1_PC2.png")
png(pca_graph, 800, 600)
plot(pca[,1], pca[,2], pch=19)
dev.off()



## ------------------------------------------------------------------------------------------------
# création de la colonne ID
Phenotype$ID <- rownames(Phenotype)

if (config$analysis$k_r) {
  # Calcul de la matrice de parenté selon VanRaden
  Kinship_Poplar_matrix <- MVP.K.VanRaden(geno_big, verbose = TRUE)
}
# GWAS rMVP (boucle multi-phénotypes)
for (trait in config$analysis$phenotypes) {

  out_dir <- file.path("results/GWAS", trait)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  phe_trait <- Phenotype[, c("ID", trait)]

  res <- MVP(
    phe = phe_trait,
    geno = geno_big,
    map = map,
    K = Kinship_Poplar_matrix,
    nPC.MLM = config$analysis$n_pcs,      # N composantes de l'ACP calculé interne
    nPC.FarmCPU = config$analysis$n_pcs,
    method = config$gwas$methods,         # "MLM" et "FarmCPU"
    vc.method = config$gwas$vc_method,    # "EMMA" ou "P3D"
    ncpus = config$settings$n_cores,
    file.output = TRUE,
    out = out_dir
  )
}


