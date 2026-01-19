library(data.table)
library(ggplot2)
library(mclust)

library(transport)




# ============================================================
# Outils internes
# ============================================================

.check_methods_available <- function(dt, pheno, m1, m2) {
  ok1 <- dt[phenotype == pheno & method_context == m1, .N] > 0
  ok2 <- dt[phenotype == pheno & method_context == m2, .N] > 0
  if (!ok1 || !ok2) {
    stop("Données manquantes pour au moins une méthode et ce phénotype.")
  }
}

.compute_window_ari <- function(score1, score2, topK) {
  if (length(score1) < 2) return(NA_real_)
  k_eff <- min(topK, length(score1) - 1)

  thr1 <- quantile(score1, probs = 1 - k_eff / length(score1), na.rm = TRUE)
  thr2 <- quantile(score2, probs = 1 - k_eff / length(score2), na.rm = TRUE)

  part1 <- as.integer(score1 >= thr1)
  part2 <- as.integer(score2 >= thr2)

  mclust::adjustedRandIndex(part1, part2)
}

# ============================================================
# ARI par fenêtres physiques (paires de bases)
# ============================================================

compute_ari_by_bp_window <- function(
    dt, phenotype, method1, method2,
    window_bp = 5e5, topK = 100
) {

  message(sprintf(
    "ARI par fenêtres physiques (%s vs %s) – phénotype %s",
    method1, method2, phenotype
  ))

  .check_methods_available(dt, phenotype, method1, method2)

  dt_work <- dt[
    phenotype == phenotype & method_context %in% c(method1, method2),
    .(score = mean(mean)),
    by = .(SNP, chrom_num, pos, method_context)
  ]

  dt_work[, window_id := floor(pos / window_bp)]

  dt_wide <- dcast(
    dt_work,
    SNP + chrom_num + pos + window_id ~ method_context,
    value.var = "score"
  )

  setnames(dt_wide, c(method1, method2), c("score1", "score2"))

  res <- dt_wide[
    ,
    .(
      ari = .compute_window_ari(score1, score2, topK),
      n_snps = .N,
      mid_pos = window_id[1] * window_bp + window_bp / 2
    ),
    by = .(chrom_num, window_id)
  ]

  res[!is.na(ari)]
}

# ============================================================
# ARI par fenêtres de SNPs (effectif fixe)
# ============================================================

compute_ari_by_snp_window <- function(
    dt, phenotype, method1, method2,
    snps_per_window = 300, topK = 100
) {

  message(sprintf(
    "ARI par fenêtres de %s SNPs (%s vs %s)",
    snps_per_window, method1, method2
  ))

  .check_methods_available(dt, phenotype, method1, method2)

  dt_work <- dt[
    phenotype == phenotype & method_context %in% c(method1, method2),
    .(score = mean(mean)),
    by = .(SNP, chrom_num, pos, method_context)
  ]

  setorder(dt_work, chrom_num, pos)

  dt_work[, snp_rank := seq_len(.N), by = chrom_num]
  dt_work[, window_id := floor((snp_rank - 1) / snps_per_window)]

  dt_wide <- dcast(
    dt_work,
    SNP + chrom_num + pos + window_id ~ method_context,
    value.var = "score"
  )

  setnames(dt_wide, c(method1, method2), c("score1", "score2"))

  res <- dt_wide[
    ,
    .(
      ari = .compute_window_ari(score1, score2, topK),
      n_snps = .N,
      mid_pos = mean(pos, na.rm = TRUE)
    ),
    by = .(chrom_num, window_id)
  ]

  res[!is.na(ari)]
}

# ============================================================
# Comparaisons multiples
# ============================================================

compute_multi_ari <- function(
    dt, phenotype, method_pairs,
    mode = c("snp", "bp"),
    snps_per_window = 300,
    window_bp = 5e5,
    topK = 100
) {

  mode <- match.arg(mode)

  rbindlist(lapply(method_pairs, function(p) {
    res <- if (mode == "snp") {
      compute_ari_by_snp_window(
        dt, phenotype, p[1], p[2],
        snps_per_window, topK
      )
    } else {
      compute_ari_by_bp_window(
        dt, phenotype, p[1], p[2],
        window_bp, topK
      )
    }
    res[, pair := sprintf("%s vs %s", p[1], p[2])]
    res
  }))
}



# ============================================================
# Visualisation
# ============================================================

plot_ari_profiles <- function(dt_ari, window_label) {

  ggplot(dt_ari, aes(x = mid_pos / 1e6, y = ari)) +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf, ymin = -0.6, ymax = 0,
      fill = "firebrick3", alpha = 0.05
    ) +
    geom_step(
      aes(group = interaction(pair, chrom_num)),
      linewidth = 0.4, color = "grey40", alpha = 0.7
    ) +
    geom_point(aes(color = ari), size = 1.2) +
    facet_grid(pair ~ chrom_num, scales = "free_x", space = "free_x") +
    scale_color_gradient2(
      low = "firebrick3",
      mid = "grey90",
      high = "darkturquoise",
      midpoint = 0,
      name = "ARI"
    ) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "darkturquoise", linewidth = 0.3) +
    theme_minimal() +
    theme(
      strip.text.y = element_text(angle = 0, size = 7, face = "bold"),
      strip.text.x = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.2, "lines"),
      axis.text.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = sprintf("Concordance génomique – ARI (%s)", window_label),
      subtitle = "Consensus élevé (bleu) vs divergence non linéaire (rouge)",
      x = "Position génomique (Mb)",
      y = "Adjusted Rand Index"
    )
}



# ============================================================
# Outils internes : Calcul EMD (Wasserstein)
# ============================================================

.check_methods_available <- function(dt, pheno, m1, m2) {
  # Vérifie si les méthodes demandées existent dans la colonne 'method_context'
  available <- unique(dt[phenotype == pheno, method_context])
  if (!all(c(m1, m2) %in% available)) {
    stop(sprintf("Méthode introuvable pour ce phénotype.\nDisponibles: %s\nDemandées: %s, %s",
                 paste(available, collapse=", "), m1, m2))
  }
}

.compute_window_emd <- function(score1, score2, positions) {

  # 1. Sécurité : On a besoin de données et de sommes > 0
  # On remplace les NA par 0 pour le calcul de masse
  s1 <- replace(score1, is.na(score1), 0)
  s2 <- replace(score2, is.na(score2), 0)

  # Si une méthode n'a aucun signal dans cette fenêtre (somme = 0),
  # la distance n'est pas définie ou est maximale. On retourne NA.
  if (sum(s1) <= 0 || sum(s2) <= 0 || length(positions) < 2) {
    return(NA_real_)
  }

  # 2. Normalisation des MASSES (Somme = 1)
  # C'est impératif pour l'EMD : on compare la répartition de probabilité
  w1 <- s1 / sum(s1)
  w2 <- s2 / sum(s2)

  # 3. Normalisation des POSITIONS (Scaling 0-1)
  # Permet de comparer des fenêtres de tailles différentes
  # Une distance de 0.1 signifie un déplacement de 10% de la taille de la fenêtre
  min_p <- min(positions)
  max_p <- max(positions)

  if (max_p == min_p) return(0)

  pos_norm <- (positions - min_p) / (max_p - min_p)

  # 4. Calcul Wasserstein
  tryCatch({
    transport::wasserstein1d(a = pos_norm, b = pos_norm, p = 1, wa = w1, wb = w2)
  }, error = function(e) return(NA_real_))
}


# ============================================================
# Calcul de l'EMD Globale (Génome Linéaire)
# ============================================================
compute_global_linear_emd <- function(dt, phenotype_target, method1, method2, score_col = "mean") {

  # 1. Nettoyage des noms (Sécurité factor/character)
  m1 <- as.character(method1)
  m2 <- as.character(method2)
  pheno <- as.character(phenotype_target)

  # 2. Extraction et conversion en positif
  # On ne garde que ce qui est nécessaire pour alléger
  dt_sub <- dt[phenotype == pheno & method_context %in% c(m1, m2)]

  if (nrow(dt_sub) == 0) return(NA_real_)

  # 3. Position Cumulative (Enchaîner les chromosomes)
  # On s'assure que chrom_num est traité de la même façon partout
  dt_sub[, chr_id := as.numeric(as.character(chrom_num))]
  chr_info <- dt_sub[, .(max_p = max(pos, na.rm = TRUE)), by = chr_id][order(chr_id)]
  chr_info[, offset := shift(cumsum(as.numeric(max_p)), fill = 0)]

  dt_work <- merge(dt_sub, chr_info[, .(chr_id, offset)], by = "chr_id")
  dt_work[, pos_global := pos + offset]

  # 4. Préparation des scores (Valeur absolue pour l'EMD)
  dt_work[, score_abs := abs(get(score_col))]
  dt_work[is.na(score_abs), score_abs := 0] # Remplacer NA par 0

  # 5. Pivot pour aligner les SNPs
  dt_wide <- dcast(dt_work, SNP + pos_global ~ method_context, value.var = "score_abs", fill = 0)

  # Vérification des colonnes après pivot
  if (!(m1 %in% colnames(dt_wide)) | !(m2 %in% colnames(dt_wide))) return(NA_real_)

  # 6. Calcul CDF manuel (Mathématiquement robuste)
  setorder(dt_wide, pos_global)

  s1 <- dt_wide[[m1]]
  s2 <- dt_wide[[m2]]

  sum1 <- sum(s1); sum2 <- sum(s2)
  if (sum1 == 0 | sum2 == 0) return(NA_real_)

  w1 <- s1 / sum1
  w2 <- s2 / sum2

  # Normalisation de la position entre 0 et 1
  p_all <- dt_wide$pos_global
  p_norm <- (p_all - min(p_all)) / (max(p_all) - min(p_all))

  # EMD 1D = Somme des différences cumulées
  cdf1 <- cumsum(w1)
  cdf2 <- cumsum(w2)
  diff_p <- diff(p_norm)

  # Calcul final
  res <- sum(abs(cdf1[-length(cdf1)] - cdf2[-length(cdf2)]) * diff_p)

  return(res)
}
# ============================================================
# EMD par Fenêtre Physique (Sliding Window)
# ============================================================

compute_emd_by_bp_window <- function(
    dt, phenotype_target, method1, method2,
    window_bp = 1e6,    # Taille fenêtre en pb (défaut 1Mb)
    score_col = "mean"  # Colonne à utiliser comme poids
) {
  message(sprintf("Calcul EMD (Fenêtres) : %s vs %s [%s]", method1, method2, phenotype_target))

  # Vérification
  .check_methods_available(dt, phenotype_target, method1, method2)

  # 1. Extraction des données
  # On filtre sur le phénotype et les méthodes
  dt_work <- dt[
    phenotype == phenotype_target & method_context %in% c(method1, method2),
    .(SNP, chrom_num, pos, method_context, score_val = get(score_col))
  ]

  # 2. Création de l'ID fenêtre
  dt_work[, window_id := floor(pos / window_bp)]

  # 3. Pivot (Format Large)
  # fill = 0 est crucial : si un SNP n'est pas sélectionné par une méthode, son poids est 0
  dt_wide <- dcast(
    dt_work,
    SNP + chrom_num + pos + window_id ~ method_context,
    value.var = "score_val",
    fill = 0
  )

  # Renommage des colonnes dynamiques vers des noms fixes pour le calcul
  setnames(dt_wide, c(method1, method2), c("score1", "score2"))

  # 4. Tri par position (Requis par transport::wasserstein1d)
  setorder(dt_wide, chrom_num, pos)

  # 5. Calcul par groupe (Chromosome + Fenêtre)
  res <- dt_wide[, .(
    emd = .compute_window_emd(score1, score2, pos),
    n_snps = .N,
    mid_pos = window_id[1] * window_bp + window_bp / 2
  ),
  by = .(chrom_num, window_id)
  ]

  # Nettoyage des NA (fenêtres vides)
  return(res[!is.na(emd)])
}

# ============================================================
# EMD par Chromosome (Global)
# ============================================================

compute_emd_by_chromosome <- function(
    dt, phenotype_target, method1, method2,
    score_col = "mean"
) {
  message(sprintf("Calcul EMD (Chromosome) : %s vs %s", method1, method2))
  .check_methods_available(dt, phenotype_target, method1, method2)

  dt_work <- dt[
    phenotype == phenotype_target & method_context %in% c(method1, method2),
    .(SNP, chrom_num, pos, method_context, score_val = get(score_col))
  ]

  dt_wide <- dcast(dt_work, SNP + chrom_num + pos ~ method_context, value.var = "score_val", fill = 0)
  setnames(dt_wide, c(method1, method2), c("score1", "score2"))
  setorder(dt_wide, chrom_num, pos)

  res <- dt_wide[, .(
    emd = .compute_window_emd(score1, score2, pos),
    n_snps = .N,
    mid_pos = median(pos) # Pour placer le point au centre du chr sur le graph
  ),
  by = .(chrom_num)
  ]

  return(res[!is.na(emd)])
}



plot_emd_profiles <- function(dt_emd, title_label = "") {

  ggplot(dt_emd, aes(x = mid_pos / 1e6, y = emd)) +
    # Zone "Accord" (EMD faible)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.05,
             fill = "green", alpha = 0.1) +

    # Trace
    geom_step(aes(group = chrom_num), linewidth = 0.6, color = "grey30") +
    geom_point(aes(color = emd), size = 1.8, alpha = 0.9) +

    # Facettes par chromosome (échelle X libre pour respecter la taille des chr)
    facet_grid(. ~ chrom_num, scales = "free_x", space = "free_x") +

    # Couleurs : Vert (Proche) -> Rouge (Divergent)
    scale_color_gradient(low = "#1a9850", high = "#d73027", name = "Divergence\n(EMD)") +

    theme_minimal() +
    theme(
      strip.text.x = element_text(face = "bold", size = 9),
      panel.spacing.x = unit(0.15, "lines"),
      axis.text.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = paste("Paysage de Divergence Génomique", title_label),
      subtitle = "Comparaison spatiale des scores d'importance (Wasserstein)",
      x = "Position sur le Chromosome",
      y = "EMD (Distance de transport)"
    )
}

plot_emd_multi_profiles <- function(dt_multi, title_label = "") {

  # On s'assure que comparaison est un facteur pour garder l'ordre de ta liste pairs
  dt_multi[, comparaison := factor(comparaison, levels = unique(comparaison))]

  ggplot(dt_multi, aes(x = mid_pos / 1e6, y = emd)) +
    # Fond coloré pour les zones de faible divergence (Accord)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.02,
             fill = "#1a9850", alpha = 0.05) +

    # Lignes de liaison pour l'aspect "profil"
    geom_line(aes(group = chrom_num), color = "grey80", linewidth = 0.3) +

    # Points colorés selon l'intensité du décrochage
    geom_point(aes(color = emd), size = 0.8, alpha = 0.8) +

    # LA CLÉ : Facettes croisées (Lignes = Comparaisons, Colonnes = Chromosomes)
    facet_grid(comparaison ~ chrom_num, scales = "free_x", space = "free_x") +

    # Échelle de couleur Divergente (Bleu = Accord, Rouge = Divergence)
    scale_color_gradientn(
      colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"),
      name = "Divergence\n(EMD)"
    ) +

    theme_minimal() +
    theme(
      strip.text.y = element_text(angle = 0, face = "bold", size = 8, hjust = 0),
      strip.text.x = element_text(face = "bold", size = 7),
      panel.spacing.x = unit(0.1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = paste("Paysages de Divergence Génomique", title_label),
      subtitle = "Comparaison de la répartition spatiale des scores d'importance",
      x = "Position sur les Chromosomes (Mb)",
      y = "EMD (Distance de transport locale)"
    )
}

plot_emd_barplot_chr <- function(dt_chr, title_label = "") {

  # On s'assure que chrom_num est bien ordonné numériquement
  dt_chr[, chrom_num := factor(chrom_num, levels = sort(as.numeric(unique(chrom_num))))]
  dt_chr[, comparaison := factor(comparaison, levels = unique(comparaison))]

  ggplot(dt_chr, aes(x = chrom_num, y = emd, fill = emd)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.2) +

    # Facettage par comparaison (une ligne par paire de méthodes)
    facet_wrap(~ comparaison, ncol = 1) +

    # Échelle de couleur pour souligner l'intensité
    scale_fill_gradientn(
      colors = c("#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
      name = "Divergence\n(EMD)"
    ) +

    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste("Divergence Structurelle par Chromosome", title_label),
      subtitle = "Comparaison de l'Earth Mover's Distance à l'échelle chromosomique entière",
      x = "Numéro du Chromosome",
      y = "EMD (Distance de Wasserstein)",
      caption = "Une EMD élevée indique que les méthodes localisent les pics d'importance à des endroits différents du chromosome."
    )
}





######################################################################
############# La Fonction "Masterclass" LD-EMD #####################
########################################################################

compute_multi_ld_emd <- function(
    geno_mat,           # Matrice individus x SNPs (0/1/2) avec colnames = SNP IDs
    scores_dt,          # Ton data.table all_scores_dt
    phenotype_target,   # Le phénotype cible (ex: "Rust")
    pairs_list,         # Liste des paires : list(c("M1","M2"), ...)
    top_k = 5000,       # Nombre de SNPs à prendre dans le Top de CHAQUE méthode
    score_col = "mean"  # La colonne d'importance
) {

  # --- 1. Initialisation ---
  results_list <- list()
  total_pairs <- length(pairs_list)

  message(sprintf("=== DÉBUT ANALYSE LD-EMD : %d Paires à traiter ===", total_pairs))
  message(sprintf("Phénotype : %s | Top K : %d SNPs par méthode", phenotype_target, top_k))

  # Barre de progression textuelle
  pb <- txtProgressBar(min = 0, max = total_pairs, style = 3)

  # Filtrage initial des scores pour le phénotype
  dt_pheno <- scores_dt[phenotype == phenotype_target]

  # --- 2. Boucle sur les paires ---
  for (i in seq_along(pairs_list)) {
    pair <- pairs_list[[i]]
    m1 <- pair[1]
    m2 <- pair[2]

    # a. Extraction des Top K SNPs pour chaque méthode
    # On trie par valeur absolue décroissante
    snps_m1 <- dt_pheno[method_context == m1][order(-abs(get(score_col)))][1:top_k, SNP]
    snps_m2 <- dt_pheno[method_context == m2][order(-abs(get(score_col)))][1:top_k, SNP]

    # b. Création de l'UNION (La zone de jeu)
    # On ne garde que les SNPs qui existent vraiment dans la matrice de génotype
    union_candidates <- unique(c(snps_m1, snps_m2))
    union_candidates <- union_candidates[!is.na(union_candidates)]

    # Intersection avec les génotypes disponibles (Sécurité anti-crash)
    valid_snps <- intersect(union_candidates, colnames(geno_mat))

    n_original <- length(union_candidates)
    n_valid <- length(valid_snps)

    if (n_valid < 10) {
      warning(sprintf("Pas assez de SNPs communs avec la matrice génotypique pour %s vs %s", m1, m2))
      setTxtProgressBar(pb, i)
      next
    }

    # c. Préparation des scores alignés
    # On extrait les scores pour ces SNPs valides
    dt_subset <- dt_pheno[SNP %in% valid_snps & method_context %in% c(m1, m2)]

    # Pivot pour avoir une colonne par méthode (alignement parfait)
    # fill = 0 car si un SNP est dans l'union mais pas dans la méthode (hors top), son score est négligeable
    dt_wide <- dcast(dt_subset, SNP ~ method_context, value.var = score_col, fill = 0)

    # Sécurité : alignement strict de l'ordre des lignes avec les colonnes de geno
    # On trie dt_wide selon l'ordre des SNPs dans valid_snps pour matcher le subset de génotype
    dt_wide <- dt_wide[match(valid_snps, SNP)]

    # d. Extraction Génotype & Calcul LD (Ground Metric)
    # C'est l'étape lourde : on extrait les colonnes correspondant aux SNPs
    geno_sub <- geno_mat[, valid_snps]

    # Calcul de r^2 (Corrélation au carré)
    # use="pairwise" gère les quelques NAs éventuels dans le génotypage
    ld_matrix <- cor(geno_sub, use = "pairwise.complete.obs")^2
    ld_matrix[is.na(ld_matrix)] <- 0 # Sécurité

    # Distance = 1 - r^2
    cost_matrix <- 1 - ld_matrix

    # e. Normalisation des masses (Poids)
    # On utilise abs() car l'EMD requiert des masses positives
    w1 <- abs(dt_wide[[m1]])
    w2 <- abs(dt_wide[[m2]])

    # Normalisation pour sommer à 1
    w1 <- w1 / sum(w1)
    w2 <- w2 / sum(w2)

    # f. Calcul EMD (Transport Optimal)
    # P = 1 pour la distance de Wasserstein standard
    emd_val <- transport::wasserstein(w1, w2, costm = cost_matrix, p = 1)

    # Stockage
    results_list[[i]] <- data.table(
      method1 = m1,
      method2 = m2,
      comparaison = paste(m1, "vs", m2),
      emd_ld = emd_val,
      n_snps_used = n_valid,
      pct_snps_found = round((n_valid / n_original) * 100, 1)
    )

    # Nettoyage mémoire immédiat (Crucial pour la RAM)
    rm(geno_sub, ld_matrix, cost_matrix, dt_wide)
    gc(verbose = FALSE)

    # Update barre
    setTxtProgressBar(pb, i)
  }

  close(pb)
  message("\n=== Calcul terminé ! ===")

  return(rbindlist(results_list))
}
