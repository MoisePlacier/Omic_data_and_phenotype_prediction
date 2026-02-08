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

.compute_window_emd <- function(s1, s2, pos) {
  # 1. Traitement des scores : Seuil à zéro (essentiel pour RF corrigée)
  # On remplace aussi les NA par 0 au passage
  w1 <- pmax(replace(s1, is.na(s1), 0), 0)
  w2 <- pmax(replace(s2, is.na(s2), 0), 0)

  sum1 <- sum(w1)
  sum2 <- sum(w2)

  # 2. Sécurité : il faut du signal positif et au moins 2 SNPs
  if (sum1 <= 0 || sum2 <= 0 || length(pos) < 2) {
    return(NA_real_)
  }

  # 3. Normalisation des MASSES (Somme = 1)
  w1 <- w1 / sum1
  w2 <- w2 / sum2

  # 4. Normalisation des POSITIONS (Scaling 0-1)
  # Permet de comparer l'EMD entre chr de tailles différentes
  min_p <- min(pos)
  max_p <- max(pos)
  if (max_p == min_p) return(0)
  pos_norm <- (pos - min_p) / (max_p - min_p)

  # 5. Calcul Wasserstein 1D
  # Note : wa/wb sont les poids, a/b sont les positions
  res <- tryCatch({
    transport::wasserstein1d(a = pos_norm, b = pos_norm, p = 1, wa = w1, wb = w2)
  }, error = function(e) NA_real_)

  return(res)
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
# EMD par Chromosome
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

  # --- Sécurité : Copie locale ---
  dt_plot <- copy(dt_chr)

  # --- Ajustement des colonnes ---
  if(!"emd" %in% names(dt_plot) & "emd_ld" %in% names(dt_plot)) {
    setnames(dt_plot, "emd_ld", "emd")
  }

  if(!"comparaison" %in% names(dt_plot)) {
    dt_plot[, comparaison := paste(method1, "vs", method2)]
  }

  # --- Préparation des facteurs ---
  dt_plot[, chrom_num := factor(chrom_num, levels = sort(as.numeric(unique(chrom_num))))]
  dt_plot[, comparaison := factor(comparaison, levels = unique(comparaison))]

  # --- Graphique ---
  ggplot(dt_plot, aes(x = chrom_num, y = emd, fill = emd)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.2) +

    # Facettage par comparaison
    facet_wrap(~ comparaison, ncol = 1) +

    # Échelle de couleur FIXÉE entre 0 et 1
    scale_fill_gradientn(
      colors = c("#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
      name = "Divergence\n(EMD)",
      limits = c(0, 0.2),            # Force l'échelle de 0 à 1
      breaks = seq(0, 1, 0.05),     # Optionnel : affiche des paliers tous les 0.2
      oob = scales::squish         # Maintient les valeurs hors-bornes dans le dégradé
    ) +

    # Fixer l'axe Y pour éviter les distorsions visuelles entre paires
    coord_cartesian(ylim = c(0, 0.15)) +

    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 10, colour = "grey20"),
      axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste("Divergence Structurelle par Chromosome", title_label),
      subtitle = "Échelle normalisée (0 = Consensus parfait, 1 = Opposition totale)",
      x = "Numéro du Chromosome",
      y = "EMD (Normalisée)",
      caption = "L'échelle de couleur est fixe (0-1) pour permettre la comparaison entre différents phénotypes."
    )
}


######################################################################
############# La Fonction "Masterclass" LD-EMD #####################
########################################################################

library(data.table)
library(transport)
library(stats)

compute_multi_ld_emd <- function(
    geno_mat,           # Matrice individus x SNPs (0/1/2) - colnames = SNP IDs
    scores_dt,          # Data.table contenant les scores (SNP, method, context, mean...)
    phenotype_target,   # Nom du phénotype cible (ex: "Rust")
    pairs_list,         # Liste des paires à comparer : list(c("M1","M2"), ...)
    top_k = 5000,       # Nombre de SNPs à considérer par méthode
    score_col = "mean"  # Colonne utilisée pour le poids (importance)
) {

  # --- 1. Initialisation & Contrôles ---
  if (!phenotype_target %in% unique(scores_dt$phenotype)) {
    stop(sprintf("Le phénotype '%s' est introuvable dans scores_dt.", phenotype_target))
  }

  results_list <- list()
  total_pairs <- length(pairs_list)

  message(sprintf("\n=== ANALYSE LD-EMD (Optimisée) ==="))
  message(sprintf("Phénotype : %s | Top K : %d | Paires : %d", phenotype_target, top_k, total_pairs))

  # Barre de progression
  pb <- txtProgressBar(min = 0, max = total_pairs, style = 3)

  # Pré-filtrage des scores pour le phénotype (Gain de temps dans la boucle)
  dt_pheno <- scores_dt[phenotype == phenotype_target]

  # --- 2. Boucle sur les paires ---
  for (i in seq_along(pairs_list)) {

    pair <- pairs_list[[i]]
    m1 <- pair[1]
    m2 <- pair[2]

    # --- A. Identification de l'Espace de Jeu (Union des Top K) ---
    # On récupère les IDs des meilleurs SNPs pour chaque méthode
    snps_m1 <- dt_pheno[method_context == m1][order(-abs(get(score_col)))][1:min(top_k, .N), SNP]
    snps_m2 <- dt_pheno[method_context == m2][order(-abs(get(score_col)))][1:min(top_k, .N), SNP]

    # Union : On veut comparer tout ce qui est important pour l'un OU l'autre
    union_candidates <- unique(c(snps_m1, snps_m2))

    # Intersection stricte avec les génotypes disponibles (Sécurité anti-crash)
    valid_snps <- intersect(union_candidates, colnames(geno_mat))

    # Si pas assez de SNPs communs, on saute
    if (length(valid_snps) < 5) {
      warning(sprintf("Paire %s vs %s ignorée : < 5 SNPs communs avec la matrice génotypique.", m1, m2))
      setTxtProgressBar(pb, i)
      next
    }

    # --- B. Construction de la Table de Poids (Alignement Robuste) ---
    # 1. Extraction des scores bruts pour l'union des SNPs
    dt_subset <- dt_pheno[SNP %in% valid_snps & method_context %in% c(m1, m2)]

    # 2. Pivot : SNP en ligne, Méthodes en colonnes
    # fill = 0 est un premier filet de sécurité, mais incomplet si le SNP manque totalement pour une méthode
    dt_wide <- dcast(dt_subset, SNP ~ method_context, value.var = score_col, fill = 0)

    # 3. Alignement FORCE sur valid_snps via Merge (Cœur du correctif)
    # Cela garantit que l'ordre des lignes correspondra exactement aux colonnes de geno_sub
    dt_final <- data.table(SNP = valid_snps)
    dt_final <- merge(dt_final, dt_wide, by = "SNP", all.x = TRUE)

    # 4. Nettoyage des NA résiduels (Transformation en 0)
    # Si un SNP est valid_snps mais n'avait pas de score dans dt_wide, merge a mis NA
    for (col in c(m1, m2)) {
      # Utilisation de set() pour la rapidité (pas de copie mémoire)
      set(dt_final, i = which(is.na(dt_final[[col]])), j = col, value = 0)
    }

    # --- C. Calcul de la Matrice de Coût (LD) ---
    # Extraction des génotypes dans le MÊME ORDRE que dt_final
    geno_sub <- geno_mat[, dt_final$SNP]

    # Calcul de r² (Corrélation au carré)
    # use = "pairwise" gère les données manquantes ponctuelles dans le génotype
    ld_matrix <- cor(geno_sub, use = "pairwise.complete.obs")^2

    # Sécurité ultime : si variance nulle quelque part, cor renvoie NA
    ld_matrix[is.na(ld_matrix)] <- 0

    # Distance = 1 - r² (Proche = 0, Distant = 1)
    cost_matrix <- 1 - ld_matrix

    # --- D. Préparation des Poids pour Wasserstein ---
    w1 <- abs(dt_final[[m1]])
    w2 <- abs(dt_final[[m2]])

    # Normalisation (Somme = 1)
    # Cas limite : si une méthode a tous ses scores à 0 (improbable mais possible), distribution uniforme
    sum_w1 <- sum(w1)
    sum_w2 <- sum(w2)

    if (sum_w1 == 0) w1 <- rep(1/length(w1), length(w1)) else w1 <- w1 / sum_w1
    if (sum_w2 == 0) w2 <- rep(1/length(w2), length(w2)) else w2 <- w2 / sum_w2

    # --- E. Calcul EMD ---
    # p = 1 correspond à la distance "Earth Mover's" standard
    tryCatch({
      emd_val <- transport::wasserstein(w1, w2, costm = cost_matrix, p = 1)

      results_list[[i]] <- data.table(
        method1 = m1,
        method2 = m2,
        comparaison = paste(m1, "vs", m2),
        emd_ld = emd_val,
        n_snps_union = length(valid_snps),
        coverage_pct = round(length(valid_snps) / length(union_candidates) * 100, 1)
      )
    }, error = function(e) {
      warning(paste("Erreur transport::wasserstein pour", m1, "vs", m2, ":", e$message))
    })

    # --- F. Nettoyage Mémoire (Optimisation) ---
    # Suppression explicite des objets lourds
    rm(geno_sub, ld_matrix, cost_matrix, dt_wide, dt_final, w1, w2)

    # Garbage Collector "intelligent" : seulement toutes les 5 itérations pour ne pas ralentir
    if (i %% 5 == 0) gc(verbose = FALSE)

    setTxtProgressBar(pb, i)
  }

  close(pb)
  message("\n=== Calcul terminé avec succès ===")

  return(rbindlist(results_list))
}



compute_multi_ld_emd_diffK <- function(
    geno_mat,           # Matrice individus x SNPs (0/1/2) - colnames = SNP IDs
    scores_dt,          # Data.table contenant les scores (SNP, phenotype, method_context, mean...)
    phenotype_target,   # Phénotype cible
    pairs_list,         # Liste de listes : list(list(m1=..., k1=..., m2=..., k2=...), ...)
    score_col = "mean"  # Colonne utilisée pour le poids
) {

  # --- 1. Initialisation & contrôles ---
  if (!phenotype_target %in% unique(scores_dt$phenotype)) {
    stop(sprintf("Le phénotype '%s' est introuvable dans scores_dt.", phenotype_target))
  }

  dt_pheno <- scores_dt[phenotype == phenotype_target]
  results_list <- vector("list", length(pairs_list))

  message("\n=== ANALYSE LD-EMD (Top-K asymétriques, Génome entier) ===")
  message(sprintf("Phénotype : %s | Comparaisons : %d", phenotype_target, length(pairs_list)))

  pb <- txtProgressBar(min = 0, max = length(pairs_list), style = 3)

  # --- 2. Boucle sur les paires ---
  for (i in seq_along(pairs_list)) {

    pair <- pairs_list[[i]]
    m1 <- pair$m1; k1 <- pair$k1
    m2 <- pair$m2; k2 <- pair$k2

    # --- A. Sélection des Top-K spécifiques à chaque méthode ---
    snps_m1 <- dt_pheno[method_context == m1][
      order(-abs(get(score_col)))][1:min(k1, .N), SNP]

    snps_m2 <- dt_pheno[method_context == m2][
      order(-abs(get(score_col)))][1:min(k2, .N), SNP]

    # Union génome entier
    union_candidates <- unique(c(snps_m1, snps_m2))
    valid_snps <- intersect(union_candidates, colnames(geno_mat))

    if (length(valid_snps) < 10) {
      warning(sprintf("Paire %s vs %s ignorée : < 10 SNPs valides.", m1, m2))
      setTxtProgressBar(pb, i)
      next
    }

    # --- B. Construction de la table de poids ---
    dt_subset <- dt_pheno[
      SNP %in% valid_snps & method_context %in% c(m1, m2)
    ]

    dt_wide <- dcast(
      dt_subset, SNP ~ method_context,
      value.var = score_col, fill = 0
    )

    dt_final <- merge(
      data.table(SNP = valid_snps),
      dt_wide, by = "SNP", all.x = TRUE
    )

    for (col in c(m1, m2)) {
      set(dt_final, i = which(is.na(dt_final[[col]])), j = col, value = 0)
    }

    # --- C. Matrice de coût LD globale ---
    geno_sub <- geno_mat[, dt_final$SNP]
    ld_matrix <- cor(geno_sub, use = "pairwise.complete.obs")^2
    ld_matrix[is.na(ld_matrix)] <- 0
    cost_matrix <- 1 - ld_matrix

    # --- D. Poids normalisés ---
    w1 <- abs(dt_final[[m1]])
    w2 <- abs(dt_final[[m2]])

    if (sum(w1) == 0) w1 <- rep(1 / length(w1), length(w1)) else w1 <- w1 / sum(w1)
    if (sum(w2) == 0) w2 <- rep(1 / length(w2), length(w2)) else w2 <- w2 / sum(w2)

    # --- E. Wasserstein ---
    tryCatch({
      emd_val <- transport::wasserstein(w1, w2, costm = cost_matrix, p = 1)

      results_list[[i]] <- data.table(
        method1 = m1,
        k1 = k1,
        method2 = m2,
        k2 = k2,
        comparaison = sprintf("%s(%d) vs %s(%d)", m1, k1, m2, k2),
        emd_ld = emd_val,
        n_snps_union = length(valid_snps),
        coverage_m1 = round(length(intersect(snps_m1, valid_snps)) / k1 * 100, 1),
        coverage_m2 = round(length(intersect(snps_m2, valid_snps)) / k2 * 100, 1)
      )
    }, error = function(e) {
      warning(sprintf("Erreur EMD pour %s vs %s : %s", m1, m2, e$message))
    })

    rm(geno_sub, ld_matrix, cost_matrix, dt_wide, dt_final, w1, w2)
    if (i %% 3 == 0) gc(verbose = FALSE)

    setTxtProgressBar(pb, i)
  }

  close(pb)
  message("\n=== Calcul terminé ===")

  return(rbindlist(results_list, fill = TRUE))
}

############################################################## emd ld chromosome ####################

compute_multi_ld_emd_by_chr <- function(
    geno_mat,           # Matrice individus x SNPs (colnames = SNP IDs)
    scores_dt,          # Data.table avec colonnes SNP, chrom_num, method_context, etc.
    phenotype_target,   # Phénotype cible
    pairs_list,         # Liste des paires : list(c("M1","M2"), ...)
    top_k_per_chr = 500,# Nombre de SNPs à prendre par CHR et par MÉTHODE
    score_col = "mean"  # Colonne d'importance
) {

  # --- 1. Initialisation ---
  if (!phenotype_target %in% unique(scores_dt$phenotype)) {
    stop(sprintf("Le phénotype '%s' est introuvable.", phenotype_target))
  }

  dt_pheno <- scores_dt[phenotype == phenotype_target]
  all_chrs <- sort(unique(dt_pheno$chrom_num))
  results_list <- list()

  message(sprintf("\n=== ANALYSE LD-EMD PAR CHROMOSOME ==="))
  message(sprintf("Phénotype : %s | Top K/Chr : %d", phenotype_target, top_k_per_chr))

  # --- 2. Boucle sur les paires ---
  for (pair in pairs_list) {
    m1 <- pair[1]
    m2 <- pair[2]
    message(sprintf("-> Comparaison : %s vs %s", m1, m2))

    # --- 3. Boucle sur les Chromosomes ---
    for (chr in all_chrs) {

      # a. Sélection du Top K spécifique au chromosome pour chaque méthode
      snps_m1 <- dt_pheno[method_context == m1 & chrom_num == chr][
        order(-abs(get(score_col)))][1:min(top_k_per_chr, .N), SNP]

      snps_m2 <- dt_pheno[method_context == m2 & chrom_num == chr][
        order(-abs(get(score_col)))][1:min(top_k_per_chr, .N), SNP]

      # b. Union et validation contre la matrice de génotype
      union_snps <- unique(c(snps_m1, snps_m2))
      valid_snps <- intersect(union_snps, colnames(geno_mat))

      if (length(valid_snps) < 10) next # Trop peu de signal sur ce chr

      # c. Préparation des poids (Alignement)
      dt_subset <- dt_pheno[SNP %in% valid_snps & method_context %in% c(m1, m2)]
      dt_wide <- dcast(dt_subset, SNP ~ method_context, value.var = score_col, fill = 0)

      # Alignement strict sur l'ordre des SNPs
      dt_final <- merge(data.table(SNP = valid_snps), dt_wide, by = "SNP", all.x = TRUE)
      for (col in c(m1, m2)) set(dt_final, i = which(is.na(dt_final[[col]])), j = col, value = 0)

      # d. Calcul du coût de transport (LD locale)
      geno_sub <- geno_mat[, dt_final$SNP]
      ld_matrix <- cor(geno_sub, use = "pairwise.complete.obs")^2
      ld_matrix[is.na(ld_matrix)] <- 0
      cost_matrix <- 1 - ld_matrix

      # e. Normalisation des masses
      w1 <- abs(dt_final[[m1]])
      w2 <- abs(dt_final[[m2]])

      if (sum(w1) == 0) w1 <- rep(1/length(w1), length(w1)) else w1 <- w1 / sum(w1)
      if (sum(w2) == 0) w2 <- rep(1/length(w2), length(w2)) else w2 <- w2 / sum(w2)

      # f. Calcul EMD
      try({
        emd_val <- transport::wasserstein(w1, w2, costm = cost_matrix, p = 1)

        results_list[[length(results_list) + 1]] <- data.table(
          method1 = m1,
          method2 = m2,
          chrom_num = chr,
          emd_ld = emd_val,
          n_snps = length(valid_snps)
        )
      }, silent = TRUE)
    }
    # Nettoyage entre paires
    gc(verbose = FALSE)
  }

  return(rbindlist(results_list))
}



compute_multi_ld_emd_by_chr_kdependant <- function(
    geno_mat,
    scores_dt,
    phenotype_target,
    pairs_list,         # Liste de listes : list(list(m1="M1", k1=50, m2="M2", k2=5000), ...)
    score_col = "mean"
) {

  # 1. Pré-filtrage global par phénotype et SNPs valides
  dt_pheno <- scores_dt[
    phenotype == phenotype_target &
      SNP %in% colnames(geno_mat),
    .(SNP, chrom_num, pos, method_context, score_val = get(score_col))
  ]

  # On s'assure que les scores sont positifs (seuil à zéro)
  dt_pheno[, score_val := pmax(score_val, 0)]

  all_results <- list()

  for (i in seq_along(pairs_list)) {
    p <- pairs_list[[i]]
    m1 <- p$m1; k1 <- p$k1
    m2 <- p$m2; k2 <- p$k2
    pair_label <- paste0(m1, " (k=", k1, ") vs ", m2, " (k=", k2, ")")

    message("Calcul : ", pair_label)

    # 2. Sélection des Top K SNPs pour chaque méthode (à l'échelle du génome entier)
    top_m1 <- dt_pheno[method_context == m1][order(-score_val)][1:min(.N, k1)]$SNP
    top_m2 <- dt_pheno[method_context == m2][order(-score_val)][1:min(.N, k2)]$SNP

    # 3. Filtrage des données pour cette paire
    dt_pair <- dt_pheno[SNP %in% union(top_m1, top_m2) & method_context %in% c(m1, m2)]

    # 4. Calcul de l'EMD par chromosome
    res_pair <- dt_pair[, {
      # Pivot local
      dt_w <- dcast(.SD, SNP + pos ~ method_context, value.var = "score_val", fill = 0)
      setnames(dt_w, c(m1, m2), c("s1", "s2"))

      # Si une méthode n'a aucun SNP dans le top K sur ce chromosome précis
      if (sum(dt_w$s1) <= 0 || sum(dt_w$s2) <= 0 || .N < 2) {
        .(emd_ld = as.numeric(NA))
      } else {
        # Masses
        w1 <- dt_w$s1 / sum(dt_w$s1)
        w2 <- dt_w$s2 / sum(dt_w$s2)

        # Matrice de coût LD (1-r²)
        geno_win <- geno_mat[, dt_w$SNP, drop = FALSE]
        cost_matrix <- 1 - (cor(geno_win, use = "pairwise.complete.obs")^2)
        cost_matrix[is.na(cost_matrix)] <- 1 # Cas sans variance

        # Transport Optimal
        val <- tryCatch({
          transport::wasserstein(w1, w2, costm = cost_matrix, p = 1)
        }, error = function(e) as.numeric(NA))

        .(emd_ld = val)
      }
    }, by = .(chrom_num)]

    res_pair[, comparaison := pair_label]
    all_results[[i]] <- res_pair
  }

  return(rbindlist(all_results)[!is.na(emd_ld)])
}



plot_ld_emd_barplot_chr <- function(dt_chr, title_label = "") {

  # --- Sécurité : Copie locale ---
  dt_plot <- copy(dt_chr)

  # --- Ajustement des colonnes ---
  if(!"emd" %in% names(dt_plot) & "emd_ld" %in% names(dt_plot)) {
    setnames(dt_plot, "emd_ld", "emd")
  }

  if(!"comparaison" %in% names(dt_plot)) {
    dt_plot[, comparaison := paste(method1, "vs", method2)]
  }

  # --- Préparation des facteurs ---
  dt_plot[, chrom_num := factor(chrom_num, levels = sort(as.numeric(unique(chrom_num))))]
  dt_plot[, comparaison := factor(comparaison, levels = unique(comparaison))]

  # --- Graphique ---
  ggplot(dt_plot, aes(x = chrom_num, y = emd, fill = emd)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.2) +

    # Facettage par comparaison (une ligne par paire de méthodes)
    facet_wrap(~ comparaison, ncol = 1) +

    # Échelle de couleur FIXÉE (0 = Haplotypes identiques, 1 = Haplotypes différents)
    scale_fill_gradientn(
      colors = c("#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
      name = "Divergence\n(LD-EMD)",
      limits = c(0, 0.75),
      oob = scales::squish
    ) +

    # On fixe l'axe Y pour une comparaison visuelle directe
    coord_cartesian(ylim = c(0, 0.75)) +

    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 10, colour = "grey20"),
      axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste("Divergence Structurelle par Chromosome (LD)", title_label),
      subtitle = "Coût basé sur le déséquilibre de liaison (1 - r²)",
      x = "Numéro du Chromosome",
      y = "EMD-LD (Normalisée)",
      caption = "Une valeur proche de 0 indique que les deux méthodes pointent vers les mêmes blocs de LD."
    )
}

####################################### LD par fenètre en pdb ###########################
compute_emd_ld_by_window <- function(
    geno_mat,           # Matrice individus x SNPs (colnames = SNP IDs)
    dt,                 # Data.table avec SNP, chrom_num, pos, method_context, etc.
    phenotype_target,
    pairs_list,         # Liste de paires : list(c("M1", "M2"), ...)
    window_bp = 1e6,    # Taille fenêtre (ex: 1 Mb)
    score_col = "mean"
) {

  # --- 1. Pré-filtrage global ---
  # On ne garde que le phénotype et les SNPs présents dans le génotype pour tout le calcul
  dt_pheno <- dt[
    phenotype == phenotype_target &
      SNP %in% colnames(geno_mat),
    .(SNP, chrom_num, pos, method_context, score_val = get(score_col))
  ]

  if (nrow(dt_pheno) == 0) stop("Aucune donnée trouvée pour ce phénotype ou SNPs absents de geno_mat.")

  all_results <- list()

  message(sprintf("\n=== ANALYSE LD-EMD MULTI-PAIRES [%s] ===", phenotype_target))

  # --- 2. Boucle sur les paires de méthodes ---
  for (i in seq_along(pairs_list)) {

    m1 <- pairs_list[[i]][1]
    m2 <- pairs_list[[i]][2]
    pair_name <- paste(m1, "vs", m2)

    message(sprintf("Calcul en cours (%d/%d) : %s", i, length(pairs_list), pair_name))

    # Sous-ensemble pour les deux méthodes de la paire
    dt_pair <- dt_pheno[method_context %in% c(m1, m2)]

    # Création de l'ID fenêtre
    dt_pair[, window_id := floor(pos / window_bp)]

    # Pivot (Format Large)
    dt_wide <- dcast(
      dt_pair,
      SNP + chrom_num + pos + window_id ~ method_context,
      value.var = "score_val",
      fill = 0
    )

    # Sécurité : vérifier que les colonnes existent bien après le dcast
    if (!m1 %in% names(dt_wide)) dt_wide[[m1]] <- 0
    if (!m2 %in% names(dt_wide)) dt_wide[[m2]] <- 0

    setnames(dt_wide, c(m1, m2), c("score1", "score2"))

    # --- 3. Calcul par groupe (Chromosome + Fenêtre) ---
    res_pair <- dt_wide[, {

      if (.N < 2) {
        .(emd_ld = as.numeric(NA), n_snps = .N, mid_pos = pos[1])
      } else {
        # Masses normalisées
        w1 <- abs(score1); if(sum(w1) > 0) w1 <- w1 / sum(w1) else w1 <- rep(1/.N, .N)
        w2 <- abs(score2); if(sum(w2) > 0) w2 <- w2 / sum(w2) else w2 <- rep(1/.N, .N)

        # Coût (LD locale)
        geno_win <- geno_mat[, SNP]
        ld_matrix <- cor(geno_win, use = "pairwise.complete.obs")^2
        ld_matrix[is.na(ld_matrix)] <- 0
        cost_matrix <- 1 - ld_matrix

        # EMD
        val <- tryCatch({
          transport::wasserstein(w1, w2, costm = cost_matrix, p = 1)
        }, error = function(e) as.numeric(NA))

        .(
          emd_ld = val,
          n_snps = .N,
          mid_pos = (window_id[1] * window_bp) + (window_bp / 2)
        )
      }
    }, by = .(chrom_num, window_id)]

    # Ajout des métadonnées de la paire
    res_pair[, `:=`(
      method1 = m1,
      method2 = m2,
      pair = pair_name
    )]

    all_results[[i]] <- res_pair[!is.na(emd_ld)]

    # Nettoyage mémoire entre les paires
    rm(dt_pair, dt_wide, res_pair); gc(verbose = FALSE)
  }

  message("=== Analyse multi-paires terminée ===")
  return(rbindlist(all_results))
}


plot_emd_ld_window_profile <- function(dt_emd, window_label = "") {

  # 1. Sécurité : Création de la colonne pair
  if (!"pair" %in% names(dt_emd) & all(c("method1", "method2") %in% names(dt_emd))) {
    dt_emd[, pair := paste(method1, "vs", method2)]
  }

  # 2. Nettoyage : On ne garde que les groupes ayant au moins 2 fenêtres
  # pour éviter l'erreur "Each group consists of only one observation"
  dt_plot <- dt_emd[, if (.N >= 2) .SD, by = .(pair, chrom_num)]

  # Si jamais il y a des points isolés qu'on veut quand même voir :
  dt_isolated <- dt_emd[, if (.N < 2) .SD, by = .(pair, chrom_num)]

  p <- ggplot(dt_plot, aes(x = mid_pos / 1e6, y = emd_ld)) +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf, ymin = 0.7, ymax = Inf,
      fill = "firebrick3", alpha = 0.05
    ) +

    # Dessine l'escalier uniquement là où il y a au moins 2 points
    geom_step(
      aes(group = interaction(pair, chrom_num)),
      linewidth = 0.4, color = "grey40", alpha = 0.6
    ) +

    # On ajoute les points de TOUTES les données (même isolées)
    geom_point(data = dt_emd, aes(color = emd_ld), size = 1.2) +

    facet_grid(pair ~ chrom_num, scales = "free_x", space = "free_x") +

    scale_color_gradientn(
      colors = c("darkturquoise", "grey90", "orange", "firebrick3"),
      values = c(0, 0.3, 0.7, 1),
      name = "Divergence\n(EMD-LD)"
    ) +

    geom_hline(yintercept = 0, linewidth = 0.3, color = "darkturquoise") +

    theme_minimal() +
    theme(
      strip.text.y = element_text(angle = 0, size = 7, face = "bold"),
      strip.text.x = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.1, "lines"),
      axis.text.x = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = sprintf("Profil de divergence génomique – EMD LD (%s)", window_label),
      subtitle = "Consensus local (bleu) vs divergence structurelle/haplotypique (rouge)",
      x = "Position génomique (Mb)",
      y = "Distance de Wasserstein (LD-based)",
      caption = "L'EMD-LD mesure le coût de transport basé sur 1-r². Les segments représentent la continuité des fenêtres."
    )

  return(p)
}
