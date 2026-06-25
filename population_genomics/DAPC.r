# ============================================================
# DAPC using metadata Source / Source1 as grouping variable
# ============================================================

# install.packages(c("adegenet", "tidyverse", "ggrepel"))
library(adegenet)
library(tidyverse)
library(ggrepel)

# ------------------------------------------------------------
# 1. Read gene / feature matrix
# ------------------------------------------------------------
# Assumption:
# rows = genes/features
# columns = isolates/samples
# first column = gene ID

vf <- read.delim(
  "data/90vf_count.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

vf_mat <- vf %>%
  column_to_rownames(colnames(vf)[1]) %>%
  as.matrix()

mode(vf_mat) <- "numeric"

# DAPC requires rows = samples, columns = genes/features
X <- t(vf_mat)

# ------------------------------------------------------------
# 2. Read metadata
# ------------------------------------------------------------
# Your metadata should contain:
# tip = sample name
# Source or Source1 = group

meta <- read.delim(
  "data/phylogeny/metadata.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# If Source1 does not exist, use Source as Source1
if (!"Source1" %in% colnames(meta)) {
  if ("Source" %in% colnames(meta)) {
    meta$Source1 <- meta$Source
  } else {
    stop("No column named Source1 or Source found in metadata.")
  }
}

# Check required sample ID column
if (!"tip" %in% colnames(meta)) {
  stop("Metadata must contain a column named 'tip' with sample IDs.")
}

# ------------------------------------------------------------
# 3. Match metadata with feature matrix
# ------------------------------------------------------------

common_samples <- intersect(rownames(X), meta$tip)

X <- X[common_samples, , drop = FALSE]

meta2 <- meta %>%
  filter(tip %in% common_samples) %>%
  arrange(match(tip, rownames(X)))

# Confirm matching order
stopifnot(all(meta2$tip == rownames(X)))

# ------------------------------------------------------------
# 4. Prepare grouping variable: Source1
# ------------------------------------------------------------

grp <- factor(meta2$Source1)

table(grp)

# Optional: remove groups with fewer than 2 isolates
# This is useful because DAPC cross-validation can fail with tiny groups.

grp_counts <- table(grp)
keep_groups <- names(grp_counts[grp_counts >= 2])

keep_samples <- grp %in% keep_groups

X <- X[keep_samples, , drop = FALSE]
meta2 <- meta2[keep_samples, ]
grp <- droplevels(factor(meta2$Source1))

table(grp)

# ------------------------------------------------------------
# 5. Remove invariant genes/features
# ------------------------------------------------------------

keep_features <- apply(X, 2, function(z) sd(z, na.rm = TRUE) > 0)

X <- X[, keep_features, drop = FALSE]

# Optional: convert copy number to presence/absence
# Recommended for gene presence/absence data.

X_use <- ifelse(X > 0, 1, 0)

# ------------------------------------------------------------
# 6. Choose number of PCs by cross-validation
# ------------------------------------------------------------

set.seed(123)

max_pc <- min(50, nrow(X_use) - 1, ncol(X_use))

xval <- xvalDapc(
  x = X_use,
  grp = grp,
  n.pca.max = max_pc,
  training.set = 0.8,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.rep = 30,
  xval.plot = TRUE
)

best_pca <- xval$`Number of PCs Achieving Highest Mean Success`

best_pca

# ------------------------------------------------------------
# Set maximum possible PCs
# ------------------------------------------------------------

max_possible_pcs <- min(
  nrow(X_use) - 1,
  ncol(X_use)
)

max_possible_pcs

# ------------------------------------------------------------
# 7. Run DAPC using Source1 as group
# ------------------------------------------------------------

xval <- xvalDapc(
  x = X_use,
  grp = grp,
  n.pca.max = max_possible_pcs,
  training.set = 0.8,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.rep = 30,
  xval.plot = TRUE
)

str(xval)

best_pca <- xval$`Number of PCs Achieving Highest Mean Success`

best_pca

best_pca <- as.numeric(best_pca[1])

if (is.na(best_pca) || best_pca < 1) {
  best_pca <- min(10, max_possible_pcs)
}

if (best_pca > max_possible_pcs) {
  best_pca <- max_possible_pcs
}

best_pca

n_da_use <- min(
  length(levels(grp)) - 1,
  nrow(X_use) - length(levels(grp))
)

n_da_use

dapc_res <- dapc(
  X_use,
  grp = grp,
  n.pca = best_pca,
  n.da = n_da_use
)


scatter(
  dapc_res,
  scree.da = TRUE,
  bg = "white",
  pch = 20,
  cell = 0,
  cstar = 0,
  solid = 0.7,
  clabel = 0,
  legend = TRUE
)

library(ggplot2)
library(ggrepel)

dapc_scores <- as.data.frame(dapc_res$ind.coord)

dapc_scores$tip <- rownames(dapc_scores)
dapc_scores$Source1 <- grp

ggplot(dapc_scores, aes(x = LD1, y = LD2, color = Source1)) +
  geom_point(size = 3, alpha = 0.85) +
  theme_bw() +
  labs(
    x = "Discriminant function 1",
    y = "Discriminant function 2",
    color = "Source1"
  )

dapc_res <- dapc(
  X_use,
  grp = grp,
  n.pca = best_pca,
  n.da = length(levels(grp)) - 1
)

# ------------------------------------------------------------
# 8. Basic DAPC plot
# ------------------------------------------------------------

scatter(
  dapc_res,
  scree.da = TRUE,
  bg = "white",
  pch = 20,
  cell = 0,
  cstar = 0,
  solid = 0.7,
  clabel = 0,
  legend = TRUE
)

# ------------------------------------------------------------
# 9. ggplot DAPC plot
# ------------------------------------------------------------

dapc_scores <- as.data.frame(dapc_res$ind.coord)

dapc_scores$tip <- rownames(dapc_scores)
dapc_scores$Source1 <- grp

ggplot(dapc_scores, aes(x = LD1, y = LD2, color = Source1)) +
  geom_point(size = 3, alpha = 0.85) +
  theme_bw() +
  labs(
    x = "Discriminant function 1",
    y = "Discriminant function 2",
    color = "Source1"
  )

# With sample labels

ggplot(dapc_scores, aes(x = LD1, y = LD2, color = Source1, label = tip)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size = 2.5, max.overlaps = 30) +
  theme_bw() +
  labs(
    x = "Discriminant function 1",
    y = "Discriminant function 2",
    color = "Source1"
  )

# ------------------------------------------------------------
# 10. Save DAPC coordinates
# ------------------------------------------------------------

write.csv(
  dapc_scores,
  "DAPC_scores_by_Source1.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Posterior membership probabilities
# ------------------------------------------------------------

posterior <- as.data.frame(dapc_res$posterior)

posterior$tip <- rownames(posterior)
posterior$Source1 <- grp

write.csv(
  posterior,
  "DAPC_posterior_probabilities_by_Source1.csv",
  row.names = FALSE
)

# Plot membership probabilities

compoplot(
  dapc_res,
  show.lab = TRUE,
  legend = TRUE,
  lab = rownames(X_use)
)

# ------------------------------------------------------------
# 12. Identify genes/features contributing most to separation
# ------------------------------------------------------------

var_contrib <- as.data.frame(dapc_res$var.contr)

var_contrib$feature <- rownames(var_contrib)

write.csv(
  var_contrib,
  "DAPC_feature_contributions_by_Source1.csv",
  row.names = FALSE
)

# Top contributors to LD1

top_LD1 <- var_contrib %>%
  arrange(desc(abs(LD1))) %>%
  slice(1:30)

ggplot(top_LD1, aes(x = reorder(feature, abs(LD1)), y = abs(LD1))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Gene / feature",
    y = "Absolute contribution to LD1"
  )

write.csv(
  top_LD1,
  "DAPC_top_LD1_features_by_Source1.csv",
  row.names = FALSE
)
