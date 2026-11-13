# PANX1 knockout and spironolactone in A375-P melanoma cells
# Data: GEO GSE348006 (Partek-normalized counts, 16 samples)
# Design: 2x2 (SCR/KO x Ethanol/SPIR). 4 clones per genotype, and each clone
#         was treated with both Ethanol and SPIR (paired samples).
# Run from the RStudio Project folder; input files go in data/.

# ---- 0. One-time setup (run once by hand, then leave commented) -------------
# install.packages(c("tidyverse", "ggrepel", "BiocManager"))
# BiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))

# ---- 1. Libraries and folders -----------------------------------------------
library(tidyverse)
library(ggrepel)
library(limma)
library(clusterProfiler)   # load after tidyverse; it masks select(), so we
library(org.Hs.eg.db)      # call dplyr::select() explicitly below
library(enrichplot)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---- 2. Load data -----------------------------------------------------------
counts_raw <- read_tsv("data/GSE348006_Partek_Penuela_16_Normalization_Normalized_counts.txt.gz")

# Columns 12-27 are the 16 samples. Note: the "Gene Symbol" column actually
# holds Ensembl IDs; the real symbols are in gene_name.
sample_cols <- colnames(counts_raw)[12:27]

expr <- counts_raw %>%
  dplyr::select(all_of(sample_cols)) %>%
  as.matrix()
rownames(expr) <- counts_raw$gene_id

gene_info <- counts_raw %>% dplyr::select(gene_id, gene_name)

# ---- 3. Sample metadata (built from the sample names, never by position) ----
metadata <- tibble(sample = sample_cols) %>%
  mutate(
    genotype  = ifelse(str_detect(sample, "^KO"), "KO", "SCR"),
    treatment = ifelse(str_detect(sample, "SPIR"), "SPIR", "Ethanol"),
    clone     = str_extract(sample, "^(KO|SCR)[0-9]")
  ) %>%
  mutate(
    genotype  = factor(genotype,  levels = c("SCR", "KO")),        # SCR = baseline
    treatment = factor(treatment, levels = c("Ethanol", "SPIR"))   # Ethanol = baseline
  )

print(table(metadata$genotype, metadata$treatment))   # expect 4 in every cell
print(table(metadata$clone))                          # expect 2 per clone
stopifnot(all(colnames(expr) == metadata$sample))     # labels match the data

# ---- 4. Filter low-expressed genes, then log-transform -----------------------
# Keep genes with normalized count >= 10 in at least 4 samples (one group's size)
keep <- rowSums(expr >= 10) >= 4
cat("Genes kept:", sum(keep), "of", nrow(expr), "\n")
log_expr <- log2(expr[keep, ] + 1)

# ---- 5. PCA (quality control) -----------------------------------------------
pca <- prcomp(t(log_expr), scale. = TRUE)
pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

pca_df <- as_tibble(pca$x[, 1:2]) %>% bind_cols(metadata)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = treatment)) +
  geom_line(aes(group = clone), alpha = 0.5) +        # joins the two samples of each clone
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  labs(x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)")) +
  theme_minimal()
p_pca
ggsave("figures/pca_clones.png", p_pca, width = 8, height = 6, dpi = 300)

# ---- 6. Model: genotype x treatment, accounting for clone pairing -----------
design <- model.matrix(~ genotype * treatment, data = metadata)
colnames(design)
# genotypeKO               = KO vs SCR (in Ethanol)
# treatmentSPIR            = SPIR vs Ethanol (in SCR)
# genotypeKO:treatmentSPIR = does SPIR act differently in KO? (interaction)

corfit <- duplicateCorrelation(log_expr, design, block = metadata$clone)
corfit$consensus.correlation

fit <- lmFit(log_expr, design,
             block = metadata$clone,
             correlation = corfit$consensus.correlation)
fit <- eBayes(fit, trend = TRUE)

get_res <- function(coef_name) {
  topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_info, by = "gene_id")
}

res_ko    <- get_res("genotypeKO")
res_spir  <- get_res("treatmentSPIR")
res_inter <- get_res("genotypeKO:treatmentSPIR")

# Number of genes with adjusted p < 0.05 in each comparison
print(sapply(list(KO = res_ko, SPIR = res_spir, interaction = res_inter),
       function(r) sum(r$adj.P.Val < 0.05)))

head(res_ko, 10)
res_ko %>% filter(gene_name == "PANX1")   # PANX1 mRNA itself is not expected to drop much

write_csv(res_ko,    "results/limma_KO_vs_SCR.csv")
write_csv(res_spir,  "results/limma_SPIR_vs_Ethanol.csv")
write_csv(res_inter, "results/limma_interaction.csv")

# ---- 7. Gene set enrichment (GSEA on GO biological processes) ---------------
make_ranks <- function(res) {
  v <- res$t                      # rank by t statistic (size and reliability of change)
  names(v) <- res$gene_id
  sort(v, decreasing = TRUE)
}

run_gsea <- function(res) {
  set.seed(1)
  gseGO(geneList = make_ranks(res),
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        minGSSize = 15,
        maxGSSize = 500,
        pvalueCutoff = 0.05)
}

cell_cycle_hits <- function(gsea) {
  as.data.frame(gsea) %>%
    dplyr::filter(str_detect(Description,
                             "cell cycle|mitotic|DNA replication|chromosome segregation")) %>%
    dplyr::select(Description, NES, p.adjust) %>%
    arrange(NES)
}

gsea_ko    <- run_gsea(res_ko)
gsea_spir  <- run_gsea(res_spir)
gsea_inter <- run_gsea(res_inter)

nrow(as.data.frame(gsea_ko))
nrow(as.data.frame(gsea_spir))
nrow(as.data.frame(gsea_inter))

# Cell-cycle terms: negative NES = lower in KO / lower with SPIR / lower in KO+SPIR
cc_ko    <- cell_cycle_hits(gsea_ko);    print(cc_ko)
cc_spir  <- cell_cycle_hits(gsea_spir);  print(cc_spir)
cc_inter <- cell_cycle_hits(gsea_inter); print(cc_inter)

write_csv(as.data.frame(gsea_ko),    "results/gsea_KO_vs_SCR.csv")
write_csv(as.data.frame(gsea_spir),  "results/gsea_SPIR_vs_Ethanol.csv")
write_csv(as.data.frame(gsea_inter), "results/gsea_interaction.csv")

# ---- 8. Figures -------------------------------------------------------------
# Merge near-duplicate GO terms so the plot is readable
gsea_simple <- simplify(gsea_ko, cutoff = 0.7)
nrow(as.data.frame(gsea_simple))

p_dot <- dotplot(gsea_simple, showCategory = 8, split = ".sign") +
  facet_grid(. ~ .sign) +
  theme(axis.text.y = element_text(size = 8))
p_dot
ggsave("figures/gsea_ko_dotplot.png", p_dot, width = 10, height = 8, dpi = 300)

# Classic GSEA curve for one cell-cycle term
p_curve <- gseaplot2(gsea_ko, geneSetID = "GO:0007076",
                     title = "Mitotic chromosome condensation")
p_curve
ggsave("figures/gsea_curve_mitotic_condensation.png", p_curve,
       width = 8, height = 6, dpi = 300)

# ---- 9. OPTIONAL: compare with the authors' own results ---------------------
# Not yet tested. Put GSE348006_gene_list.txt.gz in data/ first.
# genotypeKO here = KO vs SCR in Ethanol, matching the authors' "KO,Ethanol vs SCR,Ethanol".
# authors <- read_tsv("data/GSE348006_gene_list.txt.gz")
# authors_ko <- authors %>%
#   transmute(gene_id     = `Gene ID`,
#             fdr_authors = as.numeric(`FDR step up (KO,Ethanol vs SCR,Ethanol)`),
#             fc_authors  = as.numeric(`Fold change (KO,Ethanol vs SCR,Ethanol)`))
# cmp <- inner_join(res_ko, authors_ko, by = "gene_id")
# cor(cmp$logFC, sign(cmp$fc_authors) * log2(abs(cmp$fc_authors)),
#     method = "spearman", use = "complete.obs")
# ours    <- res_ko$gene_id[res_ko$adj.P.Val < 0.05]
# theirs  <- authors_ko$gene_id[which(authors_ko$fdr_authors < 0.05)]
# length(intersect(ours, theirs)); length(ours); length(theirs)

# ---- 10. Record package versions for reproducibility ------------------------
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
