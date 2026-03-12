## 7validation.R — 关键靶点验证（PPI Knock + 表达箱线图 + ROC）
rm(list = ls())
dir.create("output", showWarnings = FALSE)

library(TCMDATA)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pROC)
library(patchwork)
library(ggpubr)
library(grid)
library(caret)

theme_sci <- function(base_size = 9) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text  = element_text(color = "black", size = base_size),
      axis.line  = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35)
    )
}

pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
         "#8491B4", "#91D1C2", "#DC9A6C", "#7E6148", "#B09C85")

## ---- 加载数据 ----
key_targets <- readRDS("output/key_targets.rds")
ppi_scored  <- readRDS("output/ppi_scored.rds")
vsd_mat     <- readRDS("output/vsd.rds")
meta        <- readRDS("output/sample_meta.rds")

manual_remove <- c("GSM4683591", "GSM4683603")
keep <- setdiff(colnames(vsd_mat), manual_remove)
vsd_mat <- vsd_mat[, keep, drop = FALSE]
meta    <- meta[keep, , drop = FALSE]

cat("Key targets:", paste(key_targets, collapse = ", "), "\n")

## ############################################################
## PART 1 — PPI Knock
## ############################################################

knock_res <- ppi_knock(ppi_scored, targets = key_targets, n_perm = 500, seed = 42)
knock_df  <- knock_res$Summary

p7a <- ggplot(knock_df, aes(x = Metric, y = Normalized_RI, fill = Metric)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed",
       color = "grey55", linewidth = 0.4) +
  scale_fill_manual(values = pal[1:4]) +
  labs(title = paste0("PPI knock (P = ", formatC(knock_res$Total_Pvalue, format = "e", digits = 2), ")"),
    x = NULL, y = "Normalized RI (Z-score)") +
  theme_sci() +
  theme(plot.title = element_text(size = 10),
     axis.text.x = element_text(face = "bold"))

## ############################################################
## PART 1 — 关键靶点表达箱线图
## ############################################################

expr_mat  <- t(vsd_mat[key_targets, , drop = FALSE])
expr_df   <- data.frame(group = meta$condition, expr_mat, check.names = FALSE)
expr_long <- pivot_longer(expr_df, cols = all_of(key_targets),
                          names_to = "Gene", values_to = "Expression")
expr_long$Gene <- factor(expr_long$Gene, levels = key_targets)

p7b <- ggplot(expr_long, aes(x = group, y = Expression, fill = group)) +
  geom_violin(width = 0.82, alpha = 0.18, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.36, outlier.shape = NA, alpha = 0.92,
               color = "grey25", linewidth = 0.35) +
  geom_jitter(aes(color = group), width = 0.10, size = 0.8, alpha = 0.42,
              show.legend = FALSE) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Hlty" = "#5B78B0", "Seps_P" = "#D85C4A")) +
  scale_color_manual(values = c("Hlty" = "#355C9A", "Seps_P" = "#BE3E2F")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.y.npc = 0.96, size = 4.2, bracket.size = 0.35) +
  labs(x = NULL, y = "Expression", fill = "Group") +
  theme_sci() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_rect(fill = "grey96", color = NA),
        panel.spacing = unit(8, "mm"))

p7b
ggsave("output/fig7a_boxplot.pdf", p7b,
       units = "mm", dpi = 600, device = cairo_pdf)

## ############################################################
## PART 2 — 单基因诊断 ROC 曲线
## ############################################################

folds <- createFolds(meta$condition, k = 5, returnTrain = FALSE)

get_cv_scores <- function(gene_name) {
  expr <- as.numeric(vsd_mat[gene_name, ])
  pred <- rep(NA_real_, length(expr))

  for (idx in folds) {
    train_idx <- setdiff(seq_along(expr), idx)
    train_df <- data.frame(
      y = meta$condition[train_idx],
      expr = expr[train_idx]
    )
    test_df <- data.frame(expr = expr[idx])

    fit <- glm(y ~ expr, data = train_df, family = binomial())
    pred[idx] <- predict(fit, newdata = test_df, type = "response")
  }

  pred
}

cv_scores <- lapply(key_targets, get_cv_scores)
names(cv_scores) <- key_targets

roc_list <- lapply(key_targets, function(g) {
  roc(response  = meta$condition,
      predictor = cv_scores[[g]],
      levels    = c("Hlty", "Seps_P"),
      quiet     = TRUE)
})
names(roc_list) <- key_targets

roc_df <- do.call(rbind, lapply(key_targets, function(g) {
  r <- roc_list[[g]]
  data.frame(
    Gene        = g,
    Label       = paste0(g, "\nAUC = ", sprintf("%.3f", auc(r))),
    Specificity = 1 - r$specificities,
    Sensitivity = r$sensitivities
  )
}))
roc_df$Gene <- factor(roc_df$Gene, levels = key_targets)
label_df <- roc_df |>
  dplyr::distinct(.data$Gene, .data$Label)

p7c <- ggplot(roc_df, aes(x = Specificity, y = Sensitivity, color = Gene)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey65", linewidth = 0.4) +
  scale_color_manual(values = pal[seq_along(key_targets)]) +
  coord_equal() +
  facet_wrap(~ Label, ncol = 2) +
  labs(x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  theme_sci() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold.italic", size = 9),
        strip.background = element_rect(fill = "grey96", color = NA),
        panel.spacing = unit(8, "mm"))

p7c
ggsave("output/fig7c_roc.pdf", p7c,
       width = 140, height = 135, units = "mm", dpi = 600, device = cairo_pdf)

for (g in key_targets) {
  cat(sprintf("  %s  AUC = %.3f\n", g, as.numeric(auc(roc_list[[g]]))))
}

## 找单体-靶点对应关系
lingzhi <- search_herb(herb = "lingzhi", type = "Herb_pinyin_name")
molecule_target <- subset(lingzhi, lingzhi$target %in% key_targets)



