## 6ml.R — 机器学习靶点筛选（SVM-RFE + LASSO + RF，5-fold CV）
rm(list = ls())
dir.create("output", showWarnings = FALSE)

library(glmnet)
library(randomForest)
library(caret)
library(kernlab)
library(TCMDATA)
library(ggplot2)
library(dplyr)
library(patchwork)
library(reshape2)

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

## ---- 1. 加载数据 ----
vsd_mat  <- readRDS("output/vsd.rds")
meta     <- readRDS("output/sample_meta.rds")
rank_res <- readRDS("output/ppi_rank.rds")
rank_df  <- rank_res$table

top20 <- head(rank_df$name, 20)
cat("Top 20 PPI targets:", paste(top20, collapse = ", "), "\n")

## ---- 2. 准备特征矩阵 ----
manual_remove <- c("GSM4683591", "GSM4683603")
keep <- setdiff(colnames(vsd_mat), manual_remove)
vsd_mat <- vsd_mat[, keep, drop = FALSE]
meta    <- meta[keep, , drop = FALSE]

top20 <- intersect(top20, rownames(vsd_mat))
X <- t(vsd_mat[top20, , drop = FALSE])
X <- scale(X)
y <- factor(meta$condition, levels = c("Hlty", "Seps_P"))

cat("Features:", ncol(X), " Samples:", nrow(X), "\n")

## ---- 3. SVM-RFE (5-fold × 10 repeats CV) ----
set.seed(42)
svm_funcs <- caretFuncs
svm_funcs$selectSize <- function(x, metric, maximize) {
  caret::pickSizeTolerance(x, metric = metric, tol = 1.5, maximize = maximize)
}
ctrl_svm <- rfeControl(functions = svm_funcs, method = "repeatedcv",
                       number = 5, repeats = 10)
svm_res  <- rfe(X, y, sizes = seq_len(ncol(X)), rfeControl = ctrl_svm, method = "svmLinear")
svm_genes <- predictors(svm_res)
cat("SVM-RFE optimal:", svm_res$optsize, "genes →", paste(svm_genes, collapse = ", "), "\n")

## ---- 4. LASSO (5-fold CV) ----
set.seed(42)
lambda_grid <- exp(seq(log(10), log(1e-4), length.out = 200))
lasso_fit <- glmnet(X, as.numeric(y == "Seps_P"),
              alpha = 1, family = "binomial",
              lambda = lambda_grid,
              standardize = FALSE)
cv_lasso  <- cv.glmnet(X, as.numeric(y == "Seps_P"),
                       alpha = 1, family = "binomial",
                lambda = lambda_grid,
                standardize = FALSE,
                nfolds = 5, type.measure = "deviance")
lasso_coef  <- coef(cv_lasso, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
lasso_genes <- lasso_genes[lasso_genes != "(Intercept)"]
cat("LASSO:", length(lasso_genes), "genes →", paste(lasso_genes, collapse = ", "), "\n")

## ---- 5. RF-RFE (5-fold × 10 repeats CV) ----
set.seed(42)
rf_funcs <- rfFuncs
rf_funcs$selectSize <- function(x, metric, maximize) {
  caret::pickSizeTolerance(x, metric = metric, tol = 1.5, maximize = maximize)
}
ctrl_rf <- rfeControl(functions = rf_funcs, method = "repeatedcv",
                      number = 5, repeats = 10)
rf_res  <- rfe(X, y, sizes = seq_len(ncol(X)), rfeControl = ctrl_rf)
rf_genes <- predictors(rf_res)
cat("RF optimal:", rf_res$optsize, "genes →", paste(rf_genes, collapse = ", "), "\n")

## ---- 6. 交集 → 关键靶点 ----
key_targets <- Reduce(intersect, list(svm_genes, lasso_genes, rf_genes))
cat("\n=== Key targets (SVM-RFE ∩ LASSO ∩ RF):", length(key_targets), "===\n")
cat(paste(key_targets, collapse = ", "), "\n")

## ---- 7. 提取三种算法的重要性 ----
svm_fit <- train(
  x = as.data.frame(X[, svm_genes, drop = FALSE]),
  y = y,
  method = "svmLinear",
  trControl = trainControl(method = "cv", number = 5)
)
svm_imp <- varImp(svm_fit, scale = FALSE)$importance
svm_imp$gene <- rownames(svm_imp)
numeric_cols <- vapply(svm_imp, is.numeric, logical(1))
numeric_cols[match("gene", names(numeric_cols), nomatch = 0)] <- FALSE
if ("Overall" %in% colnames(svm_imp) && is.numeric(svm_imp$Overall)) {
  svm_imp$importance <- svm_imp$Overall
} else {
  svm_imp$importance <- rowMeans(as.matrix(svm_imp[, numeric_cols, drop = FALSE]))
}
svm_imp <- svm_imp |>
  dplyr::select(.data$gene, .data$importance) |>
  dplyr::arrange(dplyr::desc(abs(.data$importance))) |>
  dplyr::mutate(
    relative_importance = .data$importance / max(.data$importance) * 100,
    gene = factor(.data$gene, levels = rev(.data$gene))
  )

svm_lower <- 0.98
svm_upper <- max(svm_imp$importance, na.rm = TRUE)
svm_upper <- max(svm_lower, ceiling(svm_upper / 0.005) * 0.005)
svm_breaks <- seq(svm_lower, svm_upper, by = 0.005)
svm_imp <- svm_imp |>
  dplyr::mutate(plot_importance = .data$importance - svm_lower)

lasso_bar_df <- data.frame(
  gene = lasso_genes,
  coefficient = as.numeric(lasso_coef[lasso_genes, 1]),
  stringsAsFactors = FALSE
) |>
  dplyr::arrange(dplyr::desc(abs(.data$coefficient))) |>
  dplyr::mutate(
    abs_coefficient = abs(.data$coefficient),
    gene = factor(.data$gene, levels = rev(.data$gene))
  )

rf_fit <- randomForest(
  x = X[, rf_genes, drop = FALSE],
  y = y,
  ntree = 500,
  importance = TRUE
)
rf_imp <- data.frame(
  gene = rownames(importance(rf_fit, type = 2)),
  importance = as.numeric(importance(rf_fit, type = 2)[, 1]),
  stringsAsFactors = FALSE
) |>
  dplyr::arrange(dplyr::desc(.data$importance)) |>
  dplyr::mutate(
    relative_importance = .data$importance / max(.data$importance) * 100,
    gene = factor(.data$gene, levels = rev(.data$gene))
  )

## ---- 8. Fig 6A: SVM-RFE 性能 + 相对重要性 ----
svm_plot_df <- svm_res$results
p6a1 <- ggplot(svm_plot_df, aes(Variables, Accuracy)) +
  geom_ribbon(aes(ymin = Accuracy - AccuracySD, ymax = Accuracy + AccuracySD),
              fill = "#E64B35", alpha = 0.12) +
  geom_line(color = "#E64B35", linewidth = 0.8) +
  geom_point(color = "#E64B35", size = 2) +
  geom_vline(xintercept = svm_res$optsize, linetype = "dashed", color = "grey40") +
  annotate("text", x = svm_res$optsize + 0.3, y = min(svm_plot_df$Accuracy),
           label = paste0("n = ", svm_res$optsize), hjust = 0, size = 3, color = "grey30") +
  scale_x_continuous(breaks = seq_len(ncol(X))) +
  labs(x = "Number of features", y = "Accuracy (5-fold CV)") +
  theme_sci()

p6a2 <- ggplot(svm_imp, aes(x = .data$gene, y = .data$plot_importance, fill = .data$importance)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, svm_upper - svm_lower),
    breaks = svm_breaks - svm_lower,
    labels = function(v) sprintf("%.3f", v + svm_lower),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_gradient(low = "#F8C3BC", high = "#C4372A", name = "SVM importance") +
  labs(x = NULL, y = "SVM importance") +
  theme_sci() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.title = element_text(face = "bold"),
        legend.position = "right")

p6a <- plot_list(p6a1, p6a2)
p6a

## ---- 9. Fig 6B: LASSO 路径 + 系数柱状图 ----
coef_mat <- as.matrix(lasso_fit$beta)
coef_df  <- reshape2::melt(coef_mat)
colnames(coef_df) <- c("Gene", "Step", "Coefficient")
coef_df$log_lambda <- rep(log(lasso_fit$lambda), each = nrow(coef_mat))

p6b1 <- ggplot(coef_df, aes(log_lambda, Coefficient, group = Gene, color = Gene)) +
  geom_line(linewidth = 0.6, show.legend = FALSE) +
  geom_vline(xintercept = log(cv_lasso$lambda.min), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = log(cv_lasso$lambda.1se), linetype = "dotted", color = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  labs(x = expression(log(lambda)), y = "Coefficient") +
  theme_sci()

p6b2 <- ggplot(lasso_bar_df, aes(x = .data$gene, y = .data$coefficient, fill = .data$coefficient)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  coord_flip() +
  scale_fill_gradient(low = "#CBD7F0", high = "#3C5488", name = "LASSO coefficient") +
  labs(x = NULL, y = "Absolute coefficient") +
  theme_sci() +
  theme(axis.text.y = element_text(face = "italic"),
  legend.title = element_text(face = "bold"),
        legend.position = "right")

p6b <- plot_list(p6b1, p6b2)
p6b

## ---- 10. Fig 6C: RF-RFE 性能 + 相对重要性 ----
rf_plot_df <- rf_res$results
p6c1 <- ggplot(rf_plot_df, aes(Variables, Accuracy)) +
  geom_ribbon(aes(ymin = Accuracy - AccuracySD, ymax = Accuracy + AccuracySD),
              fill = "#00A087", alpha = 0.12) +
  geom_line(color = "#00A087", linewidth = 0.8) +
  geom_point(color = "#00A087", size = 2) +
  geom_vline(xintercept = rf_res$optsize, linetype = "dashed", color = "grey40") +
  annotate("text", x = rf_res$optsize + 0.3, y = min(rf_plot_df$Accuracy),
           label = paste0("n = ", rf_res$optsize), hjust = 0, size = 3, color = "grey30") +
  scale_x_continuous(breaks = seq_len(ncol(X))) +
  labs(x = "Number of features", y = "Accuracy (5-fold CV)") +
  theme_sci()

p6c2 <- ggplot(rf_imp, aes(x = .data$gene, y = .data$relative_importance, fill = .data$relative_importance)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  coord_flip() +
  scale_fill_gradient(low = "#B7E3D7", high = "#00896B", name = "RF importance (%)") +
  labs(x = NULL, y = "Relative importance (%)") +
  theme_sci() +
  theme(axis.text.y = element_text(face = "italic"),
  legend.title = element_text(face = "bold"),
        legend.position = "right")

p6c <- plot_list(p6c1, p6c2)
p6c

## ---- 11. Fig 6D: Venn diagram ----
venn_ml <- getvenndata(svm_genes, lasso_genes, rf_genes,
                       set_names = c("SVM-RFE", "LASSO", "RF"))
p6d <- ggvenn_plot(venn_ml,
                   set.color = c("#E64B35", "#3C5488", "#00A087"),
                   stroke.size = 0.5, text.size = 4) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

p6d

ggsave("output/fig6d_ml_venn.pdf", p6d,
       width = 120, height = 110, units = "mm", dpi = 600, device = cairo_pdf)

## ---- 12. p6abc 输出 ----
ggsave("output/fig6asvm.pdf", p6a,
  width = 330, height = 200, units = "mm", dpi = 600, device = cairo_pdf)
ggsave("output/fig6blasso.pdf", p6b,
       width = 330, height = 200, units = "mm", dpi = 600, device = cairo_pdf)
ggsave("output/fig6crf.pdf", p6c,
       width = 330, height = 200, units = "mm", dpi = 600, device = cairo_pdf)

## ---- 13. 保存 ----
saveRDS(key_targets,  "output/key_targets.rds")
saveRDS(svm_genes,    "output/svm_genes.rds")
saveRDS(lasso_genes,  "output/lasso_genes.rds")
saveRDS(rf_genes,     "output/rf_genes.rds")
saveRDS(top20,        "output/top20_ppi_genes.rds")

cat("\nDone. Key targets:", paste(key_targets, collapse = ", "), "\n")
