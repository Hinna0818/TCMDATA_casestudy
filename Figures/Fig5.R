rm(list = ls())
library(TCMDATA)
library(ggplot2)
library(patchwork)

# 1. 准备数据covid19
data("covid19", package="TCMDATA")
expr_mat <- as.matrix(covid19$expr)
group <- covid19$group_info$group

# 提取方差top 50的基因作为测试靶点池
gene_var <- apply(expr_mat, 1, var)
genes <- names(sort(gene_var, decreasing = TRUE))[1:50]

ml_data <- prepare_ml_data(expr_mat = expr_mat, 
                           group = group, 
                           positive_class = "ICU", 
                           genes = genes)

# 2. 独立运行模型，用于提取算法各自特定的特征重要性
set.seed(2025)
res_lasso <- ml_lasso(ml_data, cv_folds = 10, seed = 2025)
res_rf <- ml_rf(ml_data, n_trees = 500, seed = 2025)
res_svm <- ml_svm_rfe(ml_data, cv_folds = 5, seed = 2025)
res_xgb <- ml_xgboost(ml_data, cv_folds = 5, seed = 2025)

# 构建 tcm_ml_list 对象以支持 plot_ml_roc
ml_list <- list(lasso = res_lasso, rf = res_rf, svm = res_svm, xgboost = res_xgb)
class(ml_list) <- "tcm_ml_list" 


c_red <- "#E64B35FF"
c_blue <- "#4DBBD5FF" 
c_green <- "#00A087FF"
c_purple <- "#3C5488FF"
c_orange <- "#F39B7FFF"

my_theme <- theme_classic() +
  theme(
    plot.title = element_blank(),           # 去掉标题
    plot.subtitle = element_blank(),        # 去掉副标题(如Mode A)
    axis.text.y = element_text(face = "italic", size = 10, color = "black"), # 基因斜体
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "plain"))

# 3. 构建 8 个漂亮的子图 (A-H)

# A: LASSO CV (误差演化图) - Base R plot, 用 wrap_elements 包装以支持 patchwork 排版
pA <- wrap_elements(panel = ~{
  par(mar = c(4, 4, 1, 1)) # 设置较小的边距
  plot_enet_cv(res_lasso)
})
pA

# B: LASSO Path (系数路径图) - Base R plot, 用 wrap_elements 包装
pB <- wrap_elements(panel = ~{
  par(mar = c(4, 4, 1, 1))
  plot_enet_path(res_lasso)
})
pB

# C: LASSO Coefs (特征重要性条形图)
pC <- plot_enet_coefs(res_lasso) + 
  my_theme +
  # 兼容可能的列名大小写
  scale_fill_manual(values = c("Positive" = c_red, "Negative" = c_blue, "positive" = c_red, "negative" = c_blue)) +
  labs(x = "LASSO Coefficient", y = NULL)
pC

# D: RF Boruta 重要性分析图 (使用包内自带函数绘制Boruta)
pD <- plot_rf_boruta(res_rf, top_n = 15) + 
  my_theme  +  
  theme(legend.position = "none") +
  labs(x = "Boruta Importance", y = NULL)
pD

# E: SVM Importance (特征重量条形图)
df_svm <- res_svm$importance
df_svm <- df_svm[order(df_svm$importance, decreasing = TRUE), ]
if(nrow(df_svm) > 25) df_svm <- df_svm[1:15, ]
pE <- ggplot(df_svm, aes(x = importance, y = reorder(gene, importance))) +
  geom_col(fill = c_purple, width = 0.7) +
  my_theme +
  labs(x = "Feature Weight (SVM-RFE)", y = NULL)
pE

# F: XGBoost Importance (信息增益条形图)
df_xgb <- res_xgb$importance
df_xgb <- df_xgb[order(df_xgb$importance, decreasing = TRUE), ]
if(nrow(df_xgb) > 15) df_xgb <- df_xgb[1:15, ]
pF <- ggplot(df_xgb, aes(x = importance, y = reorder(gene, importance))) +
  geom_col(fill = c_orange, width = 0.7) +
  my_theme +
  labs(x = "Information Gain (XGBoost)", y = NULL)
pF

# G: ML Consistency Venn (原生的 ggvenn 控制配色)
library(ggvenn)
venn_list <- list(
  LASSO = ml_list$lasso$genes,
  "Random Forest" = ml_list$rf$genes,
  "SVM-RFE" = ml_list$svm$genes,
  XGBoost = ml_list$xgboost$genes
)
pG <- ggvenn(
  venn_list, 
  fill_color = c("#4DBBD599", "#E64B3599", "#00A08799", "#F39B7F99"), 
  stroke_size = 0.5, 
  set_name_size = 4, 
  text_size = 4
) + theme(plot.title = element_blank(), plot.subtitle = element_blank())
pG

# H: ROC Comparisons
pH <- plot_ml_roc(ml_list) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "plain"),
    legend.position = c(0.7, 0.25),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 9, face = "plain")
  ) +
  scale_color_manual(values = c(c_red, c_blue, c_green, c_purple, c_orange, "black"))
pH

# 4. 创建输出目录
out_dir <- "output/figure_ml"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Saving individual PDF plots...")
ggsave(file.path(out_dir, "pA_lasso_cv.pdf"), pA)
ggsave(file.path(out_dir, "pB_lasso_path.pdf"), pB)
ggsave(file.path(out_dir, "pC_lasso_imp.pdf"), pC)
ggsave(file.path(out_dir, "pD_rf_boruta.pdf"), pD)
ggsave(file.path(out_dir, "pE_svm_imp.pdf"), pE)
ggsave(file.path(out_dir, "pF_xgb_imp.pdf"), pF)
ggsave(file.path(out_dir, "pG_venn.pdf"), pG)
ggsave(file.path(out_dir, "pH_roc.pdf"), pH)

