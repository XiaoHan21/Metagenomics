### 20/08/2025
### XIAO Han Shawn
### Purpose:
### *** 1. Alpha/Beta diversity analysis 
###     2. Discard gender factor
###     3. 

# 确保已经安装并加载 Maaslin2
if(!requireNamespace("Maaslin2", quietly = TRUE)) BiocManager::install("Maaslin2")

pkgs <- c(
  "dplyr", "openxlsx", "phyloseq", "fossil", 
  "ggplot2", "vegan", "ape", "microbiome", 
  "ggpubr", "plyr", "scales", "tidyr",
  "optparse", "DESeq2", "data.table",
  "future", "future.apply","tidyverse","Maaslin2"
)

invisible(lapply(pkgs, function(pkg) {
  suppressMessages(require(pkg, character.only = TRUE))
}))

## 1. Load data
#-------------------------------------------------------------------------------
wkdir <- "/home/han_xiao/xiaohan/Others/Jiangbo/DepressionMouse/01.metagenomics/02.script/r03.Diveristy2Diff.ipynb"
indir <- "/home/han_xiao/xiaohan/Others/Jiangbo/DepressionMouse/01.metagenomics/03.result/05.metaphlan/"
outdir <- "/home/han_xiao/xiaohan/Others/Jiangbo/DepressionMouse/01.metagenomics/03.result"

# if(!dir.exists(wkdir)){dir.create(wkdir,recursive = T)}
if(!dir.exists(indir)){dir.create(indir,recursive = T)}
if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}

setwd(outdir)

if(!dir.exists("07.diversity")){dir.create("07.diversity")}
if(!dir.exists("08.diffAbund")){dir.create("08.diffAbund")}

abs.abund.fpath <- "/home/han_xiao/xiaohan/Others/Jiangbo/DepressionMouse/01.metagenomics/03.result/05.metaphlan/taxonomy.tsv"
abs.abund.df <- read_tsv(abs.abund.fpath, comment = "#")%>%
  filter(grepl("s__", ID) & !grepl("t__", ID)) %>%
  mutate(Species = sub(".*s__", "", ID)) %>%
  select(Species, everything(), -ID) %>%
  column_to_rownames("Species")
abs.abund.df <- abs.abund.df / 100
rel.abund.df <- abs.abund.df

# Prepare sample metadata

metadata.fpath <- "/home/han_xiao/xiaohan/Others/Jiangbo/DepressionMouse/01.metagenomics/01.data/metadata.xlsx"

metadata <- read.xlsx(metadata.fpath,rowNames=T)%>% select(intervention,Color)
colnames(metadata) <- c("type","color")     
metadata$sample <- rownames(metadata)
color.panel <- metadata %>% select(type,color) %>% distinct(type,color) %>% pull(color)
names(color.panel) <- metadata %>% select(type,color) %>% distinct(type,color) %>% pull(type)
                                    
# cut off of mean relative abundance
cutoff <- 0.0001

# comparision groups
groups <- names(color.panel)

# ---- 1. Set Biological Ordering of Factors (type_levels) ----
comp.groups <- list(
  # A. 与健康对照组 (WT control) 的对比：评估造模影响及治疗后是否恢复到基线水平
  c("aVNS 0.1mA 20Hz", "WT control"),
  c("aVNS 0.1mA 100Hz", "WT control"),
  c("aVNS 0.1mA 5Hz", "WT control"),
  c("aVNS 0.3mA 20Hz", "WT control"),
  c("Footshock stress", "WT control"), # 造模验证 (Model Validation)
  
  # B. 与模型组 (Footshock stress) 的对比：评估不同参数 aVNS 的治疗效应 (Treatment Efficacy)
  c("aVNS 0.1mA 20Hz", "Footshock stress"),
  c("aVNS 0.1mA 100Hz", "Footshock stress"),
  c("aVNS 0.1mA 5Hz", "Footshock stress"),
  c("aVNS 0.3mA 20Hz", "Footshock stress")
)

# ---- 3. Comparisons for Differential Analysis (e.g., DESeq2 / Limma) ----
# 注意方向：case (实验组/分子) 在前，control (对照组/分母) 在后。Fold Change = Case/Control
comparisons <- list(
  # A. 与健康对照组 (WT control) 的差异分析
  list(case = "aVNS 0.1mA 20Hz",   control = "WT control"),
  list(case = "aVNS 0.1mA 100Hz",  control = "WT control"),
  list(case = "aVNS 0.1mA 5Hz",    control = "WT control"),
  list(case = "aVNS 0.3mA 20Hz",   control = "WT control"),
  list(case = "Footshock stress",  control = "WT control"), # Stress 诱导的关键靶点
  
  # B. 与模型组 (Footshock stress) 的差异分析 (寻找被不同参数 aVNS 逆转/调节的基因)
  list(case = "aVNS 0.1mA 20Hz",   control = "Footshock stress"),
  list(case = "aVNS 0.1mA 100Hz",  control = "Footshock stress"),
  list(case = "aVNS 0.1mA 5Hz",    control = "Footshock stress"),
  list(case = "aVNS 0.3mA 20Hz",   control = "Footshock stress")
)

type_levels <- c('WT control','Footshock stress','aVNS 0.1mA 5Hz','aVNS 0.1mA 20Hz','aVNS 0.1mA 100Hz','aVNS 0.3mA 20Hz')

type_levels

head(rel.abund.df)

head(metadata)

setwd(outdir)
wkdir <- file.path(outdir,"08.diffAbund","00.Metadata")
if(!dir.exists(wkdir)){ dir.create(wkdir,recursive = T)}
setwd(wkdir)


#-------------------------------------------------------------------------------
# Function to filter species based on their mean relative abundance across samples
# Args:
#   abs.mat:      Matrix of absolute abundances (rows = species, columns = samples)
#   rel.mat:      Matrix of relative abundances (rows = species, columns = samples)
#   mean.rel.cutoff: Threshold value for filtering species by mean relative abundance (default: 0.001)
# Returns:
#   A list containing:
#     - filtered.abs: Filtered absolute abundance matrix
#     - filtered.rel: Filtered relative abundance matrix
#     - cutoff:       The mean relative abundance cutoff used
Filter_Species = function(abs.mat, rel.mat, mean.rel.cutoff = 0.0001) {
    
    # Calculate mean relative abundance for each species across all samples
    sp.mean.across.samples = rowMeans(rel.mat)
    
    # Get species names with mean relative abundance above the cutoff
    sp.filtered.samples = rownames(rel.mat)[sp.mean.across.samples > mean.rel.cutoff]
    
    # Check if any species passed the filter
    if (length(sp.filtered.samples) == 0) {
        stop("No species passed the mean.rel.cutoff filter")
    }
    
    # Subset the absolute and relative abundance matrices using filtered species
    # drop = FALSE preserves matrix structure when only 1 species remains
    filtered.abs.mat = abs.mat[sp.filtered.samples, , drop = FALSE]
    filtered.rel.mat = rel.mat[sp.filtered.samples, , drop = FALSE]
    
    # Return results as a list with cutoff value
    return(list(
        filtered.abs = filtered.abs.mat,
        filtered.rel = filtered.rel.mat,
        cutoff = mean.rel.cutoff
    ))
}

abs.mat <- abs.abund.df
rel.mat <- rel.abund.df

abs.mat.fil <- Filter_Species(abs.mat,rel.mat,mean.rel.cutoff = cutoff)[["filtered.abs"]]
rel.mat.fil <- Filter_Species(abs.mat,rel.mat,mean.rel.cutoff = cutoff)[["filtered.rel"]]

cat("Cut-off of mean relative abundace across all samples:", cutoff)

write.xlsx(abs.mat.fil, paste0("AbsoluteAbundance_Cutoff",cutoff,".xlsx"),colNames=T,rowNames=T)
write.xlsx(rel.mat.fil, paste0("RelativeAbundance_Cutoff",cutoff,".xlsx"),colNames=T,rowNames=T)

setwd(file.path(outdir,"08.diffAbund"))
if(!dir.exists("01.Diff")){dir.create("01.Diff")}
setwd("01.Diff")

#-------------------------------------------------------------------------------
# 使用 MaAsLin2 的差异分析核心函数
#-------------------------------------------------------------------------------
Caculate_DiffExp_MaAsLin2 <- function(rel.mat, metadata, 
                                      foldchange = 2, p.value = 0.05, p.adj = 0.1,
                                      case.type, control.type, outdir_base) {

    print(paste("Comparing:", case.type, "vs", control.type))
    
    # 1. 过滤并准备 metadata (行名必须是样本名)
    meta_sub <- metadata %>% filter(type %in% c(case.type, control.type))
    rownames(meta_sub) <- meta_sub$sample
    
    # 2. 过滤并转置相对丰度矩阵 (MaAsLin2要求：行是样本，列是物种)
    # 只取 meta_sub 中包含的样本
    rel_sub <- rel.mat[, meta_sub$sample] 
    rel_sub_t <- as.data.frame(t(rel_sub))
    
    # 3. 规范化输出文件夹名称 (将空格替换为下划线，避免路径报错)
    case_clean <- gsub(" ", "_", case.type)
    ctrl_clean <- gsub(" ", "_", control.type)
    comp_outdir <- file.path(outdir_base, paste0("Case-", case_clean, "_vs_Ctrl-", ctrl_clean))
    if(!dir.exists(comp_outdir)){dir.create(comp_outdir, recursive = TRUE)}
    
    # 4. 运行 MaAsLin2
    # 注意：因为输入已经是相对丰度，normalization 设为 "NONE"
    # transform 使用 "LOG" 是微生物相对丰度差异分析的标准配置
    fit_data <- Maaslin2(
        input_data     = rel_sub_t,
        input_metadata = meta_sub,
        output         = comp_outdir,
        fixed_effects  = "type",
        reference      = c("type", control.type), # 强制设定对照组作为基线
        normalization  = "NONE",
        transform      = "LOG", 
        min_abundance  = 0.0, # 假设你在外部已经做过过滤，这里不二次过滤
        min_prevalence = 0.0,
        plot_heatmap   = FALSE, # 批量跑关闭画图可加速计算
        plot_scatter   = FALSE
    )
    
    # 5. 提取 MaAsLin2 统计结果
    res <- fit_data$results
    
    # 6. 计算传统的均值和 Log2FC (保持你原来的代码习惯)
    sample_case <- meta_sub$sample[meta_sub$type == case.type]
    sample_ctrl <- meta_sub$sample[meta_sub$type == control.type]
    
    FC <- as.data.frame(matrix(nrow = nrow(rel_sub), ncol = 4, data = 0))
    colnames(FC) <- c("Species", paste0(case.type, ".AvgExp"), paste0(control.type, ".AvgExp"), "Log2FC")
    
    for (i in 1:nrow(rel_sub)) {
        sp <- rownames(rel_sub)[i]
        FC[i,1] <- sp
        FC[i,2] <- mean(as.numeric(rel_sub[i, sample_case]), na.rm=TRUE)
        FC[i,3] <- mean(as.numeric(rel_sub[i, sample_ctrl]), na.rm=TRUE)
        # 添加一个极小值避免 log2(0) 报错
        FC[i,4] <- log2((FC[i,2] + 1e-6) / (FC[i,3] + 1e-6))
    }
    
    # 由于 R 数据框的列名规则，物种名称中的特殊符号可能会被 Maaslin2 替换为小数点 '.'
    # 为了合并，我们使用 Maaslin2 内部的替换规则统一一下物种名
    FC$feature <- make.names(FC$Species)
    
    # 7. 合并均值、Log2FC与 MaAsLin2 的 p值/q值
    tmp2 <- merge(FC, res[, c("feature", "coef", "stderr", "pval", "qval")], by = "feature", all.x = TRUE)
    
    # 8. 打标签 (U = Up, D = Down)
    tmp2$type <- "_"
    # 注意 MaAsLin2 输出可能包含 NA，需要先处理
    tmp2$pval[is.na(tmp2$pval)] <- 1
    tmp2$qval[is.na(tmp2$qval)] <- 1
    
    tmp2[tmp2$Log2FC > log2(foldchange) & tmp2$pval < p.value & tmp2$qval < p.adj, "type"] <- "U"
    tmp2[tmp2$Log2FC < -log2(foldchange) & tmp2$pval < p.value & tmp2$qval < p.adj, "type"] <- "D"
    
    # 整理列顺序，去掉用于 merge 的中间变量
    tmp2 <- tmp2 %>% select(Species, everything(), -feature)
    
    # 9. 输出 Excel
    write.xlsx(tmp2, file.path(outdir_base, paste0("MaAsLin2_Case-", case_clean, ".vs.Ctrl-", ctrl_clean, ".xlsx")),
               colNames = TRUE, rowNames = FALSE)
    
    return(tmp2)
}



getwd()

#-------------------------------------------------------------------------------
# 并行执行批量对比
#-------------------------------------------------------------------------------
# 设置你想要保存差异分析表格的总目录
diffAbund_dir <- "." 

plan(multisession, workers = min(29, length(comparisons))) # workers数建议不超过比较组的数量

res.list <- future_lapply(comparisons, function(comp) {
    Caculate_DiffExp_MaAsLin2(
        rel.mat      = rel.mat.fil, # 使用你的过滤后相对丰度矩阵
        metadata     = metadata,
        p.value      = 0.05,
        p.adj        = 0.1,
        foldchange   = 2,           # 丰度差异建议使用较低阈值如1.5，因微生物方差大
        case.type    = comp$case,
        control.type = comp$control,
        outdir_base  = diffAbund_dir
    )
})

plan(sequential)

#-------------------------------------------------------------------------------
# 调整 worker 数量，建议不超过实际文件数量或 CPU 核心数
plan(multisession, workers = 20)  

# 1. 精准匹配 MaAsLin2 输出的差异分析 Excel 文件，避免抓取到日志或其他副产物
all.files <- list.files(".", pattern = "^MaAsLin2_Case-.*\\.xlsx$", recursive = TRUE, full.names = TRUE)

# 并行执行
stat.list <- future_lapply(seq_along(all.files), function(i){
  file.fullpath <- all.files[i]
  
  # 提取文件名并去掉前缀和后缀，剩下 "xxx.vs.Ctrl-yyy"
  base_name <- basename(file.fullpath)
  base_name <- gsub("^MaAsLin2_Case-|\\.xlsx$", "", base_name) 
  
  # 读取表格
  tmp <- read.xlsx(file.fullpath)
  
  # 2. 避免硬编码列号 (原来的 tmp[,7])，直接通过列名 "type" 来统计，安全且稳健
  type_col <- tmp$type
  
  # 使用 sum 统计更简洁高效
  up_count <- sum(type_col == "U", na.rm = TRUE)
  down_count <- sum(type_col == "D", na.rm = TRUE)
  nosign_count <- sum(type_col == "_", na.rm = TRUE)
  
  c(Group = base_name, Up = up_count, Down = down_count, NoSign = nosign_count)
})

# 关闭并行
plan(sequential)

# 合并结果为数据框
stat.df <- do.call(rbind, stat.list) %>% as.data.frame()

# 转换为数值型
stat.df$Up <- as.numeric(stat.df$Up)
stat.df$Down <- as.numeric(stat.df$Down)
stat.df$NoSign <- as.numeric(stat.df$NoSign)

# 3. 排序及文本清洗
stat.df2 <- stat.df %>% 
  arrange(Group) %>%
  # 按照之前设定的 ".vs.Ctrl-" 进行拆分
  tidyr::separate(Group, into = c("Case", "Control"), sep = "\\.vs\\.Ctrl-") %>%
  # 把之前为了防止路径报错加上的下划线替换回空格，让最终报告更美观
  dplyr::mutate(
    Case = gsub("_", " ", Case),
    Control = gsub("_", " ", Control)
  ) %>%
  dplyr::select(Case, Control, Up, Down, NoSign)

# 输出统计表
write.xlsx(stat.df2, "Diff_TaxonNum_Statistics.xlsx", colNames = TRUE, rowNames = FALSE)

setwd(file.path(outdir,"07.diversity"))

## 6. Calculate alpha diversity (Bypassing phyloseq wrapper for relative abundance)
#-------------------------------------------------------------------------------
# 1. 提取并转置矩阵 (vegan 包要求：行是样本，列是物种)
rel_mat_t <- t(as.matrix(rel.mat.fil))

# 2. 计算各项基础指数
# Observed: 在每个样本中，丰度大于 0 的物种数量
Observed <- rowSums(rel_mat_t > 0)
# Shannon & Simpson: vegan::diversity 完美支持相对丰度小数
Shannon <- vegan::diversity(rel_mat_t, index = "shannon")
Simpson <- vegan::diversity(rel_mat_t, index = "simpson")
InvSimpson <- vegan::diversity(rel_mat_t, index = "invsimpson")

# 3. 合并成数据框
alpha_div <- data.frame(
  Observed = Observed,
  Shannon = Shannon,
  Simpson = Simpson,
  InvSimpson = InvSimpson
)

# 4. 计算衍生指数 (Pielou 和 Dominance)
alpha_div$Pielou <- alpha_div$Shannon / log(alpha_div$Observed)
# 防止由于 Observed 为 1 导致 log(1)=0 出现的无穷大 (Inf)
alpha_div$Pielou[is.infinite(alpha_div$Pielou) | is.nan(alpha_div$Pielou)] <- 0

alpha_div$Dominance <- 1 - alpha_div$Simpson

# 5. 合并 metadata
# 确保 metadata 的行名也是样本名，这样 cbind 才不会乱序
alpha_div_meta <- cbind(alpha_div, metadata[rownames(alpha_div), ])

# 预览一下结果
head(alpha_div_meta)

alpha_long <- alpha_div_meta %>%
  # 移除了 Chao1, Fisher, GoodsCoverage，保留实际计算的指数
  select(sample, type, Observed, Shannon, Simpson, InvSimpson, Pielou, Dominance) %>%
  pivot_longer(cols = -c(sample, type), names_to = "Metric", values_to = "Value")

# 确保 type_levels 在你之前的代码中已经定义好，用来控制画图的因子顺序
alpha_long$type <- factor(alpha_long$type, levels = type_levels)

type_levels

my_theme <- 
    theme_bw() +
    theme(title = element_text(angle = 0, colour = "black", hjust = 0.5, vjust = 1, size = 6, face = "bold"),
      axis.text.x = element_text(angle = 0, colour = "black", hjust = 0.5, vjust = 1, size = 6, face = "bold"),
      axis.text.y = element_text(size = 6, colour = "black", face = "bold"),
      axis.title  = element_text(size = 6, colour = "black", face = "bold"),
      legend.title = element_text(size = 6, colour = "black", face = "bold"),
      legend.text  = element_text(size = 6, colour = "black", face = "bold"),
      panel.border = element_rect(color = "black", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 6, colour = "black", face = "bold")
    )
    

alpha_p_t_test <- ggplot(alpha_long, aes(x = type, y = Value, fill = type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = type), width = 0.2, size = 0.1, alpha = 0.8) +
  scale_fill_manual(values = color.panel) +
  scale_color_manual(values = color.panel) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +   # 分面，每个指标一个子图
  stat_compare_means(comparisons = comp.groups, method = "t.test", label = "p.format", size = 2.1) +
  my_theme

alpha_p_wilcox_test <- ggplot(alpha_long, aes(x = type, y = Value, fill = type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = type), width = 0.2, size = 0.1, alpha = 0.8) +
  scale_fill_manual(values = color.panel) +
  scale_color_manual(values = color.panel) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +   # 分面，每个指标一个子图
  stat_compare_means(comparisons = comp.groups, method = "wilcox.test", label = "p.format", size = 2.1) +
  my_theme

ggsave(file.path(outdir,"07.diversity","o1.AlphaDiv.tTest.Pvalue.pdf"),alpha_p_t_test,width = 600/72,height = 400/72,units = "in")
ggsave(file.path(outdir,"07.diversity","o1.AlphaDiv.wilcoxTest.Pvalue.pdf"),alpha_p_wilcox_test,width = 600/72,height = 400/72,units = "in")

alpha_p_t_test <- ggplot(alpha_long, aes(x = type, y = Value, fill = type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = type), width = 0.2, size = 0.1, alpha = 0.8) +
  scale_fill_manual(values = color.panel) +
  scale_color_manual(values = color.panel) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +   # 分面，每个指标一个子图
  stat_compare_means(comparisons = comp.groups, method = "t.test", label = "p.signif", size = 2.1) +
  my_theme

alpha_p_wilcox_test <- ggplot(alpha_long, aes(x = type, y = Value, fill = type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = type), width = 0.2, size = 0.1, alpha = 0.8) +
  scale_fill_manual(values = color.panel) +
  scale_color_manual(values = color.panel) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +   # 分面，每个指标一个子图
  stat_compare_means(comparisons = comp.groups, method = "wilcox.test", label = "p.signif", size = 2.1) +
  my_theme

ggsave(file.path(outdir,"07.diversity","o1.AlphaDiv.tTest.Psignif.pdf"),alpha_p_t_test,width = 600/72,height = 400/72,units = "in")
ggsave(file.path(outdir,"07.diversity","o1.AlphaDiv.wilcoxTest.Psignif.pdf"),alpha_p_wilcox_test,width = 600/72,height = 400/72,units = "in")

write.csv(alpha_div_meta, file.path(outdir,"07.diversity","o1.AlphaDiversity.tbl.csv"), row.names=T)



#-------------------------------------------------------------------------------
# 1. 构建 phyloseq 对象 (只需使用相对丰度矩阵)
OTU_rel <- otu_table(as.matrix(rel.mat.fil), taxa_are_rows = TRUE)
SAM <- sample_data(metadata)
physeq_rel <- phyloseq(OTU_rel, SAM)

# 获取保留的物种名称
kept_taxa <- taxa_names(physeq_rel)

set.seed(666)

# phyloseq 支持的 method 和 distance 列表
methods <- c("NMDS", "PCoA", "MDS", "DCA", "CCA", "RDA", "DPCoA")
distances <- c("bray", "jaccard", "euclidean", "manhattan", "wunifrac", "unifrac")

# 判断哪些方法需要 distance，哪些需要 formula
needs_formula <- c("CCA", "RDA")
ignore_distance <- c("DCA")

# 遍历所有组合
combos <- expand.grid(method = methods, dist = distances, stringsAsFactors = FALSE)

plist_all <- dlply(combos, .(method, dist), function(x) {
  meth <- x$method
  dist <- x$dist
  
  dist_use <- if (meth %in% ignore_distance) NULL else dist
  
  # CCA/RDA 需要公式
  if (meth %in% needs_formula) {
    if (!"type" %in% colnames(sample_data(physeq_rel))) return(NULL)
    ordi <- tryCatch(
      ordinate(physeq_rel, method = meth, distance = dist_use, formula = ~type),
      error = function(e) return(NULL)
    )
  } else {
    ordi <- tryCatch(
      ordinate(physeq_rel, method = meth, distance = dist_use),
      error = function(e) return(NULL)
    )
  }
  
  if (is.null(ordi)) return(NULL)
  
     p <- plot_ordination(physeq_rel, ordi, type="samples") +
              geom_point(
                aes(fill = type, color = type),
                size = 3, stroke = 0.5, alpha = 0.9
              ) +
              stat_ellipse(
                aes(color = type),   # 按组画椭圆
                type = "t",          # 常用 t 分布（也可以 "norm"）
                level = 0.95,        # 95% 置信椭圆
                linewidth = 0.6,     # 线粗细
                linetype = 1         # 实线（可改成 2 = 虚线）
              ) +
              ggtitle(paste0("Method: ", meth, "; Distance: ", dist)) +
              scale_shape_manual(values = c(21, 24)) +          
              scale_fill_manual(values = color.panel) +     
              scale_color_manual(values = color.panel) +
              my_theme
      
  return(p)
})
 

pdf(file.path(outdir,"07.diversity","BetaDiversity_All.pdf"), width = 320/72, height = 250/72)
for (p in plist_all) {
  if (!is.null(p)) print(p)
}
dev.off()


