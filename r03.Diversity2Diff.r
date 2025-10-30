### 20/08/2025
### XIAO Han Shawn
### Purpose:
### *** 1. Alpha/Beta diversity analysis 
###     2. Discard gender factor
###     3. 

pkgs <- c(
  "dplyr", "openxlsx", "phyloseq", "fossil", 
  "ggplot2", "vegan", "ape", "microbiome", 
  "ggpubr", "plyr", "scales", "tidyr",
  "optparse", "DESeq2", "data.table",
  "future", "future.apply"
)

invisible(lapply(pkgs, function(pkg) {
  suppressMessages(require(pkg, character.only = TRUE))
}))

## 1. Load data
#-------------------------------------------------------------------------------
wkdir <- "/home/han_xiao/xiaohan/Work/02-SpatialGF/02.script/Stool_Metagenomics/r03.Diveristy2Diff.ipynb"
indir <- "/home/han_xiao/xiaohan/Work/02-SpatialGF/03.result/21.Stool_metagenomics/04.bracken/"
outdir <- "/home/han_xiao/xiaohan/Work/02-SpatialGF/03.result/21.Stool_metagenomics"

# if(!dir.exists(wkdir)){dir.create(wkdir,recursive = T)}
if(!dir.exists(indir)){dir.create(indir,recursive = T)}
if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}

setwd(outdir)

if(!dir.exists("06.diversity")){dir.create("06.diversity")}
if(!dir.exists("07.diffAbund")){dir.create("07.diffAbund")}

abs.abund.fpath <- file.path(indir,"bracken.S.0.1-H.txt")
abs.abund.df <- read.delim(abs.abund.fpath,header = T,row.names = 1)
rel.abund.df <- as.data.frame(apply(abs.abund.df, 2, function(x) x / sum(x)))

# Prepare sample metadata

metadata.fpath <- "/home/han_xiao/xiaohan/Work/02-SpatialGF/01.data/Stool_metagenomics/i1.meta.csv"

metadata <- read.csv(metadata.fpath,header = T,row.names = 1)
colnames(metadata) <- c("type","color")     
metadata$sample <- rownames(metadata)
color.panel <- metadata %>% select(type,color) %>% distinct(type,color) %>% pull(color)
names(color.panel) <- metadata  %>% select(type,color) %>% distinct(type,color) %>% pull(type)
                                    
# cut off of mean relative abundance
cutoff <- 0.0001    

# comparision groups
groups <- names(color.panel)
# comp.groups <- combn(groups, 2, simplify = FALSE)
# comp.groups

metadata <- metadata %>% mutate(
    Source = sapply(strsplit(type,"_",fixed = T),"[[",1),
    Condition = sapply(strsplit(type,"_",fixed = T),"[[",2),
    Gender = gsub("[0-9]+","",sapply(strsplit(type,"_",fixed = T),"[[",3))
)
metadata

condition.panel <- c(
  `D3`    = "#add9a2",  # 淡蓝
  `1W`    = "#63c4dc",  # 粉玫瑰
  `2W`    = "#f4c098",  # 浅橙
  `4W`    = "#d093c8"   # 青绿
)

gender.panel <- c(
    "F" = "#e2d783",
    "M"   = "#519795"
)

setwd(outdir)
wkdir <- file.path(outdir,"07.diffAbund","00.Metadata")
if(!dir.exists(wkdir)){ dir.create(wkdir,recursive = T)}
setwd(wkdir)

getwd()


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

setwd(file.path(outdir,"07.diffAbund"))
if(!dir.exists("01.Diff")){dir.create("01.Diff")}
setwd("01.Diff")

#-------------------------------------------------------------------------------
Caculate_DiffExp = function(fpkm.mat, fpkm.mat.fpath=NULL,
                            count.mat, count.mat.fpath=NULL,
                            metadata, metadata.fpath=NULL,
                            foldchange = 2, p.value = 0.05, p.adj = 0.1,
                            case.type = NULL, control.type = NULL){

    print(paste0("FC: ", foldchange, "; pvalue: ", p.value, "; p.adj: ",p.adj))
    
    # 2. Load data -------------------------------------------------------------------
    # Normalize expression matrix
    if(!is.null(fpkm.mat.fpath)){
        fpkm.mat <- read.csv(fpkm.mat.fpath,header = T)
        rownames(fpkm.mat) <- fpkm.mat$ID
        fpkm.mat <- fpkm.mat[order(rownames(fpkm.mat)),]
    } else {
        fpkm.mat <- fpkm.mat
        fpkm.mat$ID <- rownames(fpkm.mat)
    }
    
    # Raw count matrix
    if(!is.null(count.mat.fpath)){
        count.mat <- read.csv(count.mat.fpath,header = T)
        rownames(count.mat) <- count.mat$ID
        count.mat <- count.mat[order(rownames(count.mat)),]
    } else {
        count.mat <- count.mat
        count.mat$ID <- rownames(count.mat)
    }

    # Sample metadata info
    if(!is.null(metadata.fpath)){
        all_info <- read.csv(metadata.fpath,header = F)
        metadata <- all_info
        metadata[,1] <- gsub("-",".",metadata[,1])
        colnames(metadata) <- c("sample","type")
    } else {
        metadata = metadata %>% dplyr::select(sample,type)
        colnames(metadata) <- c("sample","type")
    }

    # Define case/control -------------------------------------------------------------
    if (!is.null(case.type) & !is.null(control.type)) {
        type1 <- case.type
        type2 <- control.type
    } else if(!is.null(case.type)) {
        type1 <- case.type
        type2 <- setdiff(unique(metadata$type), type1)   # 默认 vs 所有其他
    } else {
        stop("Please at least denote the case.type")
    }
    
    # create outdir
    outdir <- "."
    
    # 3. log2FC + DESeq2 --------------------------------------------------------------
    sample1 <- metadata$sample[metadata$type==type1]
    sample2 <- metadata$sample[metadata$type==type2]
    
    print(paste("Comparing:", type1, "vs", type2))
    
    # 计算 log2FC
    FC <- as.data.frame(matrix(nrow = dim(fpkm.mat)[1], ncol = 4,data = 0))
    colnames(FC) <- c("ID",paste0(type1,".AvgExp"),paste0(type2,".AvgExp"),"Log2FC")
    for (i in 1:dim(fpkm.mat)[1]) {
        FC[i,1] <- fpkm.mat[i,"ID"]
        FC[i,2] <- mean(as.numeric(fpkm.mat[i,which(colnames(fpkm.mat) %in% sample1)]))
        FC[i,3] <- mean(as.numeric(fpkm.mat[i,which(colnames(fpkm.mat) %in% sample2)]))
        FC[i,4] <- log2((FC[i,2] + (1/1e6))/(FC[i,3] + (1/1e6)))
    }
    
    # DESeq2
    absData <- count.mat[,c(which(colnames(count.mat) %in% sample1),
                            which(colnames(count.mat) %in% sample2))]
    condition <- factor(c(rep(type1,length(sample1)),
                          rep(type2,length(sample2))), levels = c(type2,type1))
    colData <- data.frame(row.names=colnames(absData), condition)
    dds <- DESeqDataSetFromMatrix(absData, DataFrame(condition), design= ~ condition )
    dds2 <- DESeq(dds) 
    res <- results(dds2)
    res$ID <- rownames(res)
    
    # 合并结果
    tmp <- as.data.frame(res)[,c("ID","pvalue","padj")]
    tmp[is.na(tmp)] <- 1
    tmp2 <- merge(FC,tmp,by="ID")
    tmp2$type <- "_"
    tmp2[tmp2$Log2FC > log2(foldchange) & tmp2$pvalue < p.value & tmp2$padj < p.adj, "type"] <- "U"
    tmp2[tmp2$Log2FC < -log2(foldchange) & tmp2$pvalue < p.value & tmp2$padj < p.adj, "type"] <- "D"
    colnames(tmp2)[7] <- paste(type1,"vs",type2,sep = ".")
    
    # 输出结果 ------------------------------------------------------------------------
    write.xlsx(tmp2,paste0(outdir,"/Case-",type1,".vs.Ctrl-",type2,".xlsx"),
               colNames=T,rowNames=F)
    
    return(tmp2)
}

comparisons <- list(
  # ---- Time progression comparisons ----
  list(case = "1W", control = "D3"),
  list(case = "2W", control = "1W"),
  list(case = "4W", control = "2W"),
  
  # ---- Treatment vs Donor ----
  list(case = "D3", control = "Donor"),
  list(case = "1W", control = "Donor"),
  list(case = "2W", control = "Donor"),
  list(case = "4W", control = "Donor")
)

metadata$type <- metadata$Condition

plan(multisession, workers = 39)

# 并行执行
res.list <- future_lapply(comparisons, function(comp) {
  Caculate_DiffExp(
    fpkm.mat = rel.mat.fil,
    count.mat = abs.mat.fil,
    metadata = metadata,
    p.value = 0.05,
    p.adj = 0.1,
    foldchange = 2,
    case.type = comp$case,
    control.type = comp$control   # 这里要确保函数支持传入 control.type
  )
})

plan(sequential)

#-------------------------------------------------------------------------------

plan(multisession, workers = 20)  # 根据你的CPU核数调整

all.files <- grep("Case", list.files(".", recursive = TRUE, full.names = TRUE), value = TRUE)

# 并行执行
stat.list <- future_lapply(seq_along(all.files), function(i){
  file.fullpath <- all.files[i]
  res_name <- gsub("DiffExp_|\\.xlsx","",basename(file.fullpath))
  tmp <- read.xlsx(file.fullpath)
  
  up_count <- ifelse("U" %in% tmp[,7], table(tmp[,7])["U"], 0)
  down_count <- ifelse("D" %in% tmp[,7], table(tmp[,7])["D"], 0)
  nosign_count <- ifelse("_" %in% tmp[,7], table(tmp[,7])["_"], 0)
  
  # 添加基因symbol列
  # tmp2 <- cbind(gid2symbol.df, tmp[,-1])
  # write.xlsx(tmp2, file.fullpath, overwrite = TRUE)
  
  c(Group = res_name, Up = up_count, Down = down_count, NoSign = nosign_count)
})

# 关闭并行
plan(sequential)

# 合并结果为数据框
stat.df <- do.call(rbind, stat.list) %>% as.data.frame()
stat.df$Up <- as.numeric(stat.df$Up)
stat.df$Down <- as.numeric(stat.df$Down)
stat.df$NoSign <- as.numeric(stat.df$NoSign)

# 排序
stat.df2 <- stat.df %>% arrange(Group)

stat.df2 <- stat.df2 %>%
  tidyr::separate(Group, into = c("Case", "Control"), sep = "\\.vs\\.") %>%
  dplyr::select(Case, Control, Up, Down, NoSign)

write.xlsx(stat.df2,"Diff_TaxonNum_Statistics.xlsx",colName=T,rowNames=F)


setwd(file.path(outdir,"06.diversity"))

#-------------------------------------------------------------------------------
OTU_abs <- otu_table(as.matrix(abs.mat.fil), taxa_are_rows = TRUE)
OTU_rel <- otu_table(as.matrix(rel.mat.fil), taxa_are_rows = TRUE)
# TAX <- tax_table(taxonomy_matrix)

SAM <- sample_data(metadata)

physeq_abs <- phyloseq(OTU_abs, SAM)
physeq_rel <- phyloseq(OTU_rel, SAM)

physeq_abs

# 获取保留的物种名称
kept_taxa <- taxa_names(physeq_abs)


## 6. Calculate alpha diversity
#-------------------------------------------------------------------------------
set.seed(666)
# 3.1 For richness index ("Observed", "Chao1", "ACE") - 
physeq_rarefied <- rarefy_even_depth(
  physeq_abs,
  sample.size = min(sample_sums(physeq_abs)), # 下采样到最小测序深度
  rngseed = 666,                          # 设置随机种子保证可重复
  replace = FALSE,                        # 不放回抽样
  trimOTUs = TRUE                         # 去掉抽样后为0的OTU
)
alpha_div_abs <- estimate_richness(physeq_rarefied,
                               measures = c("Observed", "Chao1", "ACE"))

# 3.2 For evenness index ("Shannon", "Simpson", "InvSimpson", "Pielou's Evenness")
alpha_div_rel <- estimate_richness(physeq_rel,
                               measures = c("Shannon", "Simpson", "InvSimpson"))
# Add Pielou's Evenness
alpha_div_rel$Pielou <- alpha_div_rel$Shannon / log(alpha_div_abs$Observed)

# 3.3 For combinial index ("Fisher's alpha", "Dominance", "Good's Coverage")
alpha_div_com <- estimate_richness(physeq_rarefied,  # 改成稀释后的 count 数据
                               measures = c("Fisher"))

# Extended Metrics (Good’s Coverage, Dominance)
otu_abs_t <- t(as(otu_table(physeq_rarefied), "matrix"))
alpha_div_com$GoodsCoverage <- apply(otu_abs_t, 1, function(x) {
  singletons <- sum(x == 1)
  n <- sum(x)
  1 - (singletons / n)
})

# Simpson's Dominance (基于相对丰度)
alpha_div_com$Dominance <- 1 - alpha_div_rel$Simpson

## Merge Data

alpha_div <- cbind(alpha_div_abs, alpha_div_rel[,c("Shannon","Simpson","InvSimpson","Pielou")], alpha_div_com)

# Merge with metadata

alpha_div_meta <- cbind(alpha_div, data.frame(sample_data(physeq_abs)))

alpha_long <- alpha_div_meta %>%
  select(sample, type, Observed, Chao1, Shannon, Simpson, InvSimpson, Pielou, Fisher, Dominance, GoodsCoverage) %>%
  pivot_longer(cols = -c(sample, type), names_to = "Metric", values_to = "Value")

alpha_long$type <- factor(alpha_long$type,levels = c(
  "D_YSF","D_YSM","D_YCF","D_YCM","D_ASF","D_ASM","D_ACF","D_ACM",
  "R_YSF","R_YSM","R_YCF","R_YCM","R_ASF","R_ASM","R_ACF","R_ACM"
))

comp.groups <- list(
  # ---- Age (A vs Y), same diet/sex/host ----
  c("R_ACF", "R_YCF"),
  c("R_ASF", "R_YSF"),
  c("R_ACM", "R_YCM"),
  c("R_ASM", "R_YSM"),
  c("D_ACF", "D_YCF"),
  c("D_ASF", "D_YSF"),
  c("D_ACM", "D_YCM"),
  c("D_ASM", "D_YSM"),
  
  # ---- Diet (C vs S), same age/sex/host ----
  c("R_ACF", "R_ASF"),
  c("R_YCF", "R_YSF"),
  c("R_ACM", "R_ASM"),
  c("R_YCM", "R_YSM"),
  c("D_ACF", "D_ASF"),
  c("D_YCF", "D_YSF"),
  c("D_ACM", "D_ASM"),
  c("D_YCM", "D_YSM"),
  
  # ---- Sex (F vs M), same age/diet/host ----
  c("R_ACF", "R_ACM"),
  c("R_ASF", "R_ASM"),
  c("R_YCF", "R_YCM"),
  c("R_YSF", "R_YSM"),
  c("D_ACF", "D_ACM"),
  c("D_ASF", "D_ASM"),
  c("D_YCF", "D_YCM"),
  c("D_YSF", "D_YSM")
)


my_theme <- 
    theme_bw() +
    theme(plot.title = element_text(angle = 0, colour = "black", hjust = 0.5, vjust = 1, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 90, colour = "black", hjust = 1, vjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, colour = "black", face = "bold"),
      axis.title  = element_text(size = 16, colour = "black", face = "bold"),
      legend.title = element_text(size = 16, colour = "black", face = "bold"),
      legend.text  = element_text(size = 16, colour = "black", face = "bold"),
      panel.border = element_rect(color = "black", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 16, colour = "black", face = "bold")
    )

alpha_p_t_test <- ggplot(alpha_long, aes(x = type, y = Value, fill = type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = type), width = 0.2, size = 2, alpha = 0.8) +
  scale_fill_manual(values = color.panel) +
  scale_color_manual(values = color.panel) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +   # 分面，每个指标一个子图
  stat_compare_means(comparisons = comp.groups, method = "t.test", label = "p.format") +
  my_theme

alpha_p_wilcox_test <- ggplot(alpha_long, aes(x = type, y = Value, fill = type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = type), width = 0.2, size = 2, alpha = 0.8) +
  scale_fill_manual(values = color.panel) +
  scale_color_manual(values = color.panel) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +   # 分面，每个指标一个子图
  stat_compare_means(comparisons = comp.groups, method = "wilcox.test", label = "p.format") +
  my_theme

ggsave(file.path(outdir,"06.diversity","o1.AlphaDiv.tTest.Pvalue.pdf"),alpha_p_t_test,width = 25,height = 14)
ggsave(file.path(outdir,"06.diversity","o1.AlphaDiv.wilcoxTest.Pvalue.pdf"),alpha_p_wilcox_test,width = 25,height = 14)

alpha_long <- alpha_div_meta %>%
  select(sample, type, Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Pielou, Fisher, Dominance, GoodsCoverage) %>%
  pivot_longer(cols = -c(sample, type), names_to = "Metric", values_to = "Value")

alpha_p_t_test <- ggplot(alpha_long, aes(x = type, y = Value, fill = type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = type), width = 0.2, size = 2, alpha = 0.8) +
  scale_fill_manual(values = color.panel) +
  scale_color_manual(values = color.panel) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +   # 分面，每个指标一个子图
  stat_compare_means(comparisons = comp.groups, method = "t.test", label = "p.signif") +
  my_theme

alpha_p_wilcox_test <- ggplot(alpha_long, aes(x = type, y = Value, fill = type)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = type), width = 0.2, size = 2, alpha = 0.8) +
  scale_fill_manual(values = color.panel) +
  scale_color_manual(values = color.panel) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +   # 分面，每个指标一个子图
  stat_compare_means(comparisons = comp.groups, method = "wilcox.test", label = "p.signif") +
  my_theme

ggsave(file.path(outdir,"06.diversity","o1.AlphaDiversity.tTest.Psignif.pdf"),alpha_p_t_test,width = 25,height = 14)
ggsave(file.path(outdir,"06.diversity","o1.AlphaDiversity.wilcoxTest.Psignif.pdf"),alpha_p_wilcox_test,width = 25,height = 14)

write.csv(alpha_div_meta, file.path(outdir,"06.diversity","o1.AlphaDiversity.tbl.csv"), row.names=T)

#-------------------------------------------------------------------------------

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
        aes(shape = Source, fill = Condition, color = Gender),
        size = 3, stroke = 0.5, alpha = 0.9
      ) +
      ggtitle(paste0("Method: ", meth, "; Distance: ", dist)) +
      scale_shape_manual(values = c(21, 24)) +          
      scale_fill_manual(values = condition.panel, guide = guide_legend(override.aes = list(shape = 21))) +     
      scale_color_manual(values = c(
        "F" = "#C0392B",  # 深红色 (Female)
        "M" = "#2E86C1"   # 深蓝色 (Male)
      ), guide = guide_legend(override.aes = list(shape = 21))) +
      theme_bw() +
      theme(
        panel.background = element_blank(),
        axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color = "black", size = 16),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", linewidth = 1)
      )
      
  return(p)
})
 

pdf(file.path(outdir,"06.diversity","BetaDiversity_All.pdf"), width = 7, height = 6)
for (p in plist_all) {
  if (!is.null(p)) print(p)
}
dev.off()


