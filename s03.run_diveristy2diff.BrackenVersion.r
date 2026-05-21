#!/usr/bin/env Rscript

# ==============================================================================
# 0. Load Packages & Parse Arguments
# ==============================================================================
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

# Define command-line arguments
option_list = list(
  make_option(c("-a", "--abs_abund"), type="character", default=NULL, 
              help="Absolute abundance matrix file path", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="Sample metadata CSV file path", metavar="character"),
  make_option(c("-c", "--cutoff"), type="numeric", default=0.0001, 
              help="Cutoff for mean relative abundance [default= %default]", metavar="numeric"),
  make_option(c("-o", "--outdir"), type="character", default=".", 
              help="Output directory path", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$abs_abund) || is.null(opt$metadata)){
  print_help(opt_parser)
  stop("Missing required arguments: Please provide both --abs_abund and --metadata.\n", call.=FALSE)
}

# Assign arguments to variables
abs.abund.fpath <- opt$abs_abund
metadata.fpath <- opt$metadata
cutoff <- opt$cutoff
outdir <- opt$outdir

# ==============================================================================
# 1. Directory Setup & Load Data
# ==============================================================================
dir_div <- file.path(outdir, "05.diversity")
dir_diff <- file.path(outdir, "06.diffAbund", "01.Diff")
dir_meta <- file.path(outdir, "06.diffAbund", "00.Metadata")

invisible(sapply(c(dir_div, dir_diff, dir_meta), function(d) {
  if(!dir.exists(d)) dir.create(d, recursive = TRUE)
}))

# Load abundance
abs.abund.df <- read.delim(abs.abund.fpath, header = TRUE, row.names = 1)
rel.abund.df <- as.data.frame(apply(abs.abund.df, 2, function(x) x / sum(x)))

# Load metadata
metadata <- read.csv(metadata.fpath, header = TRUE, row.names = 1)
colnames(metadata) <- c("type", "color")     
metadata$sample <- rownames(metadata)

color.panel <- metadata %>% select(type, color) %>% distinct(type, color) %>% pull(color)
names(color.panel) <- metadata %>% select(type, color) %>% distinct(type, color) %>% pull(type)
                                    
groups <- names(color.panel)
type_levels <- groups

comp.groups <- list(
  c("V1_2", "V1_1"),
  c("V2_2", "V2_1"),
  c("V3_2", "V3_1"),
  c("V4_2", "V4_1")
)

comparisons <- list(
  list(case = "V1_2", control = "V1_1"),
  list(case = "V2_2", control = "V2_1"),
  list(case = "V3_2", control = "V3_1"),
  list(case = "V4_2", control = "V4_1")
)

# ==============================================================================
# 2. Filter Species
# ==============================================================================
Filter_Species = function(abs.mat, rel.mat, mean.rel.cutoff = 0.0001) {
    sp.mean.across.samples = rowMeans(rel.mat)
    sp.filtered.samples = rownames(rel.mat)[sp.mean.across.samples > mean.rel.cutoff]
    
    if (length(sp.filtered.samples) == 0) {
        stop("No species passed the mean.rel.cutoff filter")
    }
    
    filtered.abs.mat = abs.mat[sp.filtered.samples, , drop = FALSE]
    filtered.rel.mat = rel.mat[sp.filtered.samples, , drop = FALSE]
    
    return(list(
        filtered.abs = filtered.abs.mat,
        filtered.rel = filtered.rel.mat,
        cutoff = mean.rel.cutoff
    ))
}

filtered_data <- Filter_Species(abs.abund.df, rel.abund.df, mean.rel.cutoff = cutoff)
abs.mat.fil <- filtered_data[["filtered.abs"]]
rel.mat.fil <- filtered_data[["filtered.rel"]]

cat("\n[INFO] Cut-off of mean relative abundance across all samples:", cutoff, "\n")
cat("[INFO] Remaining taxa:", nrow(abs.mat.fil), "\n\n")

# Save filtered matrices in metadata folder
write.xlsx(abs.mat.fil, file.path(dir_meta, paste0("AbsoluteAbundance_Cutoff", cutoff, ".xlsx")), colNames=TRUE, rowNames=TRUE)
write.xlsx(rel.mat.fil, file.path(dir_meta, paste0("RelativeAbundance_Cutoff", cutoff, ".xlsx")), colNames=TRUE, rowNames=TRUE)


# ==============================================================================
# 3. Alpha & Beta Diversity (05.diversity)
# ==============================================================================
cat("[INFO] Running Alpha & Beta Diversity analysis...\n")

OTU_abs <- otu_table(as.matrix(abs.mat.fil), taxa_are_rows = TRUE)
OTU_rel <- otu_table(as.matrix(rel.mat.fil), taxa_are_rows = TRUE)
SAM <- sample_data(metadata)

physeq_abs <- phyloseq(OTU_abs, SAM)
physeq_rel <- phyloseq(OTU_rel, SAM)

# --- Alpha Diversity ---
set.seed(666)
physeq_rarefied <- rarefy_even_depth(physeq_abs, sample.size = min(sample_sums(physeq_abs)), rngseed = 666, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
alpha_div_abs <- estimate_richness(physeq_rarefied, measures = c("Observed", "Chao1", "ACE"))
alpha_div_rel <- estimate_richness(physeq_rel, measures = c("Shannon", "Simpson", "InvSimpson"))
alpha_div_rel$Pielou <- alpha_div_rel$Shannon / log(alpha_div_abs$Observed)
alpha_div_com <- estimate_richness(physeq_rarefied, measures = c("Fisher"))

otu_abs_t <- t(as(otu_table(physeq_rarefied), "matrix"))
alpha_div_com$GoodsCoverage <- apply(otu_abs_t, 1, function(x) {
  singletons <- sum(x == 1)
  n <- sum(x)
  1 - (singletons / n)
})
alpha_div_com$Dominance <- 1 - alpha_div_rel$Simpson

alpha_div <- cbind(alpha_div_abs, alpha_div_rel[,c("Shannon","Simpson","InvSimpson","Pielou")], alpha_div_com)
alpha_div_meta <- cbind(alpha_div, data.frame(sample_data(physeq_abs)))

write.csv(alpha_div_meta, file.path(dir_div, "o1.AlphaDiversity.tbl.csv"), row.names=TRUE)

# Plotting Alpha
alpha_long <- alpha_div_meta %>%
  select(sample, type, Observed, Chao1, Shannon, Simpson, InvSimpson, Pielou, Fisher, Dominance, GoodsCoverage) %>%
  pivot_longer(cols = -c(sample, type), names_to = "Metric", values_to = "Value")
alpha_long$type <- factor(alpha_long$type, levels = type_levels)

my_theme <- theme_bw() + theme(
  title = element_text(angle = 0, colour = "black", hjust = 0.5, vjust = 1, size = 6, face = "bold"),
  axis.text.x = element_text(angle = 45, colour = "black", hjust = 1, vjust = 1, size = 6, face = "bold"),
  axis.text.y = element_text(size = 6, colour = "black", face = "bold"),
  axis.title  = element_text(size = 6, colour = "black", face = "bold"),
  legend.title = element_text(size = 6, colour = "black", face = "bold"),
  legend.text  = element_text(size = 6, colour = "black", face = "bold"),
  panel.border = element_rect(color = "black", linewidth = 1),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  strip.text = element_text(size = 6, colour = "black", face = "bold")
)

plot_alpha <- function(method, label) {
  ggplot(alpha_long, aes(x = type, y = Value, fill = type)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(aes(color = type), width = 0.2, size = 0.1, alpha = 0.8) +
    scale_fill_manual(values = color.panel) +
    scale_color_manual(values = color.panel) +
    facet_wrap(~ Metric, scales = "free_y", nrow = 2) + 
    stat_compare_means(comparisons = comp.groups, method = method, label = label, size = 2.1) +
    my_theme
}

ggsave(file.path(dir_div, "o1.AlphaDiv.tTest.Pvalue.pdf"), plot_alpha("t.test", "p.format"), width=600/72, height=400/72, units="in")
ggsave(file.path(dir_div, "o1.AlphaDiv.wilcoxTest.Pvalue.pdf"), plot_alpha("wilcox.test", "p.format"), width=600/72, height=400/72, units="in")
ggsave(file.path(dir_div, "o1.AlphaDiv.tTest.Psignif.pdf"), plot_alpha("t.test", "p.signif"), width=600/72, height=400/72, units="in")
ggsave(file.path(dir_div, "o1.AlphaDiv.wilcoxTest.Psignif.pdf"), plot_alpha("wilcox.test", "p.signif"), width=600/72, height=400/72, units="in")

# --- Beta Diversity ---
methods <- c("NMDS", "PCoA", "MDS", "DCA", "CCA", "RDA", "DPCoA")
distances <- c("bray", "jaccard", "euclidean", "manhattan", "wunifrac", "unifrac")
needs_formula <- c("CCA", "RDA")
ignore_distance <- c("DCA")

combos <- expand.grid(method = methods, dist = distances, stringsAsFactors = FALSE)

plist_all <- dlply(combos, .(method, dist), function(x) {
  meth <- x$method
  dist <- x$dist
  dist_use <- if (meth %in% ignore_distance) NULL else dist
  
  if (meth %in% needs_formula) {
    if (!"type" %in% colnames(sample_data(physeq_rel))) return(NULL)
    ordi <- tryCatch(ordinate(physeq_rel, method = meth, distance = dist_use, formula = ~type), error = function(e) return(NULL))
  } else {
    ordi <- tryCatch(ordinate(physeq_rel, method = meth, distance = dist_use), error = function(e) return(NULL))
  }
  
  if (is.null(ordi)) return(NULL)
  
  plot_ordination(physeq_rel, ordi, type="samples") +
    geom_point(aes(fill = type, color = type), size = 3, stroke = 0.5, alpha = 0.9) +
    stat_ellipse(aes(color = type), type = "t", level = 0.95, linewidth = 0.6, linetype = 1) +
    ggtitle(paste0("Method: ", meth, "; Distance: ", dist)) +
    scale_shape_manual(values = c(21, 24)) +        
    scale_fill_manual(values = color.panel) +     
    scale_color_manual(values = color.panel) +
    my_theme
})

pdf(file.path(dir_div, "BetaDiversity_All.pdf"), width = 5, height = 4)
for (p in plist_all) {
  if (!is.null(p)) print(p)
}
invisible(dev.off())


# ==============================================================================
# 4. Differential Abundance (DESeq2) (06.diffAbund)
# ==============================================================================
cat("[INFO] Running Differential Abundance (DESeq2) analysis...\n")

Caculate_DiffExp = function(fpkm.mat, count.mat, metadata, foldchange = 2, p.value = 0.05, p.adj = 0.1, case.type = NULL, control.type = NULL, out_path = ".") {
    print(paste0("Comparing: ", case.type, " vs ", control.type, " | FC: ", foldchange, "; pvalue: ", p.value, "; p.adj: ", p.adj))
    
    fpkm.mat$ID <- rownames(fpkm.mat)
    count.mat$ID <- rownames(count.mat)
    
    type1 <- case.type
    type2 <- control.type
    
    sample1 <- metadata$sample[metadata$type==type1]
    sample2 <- metadata$sample[metadata$type==type2]
    
    FC <- as.data.frame(matrix(nrow = dim(fpkm.mat)[1], ncol = 4, data = 0))
    colnames(FC) <- c("ID", paste0(type1, ".AvgExp"), paste0(type2, ".AvgExp"), "Log2FC")
    for (i in 1:dim(fpkm.mat)[1]) {
        FC[i,1] <- fpkm.mat[i,"ID"]
        FC[i,2] <- mean(as.numeric(fpkm.mat[i, which(colnames(fpkm.mat) %in% sample1)]))
        FC[i,3] <- mean(as.numeric(fpkm.mat[i, which(colnames(fpkm.mat) %in% sample2)]))
        FC[i,4] <- log2((FC[i,2] + (1/1e6))/(FC[i,3] + (1/1e6)))
    }
    
    absData <- count.mat[, c(which(colnames(count.mat) %in% sample1), which(colnames(count.mat) %in% sample2))]
    condition <- factor(c(rep(type1,length(sample1)), rep(type2,length(sample2))), levels = c(type2,type1))
    
    dds <- DESeqDataSetFromMatrix(absData, DataFrame(condition), design= ~ condition)
    dds2 <- DESeq(dds, quiet = TRUE) 
    res <- results(dds2)
    res$ID <- rownames(res)
    
    tmp <- as.data.frame(res)[,c("ID","pvalue","padj")]
    tmp[is.na(tmp)] <- 1
    tmp2 <- merge(FC, tmp, by="ID")
    tmp2$type <- "_"
    tmp2[tmp2$Log2FC > log2(foldchange) & tmp2$pvalue < p.value & tmp2$padj < p.adj, "type"] <- "U"
    tmp2[tmp2$Log2FC < -log2(foldchange) & tmp2$pvalue < p.value & tmp2$padj < p.adj, "type"] <- "D"
    colnames(tmp2)[7] <- paste(type1, "vs", type2, sep = ".")
    
    write.xlsx(tmp2, file.path(out_path, paste0("Case-", type1, ".vs.Ctrl-", type2, ".xlsx")), colNames=TRUE, rowNames=FALSE)
    return(tmp2)
}

plan(multisession, workers = 29) # 可根据服务器配置自行调整
res.list <- future_lapply(comparisons, function(comp) {
  Caculate_DiffExp(
    fpkm.mat = rel.mat.fil,
    count.mat = abs.mat.fil,
    metadata = metadata,
    p.value = 0.05,
    p.adj = 0.1,
    foldchange = 2,
    case.type = comp$case,
    control.type = comp$control,
    out_path = dir_diff
  )
})
plan(sequential)

# Stats summary
all.files <- list.files(dir_diff, pattern = "^Case.*\\.xlsx$", full.names = TRUE)

plan(multisession, workers = 20)
stat.list <- future_lapply(seq_along(all.files), function(i){
  file.fullpath <- all.files[i]
  res_name <- gsub("Case-|\\.xlsx", "", basename(file.fullpath))
  tmp <- read.xlsx(file.fullpath)
  
  up_count <- ifelse("U" %in% tmp[,7], table(tmp[,7])["U"], 0)
  down_count <- ifelse("D" %in% tmp[,7], table(tmp[,7])["D"], 0)
  nosign_count <- ifelse("_" %in% tmp[,7], table(tmp[,7])["_"], 0)
  
  c(Group = res_name, Up = up_count, Down = down_count, NoSign = nosign_count)
})
plan(sequential)

stat.df <- do.call(rbind, stat.list) %>% as.data.frame()
stat.df$Up <- as.numeric(stat.df$Up)
stat.df$Down <- as.numeric(stat.df$Down)
stat.df$NoSign <- as.numeric(stat.df$NoSign)

stat.df2 <- stat.df %>% arrange(Group) %>%
  tidyr::separate(Group, into = c("Case", "Control"), sep = "\\.vs\\.Ctrl-") %>%
  dplyr::select(Case, Control, Up, Down, NoSign)

write.xlsx(stat.df2, file.path(outdir, "06.diffAbund", "Diff_TaxonNum_Statistics.xlsx"), colNames=TRUE, rowNames=FALSE)

cat("\n[INFO] Analysis completed successfully for cutoff:", cutoff, "\n")
