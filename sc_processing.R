library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(patchwork)

#import stuffs
mel_dat <- Read10X(data.dir="/home/sc_liftoff/dmel_testis_fb/outs/filtered_feature_bc_matrix")
write.csv(mel_dat, file="/home/sc_liftoff/bursting/dmel_testis_fb_raw.csv")
yak_dat <- Read10X(data.dir="/home/sc_liftoff/dyak_testis_fb/outs/filtered_feature_bc_matrix")
write.csv(yak_dat, file="/home/sc_liftoff/bursting/dyak_testis_fb_raw.csv")
ana_dat <- Read10X(data.dir="/home/sc_liftoff/dana_testis_fb/outs/filtered_feature_bc_matrix")
write.csv(ana_dat, file="/home/sc_liftoff/bursting/dana_testis_fb_raw.csv")

#"melanogasterize genome"
orthos_yak <- read.csv("/home/sc_liftoff/genomes_fb/dmel_orthos_yak.tsv", skip=5, header=F, sep="\t")
orthos_ana <- read.csv("/home/sc_liftoff/genomes_fb/dmel_orthos_ana.tsv", skip=5, header=F, sep="\t")

genes_yak <- unique(orthos_yak[,2])
genes_yak_uniq <- genes_yak[unlist(lapply(genes_yak, function(x) sum(orthos_yak[,2] == x) == 1))]
genes_yak_yak <- orthos_yak[which(orthos_yak[,2] %in% genes_yak_uniq), 7]
genes_yak_uniq_yak <- genes_yak_yak[unlist(lapply(genes_yak_yak, function(x) sum(orthos_yak[,7] == x) == 1))]
genes_yak_one2one <- intersect(genes_yak_uniq, orthos_yak[which(orthos_yak[,7] %in% genes_yak_uniq_yak), 2])

genes_ana <- unique(orthos_ana[,2])
genes_ana_uniq <- genes_ana[unlist(lapply(genes_ana, function(x) sum(orthos_ana[,2] == x) == 1))]
genes_ana_ana <- orthos_ana[which(orthos_ana[,2] %in% genes_ana_uniq), 7]
genes_ana_uniq_ana <- genes_ana_ana[unlist(lapply(genes_ana_ana, function(x) sum(orthos_ana[,7] == x) == 1))]
genes_ana_one2one <- intersect(genes_ana_uniq, orthos_ana[which(orthos_ana[,7] %in% genes_ana_uniq_ana), 2])

orthos_yak_uniq <- orthos_yak[which(orthos_yak[,2] %in% genes_yak_one2one),]
orthos_ana_uniq <- orthos_ana[which(orthos_ana[,2] %in% genes_ana_one2one),]

genes_shared <- intersect(genes_yak_one2one, genes_ana_one2one)

genes_shared <- intersect(genes_yak_one2one, genes_ana_one2one)
genes_shared <- genes_shared[!grepl("Rp[A-Z]", genes_shared)]
genes_shared <- genes_shared[!grepl("mt:", genes_shared)]
genes_shared <- genes_shared[!grepl("tRNA", genes_shared)]

genes_shared_fbgn_mel <- orthos_ana_uniq[which(orthos_ana_uniq[,2] %in% genes_shared), 1]
genes_shared_fbgn_yak <- orthos_yak_uniq[which(orthos_yak_uniq[,2] %in% genes_shared), 6]
genes_shared_fbgn_ana <- orthos_ana_uniq[which(orthos_ana_uniq[,2] %in% genes_shared), 6]

mel_dat_shared <- mel_dat[genes_shared_fbgn_mel,]
yak_dat_shared <- yak_dat[genes_shared_fbgn_yak,]
ana_dat_shared <- ana_dat[genes_shared_fbgn_ana,]

rownames(mel_dat_shared) <- genes_shared
rownames(yak_dat_shared) <- genes_shared
rownames(ana_dat_shared) <- genes_shared

genes_converter <- genes_shared_fbgn_mel
names(genes_converter) <- genes_shared


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col_wheel <- gg_color_hue(8)

ct_cols <- col_wheel[c(7, 3, 6, 2, 5, 1, 4, 8)]


hp6_coexprs_clust <- c("CG11929", "Elba3", "CG3251", "Taf12L", "CG15631")
mel_coexprs_fbgn <- orthos_yak[match(hp6_coexprs_clust, orthos_yak[,2]), 1]
yak_coexprs_fbgn <- orthos_yak[match(hp6_coexprs_clust, orthos_yak[,2]), 6]
ana_coexprs_fbgn <- orthos_ana[match(hp6_coexprs_clust, orthos_ana[,2]), 6]

mel_hp6_fbgn <- "FBgn0031613"
yak_hp6_fbgn <- orthos_yak[which(orthos_yak[,2] == "HP6"), 6]

#198-gene list from preprint
gene_list <- read.csv("sc_liftoff/reduced_gene_list_monocle.txt", header=F)[,1]

gene_list_mel_fbgn <- genes_converter[gene_list]
gene_list_yak_fbgn <- orthos_yak[match(gene_list, orthos_yak[,2]), 6]
gene_list_ana_fbgn <- orthos_ana[match(gene_list, orthos_ana[,2]), 6]

#just stick on HP6 because it isn't conserved
mel_dat_shared_hp6 <- rbind(mel_dat_shared, mel_dat[mel_hp6_fbgn,])
rownames(mel_dat_shared_hp6)[nrow(mel_dat_shared_hp6)] <- "HP6"

yak_dat_shared_hp6 <- rbind(yak_dat_shared, yak_dat[yak_hp6_fbgn,])
rownames(yak_dat_shared_hp6)[nrow(yak_dat_shared_hp6)] <- "HP6"

#blank ana HP6 because it doesn't exist
ana_dat_shared_hp6 <- rbind(ana_dat_shared, rep(0, 5000))
rownames(ana_dat_shared_hp6)[nrow(ana_dat_shared_hp6)] <- "HP6"

#total gene list is 198-gene list + coexpression cluster genes + HP6/Umbrea
gene_list_hp6 <- c(gene_list, hp6_coexprs_clust, "HP6")

#combining everything just to make seurat happy, but making sure to have min.cells and min.features=0 to keep all read counts
all_combined <- CreateSeuratObject(counts=cbind(cbind(mel_dat_shared_hp6, yak_dat_shared_hp6), ana_dat_shared_hp6), project="all_comb", min.cells=0, min.features=0)

#using read counts 5eva
DefaultAssay(all_combined) <- "RNA"

all_combined <- all_combined[gene_list_hp6,]

#scale otherwise it throws a fit
all_combined <- ScaleData(all_combined)

#read in cell type assignments from preprint
ct_assign <- read.csv(file="sc_liftoff/ct_assignments.csv", row.names="X")
all_combined <- AddMetaData(all_combined, factor(ct_assign$spec, levels=c("mel", "yak", "ana")), "spec")
all_combined <- AddMetaData(all_combined, factor(ct_assign$reclu_coarse, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid", "ananassae spermatid")), "reclu_coarse")
Idents(all_combined) <- "reclu_coarse"

#split it again
all_combined_list <- SplitObject(all_combined, split.by="spec")

#do the things
all_combined.list <- lapply(X = all_combined_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst")#, nfeatures = 10731)
	x <- ScaleData(x)
})

#integrate the data set to make seurat happy, but basically only retain the gene_list
all_combined_feats <- SelectIntegrationFeatures(object.list=all_combined.list, anchor.features=gene_list_hp6)
all_combined_anchors <- FindIntegrationAnchors(object.list=all_combined.list, anchor.features=gene_list_hp6)
all_combined.combined <- IntegrateData(anchorset = all_combined_anchors)

all_combined.combined <- AddMetaData(all_combined.combined, factor(ct_assign$spec, levels=c("mel", "yak", "ana")), "spec")
all_combined.combined <- AddMetaData(all_combined.combined, factor(ct_assign$reclu_coarse, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid", "ananassae spermatid")), "reclu_coarse")
Idents(all_combined.combined) <- "reclu_coarse"

#scaling over *all* reads, but then downstream PCA is on only 198-gene list
all_combined.combined <- ScaleData(all_combined.combined)
all_combined.combined <- RunPCA(all_combined.combined, weight.by.var=FALSE, features = gene_list)
#settings as used in Monocle3 for n.neighbors and min.dist, also dims based on elbow plot
all_combined.combined <- RunUMAP(all_combined.combined, dims=1:8, n.neighbors=15L, min.dist=0.1)
#forget about ananassae spermatid, because not relevant here and only just confusing
all_combined.no_ansperm <- subset(all_combined.combined, reclu_coarse %in% c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))

svg("sc_liftoff/hp6_elbow.svg", width=6, height=4, pointsize=3)
ElbowPlot(all_combined.combined, ndims = 20, reduction = "pca")
dev.off()

png("sc_liftoff/hp6_elbow.png", width=6, height=4, units="in", res=600, pointsize=3)
ElbowPlot(all_combined.combined, ndims = 20, reduction = "pca")
dev.off()


svg("sc_liftoff/hp6_pca.svg", width=6*1.25, height=4*1.25, pointsize=3)
p <- DimPlot(all_combined.no_ansperm, reduction="pca", pt.size=1e-2, shuffle=TRUE)
p
dev.off()

png("sc_liftoff/hp6_pca.png", width=6*1.25, height=4*1.25, pointsize=3, units="in", res=600)
p <- DimPlot(all_combined.no_ansperm, reduction="pca", pt.size=1e-2, shuffle=TRUE)
p
dev.off()

svg("sc_liftoff/hp6_umap.svg", width=6*1.25, height=4*1.25, pointsize=3)
p <- DimPlot(all_combined.no_ansperm, reduction="umap", pt.size=1e-2)
p
dev.off()

png("sc_liftoff/hp6_umap.png", width=6*1.25, height=4*1.25, pointsize=3, units="in", res=600)
p <- DimPlot(all_combined.no_ansperm, reduction="umap", pt.size=1e-2)
p
dev.off()

DefaultAssay(all_combined.no_ansperm) <- "RNA"

#resplit now that pca and umap have been performed for all species
all_combined.no_ansperm_spec <- SplitObject(all_combined.no_ansperm, split.by="spec")

#rest of the owl
#uses raw counts!
svg("sc_liftoff/hp6_coexprs_clust.svg", width=4*6, height=4*3, pointsize=3)

p1 <- VlnPlot(all_combined.no_ansperm_spec[[1]], features=c("HP6", hp6_coexprs_clust), pt.size=0, ncol=6, cols=ct_cols, slot="counts", log=TRUE)
p1 <- p1 & theme(text=element_text(family="Roboto"))

p2 <- VlnPlot(all_combined.no_ansperm_spec[[2]], features=c("HP6", hp6_coexprs_clust), pt.size=0, ncol=6, cols=ct_cols, slot="counts", log=TRUE)
p2 <- p2 & theme(text=element_text(family="Roboto"))

p3 <- VlnPlot(all_combined.no_ansperm_spec[[3]], features=c("HP6", hp6_coexprs_clust), pt.size=0, ncol=6, cols=ct_cols, slot="counts", log=TRUE)
p3 <- p3 & theme(text=element_text(family="Roboto"))

p1/p2/p3

dev.off()


png("sc_liftoff/hp6_coexprs_clust.png", width=4*6, height=4*3, pointsize=3, units="in", res=600)

p1 <- VlnPlot(all_combined.no_ansperm_spec[[1]], features=c("HP6", hp6_coexprs_clust), pt.size=0, ncol=6, cols=ct_cols, slot="counts", log=TRUE)
p1 <- p1 & theme(text=element_text(family="Roboto"))

p2 <- VlnPlot(all_combined.no_ansperm_spec[[2]], features=c("HP6", hp6_coexprs_clust), pt.size=0, ncol=6, cols=ct_cols, slot="counts", log=TRUE)
p2 <- p2 & theme(text=element_text(family="Roboto"))

p3 <- VlnPlot(all_combined.no_ansperm_spec[[3]], features=c("HP6", hp6_coexprs_clust), pt.size=0, ncol=6, cols=ct_cols, slot="counts", log=TRUE)
p3 <- p3 & theme(text=element_text(family="Roboto"))

p1/p2/p3

dev.off()

#idk, not really used, but makes seurat happy
DefaultAssay(all_combined.no_ansperm_spec[[1]]) <- "integrated"
DefaultAssay(all_combined.no_ansperm_spec[[2]]) <- "integrated"
DefaultAssay(all_combined.no_ansperm_spec[[3]]) <- "integrated"

#uses raw counts!
svg("sc_liftoff/hp6_coexprs_clust_umap.svg", width=4*6, height=4*3, pointsize=3)

p1 <- FeaturePlot(all_combined.no_ansperm_spec[[1]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q95", reduction="umap")
p1 <- p1 & theme(text=element_text(family="Roboto"), legend.position="none")

p2 <- FeaturePlot(all_combined.no_ansperm_spec[[2]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q95", reduction="umap")
p2 <- p2 & theme(text=element_text(family="Roboto"), legend.position="none")

p3 <- FeaturePlot(all_combined.no_ansperm_spec[[3]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q97", reduction="umap")
p3 <- p3 & theme(text=element_text(family="Roboto"), legend.position="none")

p1/p2/p3

dev.off()


png("sc_liftoff/hp6_coexprs_clust_umap.png", width=4*6, height=4*3, pointsize=3, units="in", res=600)

p1 <- FeaturePlot(all_combined.no_ansperm_spec[[1]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q95", reduction="umap")
p1 <- p1 & theme(text=element_text(family="Roboto"), legend.position="none")

p2 <- FeaturePlot(all_combined.no_ansperm_spec[[2]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q95", reduction="umap")
p2 <- p2 & theme(text=element_text(family="Roboto"), legend.position="none")

p3 <- FeaturePlot(all_combined.no_ansperm_spec[[3]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q97", reduction="umap")
p3 <- p3 & theme(text=element_text(family="Roboto"), legend.position="none")

p1/p2/p3

dev.off()



svg("sc_liftoff/hp6_coexprs_clust_pca.svg", width=4*6, height=4*3, pointsize=3)

p1 <- FeaturePlot(all_combined.no_ansperm_spec[[1]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q95", reduction="pca", combine=FALSE)
p1 <- lapply(p1, function(x){
			x <- x + theme(text=element_text(family="Roboto"), legend.position="none") + xlim(c(-0.0217, 0.0151)) + ylim(c(-0.0378, 0.0122))
			return(x)
		})
p1 <- p1[[1]] + p1[[2]] + p1[[3]] + p1[[4]] + p1[[5]] + p1[[6]] + plot_layout(ncol=6)

p2 <- FeaturePlot(all_combined.no_ansperm_spec[[2]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q95", reduction="pca", combine=FALSE)
p2 <- lapply(p2, function(x){
			x <- x + theme(text=element_text(family="Roboto"), legend.position="none") + xlim(c(-0.0217, 0.0151)) + ylim(c(-0.0378, 0.0122))
			return(x)
		})
p2 <- p2[[1]] + p2[[2]] + p2[[3]] + p2[[4]] + p2[[5]] + p2[[6]] + plot_layout(ncol=6)


p3 <- FeaturePlot(all_combined.no_ansperm_spec[[3]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q97", reduction="pca", combine=FALSE)
p3 <- lapply(p3, function(x){
			x <- x + theme(text=element_text(family="Roboto"), legend.position="none") + xlim(c(-0.0217, 0.0151)) + ylim(c(-0.0378, 0.0122))
			return(x)
		})
p3 <- p3[[1]] + p3[[2]] + p3[[3]] + p3[[4]] + p3[[5]] + p3[[6]] + plot_layout(ncol=6)

p1/p2/p3

dev.off()


png("sc_liftoff/hp6_coexprs_clust_pca.png", width=4*6, height=4*3, pointsize=3, units="in", res=600)

p1 <- FeaturePlot(all_combined.no_ansperm_spec[[1]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q95", reduction="pca", combine=FALSE)
p1 <- lapply(p1, function(x){
			x <- x + theme(text=element_text(family="Roboto"), legend.position="none") + xlim(c(-0.0217, 0.0151)) + ylim(c(-0.0378, 0.0122))
			return(x)
		})
p1 <- p1[[1]] + p1[[2]] + p1[[3]] + p1[[4]] + p1[[5]] + p1[[6]] + plot_layout(ncol=6)

p2 <- FeaturePlot(all_combined.no_ansperm_spec[[2]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q95", reduction="pca", combine=FALSE)
p2 <- lapply(p2, function(x){
			x <- x + theme(text=element_text(family="Roboto"), legend.position="none") + xlim(c(-0.0217, 0.0151)) + ylim(c(-0.0378, 0.0122))
			return(x)
		})
p2 <- p2[[1]] + p2[[2]] + p2[[3]] + p2[[4]] + p2[[5]] + p2[[6]] + plot_layout(ncol=6)


p3 <- FeaturePlot(all_combined.no_ansperm_spec[[3]], features=c("HP6", hp6_coexprs_clust), ncol=6, slot="counts", max.cutoff="q97", reduction="pca", combine=FALSE)
p3 <- lapply(p3, function(x){
			x <- x + theme(text=element_text(family="Roboto"), legend.position="none") + xlim(c(-0.0217, 0.0151)) + ylim(c(-0.0378, 0.0122))
			return(x)
		})
p3 <- p3[[1]] + p3[[2]] + p3[[3]] + p3[[4]] + p3[[5]] + p3[[6]] + plot_layout(ncol=6)

p1/p2/p3

dev.off()

