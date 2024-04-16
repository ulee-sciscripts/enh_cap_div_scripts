##load("everything.rdata")

#load("rpkm.mat.trim.rdata")
#load("rpkm.mat.tiss.rdata")

#baseline_cor <- function(dat_mat, col_idx, cor_method, size){
#	results <- matrix(nrow=size, ncol=3)

#	results <- apply(results, 1, function(x) {
#		cor_val <- NA
#		
#		while(is.na(cor_val)){
#			idxs <- sample(1:nrow(dat_mat), 2, F)
#			cor_val <- cor(dat_mat[idxs[1], col_idx], dat_mat[idxs[2], col_idx], method=cor_method)
##			print(c(idxs, cor_val))
#		}
#		return(c(idxs, cor_val))
#	})
#	
#	rownames(results) <- c("gene1", "gene2", "cor")

#	return(t(results))	
#}

#load("data_source_list.rdata")

#baseline_dist_spear_tiss <- baseline_cor(rpkm_mat_tiss, modENCODE_mRNA_Seq_tissues_list, "spearman", 10000)

#baseline_dist_ken_tiss <- baseline_cor(rpkm_mat_tiss, modENCODE_mRNA_Seq_tissues_list, "kendall", 10000)

#baseline_dist_spear_all <- baseline_cor(rpkm_mat_trim, 1:124, "spearman", 10000)

#baseline_dist_ken_all <- baseline_cor(rpkm_mat_trim, 1:124, "kendall", 10000)

#bootstrap_cor <- function(dat_mat, col_idx, cor_method, size){
#	result <- c()
#	for(i in 1:size){
#		inner_mat <- matrix(nrow=2, ncol=length(col_idx))		
#		colnames(inner_mat) <- col_idx
#		for(j in col_idx){
#			inner_mat[1:2,j] <- sample(dat_mat[,j], 2, replace=T)
#		}
#		result <- c(result, cor(inner_mat[1,], inner_mat[2,], method=cor_method))
#	}
#	return(result)
#}

#bootstrap_dist_spear_tiss <- bootstrap_cor(rpkm_mat_tiss, modENCODE_mRNA_Seq_tissues_list, "spearman", 10000)
#bootstrap_dist_ken_tiss <- bootstrap_cor(rpkm_mat_tiss, modENCODE_mRNA_Seq_tissues_list, "kendall", 10000)
#bootstrap_dist_spear_all <- bootstrap_cor(rpkm_mat_trim, 1:124, "spearman", 10000)
#bootstrap_dist_ken_all <- bootstrap_cor(rpkm_mat_trim, 1:124, "kendall", 10000)

#svg("hist_baseline_tiss.svg")
#hist(baseline_dist_spear_tiss[,3], main="Distribution of Random Pairs of Genes (Spearman)")
#dev.off()
#png("hist_baseline_tiss.png")
#hist(baseline_dist_spear_tiss[,3], main="Distribution of Random Pairs of Genes (Spearman)")
#dev.off()

#svg("hist_bootstrap_tiss.svg")
#hist(bootstrap_dist_spear_tiss, main="Distribution of Bootstrapped Pairs (Spearman)")
#dev.off()
#png("hist_bootstrap_tiss.png")
#hist(bootstrap_dist_spear_tiss, main="Distribution of Bootstrapped Pairs (Spearman)")
#dev.off()

#svg("qq_boot_vs_random.svg")
#qqplot(bootstrap_dist_spear_tiss, baseline_dist_spear_tiss[,3], main="QQ Plot of Correlation Distributions", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()
#png("qq_boot_vs_random.png")
#qqplot(bootstrap_dist_spear_tiss, baseline_dist_spear_tiss[,3], main="QQ Plot of Correlation Distributions", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()

#cor_mtx_spear <- matrix(0, nrow=length(modENCODE_mRNA_Seq_tissues_list), ncol=length(modENCODE_mRNA_Seq_tissues_list))
#cor_mtx_ken <- matrix(0, nrow=length(modENCODE_mRNA_Seq_tissues_list), ncol=length(modENCODE_mRNA_Seq_tissues_list))

#for(i in 1:length(modENCODE_mRNA_Seq_tissues_list)){
#	for(j in i:length(modENCODE_mRNA_Seq_tissues_list)){
#		cor_mtx_spear[i,j] <- cor(rpkm_mat_tiss[,modENCODE_mRNA_Seq_tissues_list[i]], rpkm_mat_tiss[,modENCODE_mRNA_Seq_tissues_list[j]], method="spearman")
#		cor_mtx_ken[i,j] <- cor(rpkm_mat_tiss[,modENCODE_mRNA_Seq_tissues_list[i]], rpkm_mat_tiss[,modENCODE_mRNA_Seq_tissues_list[j]], method="kendall")
#	
#	}
#}

#rownames(cor_mtx_spear) <- modENCODE_mRNA_Seq_tissues_list
 #rownames(cor_mtx_ken) <- modENCODE_mRNA_Seq_tissues_list
 
# cor_mtx_spear <- cor_mtx_spear + t(cor_mtx_spear)
# cor_mtx_ken <- cor_mtx_ken + t(cor_mtx_ken)
# for(i in 1:length(modENCODE_mRNA_Seq_tissues_list)){
# 	cor_mtx_spear[i,i] <- 1
# 	cor_mtx_ken[i,i] <- 1
# }


#library(gplots)
#svg("heatmap_cor_raw_spear.svg", width=14, height=14)
#heatmap.2(cor_mtx_spear, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues (Spearman)")
#dev.off()
#png("heatmap_cor_raw_spear.png", width=480*2, height=480*2)
#heatmap.2(cor_mtx_spear, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues (Spearman)")
#dev.off()

#svg("heatmap_cor_raw_ken.svg", width=14, height=14)
#heatmap.2(cor_mtx_ken, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues (Kendall)")
#dev.off()
#png("heatmap_cor_raw_ken.png", width=480*2, height=480*2)
#heatmap.2(cor_mtx_ken, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues (Kendall)")
#dev.off()

#take PCs for correlated groups, extract 1st PC, quantile-normalize to distribution of contained groups
#collapse head exprs, ovaries, 1+4 day carcass, and 20+4 day dig_sys using PCA
#conds_head <- c("mE_mRNA_A_MateF_4d_head", "mE_mRNA_A_VirF_4d_head", "mE_mRNA_A_MateM_4d_head", "mE_mRNA_A_MateF_1d_head", "mE_mRNA_A_VirF_1d_head", "mE_mRNA_A_MateM_1d_head", "mE_mRNA_A_MateF_20d_head", "mE_mRNA_A_VirF_20d_head", "mE_mRNA_A_MateM_20d_head")
#conds_ovar <- c("mE_mRNA_A_VirF_4d_ovary", "mE_mRNA_A_MateF_4d_ovary")
#conds_carc <- c("mE_mRNA_A_4d_carcass", "mE_mRNA_A_1d_carcass")
#conds_dige <- c("mE_mRNA_A_20d_dig_sys", "mE_mRNA_A_4d_dig_sys")

#pc_head <- prcomp(t(rpkm_mat_tiss[,conds_head]))$rotation[,1]
#pc_ovar <- prcomp(t(rpkm_mat_tiss[,conds_ovar]))$rotation[,1]
#pc_carc <- prcomp(t(rpkm_mat_tiss[,conds_carc]))$rotation[,1]
#pc_dige <- prcomp(t(rpkm_mat_tiss[,conds_dige]))$rotation[,1]

###calculate inverse of empricial CDF using all data from each condition
#i_eCDF_head <- quantile(as.vector(rpkm_mat_tiss[,conds_head]), type=1, probs=seq(0, 1, length.out=nrow(rpkm_mat_tiss)))
#i_eCDF_ovar <- quantile(as.vector(rpkm_mat_tiss[,conds_ovar]), type=1, probs=seq(0, 1, length.out=nrow(rpkm_mat_tiss)))
#i_eCDF_carc <- quantile(as.vector(rpkm_mat_tiss[,conds_carc]), type=1, probs=seq(0, 1, length.out=nrow(rpkm_mat_tiss)))
#i_eCDF_dige <- quantile(as.vector(rpkm_mat_tiss[,conds_dige]), type=1, probs=seq(0, 1, length.out=nrow(rpkm_mat_tiss)))

#conds_head_eCDF <- rank(pc_head)
#conds_ovar_eCDF <- rank(pc_ovar)
#conds_carc_eCDF <- rank(pc_carc)
#conds_dige_eCDF <- rank(pc_dige)

#conds_head_qnorm <- pc_head
#conds_ovar_qnorm <- pc_ovar
#conds_carc_qnorm <- pc_carc
#conds_dige_qnorm <- pc_dige

#for(i in 1:nrow(rpkm_mat_tiss)){
#	conds_head_qnorm[i] <- i_eCDF_head[conds_head_eCDF[i]]
#	conds_ovar_qnorm[i] <- i_eCDF_ovar[conds_ovar_eCDF[i]]
#	conds_carc_qnorm[i] <- i_eCDF_carc[conds_carc_eCDF[i]]
#	conds_dige_qnorm[i] <- i_eCDF_dige[conds_dige_eCDF[i]]
#}

#rpkm_mat_tiss_qnorm <- rpkm_mat_tiss[,-which(colnames(rpkm_mat_tiss) %in% c(conds_head,conds_ovar,conds_carc,conds_dige))]
#rpkm_mat_tiss_qnorm <- cbind(rpkm_mat_tiss_qnorm, conds_head_qnorm)
#rpkm_mat_tiss_qnorm <- cbind(rpkm_mat_tiss_qnorm, conds_ovar_qnorm)
#rpkm_mat_tiss_qnorm <- cbind(rpkm_mat_tiss_qnorm, conds_carc_qnorm)
#rpkm_mat_tiss_qnorm <- cbind(rpkm_mat_tiss_qnorm, conds_dige_qnorm)

#colnames(rpkm_mat_tiss_qnorm)[(ncol(rpkm_mat_tiss_qnorm)-3):ncol(rpkm_mat_tiss_qnorm)] <- c("head_pc", "ovar_pc", "carc_pc", "dige_pc")

#modENCODE_mRNA_Seq_tissues_list_pc <- modENCODE_mRNA_Seq_tissues_list[-which(modENCODE_mRNA_Seq_tissues_list %in% c(conds_head, conds_ovar, conds_carc, conds_dige))]
#modENCODE_mRNA_Seq_tissues_list_pc <- c(modENCODE_mRNA_Seq_tissues_list_pc, c("head_pc", "ovar_pc", "carc_pc", "dige_pc"))

#save(list=c("rpkm_mat_tiss_qnorm", "modENCODE_mRNA_Seq_tissues_list_pc"), file="rpkm.mat.tiss.qnorm.rdata")

#bootstrap_dist_spear_tiss_qnorm <- bootstrap_cor(rpkm_mat_tiss_qnorm, modENCODE_mRNA_Seq_tissues_list_pc, "spearman", 10000)
#bootstrap_dist_ken_tiss_qnorm <- bootstrap_cor(rpkm_mat_tiss_qnorm, modENCODE_mRNA_Seq_tissues_list_pc, "kendall", 10000)

#baseline_dist_spear_tiss_qnorm <- baseline_cor(rpkm_mat_tiss_qnorm, modENCODE_mRNA_Seq_tissues_list_pc, "spearman", 10000)

#baseline_dist_ken_tiss_qnorm <- baseline_cor(rpkm_mat_tiss_qnorm, modENCODE_mRNA_Seq_tissues_list_pc, "kendall", 10000)

#cor_mtx_spear_pc <- matrix(0, nrow=length(modENCODE_mRNA_Seq_tissues_list_pc), ncol=length(modENCODE_mRNA_Seq_tissues_list_pc))
#cor_mtx_ken_pc <- matrix(0, nrow=length(modENCODE_mRNA_Seq_tissues_list_pc), ncol=length(modENCODE_mRNA_Seq_tissues_list_pc))

#for(i in 1:length(modENCODE_mRNA_Seq_tissues_list_pc)){
#	for(j in i:length(modENCODE_mRNA_Seq_tissues_list_pc)){
#		cor_mtx_spear_pc[i,j] <- cor(rpkm_mat_tiss_qnorm[,modENCODE_mRNA_Seq_tissues_list_pc[i]], rpkm_mat_tiss_qnorm[,modENCODE_mRNA_Seq_tissues_list_pc[j]], method="spearman")
#		cor_mtx_ken_pc[i,j] <- cor(rpkm_mat_tiss_qnorm[,modENCODE_mRNA_Seq_tissues_list_pc[i]], rpkm_mat_tiss_qnorm[,modENCODE_mRNA_Seq_tissues_list_pc[j]], method="kendall")
#	
#	}
#}

#rownames(cor_mtx_spear_pc) <- modENCODE_mRNA_Seq_t#svg("hist_baseline_tiss_pc.svg")
#hist(baseline_dist_spear_tiss_qnorm[,3], main="Distribution of Random Pairs of Genes (PC, Spearman)")
#dev.off()
#png("hist_baseline_tiss_pc.png")
#hist(baseline_dist_spear_tiss_qnorm[,3], main="Distribution of Random Pairs of Genes (PC, Spearman)")
#dev.off()

#svg("hist_bootstrap_tiss_pc.svg")
#hist(bootstrap_dist_spear_tiss_qnorm, main="Distribution of Bootstrapped Pairs (PC, Spearman)")
#dev.off()
#png("hist_bootstrap_tiss_pc.png")
#hist(bootstrap_dist_spear_tiss_qnorm, main="Distribution of Bootstrapped Pairs (PC, Spearman)")
#dev.off()

#svg("hist_baseline_tiss_avg.svg")
#hist(baseline_dist_spear_tiss_avg[,3], main="Distribution of Random Pairs of Genes (Avg, Spearman)")
#dev.off()
#png("hist_baseline_tiss_avg.png")
#hist(baseline_dist_spear_tiss_avg[,3], main="Distribution of Random Pairs of Genes (Avg, Spearman)")
#dev.off()

#svg("hist_bootstrap_tiss_avg.svg")
#hist(bootstrap_dist_spear_tiss_avg, main="Distribution of Bootstrapped Pairs (Avg, Spearman)")
#dev.off()
#png("hist_bootstrap_tiss_avg.png")
#hist(bootstrap_dist_spear_tiss_avg, main="Distribution of Bootstrapped Pairs (Avg, Spearman)")
#dev.off()

#svg("qq_boot_vs_random_pc.svg")
#qqplot(bootstrap_dist_spear_tiss_qnorm, baseline_dist_spear_tiss_qnorm[,3], main="QQ Plot of Correlation Distributions (PC)", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()
#png("qq_boot_vs_random_pc.png")
#qqplot(bootstrap_dist_spear_tiss_qnorm, baseline_dist_spear_tiss_qnorm[,3], main="QQ Plot of Correlation Distributions (PC)", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()

#svg("qq_boot_vs_random_avg.svg")
#qqplot(bootstrap_dist_spear_tiss_avg, baseline_dist_spear_tiss_avg[,3], main="QQ Plot of Correlation Distributions (Avg)", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()
#png("qq_boot_vs_random_avg.png")
#qqplot(bootstrap_dist_spear_tiss_avg, baseline_dist_spear_tiss_avg[,3], main="QQ Plot of Correlation Distributions (Avg)", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()issues_list_pc
#rownames(cor_mtx_ken_pc) <- modENCODE_mRNA_Seq_tissues_list_pc
#colnames(cor_mtx_spear_pc) <- modENCODE_mRNA_Seq_tissues_list_pc
#colnames(cor_mtx_ken_pc) <- modENCODE_mRNA_Seq_tissues_list_pc
# 
# cor_mtx_spear_pc <- cor_mtx_spear_pc + t(cor_mtx_spear_pc)
# cor_mtx_ken_pc <- cor_mtx_ken_pc + t(cor_mtx_ken_pc)
# for(i in 1:length(modENCODE_mRNA_Seq_tissues_list_pc)){
# 	cor_mtx_spear_pc[i,i] <- 1
# 	cor_mtx_ken_pc[i,i] <- 1
# }

#svg("heatmap_cor_raw_spear_pc.svg", width=14, height=14)
#heatmap.2(cor_mtx_spear_pc, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues, PC-reduced (Spearman)")
#dev.off()
#png("heatmap_cor_raw_spear_pc.png", width=480*2, height=480*2)
#heatmap.2(cor_mtx_spear_pc, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues, PC-reduced (Spearman)")
#dev.off()

#svg("heatmap_cor_raw_ken_pc.svg", width=14, height=14)
#heatmap.2(cor_mtx_ken_pc, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues, PC-reduced (Kendall)")
#dev.off()
#png("heatmap_cor_raw_ken_pc.png", width=480*2, height=480*2)
#heatmap.2(cor_mtx_ken_pc, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues, PC-reduced (Kendall)")
#dev.off()

#cor_pc_conds_head <- sapply(conds_head, function(x) cor(rpkm_mat_tiss[,x], rpkm_mat_tiss_qnorm[,15], method="spearman"))
#cor_pc_conds_ovar <- sapply(conds_ovar, function(x) cor(rpkm_mat_tiss[,x], rpkm_mat_tiss_qnorm[,16], method="spearman"))
#cor_pc_conds_carc <- sapply(conds_carc, function(x) cor(rpkm_mat_tiss[,x], rpkm_mat_tiss_qnorm[,17], method="spearman"))
#cor_pc_conds_dige <- sapply(conds_dige, function(x) cor(rpkm_mat_tiss[,x], rpkm_mat_tiss_qnorm[,18], method="spearman"))

#avg_head <- rowMeans(rpkm_mat_tiss[,conds_head])
#avg_dige <- rowMeans(rpkm_mat_tiss[,conds_dige])
#avg_ovar <- rowMeans(rpkm_mat_tiss[,conds_ovar])
#avg_carc <- rowMeans(rpkm_mat_tiss[,conds_carc])

#rpkm_mat_tiss_avg <- rpkm_mat_tiss[,-which(colnames(rpkm_mat_tiss) %in% c(conds_head,conds_ovar,conds_carc,conds_dige))]
#rpkm_mat_tiss_avg <- cbind(rpkm_mat_tiss_avg, avg_head)
#rpkm_mat_tiss_avg <- cbind(rpkm_mat_tiss_avg, avg_ovar)
#rpkm_mat_tiss_avg <- cbind(rpkm_mat_tiss_avg, avg_carc)
#rpkm_mat_tiss_avg <- cbind(rpkm_mat_tiss_avg, avg_dige)

#colnames(rpkm_mat_tiss_avg)[(ncol(rpkm_mat_tiss_avg)-3):ncol(rpkm_mat_tiss_avg)] <- c("head_avg", "ovar_avg", "carc_avg", "dige_avg")

#modENCODE_mRNA_Seq_tissues_list_avg <- colnames(rpkm_mat_tiss_avg)

#save(list=c("rpkm_mat_tiss_avg", "modENCODE_mRNA_Seq_tissues_list_avg"), file="rpkm.mat.tiss.avg.rdata")

#bootstrap_dist_spear_tiss_avg <- bootstrap_cor(rpkm_mat_tiss_avg, modENCODE_mRNA_Seq_tissues_list_avg, "spearman", 10000)
#bootstrap_dist_ken_tiss_avg <- bootstrap_cor(rpkm_mat_tiss_avg, modENCODE_mRNA_Seq_tissues_list_avg, "kendall", 10000)

#baseline_dist_spear_tiss_avg <- baseline_cor(rpkm_mat_tiss_avg, modENCODE_mRNA_Seq_tissues_list_avg, "spearman", 10000)

#baseline_dist_ken_tiss_avg <- baseline_cor(rpkm_mat_tiss_avg, modENCODE_mRNA_Seq_tissues_list_avg, "kendall", 10000)

#cor_mtx_spear_avg <- matrix(0, nrow=length(modENCODE_mRNA_Seq_tissues_list_avg), ncol=length(modENCODE_mRNA_Seq_tissues_list_avg))
#cor_mtx_ken_avg <- matrix(0, nrow=length(modENCODE_mRNA_Seq_tissues_list_avg), ncol=length(modENCODE_mRNA_Seq_tissues_list_avg))

#for(i in 1:length(modENCODE_mRNA_Seq_tissues_list_avg)){
#	for(j in i:length(modENCODE_mRNA_Seq_tissues_list_avg)){
#		cor_mtx_spear_avg[i,j] <- cor(rpkm_mat_tiss_avg[,modENCODE_mRNA_Seq_tissues_list_avg[i]], rpkm_mat_tiss_avg[,modENCODE_mRNA_Seq_tissues_list_avg[j]], method="spearman")
#		cor_mtx_ken_avg[i,j] <- cor(rpkm_mat_tiss_avg[,modENCODE_mRNA_Seq_tissues_list_avg[i]], rpkm_mat_tiss_avg[,modENCODE_mRNA_Seq_tissues_list_avg[j]], method="kendall")
#	
#	}
#}

#rownames(cor_mtx_spear_avg) <- modENCODE_mRNA_Seq_tissues_list_avg
#rownames(cor_mtx_ken_avg) <- modENCODE_mRNA_Seq_tissues_list_avg
#colnames(cor_mtx_spear_avg) <- modENCODE_mRNA_Seq_tissues_list_avg
#colnames(cor_mtx_ken_avg) <- modENCODE_mRNA_Seq_tissues_list_avg
# 
# cor_mtx_spear_avg <- cor_mtx_spear_avg + t(cor_mtx_spear_avg)
# cor_mtx_ken_avg <- cor_mtx_ken_avg + t(cor_mtx_ken_avg)
# for(i in 1:length(modENCODE_mRNA_Seq_tissues_list_avg)){
# 	cor_mtx_spear_avg[i,i] <- 1
# 	cor_mtx_ken_avg[i,i] <- 1
# }

#svg("heatmap_cor_raw_spear_avg.svg", width=14, height=14)
#heatmap.2(cor_mtx_spear_avg, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues, Avg-reduced (Spearman)")
#dev.off()
#png("heatmap_cor_raw_spear_avg.png", width=480*2, height=480*2)
#heatmap.2(cor_mtx_spear_avg, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues, Avg-reduced (Spearman)")
#dev.off()

#svg("heatmap_cor_raw_ken_avg.svg", width=14, height=14)
#heatmap.2(cor_mtx_ken_avg, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues, Avg-reduced (Kendall)")
#dev.off()
#png("heatmap_cor_raw_ken_avg.png", width=480*2, height=480*2)
#heatmap.2(cor_mtx_ken_avg, trace="none", Colv=NA, margin=c(15,15), main="Correlation Matrix of Tissues, Avg-reduced (Kendall)")
#dev.off()

#cor_avg_conds_head <- sapply(conds_head, function(x) cor(rpkm_mat_tiss[,x], rpkm_mat_tiss_avg[,15], method="spearman"))
#cor_avg_conds_ovar <- sapply(conds_ovar, function(x) cor(rpkm_mat_tiss[,x], rpkm_mat_tiss_avg[,16], method="spearman"))
#cor_avg_conds_carc <- sapply(conds_carc, function(x) cor(rpkm_mat_tiss[,x], rpkm_mat_tiss_avg[,17], method="spearman"))
#cor_avg_conds_dige <- sapply(conds_dige, function(x) cor(rpkm_mat_tiss[,x], rpkm_mat_tiss_avg[,18], method="spearman"))

#svg("hist_baseline_tiss_pc.svg")
#hist(baseline_dist_spear_tiss_qnorm[,3], main="Distribution of Random Pairs of Genes (PC, Spearman)")
#dev.off()
#png("hist_baseline_tiss_pc.png")
#hist(baseline_dist_spear_tiss_qnorm[,3], main="Distribution of Random Pairs of Genes (PC, Spearman)")
#dev.off()

#svg("hist_bootstrap_tiss_pc.svg")
#hist(bootstrap_dist_spear_tiss_qnorm, main="Distribution of Bootstrapped Pairs (PC, Spearman)")
#dev.off()
#png("hist_bootstrap_tiss_pc.png")
#hist(bootstrap_dist_spear_tiss_qnorm, main="Distribution of Bootstrapped Pairs (PC, Spearman)")
#dev.off()

#svg("hist_baseline_tiss_avg.svg")
#hist(baseline_dist_spear_tiss_avg[,3], main="Distribution of Random Pairs of Genes (Avg, Spearman)")
#dev.off()
#png("hist_baseline_tiss_avg.png")
#hist(baseline_dist_spear_tiss_avg[,3], main="Distribution of Random Pairs of Genes (Avg, Spearman)")
#dev.off()

#svg("hist_bootstrap_tiss_avg.svg")
#hist(bootstrap_dist_spear_tiss_avg, main="Distribution of Bootstrapped Pairs (Avg, Spearman)")
#dev.off()
#png("hist_bootstrap_tiss_avg.png")
#hist(bootstrap_dist_spear_tiss_avg, main="Distribution of Bootstrapped Pairs (Avg, Spearman)")
#dev.off()

#svg("qq_boot_vs_random_pc.svg")
#qqplot(bootstrap_dist_spear_tiss_qnorm, baseline_dist_spear_tiss_qnorm[,3], main="QQ Plot of Correlation Distributions (PC)", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()
#png("qq_boot_vs_random_pc.png")
#qqplot(bootstrap_dist_spear_tiss_qnorm, baseline_dist_spear_tiss_qnorm[,3], main="QQ Plot of Correlation Distributions (PC)", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()

#svg("qq_boot_vs_random_avg.svg")
#qqplot(bootstrap_dist_spear_tiss_avg, baseline_dist_spear_tiss_avg[,3], main="QQ Plot of Correlation Distributions (Avg)", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()
#png("qq_boot_vs_random_avg.png")
#qqplot(bootstrap_dist_spear_tiss_avg, baseline_dist_spear_tiss_avg[,3], main="QQ Plot of Correlation Distributions (Avg)", xlab="Bootstrap", ylab="Random Gene Pairs")
#dev.off()


#load("ordered_map_table_tiss_fix.rdata")
#load("rpkm.mat.tiss.rdata")

#load("rpkm.mat.tiss.qnorm.rdata")

#load("rpkm.mat.tiss.avg.rdata")

#new_gene_list <- read.csv("new_genes.csv", sep=";")

neighbor_cor <- function(dat_mat, col_idx, cor_method, gene_list, ord_map){
	

	results <- lapply(gene_list, function(x){
		x <- as.character(x)

		res_mat <- matrix(NA, nrow=4, ncol=6)
		
		if(x %in% ord_map[,1]){
			ord_map_idx <- which(ord_map[,1] == x)
			dat_mat_idx <- which(rownames(dat_mat) == ord_map[ord_map_idx ,2])
			new_gene_chr <- ord_map[ord_map_idx, "chr"]
		

			res_mat[,1] <- rep(rownames(dat_mat)[dat_mat_idx], 4)
			res_mat[,2] <- rep(x, 4)
		
			idx <-ord_map_idx - c(2, 1, -1, -2)
#			print(c(x, ord_map_idx, new_gene_chr))
			for(i in 1:4){
				neigh_idx <- which(rownames(dat_mat) == ord_map[idx[i],2])[1]
			
				if(ord_map[idx[i], "chr"] == new_gene_chr){
#				print(x)
#				print(neigh_idx)	
					cor_val <- cor(dat_mat[dat_mat_idx, col_idx], dat_mat[neigh_idx, col_idx] , method=cor_method)
				}
				else cor_val <- NA
			
				res_mat[i, 3] <- rownames(dat_mat)[neigh_idx]
				res_mat[i, 4] <- ord_map[idx[i], 1]
				res_mat[i, 5] <- c(-2, -1, 1, 2)[i]
				res_mat[i, 6] <- cor_val
			}
		}
		
		return(res_mat)
	})

	result_matrix <- results[[1]]
	
	for(i in 2:length(results)){
		result_matrix <- rbind(result_matrix, results[[i]])
	}
	
	result_matrix <- result_matrix[!is.na(result_matrix[,1]),]
	result_matrix <- result_matrix[!is.na(result_matrix[,6]),]

	result_matrix[,5] <- as.numeric(result_matrix[,5])
	result_matrix[,6] <- as.numeric(result_matrix[,6])

	return(result_matrix)

}

#neighbor_dist_spearman_fix <- neighbor_cor(rpkm_mat_tiss, modENCODE_mRNA_Seq_tissues_list, "spearman", new_gene_list[,"new.fbgn"], gene_map_ordered_tiss_fix)

#neighbor_dist_pc_fix <- neighbor_cor(rpkm_mat_tiss_qnorm, modENCODE_mRNA_Seq_tissues_list_pc, "spearman", new_gene_list[,"new.fbgn"], gene_map_ordered_tiss_fix)

#neighbor_dist_avg_fix <- neighbor_cor(rpkm_mat_tiss_avg, modENCODE_mRNA_Seq_tissues_list_avg, "spearman", new_gene_list[,"new.fbgn"], gene_map_ordered_tiss_fix)

#par_neighbor_dist_avg_fix <- neighbor_cor(rpkm_mat_tiss_avg, modENCODE_mRNA_Seq_tissues_list_avg, "spearman", new_gene_list[,"par.fbgn"], gene_map_ordered_tiss_fix)

#top_coexprs_spearman_fix <- neighbor_dist_spearman_fix[which(as.numeric(neighbor_dist_spearman_fix[,6]) > 0.9),]

#top_coexprs_pc_fix <- neighbor_dist_pc_fix[which(as.numeric(neighbor_dist_pc_fix[,6]) > 0.9),]

#top_coexprs_avg_fix <- neighbor_dist_pc_fix[which(as.numeric(neighbor_dist_avg_fix[,6]) > 0.9),]

##CG34041 and PH4alphaNE2 are *overlapping*!!! good target for CRISPR mediated deletion!

#save(list=c("neighbor_dist_spearman_fix", "neighbor_dist_pc_fix", "top_coexprs_spearman_fix", "top_coexprs_pc_fix", "neighbor_dist_avg_fix", "top_coexprs_avg_fix"), file="neighbor_fix.rdata")

#png("cor_constituents_head.png", width=960, height=480)
#barplot(rbind(cor_pc_conds_head, cor_avg_conds_head),names.arg=c("MatF_4d", "VirF_4d", "MatM_4d", "MatF_1d", "VirF_1d", "MatM_1d", "MatF_20d", "VirF_20d", "MatM_20d"), main="Correlation to Constituents (Head, Spearman)", col=c("red", "blue"), legend=c("PC/Quantile", "Average"), beside=T, ylim=c(0, 1.2), xlab="Tissues", ylab="Corr.")
#dev.off()

#svg("cor_constituents_head.svg", width=14, height=7)
#barplot(rbind(cor_pc_conds_head, cor_avg_conds_head),names.arg=c("MatF_4d", "VirF_4d", "MatM_4d", "MatF_1d", "VirF_1d", "MatM_1d", "MatF_20d", "VirF_20d", "MatM_20d"), main="Correlation to Constituents (Head, Spearman)", col=c("red", "blue"), legend=c("PC/Quantile", "Average"), beside=T, ylim=c(0, 1.2), xlab="Tissues", ylab="Corr.")
#dev.off()

#png("cor_constituents_dige.png")
#barplot(rbind(cor_pc_conds_dige, cor_avg_conds_dige),names.arg=c("20d", "4d"), main="Correlation to Constituents (Dige, Spearman)", col=c("red", "blue"), legend=c("PC/Quantile", "Average"), beside=T, ylim=c(-0.1, 1.2), xlab="Tissues", ylab="Corr.")
#dev.off()

#svg("cor_constituents_dige.svg")
#barplot(rbind(cor_pc_conds_dige, cor_avg_conds_dige),names.arg=c("20d", "4d"), main="Correlation to Constituents (Gut, Spearman)", col=c("red", "blue"), legend=c("PC/Quantile", "Average"), beside=T, ylim=c(-0.1, 1.2), xlab="Tissues", ylab="Corr.")
#dev.off()

#png("cor_constituents_carc.png")
#barplot(rbind(cor_pc_conds_carc, cor_avg_conds_carc),names.arg=c("4d", "1d"), main="Correlation to Constituents (Carcass, Spearman)", col=c("red", "blue"), legend=c("PC/Quantile", "Average"), beside=T, ylim=c(-0.4, 1.4), xlab="Tissues", ylab="Corr.")
#dev.off()

#svg("cor_constituents_carc.svg")
#barplot(rbind(cor_pc_conds_carc, cor_avg_conds_carc),names.arg=c("4d", "1d"), main="Correlation to Constituents (Carcass, Spearman)", col=c("red", "blue"), legend=c("PC/Quantile", "Average"), beside=T, ylim=c(-0.4, 1.4), xlab="Tissues", ylab="Corr.")
#dev.off()

#png("cor_constituents_ovar.png")
#barplot(rbind(cor_pc_conds_ovar, cor_avg_conds_ovar),names.arg=c("VirF_4d", "MatF_1d"), main="Correlation to Constituents (Ovary, Spearman)", col=c("red", "blue"), legend=c("PC/Quantile", "Average"), beside=T, ylim=c(-0.6, 1.4), xlab="Tissues", ylab="Corr.")
#dev.off()

#svg("cor_constituents_ovar.svg")
#barplot(rbind(cor_pc_conds_ovar, cor_avg_conds_ovar),names.arg=c("VirF_4d", "MatF_1d"), main="Correlation to Constituents (Ovary, Spearman)", col=c("red", "blue"), legend=c("PC/Quantile", "Average"), beside=T, ylim=c(-0.6, 1.4), xlab="Tissues", ylab="Corr.")
#dev.off()


##parent analysis 
#parent_cor <- function(dat_mat, col_idx, cor_method, gene_list, par_list, ordered_map){
#	results <- matrix(nrow=length(gene_list), ncol=5)
#	colnames(results) <- c("new_gene_fbgn", "par_gene_fbgn", "new_gene_sym", "par_gene_sym", "cor")
#	
#	for(i in 1:length(gene_list)){
#		results[i,1] <- as.character(gene_list[i])
#		results[i,2] <- as.character(par_list[i])

#		valid_bool <- !is.na(gene_list[i]) && !is.na(par_list[i])
#		valid_bool <- valid_bool && gene_list[i] %in% ordered_map[,1]
#		valid_bool <- valid_bool && par_list[i] %in% ordered_map[,1]

#		if(valid_bool){
#			new_gene_idx <- which(ordered_map[,1] == gene_list[i])[1]
#			new_gene_sym <- ordered_map[new_gene_idx, 2]
#			par_gene_idx <- which(ordered_map[,1] == par_list[i])[1]
#			par_gene_sym <- ordered_map[par_gene_idx, 2]
##			print(c(i, new_gene_sym, par_gene_sym))
#			results[i,3] <- new_gene_sym
#			results[i,4] <- par_gene_sym
#			results[i,5] <- cor(dat_mat[new_gene_sym, col_idx], dat_mat[par_gene_sym, col_idx], method=cor_method)
##			print(i)
#		}
#		else results[i,3] <- NA
#	}

#	return(results)
#}

#par_dist_spearman_avg <- parent_cor(rpkm_mat_tiss_avg, modENCODE_mRNA_Seq_tissues_list_avg, "spearman", new_gene_list[,"new.fbgn"], new_gene_list[,"par.fbgn"], gene_map_ordered_tiss_fix)

#orig_mech_list_denovo <- as.character(new_gene_list[which(new_gene_list[,22] == "A"),2])
#orig_mech_list_dup <- as.character(new_gene_list[which(new_gene_list[,22] == "D"),2])
#orig_mech_list_ret <- as.character(new_gene_list[which(new_gene_list[,22] == "R"),2])

#save(list=c("orig_mech_list_denovo", "orig_mech_list_dup", "orig_mech_list_ret"), file="orig_mech_lists.rdata")

#orig_mech_list_fbgn_denovo <- as.character(new_gene_list[which(new_gene_list[,22] == "A"),4])
#orig_mech_list_fbgn_dup <- as.character(new_gene_list[which(new_gene_list[,22] == "D"),4])
#orig_mech_list_fbgn_ret <- as.character(new_gene_list[which(new_gene_list[,22] == "R"),4])

##save(list=c("orig_mech_list_fbgn_denovo", "orig_mech_list_fbgn_dup", "orig_mech_list_fbgn_ret"), file="orig_mech_lists_fbgn.rdata")

#coexprs_dist_list_fix <- list(bootstrap_dist_spear_tiss_avg, baseline_dist_spear_tiss_avg[,3], as.numeric(par_dist_spearman_avg[,5]), as.numeric(neighbor_dist_avg_fix[,6]))

#names(coexprs_dist_list_fix) <- c("boot", "baseline", "parent_to_new", "new_to_neigh")

#neighbor_dist_avg_imm_fix <- neighbor_dist_avg_fix[which(neighbor_dist_avg_fix[,5] %in% c("-1", "1")),]

#paired_cor_vals_fix <- matrix(nrow=nrow(new_gene_list), ncol=4)
#for(i in 1:length(new_gene_list[,4])){
#	this_gene <- new_gene_list[i,4]
#		
#	neighbor_cor <- mean(as.numeric(neighbor_dist_avg_fix[which(neighbor_dist_avg_fix[,2] == this_gene), 6]), na.rm=T)
#	neighbor_cor_imm <- mean(as.numeric(neighbor_dist_avg_imm_fix[which(neighbor_dist_avg_imm_fix[,2] == this_gene), 6]), na.rm=T)
#	neighbor_cor_max <- max(as.numeric(neighbor_dist_avg_fix[which(neighbor_dist_avg_fix[,2] == this_gene), 6]), na.rm=T)
#	if(is.infinite(neighbor_cor_max)) neighbor_cor_max <- NaN
#	par_cor <- as.numeric(par_dist_spearman_avg[which(par_dist_spearman_avg[,1] == this_gene), 5])
#print(i)
#	paired_cor_vals_fix[i,] <- c(neighbor_cor, neighbor_cor_imm, neighbor_cor_max, par_cor)
#}

#coexprs_dist_list_fix <- list(bootstrap_dist_spear_tiss_avg, baseline_dist_spear_tiss_avg[,3], as.numeric(par_dist_spearman_avg[,5]), as.numeric(neighbor_dist_avg_fix[,6]), paired_cor_vals_fix[,3])

#names(coexprs_dist_list_fix) <- c("boot", "baseline", "parent_to_new", "new_to_neigh", "new_to_neigh_max")

#dist_to_genes <- function(dat_mat, new_fbgn, par_fbgn){
##        print(new_fbgn)
#        
#        if(! (is.na(new_fbgn) || is.na(par_fbgn))){
#                new_idx <- which(dat_mat[,1] == new_fbgn)
#                par_idx <- which(dat_mat[,1] == par_fbgn)
#                
#                if(length(new_idx) + length(par_idx) < 2) return(NA)

#                if(dat_mat[new_idx,3] == dat_mat[par_idx,3]){
#                        dist <- as.numeric(dat_mat[new_idx,4]) - as.numeric(dat_mat[par_idx,4])
#                        return(abs(dist))
#                }
#                else return(Inf)
#        }

#        return(NA)
#}

#gene_pair_linkage <- sapply(1:42, function(x) dist_to_genes(gene_map_ordered_tiss_fix, new_gene_list[x,4], new_gene_list[x,5]))

##noticed doing dist_to_genes that tss values for gene_map_ordered_trim and gene_map_ordered_tiss are INCORRECT (2016/05/12). running fix to correct. pre 2016/05/12 tss values are inaccurate! analyses are redone from here using _fix append - UL
#shared_gene_trim <- which(gene_map_ordered[,2] %in% rownames(rpkm_mat_tiss))

#gene_map_ordered_trim_fix <- gene_map_ordered[shared_gene_trim,]

#shared_gene_tiss_sums <- rowSums(rpkm_mat_tiss_avg[gene_map_ordered_trim_fix[,2], modENCODE_mRNA_Seq_tissues_list_avg])

#shared_gene_tiss_exprd <- which(!shared_gene_tiss_sums == 0)

#gene_map_ordered_tiss_fix <- gene_map_ordered_trim_fix[shared_gene_tiss_exprd,]
#save(list="gene_map_ordered_tiss_fix", file="ordered_map_table_tiss_fix.rdata")
#save(list="gene_map_ordered_trim_fix", file="ordered_map_table_trim_fix.rdata")



#png("coexprs_classes_fix.png", width=480*2, height=480*2)
#par(mar=c(5,6,4,2) + 0.1) 
#plot(paired_cor_vals_fix[,3], paired_cor_vals_fix[,4], xlab="Max Corr. to Neighbor", ylab="Corr. to Parent", main="Classification of Regulatory Interactions (Spearman)", col="white", cex.lab=3, cex.main=3, cex.axis=2)
#for(i in 1:42){
#	if(!new_gene_list[i,22] == "A"){
#		if(new_gene_list[i,22] == "D") point_type <- 23
#		else point_type <- 22

#		if(is.na(gene_pair_linkage[i])) point_col <- "white"
#		else{
#			if(gene_pair_linkage[i] < 100000) point_col <- "orange"
#			else point_col <- "blue"
#		}
#		if(new_gene_list[i, 1] == "0~3") point_cex <- 1
#		if(new_gene_list[i, 1] == "3~6") point_cex <- 2
#		if(new_gene_list[i, 1] == "6~11") point_cex <- 3
#		if(new_gene_list[i, 1] == "11~25") point_cex <- 4
#		if(new_gene_list[i, 1] == "25~35") point_cex <- 5
#		point_lwd <- sum(new_gene_list[i, 23] <= quantile(new_gene_list[,23], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=T))
#		
#		points(paired_cor_vals_fix[i, 3], paired_cor_vals_fix[i, 4], pch=point_type, col="black", bg=point_col, cex=3*point_cex-2, lwd=3*point_lwd)
#	}
#}
#points(paired_cor_vals_fix[40, 3], paired_cor_vals_fix[40, 4], pch=23, col="black", cex=3*4-2, bg="orange", lwd=3*point_lwd)

#abline(h=0.5)
#abline(v=0.5)
#dev.off()

#svg("coexprs_classes_fix.svg", width=14, height=14)
#par(mar=c(5,6,4,2) + 0.1) 
#plot(paired_cor_vals_fix[,3], paired_cor_vals_fix[,4], xlab="Max Corr. to Neighbor", ylab="Corr. to Parent", main="Classification of Regulatory Interactions (Spearman)", col="white", cex.lab=3, cex.main=3, cex.axis=2)
#for(i in 1:42){
#	if(!new_gene_list[i,22] == "A"){
#		if(new_gene_list[i,22] == "D") point_type <- 23
#		else point_type <- 22

#		if(is.na(gene_pair_linkage[i])) point_col <- "white"
#		else{
#			if(gene_pair_linkage[i] < 100000) point_col <- "orange"
#			else point_col <- "blue"
#		}
#		if(new_gene_list[i, 1] == "0~3") point_cex <- 1
#		if(new_gene_list[i, 1] == "3~6") point_cex <- 2
#		if(new_gene_list[i, 1] == "6~11") point_cex <- 3
#		if(new_gene_list[i, 1] == "11~25") point_cex <- 4
#		if(new_gene_list[i, 1] == "25~35") point_cex <- 5
#		point_lwd <- sum(new_gene_list[i, 23] <= quantile(new_gene_list[,23], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=T))
#		
#		points(paired_cor_vals_fix[i, 3], paired_cor_vals_fix[i, 4], pch=point_type, col="black", bg=point_col, cex=3*point_cex-2, lwd=3*point_lwd)
#	}
#}
#points(paired_cor_vals_fix[40, 3], paired_cor_vals_fix[40, 4], pch=23, col="black", cex=3*4-2, bg="orange", lwd=3*point_lwd)

#abline(h=0.5)
#abline(v=0.5)
#dev.off()

#png("coexprs_classes_fix_mean.png", width=480*2, height=480*2)
#par(mar=c(5,6,4,2) + 0.1) 
#plot(paired_cor_vals_fix[,3], paired_cor_vals_fix[,4], xlab="Max Corr. to Neighbor", ylab="Corr. to Parent", main="Classification of Regulatory Interactions (Spearman)", col="white", cex.lab=3, cex.main=3, cex.axis=2)
#for(i in 1:42){
#	if(!new_gene_list[i,22] == "A"){
#		if(new_gene_list[i,22] == "D") point_type <- 23
#		else point_type <- 22

#		if(is.na(gene_pair_linkage[i])) point_col <- "white"
#		else{
#			if(gene_pair_linkage[i] < 100000) point_col <- "orange"
#			else point_col <- "blue"
#		}
#		if(new_gene_list[i, 1] == "0~3") point_cex <- 1
#		if(new_gene_list[i, 1] == "3~6") point_cex <- 2
#		if(new_gene_list[i, 1] == "6~11") point_cex <- 3
#		if(new_gene_list[i, 1] == "11~25") point_cex <- 4
#		if(new_gene_list[i, 1] == "25~35") point_cex <- 5
#		point_lwd <- sum(new_gene_list[i, 23] <= quantile(new_gene_list[,23], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=T))
#		
#		points(paired_cor_vals_fix[i, 2], paired_cor_vals_fix[i, 4], pch=point_type, col="black", bg=point_col, cex=3*point_cex-2, lwd=3*point_lwd)
#	}
#}
#points(paired_cor_vals_fix[40, 2], paired_cor_vals_fix[40, 4], pch=23, col="black", cex=3*4-2, bg="orange", lwd=3*point_lwd)

#abline(h=0.5)
#abline(v=0.5)
#dev.off()

#svg("coexprs_classes_fix_mean.svg", width=14, height=14)
#par(mar=c(5,6,4,2) + 0.1) 
#plot(paired_cor_vals_fix[,3], paired_cor_vals_fix[,4], xlab="Max Corr. to Neighbor", ylab="Corr. to Parent", main="Classification of Regulatory Interactions (Spearman)", col="white", cex.lab=3, cex.main=3, cex.axis=2)
#for(i in 1:42){
#	if(!new_gene_list[i,22] == "A"){
#		if(new_gene_list[i,22] == "D") point_type <- 23
#		else point_type <- 22

#		if(is.na(gene_pair_linkage[i])) point_col <- "white"
#		else{
#			if(gene_pair_linkage[i] < 100000) point_col <- "orange"
#			else point_col <- "blue"
#		}
#		if(new_gene_list[i, 1] == "0~3") point_cex <- 1
#		if(new_gene_list[i, 1] == "3~6") point_cex <- 2
#		if(new_gene_list[i, 1] == "6~11") point_cex <- 3
#		if(new_gene_list[i, 1] == "11~25") point_cex <- 4
#		if(new_gene_list[i, 1] == "25~35") point_cex <- 5
#		point_lwd <- sum(new_gene_list[i, 23] <= quantile(new_gene_list[,23], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=T))
#		
#		points(paired_cor_vals_fix[i, 2], paired_cor_vals_fix[i, 4], pch=point_type, col="black", bg=point_col, cex=3*point_cex-2, lwd=3*point_lwd)
#	}
#}
#points(paired_cor_vals_fix[40, 2], paired_cor_vals_fix[40, 4], pch=23, col="black", cex=3*4-2, bg="orange", lwd=3*point_lwd)

#abline(h=0.5)
#abline(v=0.5)
#dev.off()

#png("coexprs_classes_key_fix.png", width=480*2, height=480*2)
#par(mar=c(5,6,4,2) + 0.1) 
#plot(paired_cor_vals_fix[,3], paired_cor_vals_fix[,4], xlab="Max Corr. to Neighbor", ylab="Corr. to Parent", main="Classification of Regulatory Interactions (Spearman)", col="white", cex.lab=3, cex.main=3, cex.axis=2)
#for(i in 1:42){
#	if(!new_gene_list[i,22] == "A"){
#		if(new_gene_list[i,22] == "D") point_type <- 23
#		else point_type <- 22

#		if(is.na(gene_pair_linkage[i])) point_col <- "white"
#		else{
#			if(gene_pair_linkage[i] < 100000) point_col <- "orange"
#			else point_col <- "blue"
#		}
#		if(new_gene_list[i, 1] == "0~3") point_cex <- 1
#		if(new_gene_list[i, 1] == "3~6") point_cex <- 2
#		if(new_gene_list[i, 1] == "6~11") point_cex <- 3
#		if(new_gene_list[i, 1] == "11~25") point_cex <- 4
#		if(new_gene_list[i, 1] == "25~35") point_cex <- 5
#		point_lwd <- sum(new_gene_list[i, 23] <= quantile(new_gene_list[,23], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=T))
#		
#		points(paired_cor_vals_fix[i, 3], paired_cor_vals_fix[i, 4], pch=point_type, col="black", bg=point_col, cex=3*point_cex-2, lwd=3*point_lwd)
#		text(paired_cor_vals_fix[i,3], paired_cor_vals_fix[i,4], new_gene_list[i,2])
#	}
#}
#points(paired_cor_vals_fix[40, 3], paired_cor_vals_fix[40, 4], pch=23, col="black", cex=3*4-2, bg="orange", lwd=3*point_lwd)
#text(paired_cor_vals_fix[40,3], paired_cor_vals_fix[40,4], new_gene_list[40,2])

#abline(h=0.5)
#abline(v=0.5)
#dev.off()

#svg("coexprs_classes_key_fix.svg", width=480*2, height=480*2)
#par(mar=c(5,6,4,2) + 0.1) 
#plot(paired_cor_vals_fix[,3], paired_cor_vals_fix[,4], xlab="Max Corr. to Neighbor", ylab="Corr. to Parent", main="Classification of Regulatory Interactions (Spearman)", col="white", cex.lab=3, cex.main=3, cex.axis=2)
#for(i in 1:42){
#	if(!new_gene_list[i,22] == "A"){
#		if(new_gene_list[i,22] == "D") point_type <- 23
#		else point_type <- 22

#		if(is.na(gene_pair_linkage[i])) point_col <- "white"
#		else{
#			if(gene_pair_linkage[i] < 100000) point_col <- "orange"
#			else point_col <- "blue"
#		}
#		if(new_gene_list[i, 1] == "0~3") point_cex <- 1
#		if(new_gene_list[i, 1] == "3~6") point_cex <- 2
#		if(new_gene_list[i, 1] == "6~11") point_cex <- 3
#		if(new_gene_list[i, 1] == "11~25") point_cex <- 4
#		if(new_gene_list[i, 1] == "25~35") point_cex <- 5
#		point_lwd <- sum(new_gene_list[i, 23] <= quantile(new_gene_list[,23], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=T))
#		
#		points(paired_cor_vals_fix[i, 3], paired_cor_vals_fix[i, 4], pch=point_type, col="black", bg=point_col, cex=3*point_cex-2, lwd=3*point_lwd)
#		text(paired_cor_vals_fix[i,3], paired_cor_vals_fix[i,4], new_gene_list[i,2])
#	}
#}
#points(paired_cor_vals_fix[40, 3], paired_cor_vals_fix[40, 4], pch=23, col="black", cex=3*4-2, bg="orange", lwd=3*point_lwd)
#text(paired_cor_vals_fix[40,3], paired_cor_vals_fix[40,4], new_gene_list[40,2])

#abline(h=0.5)
#abline(v=0.5)
#dev.off()



#png("new_par_prot_ident.png")
#par(mar=c(5,5,4,2) + 0.1) 
#hist(new_gene_list[,23], xlab="% Protein Identity w/ Par. Gene", ylab="Counts", main="Protein Sequence Conservation in New Essential Genes", cex.lab=2, cex.main=1.25, cex.axis=2)
#dev.off()

#svg("new_par_prot_ident.svg")
#par(mar=c(5,5,4,2) + 0.1) 
#hist(new_gene_list[,23], xlab="% Protein Identity w/ Par. Gene", ylab="Counts", main="Protein Sequence Conservation in New Essential Genes", cex.lab=2, cex.main=1.25, cex.axis=2)
#dev.off()

#par_neighbor_dist_avg_fix_max <- matrix(0, nrow=length(levels(as.factor(par_neighbor_dist_avg_fix[,1]))), ncol=6)

#for(i in 1:length(levels(as.factor(par_neighbor_dist_avg_fix[,1])))){
#	goi <- levels(as.factor(par_neighbor_dist_avg_fix[,1]))[i]
#	
#	neighs <- par_neighbor_dist_avg_fix[which(par_neighbor_dist_avg_fix[,1] == goi),]
#	max_neigh_idx <- which(as.numeric(neighs[,6]) == max(as.numeric(neighs[,6])))[1]
#	
#	par_neighbor_dist_avg_fix_max[i,] <- neighs[max_neigh_idx,]
#}

#neighbor_dist_avg_fix_max <- matrix(0, nrow=length(levels(as.factor(neighbor_dist_avg_fix[,1]))), ncol=6)

#for(i in 1:length(levels(as.factor(neighbor_dist_avg_fix[,1])))){
#	goi <- levels(as.factor(neighbor_dist_avg_fix[,1]))[i]
#	
#	neighs <- neighbor_dist_avg_fix[which(neighbor_dist_avg_fix[,1] == goi),]
#	max_neigh_idx <- which(as.numeric(neighs[,6]) == max(as.numeric(neighs[,6])))[1]
#	
#	neighbor_dist_avg_fix_max[i,] <- neighs[max_neigh_idx,]
#}

#new_par_neighbors <- matrix(0, nrow=1, ncol=6)
#for(i in 1:nrow(new_gene_list)){
#	if(new_gene_list[i,5] %in% par_neighbor_dist_avg_fix_max[,2] && new_gene_list[i,4] %in% neighbor_dist_avg_fix_max[,2]){
#		neigh_par_coexp <- par_neighbor_dist_avg_fix_max[which(par_neighbor_dist_avg_fix_max[,2] == new_gene_list[i,5]), 6]
#		neigh_new_coexp <- neighbor_dist_avg_fix_max[which(neighbor_dist_avg_fix_max[,2] == new_gene_list[i,4]),6]

#		new_par_neighbors <- rbind(new_par_neighbors, c(as.character(new_gene_list[i,2]), as.character(new_gene_list[i,3]), as.character(new_gene_list[i,4]), as.character(new_gene_list[i,5]), neigh_par_coexp, neigh_new_coexp)) 
#	}
#}

#png("max_coxprs_new_par.png", width=480, height=480)
#plot(as.numeric(new_par_neighbors[,5]), as.numeric(new_par_neighbors[,6]), xlab="Max Co-expression (New)", ylab="Max Co-expression (Parent)")
#abline(0,1)
#dev.off()

#svg("max_coxprs_new_par.svg", width=7, height=7)
#plot(as.numeric(new_par_neighbors[,5]), as.numeric(new_par_neighbors[,6]), xlab="Max Co-expression (New)", ylab="Max Co-expression (Parent)")
#abline(0,1)
#dev.off()

#non_lethal <- read.csv("non_lethal.csv", sep=",")

#par_dist_non_lethal_fix <- parent_cor(rpkm_mat_tiss_avg, modENCODE_mRNA_Seq_tissues_list_avg, "spearman", non_lethal[,"new.fbgn"], non_lethal[,"par.fbgn"], gene_map_ordered_tiss_fix)

#neighbor_dist_non_lethal_fix <- neighbor_cor(rpkm_mat_tiss, modENCODE_mRNA_Seq_tissues_list, "spearman", non_lethal[,"new.fbgn"], gene_map_ordered_tiss_fix)

#coexprs_dist_non_lethal <- rep(0, 9)
#for(i in 1:nrow(non_lethal)){
#	goi <- non_lethal[i,3]

#	par_row <- par_dist_non_lethal_fix[which(par_dist_non_lethal_fix[,1] == goi),]
#	neigh_mtx <- neighbor_dist_non_lethal_fix[which(neighbor_dist_non_lethal_fix[,2] == goi),]

#	if(!is.na(par_row[5])){
#		neigh_mtx_max <- max(as.numeric(neigh_mtx[,6]))
#		neigh_mtx_max_idx <- which(as.numeric(neigh_mtx[,6]) == neigh_mtx_max)[1]

#		res <- c(unlist(lapply(non_lethal[i,1:5], as.character)), neigh_mtx[neigh_mtx_max_idx,3:4], par_row[5], neigh_mtx_max)

#		coexprs_dist_non_lethal <- rbind(coexprs_dist_non_lethal, res)
#	}
#}

#coexprs_dist_non_lethal <- coexprs_dist_non_lethal[-1,]

png("nonlethal_coexprs.png", width=480, height=480)
plot(as.numeric(coexprs_dist_non_lethal[,9]), as.numeric(coexprs_dist_non_lethal[,8]), ylim=c(-0.6, 1), xlab="Max. Co-expression with Neighbor", ylab="Co-expression with Parent", main="Non-Lethal, Drosophila")
abline(v=0.5)
abline(h=0.5)
dev.off()

svg("nonlethal_coexprs.svg", width=14, height=14)
plot(as.numeric(coexprs_dist_non_lethal[,9]), as.numeric(coexprs_dist_non_lethal[,8]), ylim=c(-0.6, 1), xlab="Max. Co-expression with Neighbor", ylab="Co-expression with Parent", main="Non-Lethal, Drosophila")
abline(v=0.5)
abline(h=0.5)
dev.off()

png("density_par.png", width=480, height=480)
d1 <- density(paired_cor_vals_fix[,4], na.rm=T)
d2 <- density(as.numeric(coexprs_dist_non_lethal[,8]), na.rm=T)
plot(range(d1$x, d2$y), range(d1$y, d2$y), type="n", xlab="Parental Co-expression", ylab="Density", main="Lethal vs Non-lethal, p < 0.0007")
lines(d1, col="red", lwd=3)
lines(d2, col="black", lwd=3)
dev.off()

svg("density_par.svg", width=14, height=14)
d1 <- density(paired_cor_vals_fix[,4], na.rm=T)
d2 <- density(as.numeric(coexprs_dist_non_lethal[,8]), na.rm=T)
plot(range(d1$x, d2$y), range(d1$y, d2$y), type="n", xlab="Parental Co-expression", ylab="Density", main="Lethal vs Non-lethal, p < 0.0007")
lines(d1, col="red", lwd=3)
lines(d2, col="black", lwd=3)
dev.off()

png("density_neigh.png", width=480, height=480)
d1 <- density(paired_cor_vals_fix[,3], na.rm=T)
d2 <- density(as.numeric(coexprs_dist_non_lethal[,9]), na.rm=T)
plot(range(d1$x, d2$y), range(d1$y, d2$y), type="n", xlab="Max Neighbor Co-expression", ylab="Density", main="Lethal vs Non-lethal, p = 0.88")
lines(d1, col="red", lwd=3)
lines(d2, col="black", lwd=3)
dev.off()

svg("density_neigh.svg", width=14, height=14)
d1 <- density(paired_cor_vals_fix[,3], na.rm=T)
d2 <- density(as.numeric(coexprs_dist_non_lethal[,9]), na.rm=T)
plot(range(d1$x, d2$y), range(d1$y, d2$y), type="n", xlab="Max Neighbor Co-expression", ylab="Density", main="Lethal vs Non-lethal, p = 0.88")
lines(d1, col="red", lwd=3)
lines(d2, col="black", lwd=3)
dev.off()

png("coexprs_classes_fix_plain.png", width=480, height=480)
par(mar=c(5,6,4,2) + 0.1) 
plot(paired_cor_vals_fix[,3], paired_cor_vals_fix[,4], xlab="Max Co-expression with Neighbor", ylab="Co-expression with Parent", main="Lethal, Drosophila")

abline(h=0.5)
abline(v=0.5)
dev.off()

##"network" analysis

##generate correlation matrix
#gen_corr_mtx <- function(dat_mat, col_idents, margin=1, cor_meth="spearman"){
#	cor_mtx <- matrix(0, nrow=length(col_idents), ncol=length(col_idents))

#	#calculate lower left values of mtx (note inner loop begins @ index i)
#	for(i in 1:length(col_idents)){
#		for(j in i:length(col_idents)){
#			if(margin==1) cor_mtx[i,j] <- cor(dat_mat[col_idents[i],], dat_mat[col_idents[j],], method=cor_meth)
#			if(margin==2) cor_mtx[i,j] <- cor(dat_mat[,col_idents[i]], dat_mat[,col_idents[j]], method=cor_meth)
#		}
#	}

#	rownames(cor_mtx) <- col_idents
#	colnames(cor_mtx) <- col_idents
#	 
#	 #make mtx symmetric and adjust for double-counting on diagonal
#	cor_mtx <- cor_mtx + t(cor_mtx)
#	for(i in 1:length(col_idents)){
#	 	cor_mtx[i,i] <- 1
#	}

#	return(cor_mtx)
#}

##generate bootstrap (or permutation) matrix, margin=1 for bootstrapping by row, permutation=T for perm test
#gen_boot_mtx <- function(dat_mat, margin=2, permutation=F){
#	result <- apply(dat_mat, margin, function(x){
#		return(sample(x, length(x), replace=!permutation))
#	})

#	rownames(result) <- rownames(dat_mat)
#	colnames(result) <- colnames(dat_mat)
#	return(result)
#}

#corr_mtx_tiss <- gen_corr_mtx(rpkm_mat_tiss_avg, rownames(rpkm_mat_tiss_avg))
#corr_boot <- gen_corr_mtx(gen_boot_mtx(rpkm_mat_tiss), rownames(rpkm_mat_tiss_avg))

##convert to per-gene ranking, margin=1 for ranking by row (shouldn't matter if using MRM
#gen_rank_mtx <- function(cor_mat, margin=1){
#	result <- apply(cor_mat, margin, rank)
#	
#	rownames(result) <- rownames(cor_mat)
#	rownames(result) <- rownames(cor_mat)
#	
#	return(result)
#}

#rank_mtx_tiss <- gen_rank_mtx(corr_mtx_tiss)
#rank_mtx_boot <- gen_rank_mtx(corr_boot)

##calculate mutual rank mean of {R_ab, R_ba}, input should be square mtx
#gen_MRM <- function(rank_mat, square=F){
#	result <- matrix(0, nrow=nrow(rank_mat), ncol=ncol(rank_mat))
#	
#	for(i in 1:nrow(rank_mat)){
#		for(j in i:nrow(rank_mat)){
#			result[i,j] <- sqrt(rank_mat[i,j] * rank_mat[j,i])
#		}
#	}

#	if(square){
#		result <- result + t(result)
#		for(i in 1:nrow(result)) result[i,i] <- 1
#	}
#	
#	rownames(result) <- rownames(rank_mat)
#	colnames(result) <- colnames(rank_mat)
#	
#	return(result)
#}

#mrm_mtx_tiss <- gen_MRM(rank_mtx_tiss)
#mrm_mtx_boot <- gen_MRM(corr_boot)

##histogram
#mrm_dist_tiss <- as.numeric(mrm_mtx_tiss)
#mrm_dist_tiss <- mrm_dist_tiss[-which(mrm_dist_tiss == 0)]

#mrm_dist_boot <- as.numeric(mrm_mtx_boot)
#mrm_dist_boot <- mrm_dist_boot[-which(mrm_dist_boot == 0)]

#png("mrm_tiss.png")
#hist(mrm_dist_tiss, breaks=25, main="Mutual Rank Mean (Tissue)")
#dev.off()

#svg("mrm_tiss.svg")
#hist(mrm_dist_tiss, breaks=25, main="Mutual Rank Mean (Tissue)")
#dev.off()

#png("mrm_boot.png")
#hist(mrm_dist_boot, breaks=25, main="Mutual Rank Mean (Bootstrap)")
#dev.off()

#svg("mrm_boot.svg")
#hist(mrm_dist_boot, breaks=25, main="Mutual Rank Mean (Bootstrap)")
#dev.off()

##take top co-correlators, mrm_mat is half zeros
#top_mrm <- function(mrm_mat, alpha){
#	alpha_prime <- alpha*alpha/2
#	flat_mrm <- as.numeric(mrm_mtx_tiss)
#	flat_mrm <- flat_mrm[-which(flat_mrm == 0)]
#	thresh <- quantile(flat_mrm, 1-alpha_prime)
#	
#	res_row <- c()
#	res_col <- c()	

#	for(i in 1:nrow(mrm_mat)){
#		for(j in i:nrow(mrm_mat)){
#			if(mrm_mat[i,j] >= thresh){
#				res_row <- c(res_row, i)
#				res_col <- c(res_col, j)
#			}
#		}
#	}
#	
#	result <- rbind(res_row, res_col)
#	rownames(result) <- c("row_idx", "col_idx")
#	
#	return(result)
#}

#mrm_network_tiss <- top_mrm(mrm_mtx_tiss, 0.05)
#mrm_network_boot <- top_mrm(mrm_mtx_boot, 0.05)


##construct network



