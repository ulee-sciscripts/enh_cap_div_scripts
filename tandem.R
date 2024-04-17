gene_exprs <- read.csv("gene_exprs_30dev.tsv", sep="\t")
orig_mech <- read.csv("orig_mech.tsv", sep="\t")

idx_match <- which(orig_mech[,"ParentalGeneID"] %in% gene_exprs[,2])

orig_mech <- orig_mech[idx_match,]


#some 'RE' labels are mislabeled as 'RT'
orig_mech[which(orig_mech[,11] == "RT"),11] <- "RE"
tdm_idx <- which(orig_mech[,11] == "TD")
nontdm_idx <- which(orig_mech[,11] %in% c("DD", "RE"))
dd_idx <- which(orig_mech[,11] == "DD")
re_idx <- which(orig_mech[,11] == "RE")

calc_par_cor <- function(in_row, dat_mat){
	dat_cols <- 4:33
	
	par_id <- in_row["ParentalGeneID"]
	par_id_idx <- which(dat_mat[,"gene_id"] == par_id)
	new_id <- in_row["New_Gene"]
	new_id_idx <- which(dat_mat[,"gene_id"] == new_id)
	
	out_cor <- cor(as.numeric(dat_mat[par_id_idx, dat_cols]), as.numeric(dat_mat[new_id_idx, dat_cols]), method="spearman")
	return(out_cor)
}

tdm_cors <- apply(orig_mech[tdm_idx, ], 1, function(x) calc_par_cor(x, gene_exprs))
nontdm_cors <- apply(orig_mech[nontdm_idx, ], 1, function(x) calc_par_cor(x, gene_exprs))
dd_cors <- apply(orig_mech[dd_idx, ], 1, function(x) calc_par_cor(x, gene_exprs))
re_cors <- apply(orig_mech[re_idx, ], 1, function(x) calc_par_cor(x, gene_exprs))

cor_dat <- list(tdm_cors, nontdm_cors, dd_cors, re_cors)
names(cor_dat) <- c("Tandem\nN=85", "Non-Tandem\nN=71", "Distal Dup.\nN=61", "Retro-Element\nN=10")

png("box_cors.png",  width=7, height=7, res=300, pointsize=10, units="in")
par(mgp=c(3, 2, 0), mar=c(5.1, 5.5, 4.1, 2.1))
boxplot(cor_dat, xlab="", main="Developmental Stage Expression Correlation", ylab="")
title(xlab="Duplication Type", line=4, cex.lab=1.5)
title(ylab="Correlation (Spearman)", line=4, cex.lab=1.5)
dev.off()

svg("box_cors.svg",  width=7, height=7, pointsize=10)
par(mgp=c(3, 2, 0), mar=c(5.1, 5.5, 4.1, 2.1))
boxplot(cor_dat, xlab="", main="Developmental Stage Expression Correlation", ylab="")
title(xlab="Duplication Type", line=4, cex.lab=1.5)
title(ylab="Correlation (Spearman)", line=4, cex.lab=1.5)
dev.off()

t.test(tdm_cors, nontdm_cors)
t.test(tdm_cors, dd_cors)
t.test(tdm_cors, re_cors)