#take coordinates, range and filename
#read in header
#figure out correct line numbers
#pull matrix
library(stringr)

pull_sub_mtx <- function(coord_chr, coord_loc, range, resolution, filename, norm=T, idx_pass=0, square=F){
	locs <- read.csv(filename, sep="\t", header=F, nrows=1)
	chr_mask <- unlist(lapply(locs, function(x) str_detect(x, coord_chr)))
	locs_chr <- locs[chr_mask]
	locs_strings <- unlist(lapply(locs_chr, function(x) gsub(".*-", "", x)))
	locs_int <- as.numeric(locs_strings)
	idx_locs_chr <- which(locs_int >= coord_loc)[1]
	
	loc_name <- as.character(locs_chr[idx_locs_chr - 1])
	idx_locs <- which(locs == loc_name)
	
	if(!idx_pass == 0) idx_locs <- idx_pass

	chunks <- floor(range/(2* resolution))
	idx_range <- (idx_locs-chunks):(idx_locs+chunks)

	num_rows <- 1
	idx_in_file <- idx_locs-3
	if(square){
		num_rows <- chunks*2 + 1
		idx_in_file <- idx_locs - chunks - 3
	}

	pull_row <- read.csv(filename, sep="\t", header=T, nrow=num_rows, skip=idx_in_file)
	colnames(pull_row) <- locs

	if(square) {
		print(nrow(pull_row))
		

		rownames(pull_row) <- locs[idx_range]
		return(pull_row[,idx_range])
	}
	
	pull_row_vec <- as.numeric(pull_row[idx_range])
	names(pull_row_vec) <- names(pull_row[idx_range])
	
	if(norm) return(normalize_track(pull_row_vec))
	return(pull_row_vec)
}	

mel_file <- "../mel/5kb_mel.txt"
mel_chr <- "chr2L"
yak_file <- "../yak/5kb_yak.txt"
yak_chr <- "chr2L"
pse_file <- "../round2/pse/7.5kb.txt"
pse_chr <- "Muller_B"
mir_file <- "../round2/mir/7.5kb.txt"
mir_chr <- "chr4"

rnge <- 400000
reso <- 5000

reso_mir <- 7500

mel_chr_par <- "chrX"
yak_chr_par <- "chrX"
pse_chr_par <- "Muller_A-AD"
mir_chr_par <- "chrXL"


#normalize
normalize_track <- function(track){
	return(track/sum(track))
}

#gen list of chrs and sizes
chrs_sizes <- function(filename){
	locs <- read.csv(filename, sep="\t", header=F, nrows=1)
	chrs_vec <- unlist(lapply(locs, function(x) gsub("-.*", "", x)))[-1:-2]
	chrs_fac <- as.factor(chrs_vec)
	chrs_list <- levels(as.factor(chrs_vec))
	
	idx_list <- lapply(chrs_list, function(x){ 
		this_which <- which(as.factor(chrs_vec) == x)
		return(c(this_which[1] + 2, this_which[length(this_which)]+2))}
	)
	
	idx_mtx <- do.call("cbind", idx_list)
	colnames(idx_mtx) <- chrs_list
	
	return(idx_mtx)
}

track_mel_hp6 <- pull_sub_mtx(mel_chr, 4570000, rnge, reso, mel_file, T) 
track_mel_nei <- pull_sub_mtx(mel_chr, 4710000, rnge, reso, mel_file, T)
track_yak_hp6 <- pull_sub_mtx(yak_chr, 4680000, rnge, reso, yak_file, T)
track_yak_nei <- pull_sub_mtx(yak_chr, 4820000, rnge, reso, yak_file, T)
track_pse_hp6 <- pull_sub_mtx(pse_chr, 29165000, rnge, reso_mir, pse_file, T)
track_pse_nei <- pull_sub_mtx(pse_chr, 29070000, rnge, reso_mir, pse_file, T)

track_mir_hp6 <- pull_sub_mtx(mir_chr, 31160000, rnge, reso_mir, mir_file, T)
track_mir_nei <- pull_sub_mtx(mir_chr, 30970000, rnge, reso_mir, mir_file, T)

track_mel_par <- pull_sub_mtx(mel_chr_par, 9080000, rnge, reso, mel_file, T)
track_yak_par <- pull_sub_mtx(yak_chr_par, 12685000, rnge, reso, yak_file, T)
track_pse_par <- pull_sub_mtx(pse_chr_par, 9690000, rnge, reso_mir, pse_file, T)
track_mir_par <- pull_sub_mtx(mir_chr_par, 15345000, rnge, reso_mir, mir_file, T)



gc()

#sample genome
mc_genome_tracks <- function(filename, range, resolution, nsamp, norm=T){
	chunks <- floor(range/(2*resolution))
	
	idx_mtx <- chrs_sizes(filename)
	idx_diff <- abs(idx_mtx[2,] - idx_mtx[1,])
	idx_weights <- idx_diff/sum(idx_diff)

	if(sum(idx_diff < (2*chunks)) > 0){
		idx_problem <- which(idx_diff < (2*chunks)+1)
		print(idx_mtx[1,idx_problem])
		idx_mtx <- idx_mtx[,-idx_problem]
		idx_weights <- idx_weights[-idx_problem]
	}
	
	idx_samp_chrs <- sample(colnames(idx_mtx), nsamp, replace=T, prob=idx_weights)

	idx_samp_center <- unlist(lapply(idx_samp_chrs, function(x){
		bnd_l <- idx_mtx[1,x] + chunks
		bnd_r <- idx_mtx[2,x] - chunks
		
		return(sample(bnd_l:bnd_r, 1))
	}))
	
	idx_samp_mtx <- rbind(idx_samp_chrs, idx_samp_center)


	tracks_mtx <- apply(idx_samp_mtx, 2, function(x){
		this_chr <- x[1]
		this_loc <- as.numeric(x[2])
		
		print(as.character(c(this_chr, this_loc)))
		
		return(pull_sub_mtx(this_chr, this_loc, range, resolution, filename, norm, this_loc))
	})
	
	return(t(tracks_mtx))	
}

baseline_mel1 <- mc_genome_tracks(mel_file, rnge, reso, 333, T)
baseline_mel2 <- mc_genome_tracks(mel_file, rnge, reso, 333, T)
baseline_mel3 <- mc_genome_tracks(mel_file, rnge, reso, 334, T)
save(list=ls(), file="wrkspce2.rdata")
gc()
baseline_yak1 <- mc_genome_tracks(yak_file, rnge, reso, 333, T)
baseline_yak2 <- mc_genome_tracks(yak_file, rnge, reso, 333, T)
baseline_yak3 <- mc_genome_tracks(yak_file, rnge, reso, 334, T)
save(list=ls(), file="wrkspce2.rdata")
gc()
baseline_pse1 <- mc_genome_tracks(pse_file, rnge, reso_mir, 333, T)
baseline_pse2 <- mc_genome_tracks(pse_file, rnge, reso_mir, 333, T)
baseline_pse3 <- mc_genome_tracks(pse_file, rnge, reso_mir, 334, T)
save(list=ls(), file="wrkspce2.rdata")
gc()

baseline_mir1 <- mc_genome_tracks(mir_file, rnge, reso_mir, 333, T)
baseline_mir2 <- mc_genome_tracks(mir_file, rnge, reso_mir, 333, T)
baseline_mir3 <- mc_genome_tracks(mir_file, rnge, reso_mir, 334, T)
save(list=ls(), file="wrkspce_mir.rdata")
gc()


#collate baselines

rbind(baseline_mel1, baseline_mel2, baseline_mel3) -> baseline_mel
rbind(baseline_yak1, baseline_yak2, baseline_yak3) -> baseline_yak
rbind(baseline_pse1, baseline_pse2, baseline_pse3) -> baseline_pse
rbind(baseline_mir1, baseline_mir2, baseline_mir3) -> baseline_mir
	
#plot function
plot_4c <- function(track, baseline, main_lab, reverse=F, filename=""){
	chr_name <- gsub("-.*", "", names(track)[1])
	locs <- as.numeric(unlist(lapply(names(track), function(x) gsub(".*-", "", x))))
	
	baseline_mean <- colMeans(baseline, na.rm=T)
	baseline_std <- apply(baseline, 2, function(x) sd(x, na.rm=T))
	
	x_ax_lim <- range(locs)
	
	if(nchar(filename) > 0) {
		svg(filename, width=7*2, height=7)
	}	
	if(reverse) {
		x_ax_lim <- rev(x_ax_lim)
		baseline_mean <- rev(baseline_mean)
	}
	
	plot(locs, track, xlab=chr_name, ylab="Contact (Arb. Units)", main=main_lab, pch=".", col="white", xlim=c(x_ax_lim[1], x_ax_lim[2]), bty="n", cex.main=2, cex.lab=1.75, cex.axis=1.5)
	axis(side=1, lwd=2.5, cex.axis=1.5)
	axis(side=2, lwd=2.5, cex.axis=1.5)

	lines(locs, track, lwd=3)
	lines(locs, baseline_mean, lty=2)
	lines(locs, baseline_mean + 0.98*baseline_std, lty=3)
	lines(locs, baseline_mean - 0.98*baseline_std, lty=3) 	
	
	if(nchar(filename) > 0) dev.off()
}

plot_4c(track_mel_hp6, baseline_mel, "New Gene/D. melanogaster (HP6/Umbrea)", filename="mel_hp6.svg")
plot_4c(track_mel_nei, baseline_mel, "Neighbor/D. melanogaster (6-Gene Cluster)", filename="mel_nei.svg")
plot_4c(track_mel_par, baseline_mel, "Parent/D. melanogaster (HP1b)", filename="mel_par.svg")

plot_4c(track_pse_hp6, baseline_pse, "New Gene/D. pseudoobscura (Pre-insertion)", reverse=T, filename="pse_hp6.svg")
plot_4c(track_pse_nei, baseline_pse, "Neighbor/D. pseudoobscura (6-Gene Cluster)", reverse=T, filename="pse_nei.svg")
plot_4c(track_pse_par, baseline_pse, "Parent/D. pseudoobscura (HP1b)", reverse=T, filename="pse_par.svg")
##save(list=ls(), file="wrkspce_pse.rdata")

plot_4c(track_yak_hp6, baseline_yak, "New Gene/D. yakuba (HP6/Umbrea)", filename="yak_hp6.svg")
plot_4c(track_yak_nei, baseline_yak, "Neighbor/D. yakuba (6-Gene Cluster)", filename="yak_nei.svg")
plot_4c(track_yak_par, baseline_yak, "Parent/D. yakuba (HP1b)", filename="yak_par.svg")

plot_4c(track_mir_hp6, baseline_mir, "New Gene/D. miranda (Pre-insertion)", reverse=T, filename="mir_hp6.svg")
plot_4c(track_mir_nei, baseline_mir, "Neighbor/D. miranda (6-Gene Cluster)", reverse=T, filename="mir_nei.svg")
plot_4c(track_mir_par, baseline_mir, "Parent/D. miranda (HP1b)", reverse=T, filename="mir_par.svg")
save(list=ls(), file="wrkspce_mir.rdata")

square_mel_hp6 <- pull_sub_mtx(mel_chr, 4570000, rnge, reso, mel_file, T, square=T) 
square_mel_nei <- pull_sub_mtx(mel_chr, 4710000, rnge, reso, mel_file, T, square=T)
square_yak_hp6 <- pull_sub_mtx(yak_chr, 4680000, rnge, reso, yak_file, T, square=T)
square_yak_nei <- pull_sub_mtx(yak_chr, 4820000, rnge, reso, yak_file, T, square=T)
square_pse_hp6 <- pull_sub_mtx(pse_chr, 29165000, rnge, reso_mir, pse_file, T, square=T)
square_pse_nei <- pull_sub_mtx(pse_chr, 29070000, rnge, reso_mir, pse_file, T, square=T)

square_mir_hp6 <- pull_sub_mtx(mir_chr, 31160000, rnge, reso_mir, mir_file, T, square=T)
square_mir_nei <- pull_sub_mtx(mir_chr, 30970000, rnge, reso_mir, mir_file, T, square=T)

plot_sq <- function(sq_mtx, line_no, delt, file_name){
	svg(file_name)
	heatmap(as.matrix(sq_mtx), Rowv=NA, Colv=NA, add.expr=abline(h=c(line_no - delt, line_no, line_no + delt), v=c(line_no - delt, line_no, line_no + delt)), lwd=rep(2, 6), lty=rep(2, 6))
	dev.off()	
}


plot_sq1 <- function(sq_mtx, file_name){
	png(file_name, width=7, height=7, units="in", res=600)
	heatmap(as.matrix(sq_mtx), Rowv=NA, Colv=NA)
	dev.off()	
}


plot_sq(square_mel_hp6, 41, 24, "sq_mel_hp6_sciadv.svg")
plot_sq(square_yak_hp6, 41, 24, "sq_yak_hp6.svg")
plot_sq(square_pse_hp6, 27, 17, "sq_pse_hp6.svg")
plot_sq(square_mir_hp6, 27, 17, "sq_mir_hp6.svg")

plot_sq1(square_mel_hp6,"sq_mel_hp6_sciadv.png")


#save(list=ls(), file="wrkspce2.rdata")

