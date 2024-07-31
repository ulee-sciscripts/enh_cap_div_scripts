#! /bin/bash

apptHomer='srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/hi_c/reanalysis/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/homer.sif'
apptHomerBM='srun -p bigmem apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/hi_c/reanalysis/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/homer.sif'

apptBowtie2='srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/hi_c/reanalysis/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/bowtie2.sif'

apptCircos='srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/hi_c/reanalysis/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/circos.sif'

#mel
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR5820094_1.fastq&
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR5820094_2.fastq&
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR11283142_1.fastq&
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR11283142_2.fastq&

#tei
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR12331760_1.fastq&
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR12331760_2.fastq&

#pse
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR23968954_1.fastq&
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR23968954_2.fastq&
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR11813284_1.fastq&
$apptHomer homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 SRR11813284_2.fastq&

wait

#index
$apptBowtie2 bowtie2-build /home/dmel-all-chromosome-r6.58.fasta dmel1&
$apptBowtie2 bowtie2-build /home/dmel-2L.fasta dmel2L&
$apptBowtie2 bowtie2-build /home/GCF_016746235.2_Prin_Dtei_1.1_genomic.fna dtei&
$apptBowtie2 bowtie2-build /home/GCF_009870125.1_UCI_Dpse_MV25_genomic.fna dpse&

wait

#align
$apptBowtie2 bowtie2 -p 20 -x dmel1 -U SRR5820094_1.fastq.trimmed > SRR5820094_1.mel1.sam&
$apptBowtie2 bowtie2 -p 20 -x dmel1 -U SRR5820094_2.fastq.trimmed > SRR5820094_2.mel1.sam&

$apptBowtie2 bowtie2 -p 20 -x dmel2L -U SRR5820094_1.fastq.trimmed > SRR5820094_1.mel2L.sam&
$apptBowtie2 bowtie2 -p 20 -x dmel2L -U SRR5820094_2.fastq.trimmed > SRR5820094_2.mel2L.sam&

$apptBowtie2 bowtie2 -p 20 -x dtei -U SRR12331760_1.fastq.trimmed > SRR12331760_1.tei.sam&
$apptBowtie2 bowtie2 -p 20 -x dtei -U SRR12331760_2.fastq.trimmed > SRR12331760_2.tei.sam&

$apptBowtie2 bowtie2 -p 20 -x dpse -U SRR23968954_1.fastq.trimmed > SRR23968954_1.pse.sam&
$apptBowtie2 bowtie2 -p 20 -x dpse -U SRR23968954_2.fastq.trimmed > SRR23968954_2.pse.sam&

$apptBowtie2 bowtie2 -p 20 -x dpse -U SRR11813284_1.fastq.trimmed > SRR11813284_1.pse.sam&
$apptBowtie2 bowtie2 -p 20 -x dpse -U SRR11813284_2.fastq.trimmed > SRR11813284_2.pse.sam&

wait

#tag directory
$apptHomer makeTagDirectory dmel2LTagDir/ SRR5820094_1.mel2L.sam,SRR5820094_2.mel2L.sam -tbp 1&
$apptHomerBM makeTagDirectory dmel2LwholeTagDir/ SRR11283142_1.mel2L.sam,SRR11283142_2.mel2L.sam -tbp 1&
$apptHomer makeTagDirectory dteiTagDir/ SRR12331760_1.tei.sam,SRR12331760_2.tei.sam -tbp 1&
$apptHomer makeTagDirectory dpseTagDir/ SRR23968954_1.pse.sam,SRR23968954_2.pse.sam -tbp 1&
$apptHomer makeTagDirectory dpseWholeTagDir/ SRR11813284_1.pse.sam,SRR11813284_2.pse.sam -tbp 1&


wait 

#viz
$apptHomer analyzeHiC dmel2LTagDir/ -pos 2L:4365000-4765000 -res 1000 -window 5000 -balance > hic_hp6_mel.txt&
$apptHomer analyzeHiC dmel2LwholeTagDir/ -pos 2L:4365000-4765000 -res 1000 -window 5000 -balance > hic_hp6_mel_whole.txt&
$apptHomer analyzeHiC dteiTagDir/ -pos NC_053029.1:4510000-4910000 -res 1000 -window 5000 -balance > hic_hp6_tei.txt&
$apptHomer analyzeHiC dpseTagDir/ -pos NC_046681.1:29011000-29411000 -res 1000 -window 5000 -balance > hic_hp6_pse.txt&
$apptHomer analyzeHiC dpseWholeTagDir/ -pos NC_046681.1:29011000-29411000 -res 1000 -window 5000 -balance > hic_hp6_pse_whole.txt&
$apptHomer analyzeHiC dpseWholeTagDir/ -pos NC_046681.1:27100000-29411000 -res 1000 -window 5000 -balance > hic_hp6_pse_whole_bigger.txt&
$apptHomer analyzeHiC dpseTagDir/ -pos NC_046681.1:27100000-29411000 -res 1000 -window 5000  -balance > hic_hp6_pse_larva_bigger.txt&

$apptHomer analyzeHiC dmel2LwholeTagDir/ -pos 2L:3565000-5565000 -res 1000 -window 5000 -balance > hic_hp6_mel_whole_bigger.txt&
$apptHomer analyzeHiC dteiTagDir/ -pos NC_053029.1:3655000-5605000 -res 1000 -window 5000 -balance > hic_hp6_tei_bigger.txt&


$apptHomer analyzeHiC dpseTagDir/ -pos NC_046681.1:27100000-29411000 -res 5000 -window 10000 -balance > hic_hp6_pse_bigger.txt&
$apptHomer analyzeHiC dpseTagDir/ -pos NC_046681.1:25000000-30500000 -res 1000 -window 5000 -balance > hic_hp6_pse_bigger.txt&

wait

#significance

$apptHomer analyzeHiC dmel2LTagDir/ -pos 2L:4365000-4765000 -res 1000 -window 5000 -center -nomatrix -interactions hic_sig_hp6_mel.txt&
$apptHomer analyzeHiC dmel2LwholeTagDir/ -pos 2L:4365000-4765000 -res 1000 -window 5000 -center -nomatrix -interactions hic_sig_hp6_mel_whole.txt&
$apptHomer analyzeHiC dteiTagDir/ -pos NC_053029.1:4510000-4910000 -res 1000 -window 5000 -center -nomatrix -interactions hic_sig_hp6_tei.txt&
$apptHomer analyzeHiC dpseTagDir/ -pos NC_046681.1:29011000-29411000 -res 1000 -window 5000 -center -nomatrix -interactions hic_sig_hp6_pse.txt&
$apptHomer analyzeHiC dpseWholeTagDir/ -pos NC_046681.1:29011000-29411000 -res 1000 -window 5000 -center -nomatrix -interactions hic_sig_hp6_pse_whole.txt&

wait

circos

$apptHomer analyzeHiC dmel2LTagDir/ -pos 2L:4365000-4765000 -res 1000 -window 5000 -minDist 7500 -center -nomatrix -pvalue 0.0005 -circos hic_cir_hp6_mel&
$apptHomer analyzeHiC dmel2LwholeTagDir/ -pos 2L:4365000-4765000 -res 1000 -window 5000 -minDist 7500 -center -nomatrix -pvalue 0.001 -circos hic_cir_hp6_mel_whole&
$apptHomer analyzeHiC dteiTagDir/ -pos NC_053029.1:4510000-4910000 -res 1000 -window 5000 -minDist 7500 -center -nomatrix -pvalue 0.05 -circos hic_cir_hp6_tei&
$apptHomer analyzeHiC dpseTagDir/ -pos NC_046681.1:29011000-29411000 -res 1000 -window 5000 -minDist 7500 -center -nomatrix -pvalue 0.000001 -circos hic_cir_hp6_pse&

$apptHomer analyzeHiC dmel2LwholeTagDir/ -pos 2L:4365000-4765000 -b hp6_mel_4c_nozero.bed -g hp6_dmel_genes_liftoff.txt -res 1000 -window 5000 -minDist 7500 -center -nomatrix -pvalue 0.0001 -circos hic_cir_hp6_mel_whole2&
$apptHomer analyzeHiC dteiTagDir/ -pos NC_053029.1:4510000-4910000 -g hp6_dtei_genes_liftoff.txt -res 5000 -window 10000 -minDist 10000 -center -nomatrix -pvalue 0.05 -circos hic_cir_hp6_tei2&
$apptHomer analyzeHiC dpseTagDir/ -pos NC_046681.1:27100000-29411000 -g hp6_dpse_genes_liftoff.txt -res 5000 -window 10000 -minDist 10000 -center -nomatrix -pvalue 0.05 -circos hic_cir_hp6_pse2&
$apptHomer analyzeHiC dpseWholeTagDir/ -pos NC_046681.1:27100000-29411000 -g hp6_dpse_genes.txt -res 1000 -window 5000 -minDist 5000 -center -nomatrix -pvalue 0.05 -circos hic_cir_hp6_pse3&
$apptHomer analyzeHiC dpseWholeTagDir/ -pos NC_046681.1:27100000-29411000 -g hp6_dpse_genes_liftoff.txt -res 1000 -window 5000 -minDist 7500 -center -nomatrix -pvalue 0.05 -circos hic_cir_hp6_pse3&

wait

manually add posn rule into circos.conf
			<rule>
				importance = 250
mel
				condition  = ! on(2L, 4576000, 4578000)
tei
				condition  = ! on(NC_053029.1, 4710000, 4712000)
pse
				condition  = ! on(NC_046681.1, 29209000, 29211000)
				show       = no
			</rule>

$apptCircos /opt/circos/bin/circos -conf hic_cir_hp6_mel_whole2.circos_4c.conf&
$apptCircos /opt/circos/bin/circos -conf hic_cir_hp6_mel_whole2.circos_4c_flee1.conf&
$apptCircos /opt/circos/bin/circos -conf hic_cir_hp6_tei2.circos_4c.conf&
$apptCircos /opt/circos/bin/circos -conf hic_cir_hp6_tei2.circos_4c_flee1.conf&
$apptCircos /opt/circos/bin/circos -conf hic_cir_hp6_pse2.circos_4c.conf&





wait

