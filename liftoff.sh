#! /bin/bash

cd /ru-auth/local/home/ulee/scratch_store/analysis/hi_c/reanalysis

apptLO='srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/hi_c/reanalysis/:/home --no-home /ru-auth/local/home/ulee/scratch_store/docka/liftoff.sif /mambaforge/bin/liftoff'
apptLOTools='srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/hi_c/reanalysis/:/home --no-home  /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif /mambaforge/bin/liftofftools'

wget https://ftp.flybase.net/genomes/dmel/dmel_r6.58_FB2024_03/gtf/dmel-all-r6.58.gtf.gz
gunzip dmel-all-r6.58.gtf.gz

echo "transcript" > liftoff_features.txt
echo "exon" >> liftoff_features.txt
echo "CDS" >> liftoff_features.txt
echo "start_codon" >> liftoff_features.txt
echo "stop_codon" >> liftoff_features.txt

echo "transcript" > liftofftools_features.txt
echo "CDS" >> liftofftools_features.txt

cat comb_gr_nochr.gtf|grep -P '2L\t' > dmel-2L.gtf

$apptLO -polish -f /home/liftoff_features.txt -g /home/dmel-2L.gtf -o /home/dm_to_tei_2L.gff -u /home/dm_to_tei_2L_unmap.txt -dir /home/dm_to_tei_2L /home/tei_hp6chr.fasta /home/dmel-2L.fasta
$apptLO -polish -f /home/liftoff_features.txt -g /home/dmel-2L.gtf -o /home/dm_to_pse_2L.gff -u /home/dm_to_pse_2L_unmap.txt -dir /home/dm_to_pse_2L /home/pse_hp6chr.fasta /home/dmel-2L.fasta


$apptLOTools -r /home/dmel-2L.fasta -t /home/tei_hp6chr.fasta -rg /home/dmel-2L.gtf -tg /home/dm_to_tei_2L.gff_polished -f /home/liftofftools_features.txt -infer-genes -dir /home/dm_to_tei_2L_lot -force all
$apptLOTools -r /home/dmel-2L.fasta -t /home/pse_hp6chr.fasta -rg /home/dmel-2L.gtf -tg /home/dm_to_pse_2L.gff_polished -f /home/liftofftools_features.txt -infer-genes -dir /home/dm_to_pse_2L_lot -force all
