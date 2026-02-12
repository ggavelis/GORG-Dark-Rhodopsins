#PBS -N hmmscan_2b 
#PBS -M ggavelis@bigelow.org
#PBS -m bea
#PBS -m bea -j oe -o status_2b.pbs
#PBS -l ncpus=2,mem=10GB,walltime=24:55:00

set -e
cd /mnt/stepanauskas_nfs/projects/gorg-dark/FArhodopsin_gene_clusters/1_hmmsearch_for_carotenoids_and_flotillins/2_crtB/

module purge
module load  hmmer/3.1b2
hmmscan --tblout 2b_gdark4.hmmscan.txt NF045686.1.HMM /mnt/stepanauskas_nfs/projects/gorg-dark/FArhodopsin_gene_clusters/0_get_DRAM_AAs_for_SAGs_w_rhodopsins/output/897_gdark_SAGs_DRAM.faa

echo "done"
