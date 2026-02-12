#PBS -N hmmscan_3 
#PBS -M ggavelis@bigelow.org
#PBS -m bea
#PBS -m bea -j oe -o status_3.pbs
#PBS -l ncpus=2,mem=10GB,walltime=24:55:00

set -e
cd /mnt/stepanauskas_nfs/projects/gorg-dark/FArhodopsin_gene_clusters/1_hmmsearch_for_carotenoids_and_flotillins/3_crtI/

module purge
module load  hmmer/3.1b2
hmmscan --tblout 3_gdark4.hmmscan.txt TIGR02734.1.HMM /mnt/stepanauskas_nfs/projects/gorg-dark/FArhodopsin_gene_clusters/0_get_DRAM_AAs_for_SAGs_w_rhodopsins/output/897_gdark_SAGs_DRAM.faa

echo "done"
