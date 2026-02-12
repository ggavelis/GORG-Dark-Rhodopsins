#PBS -N hmmscan_6
#PBS -M ggavelis@bigelow.org
#PBS -m bea
#PBS -m bea -j oe -o status_6.pbs
#PBS -l ncpus=16,mem=10GB,walltime=24:55:00

set -e
cd /mnt/stepanauskas_nfs/projects/gorg-dark/FArhodopsin_gene_clusters/1_hmmsearch_for_carotenoids_and_flotillins/6_flotillin

module purge
module load  hmmer/3.1b2
hmmscan --cpu 16 --tblout 6_gdark4.hmmscan.txt Flotillin.tagged.hmm /mnt/stepanauskas_nfs/projects/gorg-dark/FArhodopsin_gene_clusters/0_get_DRAM_AAs_for_SAGs_w_rhodopsins/output/897_gdark_SAGs_DRAM.faa

echo "done"
