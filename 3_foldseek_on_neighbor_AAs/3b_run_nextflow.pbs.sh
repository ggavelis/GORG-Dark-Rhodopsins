#PBS -N foldseek_on_rhodopsin_neighbors
#PBS -m bea -j oe -o status_of_assembly.pbs
#PBS -l ncpus=1,mem=50MB,walltime=120:00:00
module purge
module use /mnt/scgc_nfs/opt/modulefiles/common/
module load nextflow
module load singularity
NXF_VER=25.04.6 nextflow run foldseek_on_rhodopsin_neighbors.nf -profile charlie --viral --foldseek
echo "done"
