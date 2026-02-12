#!/user/bin/nextflow
nextflow.enable.dsl=2

PATTERN_input="./0_faa_per_contig/*faa"
params.outdir = "output_foldseek"


params.dev = false // Lets user testrun nextflow command (by adding flag '--dev') which will have this pipeline run on JUST ONE SAG (again, as a test)
params.num_inputs = 1

//#    FOLDSEEK
PATH_AFUniProt50_descriptions = "/mnt/scgc_nfs/ref/foldseek/final_trembl_annot_cluster_reps_Alphafold-UniProt50_v10.941cd33.tsv.gz" // Derived from trEMBL information on UniProt
PATH_alphafold_descriptions = "/mnt/scgc_nfs/ref/foldseek/Alphafold-UniProt50_h.tsv" // Derived from this file: Alphafold-Uniprot50_h file accompanying foldseek
params.DB_uniprot = "/mnt/scgc_nfs/ref/foldseek/v10.941cd33/Alphafold-UniProt50"
params.DB_uniprot_for_gpu = "/mnt/scgc_nfs/ref/foldseek/v10.941cd33/padded_for_gpu_Alphafold-UniProt50"
params.DB_prost = "/mnt/scgc_nfs/ref/foldseek/v10.941cd33/ProstT5"
// params.annot_script = "/mnt/scgc_nfs/ref/foldseek/scripts/foldseek_anno.py"

workflow {
 CH_ID_faa = Channel.fromPath( PATTERN_input, checkIfExists: true )
          .map { file -> tuple(file.getSimpleName(), file) }

 // If running nextflow with --dev flag, only input a subset of the fastas, e.g. just 1.
 if( params.dev == true ){CH_ID_faa = CH_ID_faa.take ( params.dev ? params.num_inputs : -1) }

 GPU_FOLDSEEK_v10_9(CH_ID_faa)
 FOLDSEEK_ANNOT_W_UNIPROT(GPU_FOLDSEEK_v10_9.out.filter({ it[1].size()>222 }))
 }



process GPU_FOLDSEEK_v10_9{
  tag "${ID}"
  memory = "200.GB" // "50.GB"
  queue = 'gpu'
  cpus = 126 //16 //128 is all the cores, but keep a couple open for interactive jobs
  errorStrategy = 'ignore' //'terminate'
  publishDir "${params.outdir}/raw_foldseek/", mode: 'copy'
  beforeScript 'module load anaconda; source activate /mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/foldseek_v10.9'
  conda '/mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/foldseek_v10.9'
  input: tuple val(ID), path(FAA)
  output:
    tuple val(ID), path("${ID}_AlphafoldUniProt50_hits.tsv.gz"), emit: TSVhits
  shell:
  '''
  # Run the foldseek easy-search command AND record its runtime (saving time results to time.tmp)
  foldseek easy-search \
    !{FAA} \
    !{params.DB_uniprot_for_gpu} \
    headless_results.tsv \
    tmp \
    --prostt5-model !{params.DB_prost} \
    --format-output query,target,bits,evalue,pident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,taxid,taxname,taxlineage \
    --sort-by-structure-bits 0 \
    --remove-tmp-files 1 \
    --prefilter-mode 1 \
    --threads !{task.cpus} \
    --gpu 1
  
  echo -e "query\ttarget\tbits\tevalue\tpident\talnlen\tmismatch\tgapopen\tqstart\tqend\tqlen\ttstart\ttend\ttlen\ttaxid\ttaxname\ttaxlineage" > !{ID}_AlphafoldUniProt50_hits.tsv
  cat headless_results.tsv  >> !{ID}_AlphafoldUniProt50_hits.tsv #!{ID}_AfUniprot50_hits.tsv

  # cleanup
  gzip !{ID}_AlphafoldUniProt50_hits.tsv
  rm headless_results.tsv
  '''}


process FOLDSEEK_ANNOT_W_UNIPROT{
  tag "${ID}"
  container = 'brwnj/kmernorm:v1.0.0'
  errorStrategy = 'terminate'
  publishDir "${params.outdir}/foldseek_w_uniprot_descriptions/", mode: 'copy'
  input: tuple val(ID), path(TSV_hits)
  output: tuple val(ID), path("${ID}_AfUniprot50_annotated_hits.tsv")
  script:
  """
  #!/usr/bin/env python
  import pandas as pd

  print("reading reference database")
  DF_ref = pd.read_csv("${PATH_AFUniProt50_descriptions}", sep="\t", compression="gzip")
  DF_ref.drop(columns=['Unnamed: 0'],inplace=True)

  print("reading hits")
  DF_hits=pd.read_csv("${TSV_hits}", sep="\t")

  print("removing prefix and suffix of hit AF protein model so that it matches the original Uniprot SeqID format")
  print("E.g. 'AF-A0A7S0XJI3-F1-model_v4' -> 'A0A7S0XJI3' ")
  DF_hits['derived_target_name'] = DF_hits['target'].str.replace('^AF-','').str.replace('-F1-model_v4\$','')

  print("Merging with clusterID")
  DF = pd.merge(DF_hits, DF_ref, left_on='derived_target_name', right_on='clusterID', how='left')
  DF.to_csv("${ID}_AfUniprot50_annotated_hits.tsv", sep="\t")
  """ }