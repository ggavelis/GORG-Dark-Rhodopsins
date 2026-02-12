# GORG-Dark-Rhodopsins
This github describes how rhodopsins from the dark ocean were analyzed for the study "Single-cell-resolved genome atlas of prokaryoplankton inhabiting the ocean’s interior" by Chang et al.

**PRIOR TO ANALYSIS IN THIS GITHUB**<br>
Two steps were taken to identify rhodopsin candidates. Quotes reflect methods passages from Chang et al.
1. **DRAM to identify preliminary rhodopsins**
Each Single Amplified Genome (SAG) was subjected to DRAM v1.4 to predict its CDS, proteins, and their annotations. (DRAM was run in default mode, configured with its core set of curated HMM databases (i.e., without KEGG or other customized databases), and applied the additional parameter --min_contig_size 2000 to skip annotating contigs shorter than 2 kbp.
2. **Rhodopsin validation via alignment and phylogeny**
"We further validated these via amino acid alignment against two atomic structure-resolved rhodopsins: pdb_00001c3w and pdb_00002l6x, which have leucine and glutamine as their key tuning residues, respectively. We aligned the proteins using MAFFT37 with default parameters. Each GORG-Dark rhodopsin in the alignment was inspected for the presence or absence of canonical rhodopsin residues, such as (1) the lysine (position 591 in the alignment) in the seventh transmembrane helix to which retinal is typically bound, (2) the leucine or glutamine  (position 390) associated with spectral tuning for the absorption of blue or green or blue light, and (3) the DTD/DTE or TSA motif (positions 382, 386, and 393) that indicates proton pumping or chloride activity. We discarded rhodopsin candidates that could not be aligned or that were too short to contain phylogenetic signal." 

**PERFORMED IN THIS GITHUB:**
\[Brackets refer to github folders where the step occurred\]<br>
<ins>Rich annotation of SAGs/contigs with rhodopsin</ins>:<br>
>Moving forward with only high-quality rhodopsin proteins, we sought further clues to their cellular function by examining the genes that flanked each rhodopsin. Their proteins had already been annotated by DRAM—which emphasizes metabolism—so we took an additional step to annotate the proteins for broader functions, using a structure-based search. For each unannotated protein, we used Foldseek v5.4 to query that protein against the Alphafold/UniProt50 database 49 **[3, 4]**. If protein hits were found, we took the description of the top hit (by bitscore) and used that as the annotation. Our third and final annotation step was a targeted search for functions known to flank rhodopsin genes: the carotenoid/retinal genes CrtE, CrtB, CrtI, CrtY, brp/blh, and genes encoding flotillin **[1, 2]**. We used the command hmmscan from the program HMMER3 to search reference HMMs against our GORG-Dark proteins 50. The reference HMMs were downloaded from NCBI’s Hidden Markov Models database, and we used HMMs NF041003 and NF045549 to search for the pigment biosynthesis gene CrtE, NF042419 and NF045686 for CrtB, TIGR02734 for CrtI, TIGR01789 and TIGR01790 for CrtY, and TIGR03753 for the retinal biosynthesis gene brp/blh. We obtained an HMM for flotillin from 51.
>
<ins>Characterizing "neighborhood" of proteins within a 5 gene radius of each rhodopsin to look for conserved operons</ins>

> Having annotated the genes/proteins that neighbor rhodopsins, we then looked for functional enrichments within 5 genes on either side of each rhodopsin gene. To facilitate this, each rhodopsin gene was reoriented to run 5’ to 3’ (regardless of whether it was encoded on the + or - strand) and was assigned position 0, while proteins encoded upstrand of rhodopsin were given positions -1,-2,-3,-4,-5  (-5 is furthest up) and proteins encoded downstrand of rhodopsin were given positions +1,+2,+3,+4,+5 (+5 is furthest down) **[5]**. We tabulated the annotations at positions -5 through +5, then inspected the table semi-manually in Jupyter notebook **[6]**. For example, inspecting column -1, we asked “what are the top 20 most common annotations immediately downstrand of rhodopsin?”, and received a list of human-readable descriptions from DRAM, foldseek and HMMsearch hits. For redundant annotations (e.g. “AMP-binding enzyme PF00501.31” and “predicted AMP-binding domain [PF00501]”), we collapsed them into a shared substring (“PF00501”). To find the most enriched functions in flanking genes, we examined the 20 most common substrings at the +1 and -1 positions. Then, to see if each flanking gene was part of a conserved gene block or operon, we walked further outward. For instance, if we noticed many annotations containing “[F/l]lotillin” at position +1, we would look for which annotations were most common down strand of “[F/l]lotillin” at positions +2, +3, +4, and +5. 
>

<ins>Tabulate alphafold information (i.e. whether each rhodopsin had an 8th helix)</ins>

> "Beyond their phylogenomic context, predicting protein function also requires consideration of 3D structures. While in-silico approaches are at best heuristic, we used AlphaFold v3.0.1 52 to generate a preliminary model of each GORG-Dark rhodopsin. We visually inspected these 3D models in the webserver PAE Viewer to look for structural features that could be distinctive for GORG Dark rhodopsin clades 53. If an eight helix was found at the N-terminus of a rhodopsin, we considered it as ‘present’ if it had at least two helical rotations predicted by AlphaFold **[7]**."
>

**FINAL DATA PRODUCT**<br>
[Supplementary Table 8](https://github.com/ggavelis/GORG-Dark-Rhodopsins/blob/main/Table_S8_Rhodopsin_sequences_metadata_and_gene_neighborhoods.tsv)<br>
Which provides functional and evolutionary information about GORG-Dark rhodopsins,<br>
as visualized in [Figure 4] (https://itol.embl.de/tree/12812822615690281732644278#)<br>

