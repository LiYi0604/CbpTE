# A pipeline for ResequenceTE Matrix Construction
This pipeline is the downstream step after PanGenie genotyping.

1. `0-script_to_construct_reseq_matrix.sh` <br> The pipeline for ReseqTE matrix construction.
 
2. `0-pbs.matrix.sh` <br> Script for submitting jobs to PBS.
 
3. `2-pbs.coor_switch.sh` <br> Script for running jobs on PBS.

4. `2-coor_switch.e.py` <br> Aligns PanGenie output to the PanTE matrix by reference site.

5. `3-map.*.py` <br> Maps population-level PanTE matrix sites to the global PanTE matrix to support subsequent alignment of genotyping results.

6. `4-pop2matrix.reseq.ver2.py` <br> Aligns PanGenie output to the population-level PanTE matrix to produce resequencing-based genotypes. <br>
   To maximize genotyping accuracy, based on phylogenetic relationships, each locus in a target genome is inferred from the genotype state of the reference genome of the most closely related population in the PanTE matrix.

7. `5-precision.ipynb` <br> Uses 19 genomes with both assemblies and resequencing data to evaluate PanGenie’s genotyping performance and the ReseqTE-matrix construction method.

8. `6-matrix_merge.py` <br> Concatenates all outputs to generate a species-level ReseqTE matrix.

9. `7-reseq_mat_to_vcf.py` & `7-reseq_vcf_filter` <br> Converts the matrix to a VCF and filters sites with >10% missing rate.

10. `8-reseq_all_sample_count.ipynb` & `8-reseq_mat_site_freq.ipynb` <br> Check individual missing rates for each species and generates histograms of TE counts per species.

11. `9-reseq_*.py` <br> Classifies TE sites in the ReseqTE matrix using phylogenetic relationships. <br>
    Which defined TE sites into three types:<br>
    (1) shared in all species; <br>
    (2) shared only between tetraploid subgenome and diploid ancestor genome; <br>
    (3) Unique to each genomes.
