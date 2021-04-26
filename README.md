Human and SARS-CoV-2 RNA-Seq workflow that:

* Aligns to human and SARS-CoV-2 genomes
* Quantifies gene expression using RSEM (bulk and Smart-Seq2) or cellranger count (10x single cell)
* Jointly recalibrates base quality scores using reads from human and SARS-CoV-2 genomes
* Calls variants from SARS-CoV-2 and human reads
* Assembles SARS-CoV-2 genome using Trinity
* Aligns assembled genome to reference using minimap2
* Computes UMAP embedding, clusters, finds differentially expressed genes, and provides visualization
  using [cumulus]((https://cumulus.readthedocs.io/)) (10x or Smart-Seq2 single cell data)

## Examples

Examples for 10x and Smart-Seq2 single cell gene expression and for bulk RNA-Seq data are provided below.

Both workflows are preconfigured to use a joint human and SARS-COV-2 RNA genome reference. The SARS-COV-2 genome was
generated by [Carly Ziegler](http://shaleklab.com/author/carly/). The [viral sequence](https://raw.githubusercontent.com/broadinstitute/rna_seq_sars_cov_2/main/resources/SARSCoV2.fa)
and [gtf](https://raw.githubusercontent.com/broadinstitute/rna_seq_sars_cov_2/main/resources/SARSCoV2.gtf) are as described in 
[Kim et al. Cell 2020] (https://github.com/hyeshik/sars-cov-2-transcriptome, BetaCov/South Korea/KCDC03/2020 based on
NC_045512.2). The GTF was edited to include only CDS regions, and regions were added to describe the 5’ UTR (
“SARSCoV2_5prime”), the 3’ UTR (“SARSCoV2_3prime”), and reads aligning to anywhere within the Negative Strand(
“SARSCoV2_NegStrand”). Additionally, trailing A’s at the 3’ end of the virus were excluded from the SARSCoV2 FASTA, as
these were found to drive spurious viral alignment in non-SARS-COV-2 infected samples.


### 10X Single-Cell:

FASTQs from the
paper [Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19](https://www.nature.com/articles/s41591-020-0901-9)
are available from the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).

Example Input:

```json
  {
  "rna_seq_sars_cov_2_workflow.cumulus.max_genes": 6000,
  "rna_seq_sars_cov_2_workflow.cumulus.min_genes": 200,
  "rna_seq_sars_cov_2_workflow.cumulus.min_umis": 1000,
  "rna_seq_sars_cov_2_workflow.cumulus.percent_mito": 10,
  "rna_seq_sars_cov_2_workflow.cumulus.perform_de_analysis": true,
  "rna_seq_sars_cov_2_workflow.cumulus.plot_umap": "SARS,louvain_labels",
  "rna_seq_sars_cov_2_workflow.read1": [
    "SRR11181956_1.fastq"
  ],
  "rna_seq_sars_cov_2_workflow.read2": [
    "SRR11181956_2.fastq"
  ],
  "rna_seq_sars_cov_2_workflow.sample_id": "GSM4339771",
  "rna_seq_sars_cov_2_workflow.type": "10x",
  "rna_seq_sars_cov_2_workflow.ref_dict": "GRCh37_SARSCoV2.dict",
  "rna_seq_sars_cov_2_workflow.ref_fasta": "GRCh37_SARSCoV2.fa",
  "rna_seq_sars_cov_2_workflow.ref_fasta_index": "GRCh37_SARSCoV2.fa.fai",
  "rna_seq_sars_cov_2_workflow.reference": "GRCh37_SARSCoV2",
  "rna_seq_sars_cov_2_workflow.viral_ref_fasta": "resources/SARSCoV2.fa"
}
```

#### Output:

* Web Summary.html from cellranger count
  ![cellranger summary](https://raw.githubusercontent.com/broadinstitute/rna_seq_sars_cov_2/main/images/cellranger_summary.png)
  ![cellranger analysis](https://raw.githubusercontent.com/broadinstitute/rna_seq_sars_cov_2/main/images/cellranger_analysis.png)

* BAM file containing position-sorted reads aligned to human and SARS-CoV-2 genome
* Feature-barcode matrix containing counts for both human and SARS-CoV-2 aligned reads
* Molecule info file containing information for all molecules that contain high confidence valid barcodes and UMIs
* SARS-CoV-2 aligned BAM
* Human aligned BAM
* Assembled Trinity FASTA from SARS-CoV-2 reads
* BAM of Trinity assembly aligned to reference SARS-CoV-2 genome
* IGV SARS-CoV-2 HTML report
* Human SNPs
* SARS-CoV-2 SNPs
* h5ad file containing analysis results, including low quality cell filtration, highly variable gene selection,
  dimension reduction graph-based clustering, differentially expressed genes, UMAP visualization. Generated
  using [cumulus](https://cumulus.readthedocs.io/)
  ![clusters](https://raw.githubusercontent.com/broadinstitute/rna_seq_sars_cov_2/main/images/clusters.png)
  ![SARS gene expression](https://raw.githubusercontent.com/broadinstitute/rna_seq_sars_cov_2/main/images/SARS-expression.png)

### Bulk RNA-Seq:

FASTQs for the sample CRR119894 from the
paper [Transcriptomic characteristics of bronchoalveolar lavage fluid and peripheral blood mononuclear cells in COVID-19 patients](https://www.tandfonline.com/doi/full/10.1080/22221751.2020.1747363)
are available from the [Genome Sequence Archive](https://bigd.big.ac.cn/gsa/browse/CRA002390/CRR119894).

Example Input:

```json

{
  "rna_seq_sars_cov_2_workflow.read1": [
    "CRR119894_f1.fq.gz"
  ],
  "rna_seq_sars_cov_2_workflow.read2": [
    "CRR119894_r2.fq.gz"
  ],
  "rna_seq_sars_cov_2_workflow.type": "bulk",
  "rna_seq_sars_cov_2_workflow.sample_id": "CRR119894",
  "rna_seq_sars_cov_2_workflow.ref_dict": "GRCh37_SARSCoV2.dict",
  "rna_seq_sars_cov_2_workflow.ref_fasta": "GRCh37_SARSCoV2.fa",
  "rna_seq_sars_cov_2_workflow.ref_fasta_index": "GRCh37_SARSCoV2.fa.fai",
  "rna_seq_sars_cov_2_workflow.viral_ref_fasta": "SARSCoV2.fa",
  "rna_seq_sars_cov_2_workflow.reference": "GRCh37_SARSCoV2",
  "rna_seq_sars_cov_2_workflow.db_snp_vcf": "dbsnp.vcf.gz",
  "rna_seq_sars_cov_2_workflow.gnomad_vcf_index": "gnomad-lite.vcf.gz.csi",
  "rna_seq_sars_cov_2_workflow.repeat_mask_bed": "repeats_ucsc_gb.bed.gz",
  "rna_seq_sars_cov_2_workflow.gnomad_vcf": "gnomad-lite.vcf.gz",
  "rna_seq_sars_cov_2_workflow.rna_editing_vcf": "RNAediting.library.vcf.gz",
  "rna_seq_sars_cov_2_workflow.ref_splice_adj_regions_bed": "ref_annot.splice_adj.bed.gz",
  "rna_seq_sars_cov_2_workflow.rna_editing_vcf_index": "RNAediting.library.vcf.gz.csi",
  "rna_seq_sars_cov_2_workflow.db_snp_vcf_index": "dbsnp.vcf.gz.tbi"
}
```

#### Bulk RNA-Seq Output

![SARS-CoV-2_VARIANTS](https://raw.githubusercontent.com/broadinstitute/rna_seq_sars_cov_2/main/images/bulk_example.png)
* BAM file containing position-sorted reads aligned to human and SARS-CoV-2 genome
* Gene quantification from RSEM for both human and SARS-CoV-2 aligned reads
* SARS-CoV-2 aligned BAM
* Human aligned BAM
* Assembled Trinity FASTA from SARS-CoV-2 reads
* BAM of Trinity assembly aligned to reference SARS-CoV-2 genome
* IGV SARS-CoV-2 HTML report
* Human SNPs
* SARS-CoV-2 SNPs


### Smart-Seq2:

The Smart-Seq2 pipeline is very similar to the bulk pipeline. In addition to the bulk outputs, 
[cumulus](https://cumulus.readthedocs.io/) is used for data analysis and visualization. FASTQs from the paper
[Paracrine signalling by cardiac calcitonin controls atrial fibrogenesis and arrhythmia](https://www.nature.com/articles/s41586-020-2890-8)
are available from the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra). For Smart-Seq2, we give an array of
FASTQS as input (one FASTQ per cell). Note that in this example, read2 is set to an empty array because the reads are unpaired.
In this example, we added bulk reads aligning to the SARS-CoV-2 genome to the Smart-Seq2 data to illustrate the pipeline (We internally tested the pipeline on human SARS-CoV-2 patients, but this data is not publicly available).

Example Input:

```json

{
  "rna_seq_sars_cov_2_workflow.read1": [
    "SRR11538647_1.fastq","SRR11538648_1.fastq", "..."
  ],
  "rna_seq_sars_cov_2_workflow.read2": [],
  "rna_seq_sars_cov_2_workflow.type": "smartseq2",
  "rna_seq_sars_cov_2_workflow.cumulus.run_umap": false,
  "rna_seq_sars_cov_2_workflow.cumulus.max_genes": 10000,
  "rna_seq_sars_cov_2_workflow.cumulus.min_genes": 0,
  "rna_seq_sars_cov_2_workflow.cumulus.min_umis": 0,
  "rna_seq_sars_cov_2_workflow.sample_id": "GSE148506",
  "rna_seq_sars_cov_2_workflow.ref_dict": "GRCh37_SARSCoV2.dict",
  "rna_seq_sars_cov_2_workflow.ref_fasta": "GRCh37_SARSCoV2.fa",
  "rna_seq_sars_cov_2_workflow.ref_fasta_index": "GRCh37_SARSCoV2.fa.fai",
  "rna_seq_sars_cov_2_workflow.viral_ref_fasta": "SARSCoV2.fa",
  "rna_seq_sars_cov_2_workflow.reference": "GRCh37_SARSCoV2",
  "rna_seq_sars_cov_2_workflow.db_snp_vcf": "dbsnp.vcf.gz",
  "rna_seq_sars_cov_2_workflow.gnomad_vcf_index": "gnomad-lite.vcf.gz.csi",
  "rna_seq_sars_cov_2_workflow.repeat_mask_bed": "repeats_ucsc_gb.bed.gz",
  "rna_seq_sars_cov_2_workflow.gnomad_vcf": "gnomad-lite.vcf.gz",
  "rna_seq_sars_cov_2_workflow.rna_editing_vcf": "RNAediting.library.vcf.gz",
  "rna_seq_sars_cov_2_workflow.ref_splice_adj_regions_bed": "ref_annot.splice_adj.bed.gz",
  "rna_seq_sars_cov_2_workflow.rna_editing_vcf_index": "RNAediting.library.vcf.gz.csi",
  "rna_seq_sars_cov_2_workflow.db_snp_vcf_index": "dbsnp.vcf.gz.tbi"
}
```

#### Output

* Aggregated outputs
  * TPM-normalized count matrix
  * h5ad file containing analysis results, including low quality cell filtration, highly variable gene selection,
    dimension reduction graph-based clustering, differentially expressed genes, UMAP visualization. Generated
    using [cumulus](https://cumulus.readthedocs.io/)
    
* Per cell outputs
  * BAM file containing position-sorted reads aligned to human and SARS-CoV-2 genome
  * Gene quantification from RSEM for both human and SARS-CoV-2 aligned reads
  * SARS-CoV-2 aligned BAM
  * Human aligned BAM
  * Assembled Trinity FASTA from SARS-CoV-2 reads
  * BAM of Trinity assembly aligned to reference SARS-CoV-2 genome
  * IGV SARS-CoV-2 HTML report
  * Human SNPs
  * SARS-CoV-2 SNPs

### Software Installation

- Prerequisites
    - [Java8](https://www.oracle.com/java/technologies/java8.html)
    - [Docker](https://www.docker.com/)
    - [Cromwell](https://cromwell.readthedocs.io/) workflow runner
  
- Clone this repository:

  ```
  git clone https://github.com/broadinstitute/rna_seq_sars_cov_2.git
  ```



### Terra

- Workflows are available on [Terra](https://support.terra.bio/) if you prefer running on the cloud.


### Genome Installation

- Download plug-n-play CTAT genome lib The pipeline integrates with the standard CTAT genome lib that's leveraged as
  part of the [Trinity CTAT project](https://github.com/NCIP/Trinity_CTAT/wiki). Instructions below describe how to set
  up your CTAT genome lib and add your viral genome.

  Trinity CTAT supports both human genome builds hg19/37 and hg38, so you can choose whichever version you prefer to
  operate from, or separately install both so you have flexibility to use either.

  Download your selected plug-n-play CTAT genome lib from
  the [Trinity CTAT Data Repository](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/).

  > Note, if you already have a CTAT genome lib installed for use with other Trinity CTAT project utilities, feel free to use the genome lib you already have.

  Example:

  ```
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play.tar.gz
  ```

- Download your viral genome of interest. The SARS-COV-2 genome FASTA is available 
  [here](https://raw.githubusercontent.com/broadinstitute/rna_seq_sars_cov_2/main/resources/SARSCoV2.fa) and the GTF 
  [here](https://raw.githubusercontent.com/broadinstitute/rna_seq_sars_cov_2/main/resources/SARSCoV2.gtf).
  
  
- Create JSON input file for rna_seq_sars_cov_2_create_reference workflow:

  ```json
  {
  "rna_seq_sars_cov_2_create_reference.reference_name": "GRCh37_SARSCoV2",
  "rna_seq_sars_cov_2_create_reference.fasta": [
    "ctat_genome_lib_build_dir/ref_genome.fa",
    "SARSCoV2.fa"
  ],
  "rna_seq_sars_cov_2_create_reference.gtf": [
   "ctat_genome_lib_build_dir/ref_annot.gtf",
    "SARSCoV2.gtf"
  ],
  "rna_seq_sars_cov_2_create_reference.type": "smartseq2"
  }
  ```

  Change rna_seq_sars_cov_2_create_reference.type to "10x" to create a 10x compatible reference (references for Smart-Seq2 and bulk are identical).

- Run rna_seq_sars_cov_2_create_reference workflow:

    ```
    java -jar cromwell-55.jar run -i inputs.json rna_seq_sars_cov_2_create_reference.wdl
    ```
