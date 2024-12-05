**Containerized pipeline for DNA sequence analysis**

This Snakemake, that utilizes docker, implements facilitates analysis for detection of single nucleotide variants (SNV), copy number variants (CNV) and structural variants (SV) from DNA-sequence reads.

***Usage***

**Step 1**: clone the pipeline.
Clone the *DASA* pipeline:

```Bash
  git clone [link]
```

**Step 2**: Download required reference files:
Download following genome reference data into your [ref file path] directory.

Reference genome(FASTQ).
- BWA Index.
- GATK Index

Download following GATK reference data into your [ref file path] directory.
- know_dbsnp_vcf ([dbsnp_146.hg38.vcf.gz](wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz), [dbsnp_146.hg38.vcf.gz.tbi](wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi))
- know_dbsnp_vcf (Homo_sapiens_assembly38.known_indels.vcf.gz, Homo_sapiens_assembly38.known_indels.vcf.gz.tbi)
- gold_standard_indels (Mills_and_1000G_gold_standard.indels.hg38.vcf.gz, Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi)
- 1000G_snpshigh_confidence (1000G_phase1.snps.high_confidence.hg38.vcf.gz, 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi)
- interval_list (wgs_calling_regions.hg38.interval_list)
- af_only_gnomad (af-only-gnomad.hg38.vcf.gz, af-only-gnomad.hg38.vcf.gz.tbi)
- gatk_panel_of_normal (1000g_pon.hg38.vcf.gz, 1000g_pon.hg38.vcf.gz.tbi)
- exac_common_knownsite (small_exac_common_3.hg38.vcf.gz, small_exac_common_3.hg38.vcf.gz.tbi)


**Step 3**: Download bed files for CNVkit:
Download the target and antitarget bed files for CNVkit in CNV detection. Store bed files in your [CNVkit bed path] directory and configure the file name and path in config.yaml file (hg38gene_bed: bed_files/target.bed
hg38access_bed: bed_files/antitarget.bed). Alternatively, you can directly use the default bed files provided here.

**Step 4**: Prepare the FASTQ raw data
Make sure there only one sample's FASTQ file located in your [raw FASTQ data path] directory. The FASTQ file name must follow illumina naming convention rule depicted in this website

Eg. SampleName_S1_L001_R1_001.fastq.gz, SampleName_S2_L001_R1_001.fastq.gz


**Step 5**: Build the docker image.

Switch to the docker *DASA* project-folder:

```Bash
  cd ./dasa-master
```

Build the Docker Image (named dasa):

```Bash
  docker build -t dasa .
```

**Step 6**: Finally, you can now run the pipeline:

Check if the snakemake available on docker container (set the reference path,fastq input path and pipeline path in advance).

```Bash
  docker run -v [ref file path]:/data/ref/ \
  -v [raw FASTQ data path]:/data/InputFastqDir/ \
  -v [CNVkit bed path]:/data/bed_files/ \
  -v [local snakemake pipeline path]:/work dasa snakemake -v
```

Dry run of the snakemake DASA.

```Bash
  docker run -v [ref file path]:/data/ref/ \
  -v [raw FASTQ data path]:/data/InputFastqDir/ \
  -v [CNVkit bed path]:/data/bed_files/ \
  -v [local snakemake pipeline path]:/work dasa snakemake -j all --use-conda -n
```


Run pipeline and get the bam and vcf files.

```Bash
docker run -v [ref file path]:/data/ref/ \
-v [raw fastq data path]:/data/InputFastqDir/ \
-v [CNVkit bed path]:/data/bed_files/ \
-v [local snakemake pipeline path]:/work dasa snakemake -j all --use-conda
```

An example set up for running the pipeline could be something like:
*Modify the arguments/data paths according to match your case.*

```Bash
  docker run -v /scratch/denisawany/DASA/refData/:/data/ref/ \
  -v /scratch/denisawany/DASA/NA12878/Fastq_test/:/data/InputFastqDir/ \
  -v /scratch/denisawany/DASA/dasa_docker/dasa-master/bed_files/:/data/bed_files/ \
  -v $(pwd):/work dasa snakemake -j all --use-conda
```
