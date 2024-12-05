**Generation of CNVkit target/antitarget bed file**

The CNVkit command can be used to generate sequence-accessible coordinates file.

We'll build the target/antitarget bed file based on the hg38 (GRCh38) human genome, which is the latest annotation. To do so, we need the following reference files: (1) the gene annotation flat file, (2) the exome target bed file, and (3) the excluded bed file. Download these required reference files and generate the target/antitarget bed file as follows:

**Step 1:** ***Download required reference files***
1. Download the genome annotation file for hg38. This can be obtained from the UCSC database.<br>
   *refFlat.txt for hg38  [refFlat.txt.gz](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz)*
3. Download the exome target bed file. We can can download it from Twist Bioscience.<br>
     *Twist_Exome_Core_Covered_Targets_hg38.bed [Twist_Exome_Core_Covered_Targets_hg38.bed](https://www.twistbioscience.com/sites/default/files/resources/2022-01/Twist_Exome_Core_Covered_Targets_hg38.bed)*
4. Download the excluded bed file. We can download it from ENCODE database.<br>
    *GRCh38_unified_blacklist.bed [GRCh38_unified_blacklist.bed](https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz)*

**Step 2:** ***Generate the target/antitarget bed file***
1. Generate access.hg38.bed by executing the code:
```Python
 cnvkit.py access hg38.fa -x GRCh38_unified_blacklist.bed -o access-excludes.hg38.bed
 ```
 ```Python
 cnvkit.py access /scratch/DenisAwany/DASA_pipeline/RefData/GRCh38.d1.vd1.fa -x GRCh38_unified_blacklist.bed -o access-excludes.hg38.bed
 ```

2. Generate target/antitarget bed file by running the code:
```Python
 cnvkit.py target Twist_Exome_Core_Covered_Targets_hg38.bed --annotate hg38_refFlat.txt -o target_hg38.bed --short-names
```
 ```Python
  cnvkit.py antitarget target_hg38.bed -g access-excludes.hg38.bed -o hg38_antitarget.bed
```

