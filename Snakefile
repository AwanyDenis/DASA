configfile:"config.yaml"

####   Set VARIABLES
INPUTDIR = config["input_dir"]
OUTPUTDIR = config["output_dir"]
GENOMEDIR = config["genome_dir"]

# Genome_LINK = config["genome_link"]
# BWA_INDEXLINK = config["BWA_Index"]
# GATK_INDEXLINK = config["GATK_Index"]
Genome_NAME = config["params"]["genome_name"]
Thread_NUM = config["params"]["thread_num"]
SAMPLE_NAME,SAMPLE_NUM,LANE,READ,TAIL = glob_wildcards(INPUTDIR + "/{sample_name}_{sample_num}_{lane}_{read}_{tail}.fastq.gz")

SAMPLE_NAME_SET = set(SAMPLE_NAME)
SAMPLE_NUM_SET = set(SAMPLE_NUM)
SET_LANE = set(LANE)
SET_READ = set(READ)
SET_TAIL = set(TAIL)
tail = "".join(SET_TAIL)  ### specify the tail string for later use

SAMPLE = "".join(SAMPLE_NAME_SET)+"_"+"".join(SAMPLE_NUM_SET) +"_"+"".join(SET_TAIL)
SAMPLE_OUTPUT_DIR = OUTPUTDIR+"/"+ SAMPLE

rule all:
   input:
        genome_dir = GENOMEDIR,
        
        # filtered_vcf = expand(SAMPLE_OUTPUT_DIR + "/filterd_calls/{sample_name}_{sample_num}.snv.vcf.gz",sample_name=SAMPLE_NAME_SET, sample_num = SAMPLE_NUM_SET),
        # sv_vcf = expand(SAMPLE_OUTPUT_DIR + "/lumpySV/{sample_name}_{sample_num}.sv.vcf",sample_name=SAMPLE_NAME_SET, sample_num = SAMPLE_NUM_SET),
        # cnv_vcf = expand(SAMPLE_OUTPUT_DIR + "/CNV/{sample_name}_{sample_num}.cnv.vcf",sample_name=SAMPLE_NAME_SET, sample_num = SAMPLE_NUM_SET),
        concat_vcf = expand(SAMPLE_OUTPUT_DIR + "/SNVSVCNVConcat/{sample_name}_{sample_num}.concat.vcf",sample_name=SAMPLE_NAME_SET, sample_num = SAMPLE_NUM_SET),


# ---------------------------------------
# Quality control on raw reads
# ---------------------------------------
rule trim_galore_pe:
    input:
        [INPUTDIR +"/{sample}_{lane}_R1_"+tail+".fastq.gz", INPUTDIR +"/{sample}_{lane}_R2_"+tail+".fastq.gz"]
    output:
        temp(SAMPLE_OUTPUT_DIR+"/trimmed/{sample}_{lane}_R1_"+tail+"_val_1.fq.gz"),
        temp(SAMPLE_OUTPUT_DIR+"/trimmed/{sample}_{lane}_R1_"+tail+".fastq.gz_trimming_report.txt"),
        temp(SAMPLE_OUTPUT_DIR+"/trimmed/{sample}_{lane}_R2_"+tail+"_val_2.fq.gz"),
        temp(SAMPLE_OUTPUT_DIR+"/trimmed/{sample}_{lane}_R2_"+tail+".fastq.gz_trimming_report.txt")
    params:
        TrimmingQualityReadEnds = 20,
        TrimmingReadLengthMin = 20,
        TrimmingAdaptor_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        TrimmingAdaptor_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        output_dir = SAMPLE_OUTPUT_DIR+"/trimmed/",       ### remember change the output_dir when change the output file path.
        stringency = 5,
        extra=""
        # extra="--illumina"
    conda:
        "envs/trim.yaml"
    log:
        SAMPLE_OUTPUT_DIR+"/logs/trim_galore/{sample}_{lane}.log"
    shell:
        "trim_galore {params.extra} -q {params.TrimmingQualityReadEnds} --length {params.TrimmingReadLengthMin}  --stringency {params.stringency} "
        # "-a {params.TrimmingAdaptor_R1} -a2 {params.TrimmingAdaptor_R2} "
        "--illumina "
        "--paired {input} -o {params.output_dir} --fastqc"
        # "trim_galore --qual "${TrimmingQualityReadEnds}" --gzip --length "${TrimmingReadLengthMin}" --"${TrimmingAdaptor}" --paired "${Read1}" "${Read2}" --output_dir "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/ --fastqc"


# ---------------------------------------
# Read sequence mapping
# ---------------------------------------
rule bwa_map:
    input:
        expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        # expand("trimmed/{sample_name}_{sample_num}_{lane}_R1_001_test_val_1.fq.gz",sample_name=SAMPLE_NAME_SET, 
        # sample_num = SAMPLE_NUM_SET,lane = SET_LANE),
        # expand("trimmed/{sample_name}_{sample_num}_{lane}_R2_001_test_val_2.fq.gz",sample_name=SAMPLE_NAME_SET, 
        # sample_num = SAMPLE_NUM_SET,lane = SET_LANE),
        SAMPLE_OUTPUT_DIR+"/trimmed/{sample}_{lane}_R1_"+tail+"_val_1.fq.gz",
        SAMPLE_OUTPUT_DIR+"/trimmed/{sample}_{lane}_R2_"+tail+"_val_2.fq.gz"
    output:
        protected(SAMPLE_OUTPUT_DIR+"/mapped_reads/{sample}_{lane}.bam")
    conda:
        "envs/mapping.yaml"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA",
        thersholds=0
    log:
        SAMPLE_OUTPUT_DIR+"/logs/bwa_mem/{sample}_{lane}.log"
    # benchmark:
    #     "benchmarks/{sample}.bwa.benchmark.txt"
    threads: Thread_NUM
    shell:
        "(bwa mem -t {threads} -K 100000000 -T {params.thersholds} -R '{params.rg}' {input} | "
        "samtools view -Shb - > {output}) 2> {log}"

rule samtools_merge:
    input:
        lane1 = SAMPLE_OUTPUT_DIR+"/mapped_reads/{sample}_L001.bam",
        lane2 = SAMPLE_OUTPUT_DIR+"/mapped_reads/{sample}_L002.bam",
        lane3 = SAMPLE_OUTPUT_DIR+"/mapped_reads/{sample}_L003.bam",
        lane4 = SAMPLE_OUTPUT_DIR+"/mapped_reads/{sample}_L004.bam"
    output:
        temp(SAMPLE_OUTPUT_DIR+"/merged/{sample}.bam")
    params:
        "" # optional additional parameters as string
    threads: Thread_NUM     # This value - 1 will be sent to -@
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools merge -@ {threads} {params} {output} {input.lane1} "
        "{input.lane2} "
        "{input.lane3} "
        "{input.lane4}"

# ---------------------------------------
# Sorting & indexing mapped reads, and marking duplicates
# ---------------------------------------
rule samtools_sort:
    input:
        SAMPLE_OUTPUT_DIR+"/merged/{sample}.bam"
    output:
        SAMPLE_OUTPUT_DIR+"/sorted/{sample}.bam"
    log:
        SAMPLE_OUTPUT_DIR+"/logs/samtools/sort_sam/{sample}.log"
    params:
        ""
        # "-m 4G"
    threads: Thread_NUM     # This value - 1 will be sent to -@.
    conda:
        "envs/mapping.yaml"
    shell:
       "samtools sort {params} -@ {threads} -O bam {input} > {output}"

rule samtools_index:
    input:
        SAMPLE_OUTPUT_DIR+"/sorted/{sample}.bam"
    output:
        SAMPLE_OUTPUT_DIR+"/sorted/{sample}.bam.bai"
    threads: Thread_NUM
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools index -@ {threads} {input}"

rule mark_duplicates:
    input:
        sorted_bam = SAMPLE_OUTPUT_DIR+"/sorted/{sample}.bam",
        sorted_bam_index = SAMPLE_OUTPUT_DIR+"/sorted/{sample}.bam.bai"
    output:
        bam = SAMPLE_OUTPUT_DIR+"/dedup/{sample}.bam",
        metrics = SAMPLE_OUTPUT_DIR+"/dedup/{sample}.metrics.txt"
    log:
        SAMPLE_OUTPUT_DIR+"/logs/picard/dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=STRICT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true"
    conda:
        "envs/picard.yaml"
    shell:
        "picard MarkDuplicates {params} I={input.sorted_bam} O={output.bam} M={output.metrics}"

# ---------------------------------------
# Performing base quality score recalibration
# ---------------------------------------
rule gatk_baserecalibrator:
    input:
        bam = SAMPLE_OUTPUT_DIR+"/dedup/{sample}.bam",
        ref = expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        dict =  expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".dict"),
        known_dbsnpvcf = config["know_dbsnp_vcf"],  # optional known sites
        known_indels = config["know_indels"],
        known_gold_standardindels = config["gold_standard_indels"],
        known_snpshigh_confidence = config["1000G_snpshigh_confidence"]

    output:
        recal_table= SAMPLE_OUTPUT_DIR+"/recal/{sample}.recal_data.table"
    log:
        SAMPLE_OUTPUT_DIR+"/logs/gatk/baserecalibrator/{sample}.log"
    params:
        extra=""  # optional
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk BaseRecalibrator {params.extra} -I {input.bam} -R {input.ref} --known-sites {input.known_dbsnpvcf} "
        "--known-sites {input.known_indels} "
        "--known-sites {input.known_gold_standardindels} "
        "--known-sites {input.known_snpshigh_confidence} "
        "-O {output.recal_table}"

rule gatk_applybqsr:
    input:
        bam= SAMPLE_OUTPUT_DIR+"/dedup/{sample}.bam",
        ref=expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        dict=expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".dict"),
        recal_table= SAMPLE_OUTPUT_DIR+"/recal/{sample}.recal_data.table"
    output:
        bam= SAMPLE_OUTPUT_DIR+"/recal/{sample}.bam"
    log:
        SAMPLE_OUTPUT_DIR+"/logs/gatk/gatk_applybqsr/{sample}.log"
    params:
        extra=""  # optional
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk ApplyBQSR {params.extra} -R {input.ref} -I {input.bam} --bqsr-recal-file {input.recal_table} -O {output.bam}"

rule samtools_index_bqsr:
    input:
        SAMPLE_OUTPUT_DIR + "/recal/{sample}.bam"
    output:
        SAMPLE_OUTPUT_DIR+"/recal/{sample}.bam.bai"
    threads: Thread_NUM
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools index -@ {threads} {input}"



# ---------------------------------------
# Performing variant discovery (SNV detection)
# ---------------------------------------
# ---------------
# Mutect2 SNV caller
# ---------------
rule Mutect2:
    input:
        ref = expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        recal_bam = SAMPLE_OUTPUT_DIR+ "/recal/{sample}.bam",
        germline =  config["af_only_gnomad"],
        PON = config["gatk_panel_of_normal"],
        intervals = config["interval_list"]
    params:
        java_params = "-Xmx12g",
    output:      
        unfilterd_vcf = SAMPLE_OUTPUT_DIR+"/calls/{sample}.vcf.gz",
        filter_file = SAMPLE_OUTPUT_DIR+"/calls/{sample}.f1r2.tar.gz"
    params:
              ###############
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk Mutect2 --java-options {params.java_params} -L {input.intervals} -R {input.ref} -I {input.recal_bam} "
        "--germline-resource {input.germline} --panel-of-normals {input.PON} "
        "--f1r2-tar-gz {output.filter_file} -O {output.unfilterd_vcf}"

rule LearnReadOrientationModel:
    input:
        filter_file = SAMPLE_OUTPUT_DIR+"/calls/{sample}.f1r2.tar.gz"
    output:
        ROM_artifact = SAMPLE_OUTPUT_DIR+"/calls/{sample}.read-orientation-model.tar.gz"
    params:
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk LearnReadOrientationModel -I {input.filter_file} -O {output.ROM_artifact}"

rule GetPileupSummaries:
    input:
        tumor_bam = SAMPLE_OUTPUT_DIR+"/recal/{sample}.bam",
        germline = config["exac_common_knownsite"],
        intervals = config["exac_common_knownsite"]
    output:
        pileupsummaries = SAMPLE_OUTPUT_DIR+"/calls/{sample}.pileupsummaries.table"
    params:
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk GetPileupSummaries -I {input.tumor_bam} -V {input.germline} -L {input.intervals} -O {output}"

rule CalculateContamination:
    input:
        pileupsummaries = SAMPLE_OUTPUT_DIR+"/calls/{sample}.pileupsummaries.table"
    output:
        segments_table = SAMPLE_OUTPUT_DIR+"/calls/{sample}.segments.table",
        contamination_table = SAMPLE_OUTPUT_DIR+"/calls/{sample}.contamination.table"
    params:
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk CalculateContamination -I {input.pileupsummaries} -tumor-segmentation {output.segments_table} -O {output.contamination_table}"


rule FilterMutectCalls:
    input:
        ref = expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        unfilterd_vcf = SAMPLE_OUTPUT_DIR+"/calls/{sample}.vcf.gz",
        segments_table = SAMPLE_OUTPUT_DIR+"/calls/{sample}.segments.table",
        contamination_table = SAMPLE_OUTPUT_DIR+"/calls/{sample}.contamination.table",
        ROM_artifact = SAMPLE_OUTPUT_DIR+"/calls/{sample}.read-orientation-model.tar.gz"
    output:
        filtered_vcf = SAMPLE_OUTPUT_DIR+"/filterd_calls/{sample}.snv.vcf.gz"
    params:
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk FilterMutectCalls -V {input.unfilterd_vcf} -R {input.ref} "
        "--tumor-segmentation {input.segments_table} --contamination-table {input.contamination_table} "
        "--ob-priors {input.ROM_artifact} "      #### 
        "-O {output.filtered_vcf}"
        # " gatk FilterMutectCalls -V somatic.vcf.gz --contamination-table contamination.table -O filtered.vcf.gz"

# ---------------
# Germline HaplotypeCaller
# ---------------
rule HaplotypeCaller:
    input:
        ref = expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        recal_bam = SAMPLE_OUTPUT_DIR+ "/recal/{sample}.bam",
        germline =  config["af_only_gnomad"],
        PON = config["gatk_panel_of_normal"],
        #PON = "PON_VCF/gatk4_mutect2_pon_generated.vcf.gz",
        intervals = config["interval_list"]
    output:      
        g_vcf = SAMPLE_OUTPUT_DIR+"/snp_calls/{sample}.g.vcf.gz",
    params:
        java_params = "-Xmx12g"
              ###############
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk --java-options {params.java_params} HaplotypeCaller -L {input.intervals} -R {input.ref} -I {input.recal_bam} "
        "-ERC GVCF "
        "-O {output.g_vcf}"

rule GenotypeGVCFs:
    input:
        ref = expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        recal_bam = SAMPLE_OUTPUT_DIR+ "/recal/{sample}.bam",
        g_vcf = SAMPLE_OUTPUT_DIR+"/snp_calls/{sample}.g.vcf.gz",
        known_dbsnpvcf = config["know_dbsnp_vcf"],
        intervals = config["interval_list"]
    output:      
        snp_vcf = SAMPLE_OUTPUT_DIR+"/snp_calls/{sample}.vcf.gz",
    params:
        java_params = "-Xmx12g"
              ###############
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk --java-options {params.java_params} GenotypeGVCFs -L {input.intervals} -R {input.ref} -V {input.g_vcf} "
        "--dbsnp {input.known_dbsnpvcf} "
        "-O {output.snp_vcf}"




# ---------------------------------------
# Performing variant discovery (SV detection)
# ---------------------------------------
# Extract the discordant paired-end alignments.
# Extract the split-read alignments
rule Ext_discordant_align:
    input:
        SAMPLE_OUTPUT_DIR + "/recal/{sample}.bam"
    output:
        temp(SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.discordants.unsorted.bam")
    conda:
        "envs/svDetect.yaml"
    shell:
        "samtools view -b -F 1294 {input}> {output}"

rule Ext_splitread_align:
    input:
        SAMPLE_OUTPUT_DIR + "/recal/{sample}.bam"
    output:
        temp(SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.splitters.unsorted.bam")
    conda:
        "envs/svDetect.yaml"
    shell:
        "samtools view -h {input} | "
        "extractSplitReads_BwaMem -i stdin | "
        "samtools view -Sb - > {output}" 
# samtools view -h DNA-HD753-50ng_S8_L001.bam | ~/harddisk/BIO_TOOl/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > DNA-HD753-50ng_S8_L001.splitters.unsorted.bam

rule sort_align:
    input:
        unsorted_discordants_bam = SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.discordants.unsorted.bam",
        unsorted_splitread_bam = SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.splitters.unsorted.bam"
    output:
        discordants_bam = SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.discordants.bam",
        splitread_bam = SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.splitters.bam"
    conda:
        "envs/svDetect.yaml"
    shell:
        """
        samtools sort {input.unsorted_discordants_bam} -o {output.discordants_bam}
        samtools index {output.discordants_bam}
        samtools sort {input.unsorted_splitread_bam} -o {output.splitread_bam}
        samtools index {output.splitread_bam}
        """

rule lumpy_SV:
    input:
        bam = SAMPLE_OUTPUT_DIR + "/recal/{sample}.bam",
        discordants_bam = SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.discordants.bam",
        splitread_bam = SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.splitters.bam"
    output:
        sv_vcf = SAMPLE_OUTPUT_DIR + "/lumpySV/{sample}.sv.vcf"
    conda:
        "envs/svDetect.yaml"
    shell:
        "lumpyexpress -B {input.bam} "
        "-S {input.splitread_bam} "
        "-D {input.discordants_bam} "
        "-o {output.sv_vcf}"



# #######################################
# ---------------------------------------
# CNVdetect CNVkit
# ---------------------------------------
rule CNVkit:
    input:
        ref = expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        input_bam = SAMPLE_OUTPUT_DIR + "/recal/{sample}.bam",
        target_bed = config["hg38gene_bed"],
        access_bed = config["hg38access_bed"]
      
    params:
        outdir = SAMPLE_OUTPUT_DIR + "/CNV/",
        inputCNS = SAMPLE_OUTPUT_DIR + "/CNV/{sample}.call.cns",
        sampleID = "{sample}",
        processNum = 8
    output:
        cnvkitRef_output = SAMPLE_OUTPUT_DIR + "/CNV/{sample}_flat_reference.cnn",
        cnv_vcf = SAMPLE_OUTPUT_DIR + "/CNV/{sample}.cnv.vcf"
    conda:
        "envs/cnvkit.yaml"
    shell:
        "cnvkit.py batch {input.input_bam} "
        "-n -t {input.target_bed} "
        "-f {input.ref} --access {input.access_bed} "
        "--output-reference {output.cnvkitRef_output} -d {params.outdir} && "
        "cnvkit.py export vcf {params.inputCNS} -y  -i {params.sampleID} -o {output.cnv_vcf}"



# #######################################
# ---------------------------------------
# SNV SV CNV concat
# ---------------------------------------
rule sort_SVvcf:
    input:
        sv_vcf = SAMPLE_OUTPUT_DIR + "/lumpySV/{sample}.sv.vcf",
    output:
        sv_sort_vcf = SAMPLE_OUTPUT_DIR + "/SNVSVCNVConcat/{sample}.sv.sorted.vcf",
        sv_sort_vcfgz = SAMPLE_OUTPUT_DIR + "/SNVSVCNVConcat/{sample}.sv.sorted.vcf.gz"
    conda:
        "envs/sortSVCNV.yaml"
    shell:
        """
        cat {input.sv_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' >  {output.sv_sort_vcf}
        bgzip -c {output.sv_sort_vcf} > {output.sv_sort_vcfgz}
        tabix -p vcf {output.sv_sort_vcfgz}
        """

rule sort_CNVvcf:
    input:
        cnv_vcf = SAMPLE_OUTPUT_DIR + "/CNV/{sample}.cnv.vcf",
        #ref = expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        # splitread_bam = SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.splitters.bam"
    output:
        cnv_sort_vcf = SAMPLE_OUTPUT_DIR + "/SNVSVCNVConcat/{sample}.cnv.sorted.vcf",
        cnv_sort_vcfgz = SAMPLE_OUTPUT_DIR + "/SNVSVCNVConcat/{sample}.cnv.sorted.vcf.gz"
    conda:
        "envs/sortSVCNV.yaml"
    shell:
        """
        cat {input.cnv_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k2,2n"}}' >  {output.cnv_sort_vcf}
        bgzip -c {output.cnv_sort_vcf} > {output.cnv_sort_vcfgz}
        tabix -p vcf {output.cnv_sort_vcfgz}
        """

rule concat_SNVCNVSV:
    input:
        filtered_vcf = SAMPLE_OUTPUT_DIR+"/filterd_calls/{sample}.snv.vcf.gz",
        cnv_sort_vcfgz = SAMPLE_OUTPUT_DIR + "/SNVSVCNVConcat/{sample}.cnv.sorted.vcf.gz",
        sv_sort_vcfgz = SAMPLE_OUTPUT_DIR + "/SNVSVCNVConcat/{sample}.sv.sorted.vcf.gz",      
        #ref = expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
        # splitread_bam = SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.splitters.bam"
    output:
        Concat_vcf = SAMPLE_OUTPUT_DIR + "/SNVSVCNVConcat/{sample}.concat.vcf",
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools concat -a {input.filtered_vcf} {input.cnv_sort_vcfgz} {input.sv_sort_vcfgz} -o {output.Concat_vcf}
        """




# #######################################
# ---------------------------------------
# Evaluate final results
# ---------------------------------------
# rule unzip:
#     input:
#         # snp_vcf = SAMPLE_OUTPUT_DIR+"/snp_calls/{sample}.vcf.gz",
#         #filtered_vcf = SAMPLE_OUTPUT_DIR+"/filterd_calls/{sample}.snv.vcf.gz",
#         indel_vcf = SAMPLE_OUTPUT_DIR+"/indel_vcf/{sample}.vcf.gz",
#         # vcftool_indel_vcf = SAMPLE_OUTPUT_DIR+"/indel_vcf/{sample}.vcftool.vcf.gz",
#         # filtered_vcf = SAMPLE_OUTPUT_DIR+"/calls/{sample}.snv.vcf.gz",
#     output:
#         #  snv_vcf = SAMPLE_OUTPUT_DIR+"/snp_calls/{sample}.vcf"
#         #  snv_vcf = SAMPLE_OUTPUT_DIR+"/filterd_calls/{sample}.snv.vcf"
#         indel_vcfunzip = SAMPLE_OUTPUT_DIR+"/indel_vcf/{sample}.vcf",
#         # vcftool_indel_vcfunzip = SAMPLE_OUTPUT_DIR+"/indel_vcf/{sample}.vcftool.vcf",
#         # snv_vcf = SAMPLE_OUTPUT_DIR+"/calls/{sample}.snv.vcf"
#     shell:
#         """
#         gunzip {input.indel_vcf} 
#         """
#         # gunzip {input.vcftool_indel_vcf} 

# rule vcf_compare:
#     input:
#         # snv_vcf = "/sybig/home/jiy/Downloads/Snakemake_pipeline/EvaluationData/NA12878/TruthBam/{sample}.vcf",
#         sv_vcf = SAMPLE_OUTPUT_DIR + "/lumpySV/{sample}.sv.vcf",
#         # snv_vcf = SAMPLE_OUTPUT_DIR+"/snp_calls/{sample}.vcf",
#         # snv_vcf = SAMPLE_OUTPUT_DIR+"/filterd_calls/{sample}.snv.vcf",
#         #snv_vcf = SAMPLE_OUTPUT_DIR+"/calls/{sample}.snv.vcf",
#         gold_vcf = config["Evaluate"]["gold_vcf"],
#         # "/sybig/home/jiy/Downloads/Snakemake_pipeline/EvaluationData/NA12878/NA12878_S1.vcf",
#         highConfidentRegion = config["Evaluate"]["highConfidentRegion"],
#         # highConfidentRegion = "/sybig/home/jiy/Downloads/Snakemake_pipeline/EvaluationData/ConfidentRegions/ConfidentRegions.bed"   
#         ref = expand(GENOMEDIR+"/{genome}",genome = Genome_NAME+".fa"),
#         # ref = config["Evaluate"]["eval_ref"],
#         # splitread_bam = SAMPLE_OUTPUT_DIR+"/lumpySV/{sample}.splitters.bam"
#         vcf_compare = config["Evaluate"]["gold_vcf"],
#     output:
#         compare_vcf = SAMPLE_OUTPUT_DIR + "/compareVCF/{sample}",
#     params:
#         out_para = SAMPLE_OUTPUT_DIR + "/compareVCF/{sample}",
#     conda:
#         "envs/VCF_compare.yaml"
#     shell:
#         """
#         python2 /sybig/home/jiy/Downloads/Snakemake_pipeline/EvaluationData/genReads1-master/extra_utilities/vcf_compare.py \
#         -r {input.ref} \
#         -g {input.gold_vcf} \
#         -w {input.sv_vcf} \
#         --vcf-out \
#         -o {params.out_para} \
#         -t {input.highConfidentRegion} -T 0 --incl-homs --incl-fail --no-plot
#         """ ####### care about the vcf_compare path.
#         ##############   Liftover for PrecisionFDA
# rule Liftover:
#     input:
#         hg19_ref = config["hg19_ref"],
#         # snv_vcf = SAMPLE_OUTPUT_DIR+"/filterd_calls/{sample}.snv.vcf",
#         indel_vcfunzip = SAMPLE_OUTPUT_DIR+"/indel_vcf/{sample}.vcf",
#         # vcftool_indel_vcfunzip = SAMPLE_OUTPUT_DIR+"/indel_vcf/{sample}.vcftool.vcf",
#         liftover_chain = config["liftover_chain"],
#     output:
#         indel_hg19_vcf = SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.GATKhg19indel.vcf",
#         reject_vcf = SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.GATKhg19indelReject.vcf",
#         # indel_hg19vcftool_vcf = SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.vcftoolhg19indel.vcf",
#         # reject_vcf2 = SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.vcftoolhg19indelReject.vcf",
#     params:
#         picard_mem = "-Xmx80g"
#     conda:
#         "envs/picard.yaml"
#     shell:
#         """
#         picard {params.picard_mem} LiftoverVcf I={input.indel_vcfunzip} \
#             O={output.indel_hg19_vcf} \
#             CHAIN={input.liftover_chain} \
#                 REJECT={output.reject_vcf} \
#                 R={input.hg19_ref}
#         """

#         # picard {params.picard_mem} LiftoverVcf I={input.vcftool_indel_vcfunzip} \
#         # O={output.indel_hg19vcftool_vcf} \
#         # CHAIN={input.liftover_chain} \
#         #     REJECT={output.reject_vcf2} \
#         #     R={input.hg19_ref}
# rule gzip:
#     input:
#         # snp_vcf = SAMPLE_OUTPUT_DIR+"/snp_calls/{sample}.vcf.gz",
#         #snv_hg19_vcf = SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.hg19snv.vcf",
#         indel_hg19_vcf =  SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.GATKhg19indel.vcf",
#         # indel_hg19vcftool_vcf = SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.vcftoolhg19indel.vcf",
#         # filtered_vcf = SAMPLE_OUTPUT_DIR+"/calls/{sample}.snv.vcf.gz",
#     output:
#         # snv_vcf = SAMPLE_OUTPUT_DIR+"/snp_calls/{sample}.vcf"
#         #snv_hg19_vcfgz = SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.hg19snv.vcf.gz",
#         indel_hg19_vcfgz =  SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.GATKhg19indel.vcf.gz",
#         # indel_hg19vcftool_vcfgz = SAMPLE_OUTPUT_DIR+"/hg19_calls/{sample}.vcftoolhg19indel.vcf.gz",
#         # snv_vcf = SAMPLE_OUTPUT_DIR+"/calls/{sample}.snv.vcf"
#     shell:
#         """
#         gzip {input.indel_hg19_vcf} 
#         """
#          # gzip {input.indel_hg19vcftool_vcf} 
