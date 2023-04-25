SAMPLES, = glob_wildcards("data/samples/{sample}_R1_001.fastq.gz")

adaptors = "data/adapters.fa"

rule all:
    input:
        "results/multiqc_report.html",
        "results/multiqc_report_trimmed.html",
        directory('chr19_STAR'),
        expand("snake/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("snake/{sample}/counts.txt", sample=SAMPLES),
        r1_trimmed = expand("trimmed/{sample}_trimmed_L001_R1_001.fastq.gz", sample=SAMPLES),
        r2_trimmed = expand("trimmed/{sample}_trimmed_L001_R2_001.fastq.gz", sample=SAMPLES),
        r1_trimmed_fastqc_html = expand("results/fastqc/trimmed/{sample}_trimmed_L001_R1_001_fastqc.html", sample=SAMPLES),
        r2_trimmed_fastqc_html = expand("results/fastqc/trimmed/{sample}_trimmed_L001_R2_001_fastqc.html", sample=SAMPLES),
        r1_trimmed_fastqc_zip = expand("results/fastqc/trimmed/{sample}_trimmed_L001_R1_001_fastqc.zip", sample=SAMPLES),
        r2_trimmed_fastqc_zip = expand("results/fastqc/trimmed/{sample}_trimmed_L001_R2_001_fastqc.zip", sample=SAMPLES)

rule fastqc:
    input:
        R1 = "data/samples/{sample}_R1_001.fastq.gz",
        R2 = "data/samples/{sample}_R2_001.fastq.gz"
    output:
        html = ["results/fastqc/{sample}_R1_001_fastqc.html", "results/fastqc/{sample}_R2_001_fastqc.html"],
        zip = ["results/fastqc/{sample}_R1_001_fastqc.zip", "results/fastqc/{sample}_R2_001_fastqc.zip"]
    params: "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: 2
    shell:
        "fastqc {input.R1} {input.R2} -o results/fastqc -t {threads} &> {log}"

rule fastqc_trimmed:
    input:
        R1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        R2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz"
    output:
        html = ["results/fastqc/trimmed/{sample}_trimmed_L001_R1_001_fastqc.html", "results/fastqc/trimmed/{sample}_trimmed_L001_R2_001_fastqc.html"],
        zip = ["results/fastqc/trimmed/{sample}_trimmed_L001_R1_001_fastqc.zip", "results/fastqc/trimmed/{sample}_trimmed_L001_R2_001_fastqc.zip"]
    params: "--quiet"
    log:
        "logs/fastqc/trimmed/{sample}.log"
    threads: 2
    shell:
        "fastqc {input.R1} {input.R2} -o results/fastqc/trimmed -t {threads} &> {log}"

rule multiqc:
    input:
        expand(["results/fastqc/{sample}_R1_001_fastqc.zip", "results/fastqc/{sample}_R2_001_fastqc.zip"], sample = SAMPLES)
    output:
        "results/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    envmodules:
        "MultiQC/1.9-gimkl-2020a-Python-3.8.2"
    shell:
        "multiqc {input} -o results/ &> {log}"

rule multiqc_trimmed:
    input:
        expand(["results/fastqc/trimmed/{sample}_trimmed_L001_R1_001_fastqc.zip", "results/fastqc/trimmed/{sample}_trimmed_L001_R2_001_fastqc.zip"], sample = SAMPLES)
    output:
        "results/multiqc_report_trimmed.html"
    log:
        "logs/multiqc/trimmed/multiqc.log"
    envmodules:
        "MultiQC/1.9-gimkl-2020a-Python-3.8.2"
    shell:
        "multiqc {input} -o results/ -n multiqc_report_trimmed.html &> {log}"

rule bbduk:
    input:
        r1 = 'data/samples/{sample}_R1_001.fastq.gz',
        r2 = 'data/samples/{sample}_R2_001.fastq.gz'
    output:
        out1 = "trimmed/{sample}_trimmed_L001_R1_001.fastq.gz",
        out2 = "trimmed/{sample}_trimmed_L001_R2_001.fastq.gz",
    log: "logs/bbduk.{sample}.log"
    conda: "environment.yaml"
    shell: "bbduk.sh in1={input.r1} out1={output.out1} in2={input.r2} out2={output.out2} ref={adaptors} qtrim=r trimq=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo &>{log}; touch {output.out1} {output.out2}"

####
# 6. Star
####
rule index:
    input:
        fa = 'data/chr19_20Mb.fa',
        gtf = 'data/chr19_20Mb.gtf'
    output:
        directory('chr19_STAR')
    threads: 8
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fa} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 100'

rule snakemapping:
    input:
        R1L1 = "data/samples/{sample}_R1_001.fastq.gz",
        R2L1 = "data/samples/{sample}_R2_001.fastq.gz",
        refdir = "chr19_STAR"
    params:
        outdir = "snake/{sample}",
    output:
        "snake/{sample}/SJ.out.tab",
        "snake/{sample}/Aligned.sortedByCoord.out.bam"
    threads: 8
    shell:
        'rm -rf {params.outdir} && '
        'mkdir {params.outdir} && '
        'cd {params.outdir} && '
        'STAR --runThreadN {threads} '
        '--genomeDir ../../{input.refdir} '
        '--readFilesIn ../../{input.R1L1} ../../{input.R2L1} '
        '--readFilesCommand gunzip -c '
        '--outSAMtype BAM SortedByCoordinate && cd ..'

####
# 7. featureCounts
####
rule featurecounts:
    input:
        samples = "snake/{sample}",
        gtf = "data/chr19_20Mb.gtf",
    output:
        counts = "snake/{sample}/counts.txt",
    shell:
        'featureCounts -p -t exon -g gene_id -a {input.gtf} '
        '-o {input.samples}/counts.txt '
        '{input.samples}/Aligned.sortedByCoord.out.bam '
        '-s 1'# {1(Collibri) or 2(Kappa)}
        ## TODO think of how can 1/2 can be parameterized