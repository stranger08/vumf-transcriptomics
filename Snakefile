adaptors = "data/adapters.fa"

rule all:
    input:
        expand("results/fastqc/{dataset}/{sample}_L001_{read}_001_fastqc.html", sample=config['samples'], dataset=config['dataset'], read=config['reads']),
        "results/multiqc_report.html",
        expand("results/trimmed/{dataset}/{sample}_L001_{read}_001.fastq.gz", sample=config['samples'], dataset=config['dataset'], read=config['reads']),
        expand("results/fastqc_trimmed/{dataset}/{sample}_L001_{read}_001_fastqc.html", sample=config['samples'], dataset=config['dataset'], read=config['reads']),
        "results/multiqc_report_trimmed.html",
        directory('chr19_STAR'),
        expand("results/counts/{sample}/Aligned.sortedByCoord.out.bam", sample=config['samples']),
        expand("results/counts/{sample}/counts.txt", sample=config['samples']),

####
# 1. FastQC
####
rule fastqc:
    input:
        R1 = expand("data/{dataset}/{sample}_L001_R1_001.fastq.gz", sample=config['samples'], dataset=config['dataset']),
        R2 = expand("data/{dataset}/{sample}_L001_R2_001.fastq.gz", sample=config['samples'], dataset=config['dataset'])
    output:
        html = ["results/fastqc/{dataset}/{sample}_L001_R1_001_fastqc.html", "results/fastqc/{dataset}/{sample}_L001_R2_001_fastqc.html"],
        zip = ["results/fastqc/{dataset}/{sample}_L001_R1_001_fastqc.zip", "results/fastqc/{dataset}/{sample}_L001_R2_001_fastqc.zip"]
    params: "--quiet"
    threads: 2
    shell:
        "fastqc {input.R1} {input.R2} -o results/fastqc/{config[dataset]} -t {threads}"

###
# 2. MultiQC
###
rule multiqc:
    input:
        expand(
            [
                "results/fastqc/{dataset}/{sample}_L001_R1_001_fastqc.zip", 
                "results/fastqc/{dataset}/{sample}_L001_R2_001_fastqc.zip"
            ], 
            sample=config['samples'],
            dataset=config['dataset'])
    output:
        "results/multiqc_report.html"
    envmodules:
        "MultiQC/1.9-gimkl-2020a-Python-3.8.2"
    shell:
        "multiqc {input} -o results/"

####
# 3. bbduk
####
rule bbduk:
    input:
        r1 = ["data/{dataset}/{sample}_L001_R1_001.fastq.gz".format(sample=sample, dataset=config['dataset']) for sample in config['samples']],
        r2 = ["data/{dataset}/{sample}_L001_R2_001.fastq.gz".format(sample=sample, dataset=config['dataset']) for sample in config['samples']]
    params:
        r1=lambda wildcards, input: ','.join(input.r1),
        r2=lambda wildcards, input: ','.join(input.r2),
        out1=lambda wildcards, output: ','.join(output.out1),
        out2=lambda wildcards, output: ','.join(output.out2)
    output:
        out1 = expand("results/trimmed/{dataset}/{sample}_L001_R1_001.fastq.gz", sample=config['samples'], dataset=config['dataset']),
        out2 = expand("results/trimmed/{dataset}/{sample}_L001_R2_001.fastq.gz", sample=config['samples'], dataset=config['dataset'])
    conda: "environment.yaml"
    shell:
        """
            for sample in {config[samples]}
            do
                bbduk.sh \
                in1=data/{config[dataset]}/${{sample}}_L001_R1_001.fastq.gz \
                out1=results/trimmed/{config[dataset]}/${{sample}}_L001_R1_001.fastq.gz \
                in2=data/{config[dataset]}/${{sample}}_L001_R2_001.fastq.gz \
                out2=results/trimmed/{config[dataset]}/${{sample}}_L001_R2_001.fastq.gz \
                ref={adaptors} \
                qtrim=r \
                trimq=10 \
                ktrim=r \
                k=23 \
                mink=11 \
                hdist=1 \
                tpe \
                tbo
            done
        """

####
# 4. FastQC on trimmed data
####
rule fastqc_trimmed:
    input:
        R1 = expand("results/trimmed/{dataset}/{sample}_L001_R1_001.fastq.gz", sample=config['samples'], dataset=config['dataset']),
        R2 = expand("results/trimmed/{dataset}/{sample}_L001_R2_001.fastq.gz", sample=config['samples'], dataset=config['dataset'])
    output:
        html = ["results/fastqc_trimmed/{dataset}/{sample}_L001_R1_001_fastqc.html", "results/fastqc_trimmed/{dataset}/{sample}_L001_R2_001_fastqc.html"],
        zip = ["results/fastqc_trimmed/{dataset}/{sample}_L001_R1_001_fastqc.zip", "results/fastqc_trimmed/{dataset}/{sample}_L001_R2_001_fastqc.zip"]
    params: "--quiet"
    threads: 2
    shell:
        """
            mkdir -p results/fastqc_trimmed/{config[dataset]} &&
            fastqc {input.R1} {input.R2} -o results/fastqc_trimmed/{config[dataset]} -t {threads}
        """

###
# 5. MultiQC on trimmed data
###
rule multiqc_trimmed:
    input:
        expand(
            [
                "results/fastqc_trimmed/{dataset}/{sample}_L001_R1_001_fastqc.zip", 
                "results/fastqc_trimmed/{dataset}/{sample}_L001_R2_001_fastqc.zip"
            ], 
            sample=config['samples'],
            dataset=config['dataset'])
    output:
        "results/multiqc_report_trimmed.html"
    envmodules:
        "MultiQC/1.9-gimkl-2020a-Python-3.8.2"
    shell:
        "multiqc {input} -o results/ -n multiqc_report_trimmed.html"

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
        R1 = expand("data/{dataset}/{sample}_L001_R1_001.fastq.gz", sample=config['samples'], dataset=config['dataset']),
        R2 = expand("data/{dataset}/{sample}_L001_R2_001.fastq.gz", sample=config['samples'], dataset=config['dataset']),
        refdir = "chr19_STAR"
    params:
        outdir = "results/counts/{sample}",
        readFilesInR1 = lambda wildcards, input: ",".join(["../../../" + i for i in input.R1]),
        readFilesInR2 = lambda wildcards, input: ",".join(["../../../" + i for i in input.R2])
    output:
        "results/counts/{sample}/SJ.out.tab",
        "results/counts/{sample}/Aligned.sortedByCoord.out.bam",
    threads: 8
    shell:
        'mkdir -p {params.outdir} && '
        'cd {params.outdir} && '
        'STAR --runThreadN {threads} '
        '--genomeDir ../../../{input.refdir} '
        '--readFilesIn {params.readFilesInR1} {params.readFilesInR2} '
        '--readFilesCommand gunzip -c '
        '--outSAMtype BAM SortedByCoordinate && cd ../..'

####
# 7. featureCounts
####
rule featurecounts:
    input:
        gtf = "data/chr19_20Mb.gtf",
    output:
        counts = "results/counts/{sample}/counts.txt",
    shell:
        """
            for sample in {config[samples]}
            do
              featureCounts -p -t exon -g gene_id -s {config[strand]} -a {input.gtf} -o results/counts/${{sample}}/counts.txt results/counts/${{sample}}/Aligned.sortedByCoord.out.bam
            done
        """
