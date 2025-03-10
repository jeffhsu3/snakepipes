allelic_suffix = ['genome1', 'genome2']

## Salmon Index

rule SalmonIndex:
    input:
        "Annotation/genes.filtered.fa",
        genome_fasta
    output:
        "Salmon/SalmonIndex/decoys.txt",
        temp("Salmon/SalmonIndex/seq.fa"),
        "Salmon/SalmonIndex/seq.bin"
    benchmark:
        "Salmon/.benchmark/Salmon.index.benchmark"
    params:
        salmonIndexOptions = salmonIndexOptions
    log:
        out = "Salmon/SalmonIndex/SalmonIndex.out",
        err = "Salmon/SalmonIndex/SalmonIndex.err",
    threads: lambda wildcards: 16 if 16<max_thread else max_thread
    conda: CONDA_RNASEQ_ENV
    shell: """
        grep "^>" {input[1]} | cut -d " " -f 1 | tr -d ">" > {output[0]}
        cat {input[0]} {input[1]} > {output[1]}
        salmon index -p {threads} -t {output[1]} -d {output[0]} -i Salmon/SalmonIndex {params.salmonIndexOptions} > {log.out} 2> {log.err}
        """


def getSalmon_libtype(pairedEnd, libraryType):
    """
    Convert from a featureCounts library type to a HISAT2 option string
    """
    if pairedEnd:
        if libraryType == 1:
            return "ISF"
        elif libraryType == 2:
            return "ISR"
        else:
            return "IU"
    else:
        if libraryType == 1:
            return "SF"
        elif libraryType == 2:
            return "SR"
        else:
            return "U"

## Salmon quant, the bootstraps are needed in Sleuth
if pairedEnd:
    rule getAllelicFQ:
        input:
            allelic_bam="allelic_bams/{sample}.{allelic_suffix}.sorted.bam",
            allelic_bai="allelic_bams/{sample}.{allelic_suffix}.sorted.bam.bai"
        output:
            r1="allelicFASTQ/{sample}.{allelic_suffix}_R1.fastq.gz",
            r2="allelicFASTQ/{sample}.{allelic_suffix}_R2.fastq.gz"
        log: "allelicFASTQ/logs/bam2fq.{sample}.{allelic_suffix}.log"
        benchmark: "allelicFASTQ/.benchmark/bam2fq.{sample}.{allelic_suffix}.benchmark"
        threads: 4
        conda: CONDA_SHARED_ENV
        shell: """
               samtools collate -@ {threads} -u -O {input.allelic_bam} | \\
samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n  \\
                2> {log}
                """
    rule SalmonQuant:
        input:
            r1 = "allelicFASTQ/{sample}.{allelic_suffix}_R1.fastq.gz",
            r2 = "allelicFASTQ/{sample}.{allelic_suffix}_R2.fastq.gz",
            bin = "Salmon/SalmonIndex/seq.bin"
        output:
            quant = "SalmonAllelic/{sample}.{allelic_suffix}/quant.sf"
        log: "SalmonAllelic/logs/{sample}.{allelic_suffix}.quant.log"
        benchmark:
            "SalmonAllelic/.benchmark/SalmonQuant.{sample}.{allelic_suffix}.benchmark"
        params:
            outdir = "SalmonAllelic/{sample}.{allelic_suffix}",
            gtf = genes_gtf,
            lib_type = getSalmon_libtype(pairedEnd, libraryType)
        threads: 8
        conda: CONDA_RNASEQ_ENV
        shell: """
            salmon quant -p {threads} --softclipOverhangs --validateMappings --numBootstraps 50 -i Salmon/SalmonIndex -l {params.lib_type} -1 {input.r1} -2 {input.r2} -o {params.outdir} 2> {log}
            """
else:
    rule getAllelicFQ:
        input:
            allelic_bam="allelic_bams/{sample}.{allelic_suffix}.sorted.bam",
            allelic_bai="allelic_bams/{sample}.{allelic_suffix}.sorted.bam.bai"
        output:
            r1="allelicFASTQ/{sample}.{allelic_suffix}_R1.fastq.gz"
        log: "allelicFASTQ/logs/bam2fq.{sample}.{allelic_suffix}.log"
        benchmark: "allelicFASTQ/.benchmark/bam2fq.{sample}.{allelic_suffix}.benchmark"
        threads: 4
        conda: CONDA_SHARED_ENV
        shell: """
            samtools collate -@ {threads} -u -O {input.allelic_bam} | \\
samtools fastq -1 /dev/null -2 /dev/null -0 /dev/null -s {output.r1} -n  \\
                2> {log}
                """
    rule SalmonQuant:
        input:
            fastq = "allelicFASTQ/{sample}.{allelic_suffix}_R1.fastq.gz",
            bin = "SalmonAllelic/SalmonIndex/seq.bin",
        output:
            quant = "SalmonAllelic/{sample}.{allelic_suffix}/quant.sf"
        log: "SalmonAllelic/logs/{sample}.{allelic_suffix}.quant.log"
        benchmark:
            "SalmonAllelic/.benchmark/SalmonQuant.{sample}.{allelic_suffix}.benchmark"
        params:
            outdir = "SalmonAllelic/{sample}.{allelic_suffix}",
            gtf = genes_gtf,
            lib_type = getSalmon_libtype(pairedEnd, libraryType)
        threads: 8
        conda: CONDA_RNASEQ_ENV
        shell: """
            salmon quant -p {threads} --softclipOverhangs --validateMappings --numBootstraps 50 -i Salmon/SalmonIndex -l {params.lib_type} -r {input.fastq} -o {params.outdir} 2> {log}
            """


rule Salmon_symlinks:
    input:
        quant = "SalmonAllelic/{sample}.{allelic_suffix}/quant.sf"
    output:
        quant = "SalmonAllelic/{sample}.{allelic_suffix}.quant.sf"
    shell: """
         ln -s ../{input.quant} {output.quant}
            """


rule Salmon_TPM:
    input:
        expand("SalmonAllelic/{sample}.{allelic_suffix}.quant.sf", sample=samples,allelic_suffix=allelic_suffix)
    output:
        "SalmonAllelic/TPM.transcripts.tsv"
    benchmark:
        "SalmonAllelic/.benchmark/Salmon_TPM.benchmark"
    log:
        "SalmonAllelic/logs/Salmon_TPM.log"
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "merge_count_tables.R")+" Name TPM {output} {input} "


rule Salmon_counts:
    input:
        expand("SalmonAllelic/{sample}.{allelic_suffix}.quant.sf", sample=samples,allelic_suffix=allelic_suffix)
    output:
        "SalmonAllelic/counts.transcripts.tsv"
    benchmark:
        "SalmonAllelic/.benchmark/Salmon_counts.benchmark"
    log:
        "SalmonAllelic/logs/Salmon_counts.log"
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "merge_count_tables.R")+" Name NumReads {output} {input} "


## Prepare Salmon output for Sleuth
rule Salmon_wasabi:
    input:
        "SalmonAllelic/{sample}.{allelic_suffix}.quant.sf"
    output:
        "SalmonAllelic/{sample}.{allelic_suffix}/abundance.h5"
    log: "SalmonAllelic/logs/{sample}.{allelic_suffix}.wasabi.log"
    params:
        "SalmonAllelic/{sample}.{allelic_suffix}/"
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "wasabi.R")+" {params} 2> {log}"
