# This rule computes detailed BAM statistics using Samtools.
# It provides comprehensive alignment metrics, including insert size distribution, coverage depth, and quality scores.

rule samtools_stats:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{reads}}.sorted.bam"),
        bai = f"{path_bam}{name_genome}/mapping/{{reads}}.sorted.bam.bai"

    output:
        stats = f"{path_qc}/BAM/{name_genome}/samtools/{{reads}}.stats"

    params:
        samtools = samtools

    threads:
        config["QUALITY_CONTROL"]["THREADS"]

    log:
        stderr = f"{working_directory}/logs/samtools/{{reads}}_stats.err"

    shell:
        "{params.samtools} stats "
        "-@ {threads} {input.bam} > {output.stats} 2> {log.stderr}"


# This rule runs samtools flagstat to summarize alignment quality.
# It provides key statistics such as mapped read percentages, duplicate rate, and unmapped reads.

rule samtools_flagstat:
    input:
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{reads}}.sorted.bam"),
        bai = f"{path_bam}{name_genome}/mapping/{{reads}}.sorted.bam.bai"

    output:
        flagstat = f"{path_qc}/BAM/{name_genome}/samtools/{{reads}}.flagstat"

    params:
        samtools = samtools

    threads:
        config["QUALITY_CONTROL"]["THREADS"]

    log:
        stderr = f"{working_directory}/logs/samtools/{{reads}}_flagstat.err"

    shell:
        "{params.samtools} flagstat "
        "-@ {threads} {input.bam} > {output.flagstat} 2> {log.stderr}"


# This rule consolidates all BAM QC reports into a MultiQC summary.
# It combines alignment statistics, flagstat results, and logs from the mapping process.
# The final output is a single interactive HTML report.

rule multiqc_bam:
    input:
        samtools_stats = expand(f"{path_qc}/BAM/{name_genome}/samtools/{{reads}}.stats", reads=all_samples),
        samtools_flagstats = expand(f"{path_qc}/BAM/{name_genome}/samtools/{{reads}}.flagstat", reads=all_samples),
        log_final = expand(f"{path_bam}{name_genome}/mapping/log_star/{{reads}}_Log.final.out", reads=all_samples),
        log = expand(f"{path_bam}{name_genome}/mapping/log_star/{{reads}}_Log.out", reads=all_samples),
        log_progress = expand(f"{path_bam}{name_genome}/mapping/log_star/{{reads}}_Log.progress.out", reads=all_samples),
        tab = expand(f"{path_bam}{name_genome}/mapping/log_star/{{reads}}_SJ.out.tab", reads=all_samples)

    output:
        directory_data = directory(f"{path_qc}/multiqc/BAM/{name_genome}/{prefix}_{unique_id}_data/"),
        html = f"{path_qc}/multiqc/BAM/{name_genome}/{prefix}_{unique_id}.html"

    params:
        name = f"{prefix}_{unique_id}",
        multiqc = multiqc,
        path = f"{path_qc}/multiqc/BAM/{name_genome}/"

    threads:
        config["QUALITY_CONTROL"]["THREADS"]

    log:
        stdout = f"{working_directory}/logs/multiqc/BAM/{prefix}_{unique_id}.out",
        stderr = f"{working_directory}/logs/multiqc/BAM/{prefix}_{unique_id}.err"

    shell:
        "{params.multiqc} "
        "{input.samtools_stats} {input.samtools_flagstats} {input.log_final} {input.log} {input.log_progress} {input.tab} "
        "--filename {params.name} "
        "-o {params.path} 1> {log.stdout} 2> {log.stderr}"