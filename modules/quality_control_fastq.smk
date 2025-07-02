# This rule performs quality control on raw FASTQ files using FastQC.
# It generates an HTML report and a compressed ZIP file containing detailed QC metrics.
# Output and log directories are created if they do not already exist.

rule fastqc_raw:
    input:
        read = os.path.abspath(f"{working_directory}/1-raw_data/fastq/{{reads}}.fastq.gz")

    output:
        html = f"{path_qc}/fastqc_raw/{{reads}}_fastqc.html",
        zip = f"{path_qc}/fastqc_raw/{{reads}}_fastqc.zip"

    params:
        fastqc = fastqc,
        directory = f"{path_qc}/fastqc_raw/"

    threads:
        config["QUALITY_CONTROL"]["THREADS"]
    
    log:
        stdout = f"{working_directory}/logs/fastqc/fastq_raw/{{reads}}.out",
        stderr = f"{working_directory}/logs/fastqc/fastq_raw/{{reads}}.err"

    run:
        create_directory_if_not_exists(params.directory)
        shell("{params.fastqc} "
        "-o {params.directory} "
        "-t {threads} "
        "{input.read} > {log.stdout} 2> {log.stderr}")


# This rule aggregates the FastQC reports (HTML and ZIP files) from all raw FASTQ samples into a single MultiQC report.
# The report provides a summary of quality metrics for all samples.
# The output includes a consolidated HTML report and an optional directory containing MultiQC data.
# This step helps streamline quality control by centralizing information from multiple files.

rule multiqc_fastq_raw:
    input:
        html1 = expand(f"{path_qc}/fastqc_raw/{{reads}}.1_fastqc.html", reads=all_samples),
        zip1 = expand(f"{path_qc}/fastqc_raw/{{reads}}.1_fastqc.zip", reads=all_samples),
        html2 = expand(f"{path_qc}/fastqc_raw/{{reads}}.2_fastqc.html", reads=all_samples),
        zip2 = expand(f"{path_qc}/fastqc_raw/{{reads}}.2_fastqc.zip", reads=all_samples)

    output:
        directory_data = directory(f"{path_qc}/multiqc/fastq_raw/{prefix}_{unique_id}_data/"),
        html = f"{path_qc}/multiqc/fastq_raw/{prefix}_{unique_id}.html"
        
    params:
        name = f"{prefix}_{unique_id}",
        multiqc = multiqc,
        path = f"{path_qc}/multiqc/fastq_raw/"

    threads:
        config["QUALITY_CONTROL"]["THREADS"]
    
    log:
        stdout = f"{working_directory}/logs/multiqc/fastq_raw/{prefix}_{unique_id}.out",
        stderr = f"{working_directory}/logs/multiqc/fastq_raw/{prefix}_{unique_id}.err"
        
    shell:
        "{params.multiqc} "
        "{input.html1} {input.zip1} {input.html2} {input.zip2} "
        "--filename {params.name} "
        "-o {params.path} "
        "1> {log.stdout} 2> {log.stderr}"


# This rule performs quality control on trimmed FASTQ files using FastQC.
# It generates an HTML report and a ZIP archive containing detailed QC metrics.
# Output directories and log files are created to ensure proper organization and debugging.

rule fastqc_trimmed:
    input:
        trimmed_read = os.path.abspath(f"{path_fastq}/{{reads}}.fastq.gz")

    output:
        html = f"{path_qc}/fastqc_trimmed/{{reads}}_fastqc.html",
        zip = f"{path_qc}/fastqc_trimmed/{{reads}}_fastqc.zip"

    params:
        fastqc = fastqc,
        directory = f"{path_qc}/fastqc_trimmed/"

    threads:
        config["QUALITY_CONTROL"]["THREADS"]
    
    log:
        stdout = f"{working_directory}/logs/fastqc/fastq_trimmed/{{reads}}.out",
        stderr = f"{working_directory}/logs/fastqc/fastq_trimmed/{{reads}}.err"

    run:
        create_directory_if_not_exists(params.directory)
        shell("{params.fastqc} "
        "-o {params.directory} "
        "-t {threads} "
        "{input.trimmed_read} > {log.stdout} 2> {log.stderr}")


# This rule aggregates the FastQC reports (HTML and ZIP files) from all trimmed FASTQ samples into a single MultiQC report. 
# The report provides a summary of quality metrics for all samples.
# The output includes a consolidated HTML report and an optional directory containing MultiQC data.
# This step helps streamline quality control by centralizing information from multiple files.

rule multiqc_fastq_trimmed:
    input: 
        html1 = expand(f"{path_qc}/fastqc_trimmed/{{reads}}.1_fastqc.html", reads=all_samples),
        zip1 = expand(f"{path_qc}/fastqc_trimmed/{{reads}}.1_fastqc.zip", reads=all_samples),
        html2 = expand(f"{path_qc}/fastqc_trimmed/{{reads}}.2_fastqc.html", reads=all_samples),
        zip2 = expand(f"{path_qc}/fastqc_trimmed/{{reads}}.2_fastqc.zip", reads=all_samples),
        html = expand(f"{path_qc}/fastp_trimming/{{reads}}.html", reads=all_samples),
        json = expand(f"{path_qc}/fastp_trimming/{{reads}}.json", reads=all_samples)

    output:
        directory_data = directory(f"{path_qc}/multiqc/fastq_trimmed/{prefix}_{unique_id}_data/"),
        html = f"{path_qc}/multiqc/fastq_trimmed/{prefix}_{unique_id}.html"

    params:
        name = f"{prefix}_{unique_id}",
        multiqc = multiqc,
        path = f"{path_qc}/multiqc/fastq_trimmed/"

    threads:
        config["QUALITY_CONTROL"]["THREADS"]
    
    log:
        stdout = f"{working_directory}/logs/multiqc/fastq_trimmed/{prefix}_{unique_id}.out",
        stderr = f"{working_directory}/logs/multiqc/fastq_trimmed/{prefix}_{unique_id}.err"
    
    shell:
        "{params.multiqc} "
        "{input.html1} {input.zip1} {input.html2} {input.zip2} {input.html} {input.json} "
        "--filename {params.name} "
        "-o {params.path} "
        "1> {log.stdout} 2> {log.stderr}"