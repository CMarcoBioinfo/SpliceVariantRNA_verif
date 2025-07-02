# This rule generates reference files necessary for splicing and mapping workflows using the SpliceLauncherDB script.
# It processes a GFF3 file and a MANE annotation.
# If a MANE file is not provided, the default MANE annotation from SpliceLauncher will be used (as managed in the Snakefile).
# Outputs include BED, SJDB, and SpliceLauncher annotation files.

rule SpliceLauncherDB:
    input:
        gff3 =  gff3,
        mane = mane

    output:
        directory = directory(f"{working_directory}/2-processed_data/references/{name_genome}"),
        bed = f"{working_directory}/2-processed_data/references/{name_genome}/BEDannotation.bed",
        sjdb = f"{working_directory}/2-processed_data/references/{name_genome}/SJDBannotation.sjdb",
        annot = f"{working_directory}/2-processed_data/references/{name_genome}/SpliceLauncherAnnot.txt"

    params:
        Rscript = Rscript,
        SpliceLauncher = SpliceLauncher

    log:
        stdout = f"{working_directory}/logs/mapping/SpliceLauncherDB/{name_genome}.out",
        stderr = f"{working_directory}/logs/mapping/SpliceLauncherDB/{name_genome}.err"

    shell:
        "{params.Rscript} {params.SpliceLauncher}/scripts/generateSpliceLauncherDB.r "
        "-i {input.gff3} "
        "-o {output.directory} "
        "--mane {input.mane} 1> {log.stdout} 2> {log.stderr}"


# This function estimates the amount of memory (in megabytes) needed for genome indexing based on the genome length and the average number of bytes per base. 
# By default, it uses 10 bytes per base but can be adjusted for specific needs.
# The calculated memory is returned in megabytes, facilitating efficient resource allocation.

def estimate_memory(genome_length, bytes_per_base = 10):
    memory_required_bytes = genome_length * bytes_per_base

    memory_required_mb = memory_required_bytes / (1024 * 1024)
    return memory_required_mb

sjdbOverhang = get_config_value(config["MAPPING"]["INDEX"].get("SJDB_OVERHANG"), int(99))
ram_byte = get_config_value(config["MAPPING"]["INDEX"].get("RAM"), int(total_ram * 0.5))
genomeSAsparseD = get_config_value(config["MAPPING"]["INDEX"].get("GENOME_SA_SPARSE_D"), int(31000000000 / ram_byte) + 1)
genome_length = calculate_genome_length(genome)
genomeSAindexNbases = get_config_value(config["MAPPING"]["INDEX"].get("GENOME_SA_INDEX_NBASES"), min(14, math.log2(genome_length) / 2 - 1))

mem_mb = int(int(ram_byte) / (1024 * 1024))
memory_needed_mb = estimate_memory(genome_length, 10)

# This rule generates a genome index required for mapping reads using the STAR aligner.
# It creates various output files essential for efficient and accurate read alignment, including:
# SA, SAindex, chromosome information, and transcript-related files.
# Parameters like genome size, memory allocation, and STAR-specific settings are dynamically calculated using the provided genome length and RAM availability. 
# If certain parameters are not specified in the configuration, default values are estimated to optimize performance.
# The resources section ensures that appropriate memory is allocated during index generation.

rule star_index:
    input:
        reference = genome,
        genomeDir = f"{working_directory}/2-processed_data/references/{name_genome}",
        sjdb = f"{working_directory}/2-processed_data/references/{name_genome}/SJDBannotation.sjdb",
        gff3 = gff3
        
    output:
        reference = f"{working_directory}/2-processed_data/references/{name_genome}/{name_genome}",
        SA = f"{working_directory}/2-processed_data/references/{name_genome}/SA",
        SAindex = f"{working_directory}/2-processed_data/references/{name_genome}/SAindex",
        chrLength = f"{working_directory}/2-processed_data/references/{name_genome}/chrLength.txt",
        chrNameLength = f"{working_directory}/2-processed_data/references/{name_genome}/chrNameLength.txt",
        chrName = f"{working_directory}/2-processed_data/references/{name_genome}/chrName.txt",
        chrStart = f"{working_directory}/2-processed_data/references/{name_genome}/chrStart.txt",
        exonGeTrInfo = f"{working_directory}/2-processed_data/references/{name_genome}/exonGeTrInfo.tab",
        exonInfo = f"{working_directory}/2-processed_data/references/{name_genome}/exonInfo.tab",
        geneInfo = f"{working_directory}/2-processed_data/references/{name_genome}/geneInfo.tab",
        genome = f"{working_directory}/2-processed_data/references/{name_genome}/Genome",
        genomeParameters = f"{working_directory}/2-processed_data/references/{name_genome}/genomeParameters.txt",
        Log = f"{working_directory}/2-processed_data/references/{name_genome}/Log.out",
        sjdbInfo = f"{working_directory}/2-processed_data/references/{name_genome}/sjdbInfo.txt",
        sjdbListFromGTF = f"{working_directory}/2-processed_data/references/{name_genome}/sjdbList.fromGTF.out.tab",
        sjdblist = f"{working_directory}/2-processed_data/references/{name_genome}/sjdbList.out.tab",
        transcriptInfo = f"{working_directory}/2-processed_data/references/{name_genome}/transcriptInfo.tab"

    params:
        star = STAR,
        sjdbOverhang = sjdbOverhang,
        genomeSAsparseD = genomeSAsparseD,
        limitGenomeGenerateRAM = ram_byte,
        genomeSAindexNbases = genomeSAindexNbases

    threads:
        config["MAPPING"]["INDEX"]["THREADS"]
    
    log:
        stdout = f"{working_directory}/logs/mapping/star_index/{name_genome}.out",
        stderr = f"{working_directory}/logs/mapping/star_index/{name_genome}.err"
    
    resources:
        mem_mb = mem_mb

    shell:
        "ln -sfn {input.reference} {output.reference} && "
        "{params.star} "
        "--runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {input.genomeDir} "
        "--genomeFastaFiles {output.reference} "
        "--sjdbFileChrStartEnd {input.sjdb} "
        "--sjdbGTFfile {input.gff3} "
        "--sjdbOverhang {params.sjdbOverhang} "
        "--genomeSAsparseD {params.genomeSAsparseD} "
        "--limitGenomeGenerateRAM {params.limitGenomeGenerateRAM} "
        "--genomeSAindexNbases {params.genomeSAindexNbases} "
        "1> {log.stdout} 2> {log.stderr}"


# Calculates the optimal number of threads for STAR alignment based on available memory,
# memory needed per core, maximum cores, and requested threads, ensuring efficient resource usage.

def calculate_optimal_threads(mem_mb, memory_needed_mb, max_cores, requested_threads):
    total_mem_needed = memory_needed_mb * (max_cores // requested_threads)
    if total_mem_needed > mem_mb :
        max_parallel_rules = mem_mb // memory_needed_mb
        optimal_threads = max(1, int(max_cores / max_parallel_rules))
    else:
        optimal_threads = requested_threads
    
    optimal_threads = min(optimal_threads, max_cores)
    return optimal_threads


# Calculates the memory allocation per rule to balance parallelism and resource consumption.

def calculate_mem_per_rule(mem_mb, memory_needed_mb, optimal_threads):
    max_parallel_rules = mem_mb // memory_needed_mb
    mem_per_process = int(mem_mb / max_parallel_rules)
    return int(mem_per_process)

outFilterMismatchNmax = get_config_value(config["MAPPING"]["ALIGN"].get("OUT_FILTER_MISMATCH_NMAX"),int(2))
outFilterMultimapNmax = get_config_value(config["MAPPING"]["ALIGN"].get("OUT_FILTER_MULTIMAP_NMAX"),int(1))
outSJfilterIntronMaxVsReadN = get_config_value(config["MAPPING"]["ALIGN"].get("OUT_SJ_FILTER_INTRON_MAXVSREADN"),int(500000))

# This rule aligns paired-end reads to the reference genome using STAR. 
# Outputs include a sorted BAM file, indices, splice junction files, and log files. 
# Resource management dynamically adapts memory and threads for efficient performance.

rule star_align:
    input:
        read1 = os.path.abspath(f"{path_fastq}{{reads}}.1.fastq.gz"),
        read2 = os.path.abspath(f"{path_fastq}{{reads}}.2.fastq.gz"),
        genomeDir = f"{working_directory}/2-processed_data/references/{name_genome}",
        SA = f"{working_directory}/2-processed_data/references/{name_genome}/SA",
        SAindex = f"{working_directory}/2-processed_data/references/{name_genome}/SAindex",

    output:
        sort_bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{reads}}.sorted.bam"),
        bai = f"{path_bam}{name_genome}/mapping/{{reads}}.sorted.bam.bai",
        csi = f"{path_bam}{name_genome}/mapping/{{reads}}.sorted.bam.csi",
        log_final = f"{path_bam}{name_genome}/mapping/log_star/{{reads}}_Log.final.out",
        log = f"{path_bam}{name_genome}/mapping/log_star/{{reads}}_Log.out",
        log_progress = f"{path_bam}{name_genome}/mapping/log_star/{{reads}}_Log.progress.out",
        tab = f"{path_bam}{name_genome}/mapping/log_star/{{reads}}_SJ.out.tab",
        STARpass1 = directory(f"{path_bam}{name_genome}/mapping/log_star/{{reads}}__STARgenome"),
        STARgenome = directory(f"{path_bam}{name_genome}/mapping/log_star/{{reads}}__STARpass1")

    params:
        star = STAR,
        samtools = samtools,
        bam = os.path.abspath(f"{path_bam}{name_genome}/mapping/{{reads}}_Aligned.out.bam"),
        outFilterMismatchNmax = outFilterMismatchNmax,
        outFilterMultimapNmax = outFilterMultimapNmax,
        outSJfilterIntronMaxVsReadN = outFilterMultimapNmax,
        outFileNamePrefix = f"{path_bam}{name_genome}/mapping/{{reads}}_",
        outTmpDir = f"{path_bam}{name_genome}/mapping/{{reads}}_tmp",
        directory_log = f"{path_bam}{name_genome}/mapping/log_star/",
        log_final = f"{path_bam}{name_genome}/mapping/{{reads}}_Log.final.out",
        log = f"{path_bam}{name_genome}/mapping/{{reads}}_Log.out",
        log_progress = f"{path_bam}{name_genome}/mapping/{{reads}}_Log.progress.out",
        tab = f"{path_bam}{name_genome}/mapping/{{reads}}_SJ.out.tab",
        STARgenome = f"{path_bam}{name_genome}/mapping/{{reads}}__STARgenome",
        STARpass1 = f"{path_bam}{name_genome}/mapping/{{reads}}__STARpass1",

    threads:
        lambda wildcards: calculate_optimal_threads(mem_total_mb, memory_needed_mb, max_cores, int(config["MAPPING"]["ALIGN"]["THREADS"]))
    
    log:
        stdout_star = f"{working_directory}/logs/mapping/star_align/{{reads}}.star.out",
        stderr_star = f"{working_directory}/logs/mapping/star_align/{{reads}}.star.err",
        stdout_sort = f"{working_directory}/logs/mapping/star_align/{{reads}}.samtools_sort.out",
        stderr_sort = f"{working_directory}/logs/mapping/star_align/{{reads}}.samtools_sort.err",
        stdout_index = f"{working_directory}/logs/mapping/star_align/{{reads}}.samtools_index.out",
        stderr_index = f"{working_directory}/logs/mapping/star_align/{{reads}}.samtools_index.err"

    resources:
        mem_mb = lambda wildcards,threads: calculate_mem_per_rule(mem_total_mb, memory_needed_mb, threads)

    run:
        create_directory_if_not_exists(params["directory_log"])
        try:
            shell("{params.star} "
            "--runThreadN {threads} "
            "--outSAMstrandField intronMotif "
            "--outFilterMismatchNmax {params.outFilterMismatchNmax} "
            "--outFilterMultimapNmax {params.outFilterMultimapNmax} "
            "--genomeDir {input.genomeDir} "
            "--readFilesIn {input.read1} {input.read2} "
            "--readFilesCommand zcat "
            "--outSAMunmapped Within "
            "--outSAMtype BAM Unsorted "
            "--outSJfilterOverhangMin -1 8 8 8 "
            "--outSJfilterCountUniqueMin -1 1 1 1 "
            "--outSJfilterDistToOtherSJmin 0 0 0 0 "
            "--alignSJstitchMismatchNmax 0 -1 -1 -1 "
            "--outSJfilterIntronMaxVsReadN {params.outSJfilterIntronMaxVsReadN} "
            "--twopassMode Basic "
            "--outTmpDir {params.outTmpDir} "
            "--outSAMheaderHD \\@HD VN:1.4 SO:Unsorted "
            "--outFileNamePrefix {params.outFileNamePrefix} "
            "--genomeLoad NoSharedMemory "
            "1> {log.stdout_star} 2> {log.stderr_star}")
            shell("{params.samtools} sort -@ {threads} -o {output.sort_bam} --output-fmt BAM --write-index {params.bam} 1> {log.stdout_sort} 2> {log.stderr_sort}")
            shell("{params.samtools} index -@ {threads} -b {output.sort_bam} 1> {log.stdout_index} 2> {log.stderr_index}")
            shell("mv {params.log} {params.log_final} {params.log_progress} {params.tab} {params.STARpass1} {params.STARgenome} {params.directory_log}")
            shell("rm -rf {params.bam}")
        except Exception as e:
            shell("rm -rf {params.log} {params.log_final} {params.log_progress} {params.tab} {params.STARpass1} {params.STARgenome} {params.outTmpDir} {params.bam}")
            raise e