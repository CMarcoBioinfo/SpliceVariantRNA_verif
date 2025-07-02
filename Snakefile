#Import dependences python
import os
import sys
import shutil
from pathlib import Path
import time
import psutil
import gc
import math
import pandas as pd

# Memory management
total_ram = psutil.virtual_memory().total
mem_total_mb = int(total_ram / (1024 * 1024) * 0.8)

# Absolute path to the folder containing the Snakefile
ROOTDIR = workflow.basedir
MOD = os.path.join(ROOTDIR, "modules")
SCR = os.path.join(ROOTDIR, "scripts")
SPL = os.path.join(ROOTDIR, "SpliceLauncher")


def calculate_genome_length(fasta_file):
    genome_length = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome_length += len(line.strip())
    return genome_length

# Check if a configuration value is defined; otherwise, apply a default value.
# Convert the value to the specified type if applicable.
def get_config_value(key, default, value_type=int):
    if key is None  or (isinstance(key, str) and key.strip() == ""):
        return default  # If the value is missing, return the default.

    if value_type == str:
        return key  # If the expected type is a string, return the value as is.

    try:
        return value_type(key)  # Convert to the specified type (int, float, etc.).
    except ValueError:
        return default  # If conversion fails, return the default.

# Checking the dag and rulegraph options
args = sys.argv
dag_mode, rulegraph_mode, filegraph_mode = "--dag" in args, "--rulegraph" in args, "--filegraph" in args

if "SNAKEMAKE_PRINT" not in os.environ:
    os.environ["SNAKEMAKE_PRINT"] = "false"

def print_once(message):
    if not dag_mode and not rulegraph_mode and not filegraph_mode and os.getenv("SNAKEMAKE_PRINT") == "false":
        print(message)

# Number of cores management
if "MAX_CORES" not in os.environ:
    for arg in ["--cores", "-c", "--jobs", "-j"]:
        if arg in args:
            os.environ["MAX_CORES"] = args[args.index(arg) + 1]
            break
    else:
        os.environ["MAX_CORES"] = "1"

max_cores = int(os.environ["MAX_CORES"])
print_once(f"Valeur de max_cores : {max_cores}")

# Unique identifier based on date and time
if os.getenv("SNAKEMAKE_PRINT") == "false":
    current_time = time.localtime()
    unique_id = time.strftime("%Y%m%d%H%M%S", current_time)
    os.environ["UNIQUE_ID"] = unique_id

unique_id = os.environ["UNIQUE_ID"]

# Definition of main paths
print_once("Vérification des paramètres")
working_directory = get_config_value(config["GENERAL"].get("WORKING_DIRECTORY"), os.path.join(ROOTDIR,"data"), str)
prefix = get_config_value(config["GENERAL"].get("PREFIX"), f"run_{unique_id}", str)
path_qc = os.path.join(working_directory, "3-Quality_control")
path_results = os.path.join(working_directory, "4-Results", prefix)

print_once(f"Working directory : {working_directory} ...... OK")
print_once(f"Prefix: {prefix} ...... OK")

# Check if the specified directory exists; create it if it doesn't
def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

#Check if variable is None or empty
def check_value(var):
    if isinstance(var, int):
        return True
    if isinstance(var, str):
        return bool(var.strip())
    return var is not None

# Extracts file ID by removing the extension
# If compressed (.gz), removes two extensions, otherwise one
def get_id(file):
    return file.rsplit(".", 2)[0] if file.endswith(".gz") else file.rsplit(".", 1)[0]

# Check if file is compressed (.gz)
def is_gz(file):
    return file.endswith(".gz")

# Check if an executable exists
def check_executable(executable):
    if shutil.which(executable):
        print_once(f"{executable} ...... OK")
    else:
        print_once(f"Erreur : L'exécutable {executable} est manquant.")
        sys.exit(1)

# Get input file path based on mode (0, 1, 2)
def get_inputs(wildcards, mode):
    id = wildcards.id.rsplit("/", 1)[-1]

    if mode == 0:
        if "reads1" in wildcards.id:
            file = str(fastq2[id][0])
        elif "reads2" in wildcards.id:
            file = str(fastq2[id][1])
        else:
            file = str(fastq[id])

    elif mode == 1:
        file = str(fastq2[id][0])

    elif mode == 2:
        file = str(fastq2[id][1])

    # Handle .gz extension using is_gz()
    if is_gz(file) or os.path.exists(file):
        return file
    else:
        return f"{file}.gz"

# Define temp directories
tmp_directory = f"{working_directory}/.tmp/"
tmp_samples = f"{tmp_directory}samples/"
tmp_reads1 = f"{tmp_samples}reads1/"
tmp_reads2 = f"{tmp_samples}reads2/"

# Create only necessary directories
create_directory_if_not_exists(tmp_directory)
create_directory_if_not_exists(tmp_samples)
create_directory_if_not_exists(tmp_reads1)
create_directory_if_not_exists(tmp_reads2)

# Process samples from CSV
def process_sample(dict_csv, length_fastp=None):
    all_samples = set()

    for _, row in dict_csv.iterrows():
        id_sample, read1, read2 = row['id'], row['path_read1'], row['path_read2']

        # Check if both reads exist (compressed or uncompressed)
        if any(os.path.exists(read) or os.path.exists(f"{read}.gz") for read in [read1, read2]):
            if any(ext in read1 for ext in [".fastq", ".fq"]) and any(ext in read2 for ext in [".fastq", ".fq"]):
                name_key = f"{id_sample}.{length_fastp}bp" if length_fastp else id_sample

                read1_path = os.path.join(tmp_reads1, f"{name_key}.pre")
                read2_path = os.path.join(tmp_reads2, f"{name_key}.pre")
                fastq2[name_key] = [read1, read2]

                # Create placeholder files if they don't exist
                for path in [read1_path, read2_path]:
                    if not os.path.exists(path):
                        os.system(f"touch {path}")

                all_samples.add(name_key)

    return all_samples, fastq2

# Check bash executable
bash = get_config_value(config["DEPENDANCES"]["GENERAL"].get("BASH"), "bash", str)

# Check if trimming is enabled
use_trimming = get_config_value(config["USAGE"].get("TRIMMING"), int(0))

# Set trimming length or default value
length_fastp = int(get_config_value(config["TRIMMING"].get("LENGTH"), 100)) if use_trimming else None

# Load sample file and process if valid
csv = get_config_value(config["GENERAL"].get("SAMPLES_FILE"), "", str)
if not os.path.isfile(csv):
    print_once(f"Erreur : Le fichier {csv} n'existe pas.")
else:
    print_once(f"Fichier {csv} ...... OK")
    
all_samples, fastq2 = {}, {}

if check_value(csv):
    dict_csv = pd.read_csv(csv, sep="\t", comment="#", skip_blank_lines=True)
    all_samples, fastq2 = process_sample(dict_csv, length_fastp)

# Check config.yaml parameters
list_inputs = []

# Verify executables and general parameters
if not all_samples:
    print_once("Erreur : Aucun fichier fourni")
    sys.exit(1)
else:
    print_once("Patients ...... OK")

# Verify genome file path
genome = get_config_value(config["GENERAL"].get("GENOME"), "", str)
genome_path = os.path.join(working_directory, "1-raw_data/references", genome) if not os.path.isfile(genome) else genome

if not os.path.isfile(genome_path):
    print_once(f"Erreur : Le fichier génome {genome_path} n'existe pas.")
    sys.exit(1)
else:
    print_once(f"Fichier {genome_path} ...... OK")
    genome = os.path.abspath(genome_path)
    name_genome = os.path.splitext(os.path.basename(genome))[0]

# Check PIGZ executable
pigz = get_config_value(config["DEPENDANCES"]["FORMATING"].get("PIGZ"), "pigz", str)
check_executable(pigz)

include: os.path.join(MOD, "compress_fastq.smk")
print_once("Module compress_fastq ...... OK")

# Check trimming path
path_fastq = f"{working_directory}/1-raw_data/fastq/"
path_bam = f"{working_directory}/2-processed_data/BAM/"

if use_trimming:
    fastp = get_config_value(config["DEPENDANCES"]["FORMATING"].get("FASTP"), "fastp", str)
    check_executable(fastp)

    path_fastq = f"{working_directory}/2-processed_data/fastq_trimmed/"

    include: os.path.join(MOD, "trimming_fastq.smk")
    print_once("Module trimming_fastq ...... OK")

# Check quality control executables
use_qc = get_config_value(config["USAGE"].get("QC"), int(0))

if use_qc:
    fastqc = get_config_value(config["DEPENDANCES"]["QC"].get("FASTQC"), "fastqc", str)
    multiqc = get_config_value(config["DEPENDANCES"]["QC"].get("MULTIQC"), "multiqc", str)

    check_executable(fastqc)
    check_executable(multiqc)

    include: os.path.join(MOD,"quality_control_fastq.smk")
    print_once("Module quality_control_fastq ...... OK")

    # Define QC output directories
    directory_data_raw = f"{path_qc}/multiqc/fastq_raw/{prefix}_{unique_id}_data/"
    html_raw = f"{path_qc}/multiqc/fastq_raw/{prefix}_{unique_id}.html"

    list_inputs.append(directory_data_raw)
    list_inputs.append(html_raw)

    if use_trimming:
        directory_data_trimmed = f"{path_qc}/multiqc/fastq_trimmed/{prefix}_{unique_id}_data/"
        html_trimmed = f"{path_qc}/multiqc/fastq_trimmed/{prefix}_{unique_id}.html"

        list_inputs.append(directory_data_trimmed)
        list_inputs.append(html_trimmed)

# Check mapping executables
use_mapping = get_config_value(config["USAGE"].get("MAPPING"), int(0))

if use_mapping:
    STAR = get_config_value(config["DEPENDANCES"]["MAPPING"].get("STAR"), "STAR", str)
    samtools = get_config_value(config["DEPENDANCES"]["MAPPING"].get("SAMTOOLS"), "samtools", str)
    Rscript = get_config_value(config["DEPENDANCES"]["GENERAL"].get("RSCRIPT"), "Rscript", str)

    check_executable(STAR)
    check_executable(samtools)
    check_executable(Rscript)
    
    # Check SpliceLauncher path
    SpliceLauncher = get_config_value(config["DEPENDANCES"]["ANALYSES"].get("SPLICELAUNCHER"), os.path.join(os.getcwd(), os.path.join(ROOTDIR,"SpliceLauncher")), str)
    if not os.path.isdir(SpliceLauncher):
        print_once(f"Erreur : Le dossier {SpliceLauncher} n'existe pas.")
        sys.exit(1)
    else:
        print_once(f"Dossier {SpliceLauncher} ...... OK")

    # Validate GFF3 file
    gff3 = get_config_value(config["GENERAL"].get("GFF3"), "", str)
    gff3_path = gff3 if os.path.isfile(gff3) else f"{working_directory}/1-raw_data/annotations/{gff3}"

    if not os.path.isfile(gff3_path):
        print_once(f"Erreur : Le fichier GFF3 {gff3_path} n'existe pas.")
        sys.exit(1)
    else:
        print_once(f"Fichier {gff3_path} ...... OK")

    # Validate MANE file
    mane = get_config_value(config["GENERAL"].get("MANE"), "MANE.to.SpliceLauncher.txt", str)
    mane_path = mane if os.path.isfile(mane) else f"{working_directory}/1-raw_data/annotations/{mane}"

    if not os.path.isfile(mane_path):
        print_once(f"Erreur : Le fichier MANE {mane_path} n'existe pas.")
        sys.exit(1)
    else:
        print_once(f"Fichier {mane_path} ...... OK")

    include: os.path.join(MOD, "mapping.smk")
    print_once("Module mapping.smk ...... OK")

    if use_qc:
        include: os.path.join(MOD, "quality_control_bam.smk")
        print_once("Module quality_control_bam.smk ...... OK")

        # Define BAM QC output directories
        directory_data_bam = f"{path_qc}/multiqc/BAM/{name_genome}/{prefix}_{unique_id}_data/"
        html_bam = f"{path_qc}/multiqc/BAM/{name_genome}/{prefix}_{unique_id}.html"

        list_inputs.append(directory_data_bam)
        list_inputs.append(html_bam)

    # Check SpliceLauncher dependencies
    use_SpliceLauncher = get_config_value(config["USAGE"].get("SPLICELAUNCHER"), int(0))

    if use_SpliceLauncher:
        # Check general dependencies
        perl = get_config_value(config["DEPENDANCES"]["GENERAL"].get("PERL"), "perl", str)
        bedtools = get_config_value(config["DEPENDANCES"]["ANALYSES"].get("BEDTOOLS"), "bedtools", str)
        awk = get_config_value(config["DEPENDANCES"]["GENERAL"].get("AWK"), "awk", str)
        
        check_executable(perl)
        check_executable(bedtools)
        check_executable(awk)

        # Check sashimi plot dependencies
        use_sashimi = get_config_value(config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("USE"), int(0))

        if use_sashimi:

            python = get_config_value(config["DEPENDANCES"]["GENERAL"].get("PYTHON"), "python", str)
            ggsashimi = get_config_value(config["DEPENDANCES"]["VISUALISATION"].get("GGSASHIMI"),"ggsashimi.py", str)

            check_executable(python)
            check_executable(ggsashimi)

            # Define sashimi plot directories
            directories_non_statistical_junctions = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/sashimi_plot/non_statistical_junctions/", reads=all_samples)
            directories_statistical_junctions = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/sashimi_plot/statistical_junctions/", reads=all_samples)

            list_inputs.append(directories_non_statistical_junctions)
            list_inputs.append(directories_statistical_junctions)

        include: os.path.join(MOD, "SpliceLauncher.smk")
        print_once("Module SpliceLauncher ...... OK")

        # Define output report paths
        count_report = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_report_{date}.txt"
        directory_count_results = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results"
        
        RecapFile = expand(f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results/samples_results/{{reads}}/{{reads}}.recap{extension}",reads= all_samples)

        list_inputs.append(RecapFile)

    if not list_inputs:
        bai = expand(f"{path_bam}{name_genome}/mapping/{{reads}}.sorted.bam.bai",reads= all_samples)
        list_inputs.append(bai)

else:
    print_once("Nothing mapped ...... nothing analyzed")

if not list_inputs:
    read1 = expand(os.path.abspath(f"{path_fastq}{{reads}}.1.fastq.gz"),reads= all_samples)
    read2 = expand(os.path.abspath(f"{path_fastq}{{reads}}.2.fastq.gz"),reads= all_samples)
    list_inputs.append(read1)
    list_inputs.append(read2)

    if use_trimming:
        html = expand(f"{path_qc}/fastp_trimming/{{reads}}.html",reads= all_samples)
        json = expand(f"{path_qc}/fastp_trimming/{{reads}}.json",reads= all_samples)
        list_inputs.append(html)
        list_inputs.append(json)

os.environ["SNAKEMAKE_PRINT"] = "true"

# Define the main rule for the Snakemake workflow
rule all:
    input:
        # List of required files for a complete execution of the pipeline
        list_inputs
    
    params:
        bash = bash,
        script = os.path.join(SCR,"clear_cache.sh") 

    shell:
        # Cache cleanup after pipeline execution to free up resources
        "{params.bash} {params.script}"

# Hook `onsuccess`: Actions to execute after a successful pipeline run
onsuccess:
    # Free memory after a successful execution to avoid unnecessary load
    gc.collect()

# Hook `onerror`: Error handling and cleanup of files in case of failure
onerror:
    # Define generated files that need to be removed upon errors
    merged_files = f"{path_results}/SpliceLauncher/merged_files/{prefix}_{unique_id}.txt"
    count_report = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_report_{date}.txt"
    directory_count_results = f"{path_results}/SpliceLauncher/{prefix}_{unique_id}_results"

    # Check file existence before deletion to prevent errors
    files_to_remove = [merged_files, count_report, directory_count_results]
    existing_files = [file for file in files_to_remove if os.path.exists(file)]

    # Delete only if files exist to avoid unnecessary deletion attempts
    clear_cache = os.path.join(SCR,"clear_cache.sh")
    if existing_files:
        command = f"rm -rf {' '.join(existing_files)} && {bash} {clear_cache}"
        shell(command)

    # Free memory after error handling to optimize system performance
    gc.collect()