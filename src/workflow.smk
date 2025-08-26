'''Run NLRtracker pipeline on TSL HPC

Usage: bash scripts/do_nlrtrack.sh

'''

import os
import csv
import shutil
import gzip

from pathlib import Path
from snakemake.io import glob_wildcards
from snakemake.utils import min_version

# Specify minimum version of snakemake required
min_version("9.3.0")

configfile: "lib/config.yaml"

# Specify the container
container: config["singularity_image"]

# Specify and create path for working directory
OUTPUT_DIR = os.path.join(config["scratch"], config["workdir"])

# Create a log directory
LOG_DIR = os.path.join(OUTPUT_DIR, "logs")

# Define the final output
FINAL_OUTPUT = "done.txt"

# Function to get samples and fasta paths
def get_fa_samples(file):
    samples = []
    fa = []
    with open(file, 'r') as input:
        for line in input:
            line = line.rstrip()
            items = line.split(",")
            # Adding check to ensure samples_to_Fasta file is formatted correctly
            if len(items) < 2:
                raise ValueError(f"Malformed line in {file}: {line}")
            samples.append(items[0])
            fa.append(items[1])
    return fa, samples

# Function to map sample to fasta path
def sample_to_fa(sample, fastas, samples):
    return fastas[samples.index(sample)]

FASTAS, SAMPLES = get_fa_samples(config['sample_fasta_file'])

# Main rule to ensure correct order of execution
rule all:
    input: FINAL_OUTPUT

# Rule to run Helixer
rule run_helixer:
    input:
        # Lambda function to assign fasta paths from config .csv
        fasta=lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES)
    output:
        helixer=os.path.join(OUTPUT_DIR, "{sample}", "helixer", "{sample}_helixer.gff")
    log:
        os.path.join(LOG_DIR, "run_helixer_{sample}.log")
    resources:
        mem_mb=32000,
        partition="tsl-medium",
        slurm_extra=f"--wckey={config['wckey']}"
    params:
        lineage="land_plant",
        species=config['species'],
        model_path="/tsl/data/helixer/models/land_plant/land_plant_v0.3_a_0080.h5",
        subsequence_length=64152,
        additional_options=config.get('helixer_options', ''),
    threads: 16
    shell:
        """
        Helixer.py --lineage {params.lineage} \
        --fasta-path {input.fasta} \
        --species {params.species} \
        --gff-output-path {output.helixer} \
        --model-filepath {params.model_path} \
        --subsequence-length {params.subsequence_length} \
        {params.additional_options}
        """

# Rule to run gffread
rule run_gffread:
    input:
        # Lambda function to assign fasta paths from config .csv
        fasta=lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES),
        helixer=os.path.join(OUTPUT_DIR, "{sample}", "helixer", "{sample}_helixer.gff")
    output:
        gffread=os.path.join(OUTPUT_DIR,"{sample}", "gffread", "{sample}_gffread.fasta")
    log:
        os.path.join(LOG_DIR, "run_gffread_{sample}.log")
    resources:
        mem_mb=32000,
        partition="tsl-short",
        slurm_extra=f"--wckey={config['wckey']}"
    params:
        temp_dir=directory(os.path.join(OUTPUT_DIR, "{sample}", "gffread"))
    threads: 16
    run:
        # Use python to create temp directory
        os.makedirs(params.temp_dir, exist_ok=True)
        
        # Create a temporary uncompressed file
        temp_fasta = os.path.join(params.temp_dir, "temp.fasta")
        
        # Check if input fasta is compressed using python and without wildcards
        # Improves compatability across environments
        if input.fasta.endswith(".gz"):
            with gzip.open(input.fasta, 'rb') as fasta_in, open(temp_fasta, 'wb') as fasta_out:
                shutil.copyfileobj(fasta_in, fasta_out)
        else:
            shutil.copy(input.fasta, temp_fasta)
        
        # Run gffread
        shell("gffread -g {temp_fasta} -y {output.gffread} {input.helixer}")
        
        # Remove temporary file
        os.remove(temp_fasta)

# Rule to strip asterisks
rule strip_asterisks:
    input:
        fasta=os.path.join(OUTPUT_DIR, "{sample}", "gffread", "{sample}_gffread.fasta")
    output:
        no_asterisk=os.path.join(OUTPUT_DIR, "{sample}", "no_asterisk.fa")
    log:
        os.path.join(LOG_DIR, "strip_asterisks_{sample}.log")
    resources:
        mem_mb=4000,
        partition="tsl-short",
        slurm_extra=f"--wckey={config['wckey']}"
    shell:
        "cat {input.fasta} | sed 's/*//g' > {output.no_asterisk}"

# Checkpoint for interpro
checkpoint interpro:
    input:
        config['scratch'] + "/{sample}/no_asterisk.fa"
    output:
        fa=directory(os.path.join(OUTPUT_DIR, "{sample}", "ipro_fa"))
    resources:
        mem_mb=4000,
        partition="tsl-short",
        slurm_extra=f"--wckey={config['wckey']}"
    params:
        seqs_per_file=config['seqs_per_file'],
        hidden_done=os.path.join(OUTPUT_DIR, "{sample}", "ipro_fa")
    threads: 1
    shell: 
        """bash scripts/split_fa.sh {input} {params.seqs_per_file} {output}
        # Add hidden file for tracking if pipeline needs to be started again
        echo done > {params.hidden_done}/.done
        """

# Rule to run interpro on split files
rule run_interpro:
    input:
        fa=os.path.join(OUTPUT_DIR, "{sample}", "ipro_fa", "{n}.fa")
    output:
        gff=os.path.join(OUTPUT_DIR, "{sample}", "ipro_gff", "{n}.gff")
    log:
        os.path.join(LOG_DIR, "run_interpro_{sample}.log")
    resources:
        mem_mb=16000,
        partition="tsl-short",
        slurm_extra=f"--wckey={config['wckey']}"
    params:
        tmp=os.path.join(OUTPUT_DIR, "{sample}", "ipro_tmp")
    threads: 16
    shell:
        """
        interproscan.sh -i {input.fa} \
        -f gff3 \
        -T {params.tmp} \
        -cpu {threads} \
        -o {output.gff} \
        -appl Pfam,Gene3D,SUPERFAMILY,PRINTS,SMART,CDD,ProSiteProfiles \
        -dp
  """

# Rule to run FIMO
rule run_fimo:
    input:
        # Pass gffread output fasta directly
        fasta=config['scratch'] + "/{sample}/gffread/{sample}_gffread.fasta"
    output:
        fimo="{sample}/results/fimo_out/fimo.gff"
    params:
        mem="32G",
        queue="tsl-short",
        fimo_dir="{sample}/results/fimo_out"
    threads: 16
    run:
        # Use python to create fimo output directory
        os.makedirs(params.fimo_dir, exist_ok=True)

        shell("bash scripts/fimo.sh {params.fimo_dir} lib/meme.xml {input.fasta}")

# Rule to run HMMER
rule run_hmmer:
    input:
        # Pass gffread output fasta directly
        fasta=config['scratch'] + "/{sample}/gffread/{sample}_gffread.fasta
    output:
        hmmer="{sample}/results/CJID.txt"
    params:
        mem="32G",
        queue="tsl-short"
    threads: 16
    shell:
        "bash scripts/hmmer.sh {output.hmmer} lib/abe3069_Data_S1.hmm {input.fasta}"

# Update to ensure files within checkpoint_output are tracked
def aggregate_ipro(wildcards):
    checkpoint_output = checkpoints.interpro.get(**wildcards).output[0]

    # Ensure checkpoint completed by checking `.done` file
    done_file = Path(checkpoint_output) / ".done"
    if not done_file.exists():
        raise ValueError(f"Checkpoint output for {wildcards.sample} is incomplete!")

    # Get list of `.fa` files
    fa_files = list(Path(checkpoint_output).glob("*.fa"))
    n_values = [f.stem for f in fa_files]

    # Expand expected gff output
    return expand(config['scratch'] + "/{sample}/ipro_gff/{n}.gff",
                  sample=wildcards.sample,
                  n=n_values)

rule aggregate:
    input:
        aggregate_ipro
    output:
        "{sample}/results/interpro_result.gff"
    params:
        mem="8G",
        queue="tsl-short",
    threads: 1
    shell:
        "cat {input} > {output}"

# Rule to run NLRtracker
rule run_nlrtracker:
    input:
        interpro="{sample}/results/interpro_result.gff",
        fimo="{sample}/results/fimo_out/fimo.gff",
        hmmer="{sample}/results/CJID.txt",
        fasta=config['scratch'] + "/{sample}/gffread/{sample}_gffread.fasta"
    output:
        done="{sample}/results/done.txt"
    params:
        mem="32G",
        queue="tsl-short",
        sample_name="{sample}",
        nlrtracker_dir="{sample}/nlrtracker_out"
    threads: 16
    run:
        # Use python to create nlrtracker output directory
        os.makedirs(params.nlrtracker_dir, exist_ok=True)

        shell("""
        bash scripts/run_tracker.sh lib/interproscan_v5.71-102.0_entry.list {input.interpro} {input.fimo} \
        {input.fasta} {params.nlrtracker_dir} {params.sample_name} p {input.hmmer} lib/iTOL_NLR_template.txt
        echo "done" > {output.done}
        """)

# Rule to finalize the pipeline
rule finalize:
    input:
        expand("{sample}/results/done.txt", sample=SAMPLES)
    output:
        touch(FINAL_OUTPUT)
    params:
        mem="4G",
        queue="tsl-short"
    threads: 1
    shell:
        "echo 'Pipeline completed successfully' > {output}"
