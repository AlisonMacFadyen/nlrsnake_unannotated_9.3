#!/bin/python

# This script checks that a user has access to slurm and have installed snakemake correctly

# Requires check that snakemake is present in the directory via path snakemake_env/bin/snakemake
# Also requires check if sbatch is available within the path
# Could use running the "help" as a way to visually show the user e.g. sbatch -h and ./snakemake_env/bin/snakemake -h

import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='This script checks for correct snakemake enivronment installation and access to slurm.',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This command does not need any arguments by default as it assumes the snakemake_env is installed in the project directory at './snakemake_env'.  If the installation is located elsewhere, use the '-e/--enivronment_location' flag to specify where.
""")

parser.add_argument("-e", "--environment_location", help="path to installed snakemake environment.")

args = parser.parse_args()

if args.environment_location:
    # use that input to search
    snakemake_env = args.environment_location
    print(f"User provided the following path for snakemake_env: {snakemake_env}")
else:
    # use default
    snakemake_env = "./snakemake_env"
    print(f"\nUsing default path for snakemake_env: {snakemake_env}") 

slurm_status, slurm_result = subprocess.getstatusoutput("sbatch -h")
snake_status, snake_result = subprocess.getstatusoutput(f"source activate {snakemake_env}; snakemake -h")

missing = "command not found"
activate_help = "ActivateHelp: usage:"

slurm_install = False
snakemake_install = False


if missing in snake_result or activate_help in snake_result:
    print(f"""
          {snake_result}
          
          snakemake command is missing from your path, snakemake is not installed correctly?
          """)
else:
    print(f"""
          {snake_result}
          
          Help message printed as expected, snakemake is installed correctly.
          """)
    snakemake_install = True


if missing in slurm_result:
    print(f"""
          {slurm_result}
          
          sbatch is missing from your path, do you have access to SLURM?
          """)
else:
    print(f"""
          {slurm_result}
          
          Help message printed as expected, you have access to SLURM.
          """)
    slurm_install = True

if slurm_install == True and snakemake_install == True:
    print("Your installation has completed successfully.")
elif slurm_install == True and snakemake_install == False:
    print("You have access to SLURM but snakemake is not installed correctly.")
elif slurm_install == False and snakemake_install == True:
    print("You do not have access to SLURM but snakemake installed correctly.")
elif slurm_install == False and snakemake_install == False:
    print("You do not have access to SLURM and snakemake is not install correctly.")
else:
    print("Shouldn't see this message.")

subprocess.getstatusoutput(f"source deactivate")
