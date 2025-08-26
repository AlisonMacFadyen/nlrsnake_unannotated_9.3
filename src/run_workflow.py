#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import yaml
from pathlib import Path



with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "VERSION"), "r") as f:
    __version__ = f.read().strip()


def parse_args():
    parser = argparse.ArgumentParser(description='Run the UNNAMED DRAFT pipeline',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  run_workflow --config myproject/config.yaml              # Run pipeline with custom config
  run_workflow --dry-run                                   # Perform a dry run
  run_workflow --rule the_second_rule                      # Run until a specific rule
  run_workflow --unlock                                    # Unlock a snakemake directory
  run_workflow --version                                   # Display version information
  run_workflow --dag                                       # Generate DAG visualization

Configuration:
  The config file must contain all required parameters defined in this script.

Documentation:
  For full documentation, visit: https://github.com/publish_your_repo
""")
    parser.add_argument('--config', default='lib/config.yaml', help='Path to config file')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run')
    parser.add_argument('--unlock', action='store_true', help='Unlock a locked directory')
    parser.add_argument('--rule', help='Run a specific rule (e.g., generate_report)')
    parser.add_argument('--force', action='store_true', help='Force run the rule even if outputs exist')
    parser.add_argument('--version', action='store_true', help='Display version information')
    parser.add_argument('--dag', action='store_true', help='Generate a PDF of the DAG (workflow_dag.pdf)')

    return parser.parse_args()

class ConfigValidationError(Exception):
    """Exception raised for errors in the configuration validation."""
    pass

def validate_config(config_data, expected_keys=None):
    """
    Validate configuration file to ensure it has all required keys.
    Additional keys are allowed for user customization.
    
    Args:
        config_data (dict): The loaded configuration data
        expected_keys (list): List of required keys in the config file
        
    Raises:
        ConfigValidationError: If the configuration is invalid
    """
    if expected_keys is None:
        # Define all required keys based on your config.yaml
        expected_keys = [
            "scratch",
            "workdir",
            "singularity_image",
            "wckey",
            "main_job_partition",
            "profile_dir"
        ]
    
    # Check for missing keys
    missing_keys = [key for key in expected_keys if key not in config_data]
    if missing_keys:
        raise ConfigValidationError(f"Missing required configuration keys: {', '.join(missing_keys)}")
    
    # Check for extra keys and inform user it will be applied to workflow
    extra_keys = [key for key in config_data if key not in expected_keys]
    if extra_keys:
        print(f"Info: Additional configuration keys found (will be passed to workflow): {', '.join(extra_keys)}")
    
    # Check that required values are not None or empty strings
    empty_keys = [key for key in expected_keys if key in config_data and 
                  (config_data[key] is None or (isinstance(config_data[key], str) and config_data[key].strip() == ""))]
    if empty_keys:
        raise ConfigValidationError(f"Empty values for required configuration keys: {', '.join(empty_keys)}")
        
    return True

def load_config(config):
    """Load configuration from YAML file."""
    if not os.path.exists(config):
        print(f"Error: Configuration file not found: {config}", file=sys.stderr)
        sys.exit(1)
        
    try:
        with open(config, 'r') as config_file:
            config_data = yaml.safe_load(config_file)
        
        # Validate the configuration
        try:
            validate_config(config_data)
        except ConfigValidationError as e:
            print(f"Error validating configuration: {e}", file=sys.stderr)
            sys.exit(1)
            
        return config_data
    except yaml.YAMLError as e:
        print(f"Error parsing YAML configuration: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading configuration: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    args = parse_args()

    if args.version:
        print(f"Unnamed Draft Pipeline v{__version__}")
        sys.exit(0)

    config = load_config(args.config)
    
    # Set working directory to project root
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(project_root)

    
    OUTPUT_DIR = os.path.join(config["scratch"], config["workdir"])


    if args.unlock:
        unlock_cmd = [
            f'source activate  ./snakemake_env;',
            'snakemake -s src/workflow.smk',
            f'--configfile {args.config}',
            '--unlock'
        ]
        unlock_cmd_str = ' '.join(unlock_cmd)
        print(f"Unlocking workflow: {unlock_cmd_str}")
        subprocess.run(unlock_cmd_str, shell=True)
        sys.exit(0)

        # If DAG generation is requested, do that and exit
    if args.dag:
        dag_cmd = [
            "source activate ./snakemake_env;",
            'snakemake', '-s', 'src/workflow.smk',
            '--configfile', args.config,
            '--dag',
            '|',
            'dot',
            '-Tpdf',
            '>',
            f"{os.path.join(OUTPUT_DIR, 'workflow_dag.pdf')}"
        ]
        dag_cmd_str = ' '.join(dag_cmd)
        print("Generating DAG visualization...")
        subprocess.run(dag_cmd_str, shell=True, check=True)
        print(f"DAG visualization saved to {os.path.join(OUTPUT_DIR, 'workflow_dag.pdf')}" )
        sys.exit(0)
    
    
    # Build the command

    #make sure job goes first to sbatch
    cmd = ['sbatch', '--partition', f'{config["main_job_partition"]}', '--wckey', f'{config["wckey"]}', '-J', f'{config["workdir"]}', '--wrap="' ]

    
    # call snakemake
    cmd.extend(['source', 'activate', './snakemake_env', ';', 'snakemake', '-s', 'src/workflow.smk'])
    
    # Add config
    cmd.extend(['--configfile', args.config])

    if args.rule:
        cmd.extend([f'--until {args.rule}'])
        
        # Add force flag if specified
        if args.force:
            cmd.extend(['--forcerun', args.rule])

    # Add snakemake process params
    cmd.extend([ 
        f'--workflow-profile {config["profile_dir"]}',
        '--executor slurm'
        ])
    
    
    # Dry run
    if args.dry_run:
        cmd.append('-n')
    
    
    # Run the command
    cmd.extend('"') #conclude the wrap part
    cmd_str = ' '.join(cmd)
    print(f"Running command: {cmd_str}")
    subprocess.run(cmd_str, shell=True)

if __name__ == '__main__':
    main()
