# Genome-Repair-Tools
Genome Repair Tools (GRT) is a genomic repair toolkit primarily designed for filling missing regions (Gaps) in genome FASTA files, with additional telomere recovery functionality
Usage

Place all scripts in a directory, then run env.sh to set up the environment dependencies. The env.sh script will use Mamba to quickly install required software and create a conda environment named GRT. The genome_repair_tools.py software will run within this environment.

bash
### Run environment setup script
bash env.sh

### Add the software to your PATH (add current directory permanently to ~/.bashrc)
echo "export PATH=\"\$PATH:$PWD\"" >> ~/.bashrc

### Refresh environment variables
source ~/.bashrc

### Activate the environment
mamba activate GRT

### Make the software executable
chmod +x genome_repair_tools.py

### Test if the software runs correctly (from any directory)
genome_repair_tools.py -h
Genome Repair Tools Software Usage

bash
### Basic usage
genome_repair_tools.py -q genome_to_repair.fasta -c reference_contigs.fasta -t threads -o output_directory

### When providing both reference contigs and raw reads, you can choose the assembly software.
### The default will run these 5 assemblers: hifiasm, verkko, nextDenovo, flye, and shasta.
genome_repair_tools.py -q genome_to_repair.fasta -c reference_contigs.fasta --hifi HiFi_data.fastq.gz -t threads -o output_directory --run-hifiasm --run-verkko

## Genome repair tools software flowchart
<img src="https://github.com/user-attachments/assets/3e6c5d49-46dd-4ba3-b9cd-dd575bd08832" width="80%" />
