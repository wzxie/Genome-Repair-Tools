# Genome-Repair-Tools
Genome Repair Tools (GRT) is a genomic repair toolkit primarily designed for filling missing regions (Gaps) in genome FASTA files, with additional telomere recovery functionality
Usage


## Software installation method
### download software
git clone  https://github.com/wzxie/Genome-Repair-Tools.git
### Enter the software directory
cd  Genome-Repair-Tools/scripts/
### Run the env.sh script to install the environment dependencies. The env.sh script will call Mamba to quickly install the dependent software and create an environment called GRT.
bash env.sh

### Add the software to your PATH (add current directory permanently to ~/.bashrc)
echo "export PATH=\"\$PATH:$PWD\"" >> ~/.bashrc

### Refresh environment variables
source ~/.bashrc

### Activate the environment
mamba activate GRT

### Make the software executable
chmod +x genome_repair_tools.py

### Test whether the software can run and switch to other directories
cd ~
genome_repair_tools.py -h

## Genome Repair Tools Software Usage
### Basic usage
genome_repair_tools.py -q genome_to_repair.fasta -c reference_contigs.fasta -t threads -o output_directory
### When providing both reference contigs and raw reads, you can choose the assembly software.
### The default will run these 5 assemblers: hifiasm, verkko, nextDenovo, flye, and shasta.
genome_repair_tools.py -q genome_to_repair.fasta -c reference_contigs.fasta --hifi HiFi_data.fastq.gz -t threads -o output_directory --run-hifiasm --run-verkko

## Genome repair tools software flowchart
<p align="center">
  <img src="https://github.com/user-attachments/assets/3e6c5d49-46dd-4ba3-b9cd-dd575bd08832" width="90%" alt="GRT Software Flowchart" />
</p>

## Detailed usage instructions 
https://share.note.sx/iun7yhmf#TG3w7haQi1Xj3BjMcNYvq2AatNUlb7InqQENd0o1AUA

# Contact
We hope this tools could be helpful for the groups which focused on plants genome assembly, you can use the GitHub page to report issues or email us with any suggestions.
- Xie wenzhao: wzxie@hebtu.edu.cn
- Zhu Yirui: zyr18803205362@163.com
- Zhao rupeng: 2247290650@qq.com

