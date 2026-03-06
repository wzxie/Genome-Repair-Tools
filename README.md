# Genome-Repair-Tools
Genome Repair Tools (GRT) is a genomic repair toolkit primarily designed for filling missing regions (Gaps) in genome FASTA files, with additional telomere recovery functionality.

![流程图1 1](https://github.com/user-attachments/assets/62ce2b00-5166-4b34-9a81-fad947eaf2d4)

## Installation
### Download software
```bash
git clone  https://github.com/wzxie/Genome-Repair-Tools.git
```
### Enter the software directory
```bash
cd  Genome-Repair-Tools/scripts/
```
### Run the env.sh script to install the environment dependencies. The env.sh script will call Mamba to quickly install the dependent software and create an environment called GRT.
```bash
bash env.sh
```
### Add the software to your PATH (add current directory permanently to ~/.bashrc)
```bash
echo "export PATH=\"\$PATH:$PWD\"" >> ~/.bashrc
```
### Refresh environment variables
```bash
source ~/.bashrc
```
### Activate the environment
```bash
mamba activate GRT
```
### Make the software executable
```bash
chmod +x genome_repair_tools.py
```
### Test whether the software can run and switch to other directories
```bash
cd ~
genome_repair_tools.py -h
```
## Or you can install it manually
### Installation dependency required
mamba install -n $ENV_NAME \
```bash
    python=3.11.14 \
    tgsgapcloser=1.2.1 \
    flye=2.9.6 \
    nextdenovo=2.5.2 \
    verkko=2.2.1 \
    mummer4 \
    hifiasm=0.25.0 \
    bioconda::shasta=0.14.0 \
    bioconda::craq \
    bioconda::merqury
```
### Download the script in the same way as above

## Genome Repair Tools Software Usage
### Basic usage
```bash
genome_repair_tools.py -q genome_to_repair.fasta -c reference_contigs.fasta -t threads -o output_directory
```
### Full usage
Providing raw sequencing data. The default will run these 5 assemblers: hifiasm, verkko, nextDenovo, flye, and shasta.You can choose the assembly software,such as ```--run-hifiasm --run-verkko```
```bash
genome_repair_tools.py -q genome_to_repair.fasta -c reference_contigs.fasta --hifi HiFi_data.fastq.gz -t threads -o output_directory 
```
## Sample data download
### Download sample genome data (5MB)
```bash
wget https://github.com/wzxie/Genome-Repair-Tools/releases/download/v1.0-data/example_genome.fa
```
### Download sample sequencing data (71MB)
```bash
wget https://github.com/wzxie/Genome-Repair-Tools/releases/download/v1.0-data/example_HiFi.fastq.gz
```
### Test data usage: 
Run the following command to generate five output directories: ```merged_contigs```, ```patch_repair```, ```correct_refill```, ```telomere_recover```, and ```assemble_fill```. The gap-filled test genome result, assembled_filled_renamed.fasta, can be found in the assemble_fill folder. The entire process takes approximately 5 minutes.
```bash
genome_repair_tools.py -q PATH/example_genome.fa --hifi PATH/example_HiFi.fastq.gz -t 32 -o output_directory
```
## Detailed usage instructions 
- Chinese version: https://share.note.sx/iun7yhmf#TG3w7haQi1Xj3BjMcNYvq2AatNUlb7InqQENd0o1AUA
- English version: https://share.note.sx/6ga80su1#G3Yetb5y7mtHswv4QqC58Q9fnXQ0wfNgnxU+xIpMOwM
# Contact
We hope this tools could be helpful for the groups which focused on plants genome assembly, you can use the GitHub page to report issues or email us with any suggestions.

- Yi-Rui Zhu: zyr18803205362@163.com
- Zu-Wen Zhou: 784012725@qq.com
- Ru-Peng Zhao: 2247290650@qq.com
- Wen-Zhao Xie: wzxie@hebtu.edu.cn
