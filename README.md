# Genome-Repair-Tools
Genome Repair Tools (GRT) is an automated pipeline for closing gaps in genome assemblies, with the primary goal of improving assembly contiguity. It is specifically designed to fill missing regions (gaps) in scaffolded or chromosome‑level draft genomes using long‑read sequencing data. Users only need to provide a draft genome together with third‑generation sequencing reads (e.g., PacBio HiFi or ONT); GRT will then automatically locate and repair gaps without manual intervention. For advanced use, individual assembly modules can be selectively enabled or run independently, offering flexibility for custom workflows. In addition, GRT provides a lightweight telomere recovery function that restores terminal repeats when feasible. A schematic overview of the entire workflow is shown in the figure below.

![流程图1 1](https://github.com/user-attachments/assets/62ce2b00-5166-4b34-9a81-fad947eaf2d4)

## Dependencies
1. python (3.11.14)
2. tgsgapcloser (v1.2.1)
3. flye (v2.9.6)
4. nextdenovo (v2.5.2)
5. verkko (v2.2.1)
6. mummer4 (4.0.0rc1)
7. hifiasm (0.25.0)
8. craq (latest)
9. merqury (latest)

## Installation
Run the following command to install GRT and its dependencies.
```
# 1. Clone the repository and install the environment dependencies
git clone https://github.com/wzxie/Genome-Repair-Tools.git && cd Genome-Repair-Tools/scripts/
bash env.sh

# 2. Activate the GRT environment, set executable permissions, and add to PATH
mamba activate GRT
chmod +x genome_repair_tools.py
echo "export PATH=\"\$PATH:$PWD\"" >> ~/.bashrc

# 3. Verify the installation
genome_repair_tools.py -h
```


## Genome Repair Tools Software Usage
### Run the following command to test the pipeline on the example data (approx. 5 min):
```
genome_repair_tools.py -q ./example.fa \
    --hifi ./example_HiFi.fastq.gz \
    -t 32 \
    -o ./output_directory
```
### Expected output
| Directory | Description |
|-----------|-------------|
| `merged_contigs` | Contigs merged from multiple assemblies (intermediate) |
| `patch_repair` | Locally patched and repaired regions (intermediate) |
| `correct_refill` | Corrected and refilled gaps (intermediate) |
| `assemble_fill` | Initially filled genomes (intermediate) |
| `telomere_recover` | input/complete_genome.fasta ；output/final_genome.fasta |

## Detailed usage instructions 
- Chinese version: https://share.note.sx/iun7yhmf#TG3w7haQi1Xj3BjMcNYvq2AatNUlb7InqQENd0o1AUA
- English version: https://share.note.sx/6ga80su1#G3Yetb5y7mtHswv4QqC58Q9fnXQ0wfNgnxU+xIpMOwM
# Contact
We hope this tools could be helpful for the groups which focused on plants genome assembly, you can use the GitHub page to report issues or email us with any suggestions.

- Yi-Rui Zhu: zyr18803205362@163.com
- Zu-Wen Zhou: 784012725@qq.com
- Ru-Peng Zhao: 2247290650@qq.com
- Wen-Zhao Xie: wzxie@hebtu.edu.cn
