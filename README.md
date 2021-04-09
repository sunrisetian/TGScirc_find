# TGScirc_find
  A software that can predict circRNA from third generation sequencing data

---

## Requirement
### Data:

- Genome sequence (fasta format)
- RNA_seq data (fastq format)
### Software:

- blast(v.2.6.0+)
- blat
- gamp(v.2015-09-21)
- minimap2(v2.1)
- Python3 (v.3.7.7): (https://www.python.org/)

### Python3 package:

- Biopython (v.1.76): (https://pypi.org/project/biopython/)
- pandas
- numpy

---
## Download
  Open the terminal and input:
  ```bash
  git clone https://github.com/sunrisetian/TGScirc_find.git
  ```
---
## Usage

### - You can run TGScirc_find step by step use command line

   Circular RNA predicting
      You need to provide the following files:
     - Genome sequence (fasta format)
     - query data file (fasta format)
      
  ```bash
   python3 TGScirc_find.py -g genome_filename  -q query_filename -o output_file -t 4
  ```
  You can get the forecast in ***output_dfile***.You can choose the number of threads through the *-t*.

## Contact us

If you encounter any problems while using Pcirc, please send an email (glli@snnu.edu.cn) or submit the issues on GitHub (https://github.com/Lilab-SNNU/TGScirc_find/issues) and we will resolve it as soon as possible.
