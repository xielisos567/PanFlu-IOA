# PanFlu-IOA (Pan-Influenza Iteratively Optimized Antigen)

PanFlu-IOA is a bioinformatics toolkit designed to overcome systemic bias and extract universal antigenic features from influenza sequence data through iterative multi‑sequence alignment and consensus generation. It can be applied to any viral protein sequence dataset to derive representative antigenic signatures.

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
  - [1. Install MAFFT](#1-install-mafft)
  - [2. Clone the repository and install Python dependencies](#2-clone-the-repository-and-install-python-dependencies)
- [Demo](#demo)
  - [Run with provided test_data](#run-with-provided-test_data)
  - [Expected output & run time](#expected-output--run-time)
- [Instructions for Use](#instructions-for-use)
  - [Input data format](#input-data-format)
  - [Output files](#output-files)
- [License](#license)
- [Contact](#contact)

---

## Overview

PanFlu-IOA provides two main operational modes to extract robust antigenic signatures from large sequence sets:

1. **Mode 1 (`PanFlu-IOA_code1.py`)** – Iterative random subsampling, alignment, and consensus building over multiple rounds.
2. **Mode 2 (`PanFlu-IOA_code2.py`)** – Two‑stage random subsampling and alignment, followed by final consensus merging – suitable for larger datasets.

## System Requirements

### Hardware requirements
- A standard computer with sufficient RAM to hold the working data in memory. For medium‑sized datasets (thousands of sequences), at least 8 GB RAM is recommended.

### Software requirements
#### Operating Systems
- **Linux** (Ubuntu 16.04+, CentOS 7+) and **macOS** (10.15+) are supported.
- Windows may work via WSL2.

The package has been tested on:
- Linux: Ubuntu 18.04, 20.04
- macOS: Big Sur (11.6), Monterey (12.3)

#### Dependencies
- **Python** >= 3.7
- **External software** (must be installed separately and be available in your system `PATH`):
  - **MAFFT** >= 7.450 (used for multiple sequence alignment)
- **Python packages** (to be installed via `pip`):
  - `pandas` >= 1.0.0
  - `numpy` >= 1.18.0
  - `biopython` >= 1.78

## Installation Guide

### 1. Install MAFFT
MAFFT is **required** before running PanFlu-IOA.

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install mafft

# macOS (using Homebrew)
brew install mafft

# Verify installation
mafft --version

2. Clone the repository and install Python dependencies
# Clone the repository
git clone https://github.com/[your-username]/PanFlu-IOA.git
cd PanFlu-IOA

# (Recommended) Create and activate a Python virtual environment
python3 -m venv panflu_env
source panflu_env/bin/activate  # Linux/macOS
# or using conda: conda create -n panflu python=3.8 && conda activate panflu

# Install required Python packages
pip install pandas numpy biopython

# Verify installation
python -c "import pandas, numpy, Bio; print('All dependencies installed.')"

Demo
A test_data folder is provided in the repository, containing sample FASTA sequence files for demonstration.

Run with provided test_data
Run Mode 1 (iterative optimisation)
# -t 5: 5 iterations, -p 2: use 2 parallel processes
python PanFlu-IOA_code1.py -i ./test_data -o ./result_mode1 -t 5 -p 2
Run Mode 2 (two‑stage subsampling)
# -n 30,150: 30 sequences in stage 1, 150 in stage 2
python PanFlu-IOA_code2.py -i ./test_data -o ./result_mode2 -p 2 -n 30,150
Expected output & run time
Expected output:

Intermediate files: .extract.fasta, .aln.fasta

Consensus files: .cons.fasta

Statistics files: .stat.xls (tab‑separated, with columns Cons_seq, Cons_count, Cons_ratio)

Final results are placed in the final subdirectory.

Expected run time on a normal desktop computer (4‑core CPU, using provided test_data with ~10 files × 100 sequences each):

Mode 1 (5 iterations): ~1–2 minutes

Mode 2: ~2–3 minutes
(Actual time depends on MAFFT alignment speed.)

Instructions for Use
Basic usage (Mode 1)
python PanFlu-IOA_code1.py -i <input_dir> -o <output_dir> -t <iterations> -p <parallel_jobs>

Basic usage (Mode 2)
python PanFlu-IOA_code2.py -i <input_dir> -o <output_dir> -p <parallel_jobs> -n <stage1_count,stage2_count>

Input data format
The input directory must contain FASTA files with extension .fa or .fasta.

Sequences can be nucleotide or protein.

Output files
*.cons.fasta-final consensus sequence (gap‑rich positions removed).

*.stat.xls-tab‑separated file with columns: Cons_seq (consensus residue), Cons_count (number of sequences supporting that residue), and Cons_ratio (percentage support).

License
This project is distributed under the MIT License. See the LICENSE file for details.

Contact
Hongyang Shi (13918445114@163.com)
Shulei Jia (jiashu320lei@126.com)
School of Basic Medical Sciences, Tianjin Medical University, Tianjin, 300070, China
