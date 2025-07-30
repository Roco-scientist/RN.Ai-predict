# RN.Ai-Predict

A command-line tool to predict the efficacy of small interfering RNAs (siRNAs). This tool is the official implementation of the model described in the paper:

> Coffey, R. (2024). *Systematic feature and architecture evaluation reveals tokenized learned embeddings enhance siRNA efficacy prediction*. (Journal pending publication).

RN.Ai-Predict can generate and evaluate all possible siRNAs for a target in two ways:
1.  **From a FASTA file:** Provide your own mRNA sequence.
2.  **By gene name:** Provide a gene name (e.g., "HIF1A"), and the tool will automatically download the necessary transcript data from Ensembl.

---

## Overview

This tool takes either an mRNA sequence or a gene name and generates a list of all possible siRNAs of a specified length (19, 20, or 21 nt). For each potential siRNA, it predicts the silencing efficacy and provides a percentile rank relative to all other siRNAs for that gene.

The core of the prediction is a model trained on curated public datasets, using a gene-based cross-validation strategy to ensure robust performance on unseen gene targets.

---

## Installation

### 1. Prerequisites

- Python 3.12+ (not tested on older versions)

### 2. Clone the Repository

```bash
git clone https://github.com/Roco-scientist/RN.Ai-Predict.git
cd RN.Ai-Predict
```

### 3. Create a Virtual Environment and Install Dependencies

Using a virtual environment is highly recommended.

```bash
# Create and activate the virtual environment
python3.12 -m venv venv
source venv/bin/activate

# Install required packages
pip install -r requirements.txt
```

---

## Usage

The script `predict.py` requires exactly one of two input modes: predicting from a local FASTA file (`--mrna_fasta`) or predicting directly from a gene name (`--gene`).

### Command-Line Arguments

*   **Required (choose one):**
    *   `--mrna_fasta <path>`: Path to the input FASTA file containing one or more mRNA sequences.
    *   `--gene <GENE_NAME>`: The official gene symbol (e.g., `HIF1A`, `GAPDH`) for which to predict siRNAs. Currently supports human genes.

*   **Optional:**
    *   `--size <int>`: The length of the siRNAs to generate. Choices are `19`, `20`, or `21`. Defaults to `21`.
    *   `--model <path>`: Path to the pre-trained `.keras` model file. Defaults to `./RN.Ai-predict.model.keras`.
    *   `--db_dir <path>`: Directory to store large, downloaded Ensembl data files. Defaults to `./db`.

### Mode 1: Predict from a FASTA file

This mode is useful when you have a specific mRNA sequence you want to analyze.

**Example:**
```bash
python predict.py --mrna_fasta ./example_data/hif1a.fasta --size 21
```

### Mode 2: Predict from a Gene Name (Human)

This is the easiest way to analyze a human gene. It automatically fetches the latest transcript data from Ensembl.

**First-Time Setup:**
The first time you run a prediction using `--gene`, the script will download the necessary genome (cDNA) and annotation files from Ensembl.
- By default, these files will be stored in a `./db` directory. You can change this location with the `--db_dir` argument.
- This download can take several minutes and may require a few gigabytes of disk space. This is a one-time process for each species/dataset.

**Interactive Transcript Selection:**
After the data is available, you will be presented with an interactive menu in your terminal. This allows you to select the specific transcript(s) you want to target for siRNA design. Use the arrow keys and Enter/Space to make your selections.

**Example:**
```bash
# Predict for the HIF1A gene using default settings
python predict.py --gene HIF1A

# Predict for the GAPDH gene and store Ensembl data in a different directory
python predict.py --gene GAPDH --db_dir /path/to/my/databases/
```

---

## Output

The script generates a CSV file with prediction results.
- For FASTA input, the output is `<your_fasta_filename>.sirna_prediction.csv`.
- For gene input, the output is `<GENE_NAME>.sirna_prediction.csv`.

The CSV file contains the following columns:

| Column Name           | Description                                                                 |
|-----------------------|-----------------------------------------------------------------------------|
| `siRNA_ID`            | A unique identifier for the siRNA, based on its position in the mRNA.       |
| `siRNA_Sense`         | The sequence of the siRNA sense strand.                                     |
| `siRNA_Antisense`     | The sequence of the siRNA antisense strand (complementary to the mRNA).     |
| `Efficacy_Prediction` | The predicted knockdown efficacy percentage (0-100).                        |
| `Rank_Order`          | The rank of the siRNA's efficacy among all siRNAs for that gene.            |

---

## Citation

If you use this tool or the RN.Ai-Predict model in your research, please cite my paper:

```bibtex
@article{Coffey2024,
  author    = {Coffey, Rory},
  title     = {Systematic feature and architecture evaluation reveals tokenized learned embeddings enhance siRNA efficacy prediction},
  journal   = {BioRxiv},
  year      = {2025},
  note      = {(pending publication)}
}
```

---

## License

This project is licensed under the terms of the [LICENSE](LICENSE) file.
