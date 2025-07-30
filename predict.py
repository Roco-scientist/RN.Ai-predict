import argparse
import genes
import numpy as np
import tensorflow as tf

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from features import Feature, mRNA, siRNA
from pandas import DataFrame
from pathlib import Path
from typing import List, Dict, Union


def arguments():
    """
    Parses command-line arguments for the prediction script.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.
    """
    args = argparse.ArgumentParser(description="Predicts siRNA efficacy")
    group = args.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene", dest="gene", type=str, help="Target mRNA Gene Name")
    group.add_argument("--mrna_fasta", dest="mrna_fasta", type=str, help="mRNA Fasta file")
    args.add_argument(
        "--size",
        dest="sirna_size",
        type=int,
        default=21,
        choices=(19, 20, 21),
        help="Desired siRNA size",
    )
    args.add_argument(
        "--model",
        type=str,
        default="./RN.Ai-predict.model.keras",
        help="Model to use",
    )
    args.add_argument(
        "--db_dir",
        type=str,
        default="./db",
        help="Directory to store downloaded Ensembl data.",
    )
    return args.parse_args()

def load_and_validate_fasta(fasta_path: Path) -> List[SeqRecord]:
    """
    Loads and validates a FASTA file.
    Checks for file existence, valid FASTA format, and that the file is not empty.

    Args:
        fasta_path (Path): The path to the FASTA file.
    Returns:
        List[SeqRecord]: A list of sequence records from the file.
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is not a valid FASTA format or is empty.
    """
    if not fasta_path.is_file():
        raise FileNotFoundError(f"Input file not found at {fasta_path}")

    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            raise ValueError(f"The file '{fasta_path.name}' is empty or contains no valid FASTA records.")
        return records
    except Exception:
        # Catch any other parsing errors from BioPython
        raise ValueError(f"The file '{fasta_path.name}' does not appear to be in a valid FASTA format.")

def get_predictions(x: Dict[str, np.ndarray], model_file: Union[str, Path]) -> np.ndarray:
    """
    Loads a Keras model and predicts siRNA efficacy from feature vectors.

    Args:
        x (Dict[str, np.ndarray]): A dictionary of feature matrices, where each key
            is a feature name and each value is a NumPy array of tokenized sequences.
        model_file (Union[str, Path]): The file path to the Keras model.

    Returns:
        np.ndarray: A 1D NumPy array of predicted efficacy scores (0-100), rounded.
    """
    model = tf.keras.models.load_model(model_file)
    for key in x.keys():
        x[key] = np.array(x[key])
    predictions = model.predict(x).flatten()
    return np.array([round(eff * 100, 2) for eff in predictions])


def get_nucleotide_features(mrna: mRNA, sirnas: List[siRNA]) -> Dict[str, List[List[int]]]:
    """
    Generates tokenized nucleotide features for each siRNA and its mRNA context.

    Args:
        mrna (mRNA): The parent mRNA sequence object.
        sirnas (List[siRNA]): A list of siRNA objects to be featurized.

    Returns:
        Dict[str, List[List[int]]]: A dictionary containing lists of tokenized sequences
            for the siRNA, 5' flank, and 3' flank, ready for model input.
    """
    x = {
        str(Feature.SIRNA_TOKEN): [],
        str(Feature.MRNA_5_FLANK_TOKEN): [],
        str(Feature.MRNA_3_FLANK_TOKEN): [],
    }
    for sirna in sirnas:
        mrna_parent = mrna
        sirna_token = sirna.nucleotide_tokenization()
        five_mrna_token, three_mrna_token = mrna_parent.nucleotide_tokenization(
            sirna, split=True
        )
        x[str(Feature.SIRNA_TOKEN)].append(sirna_token)
        x[str(Feature.MRNA_5_FLANK_TOKEN)].append(five_mrna_token)
        x[str(Feature.MRNA_3_FLANK_TOKEN)].append(three_mrna_token)
    return x


def split_sirnas(mrna: mRNA, size: int) -> List[siRNA]:
    """
    Generates all possible siRNA candidates of a given size from an mRNA sequence.

    This function slides a window across the mRNA to define target sites,
    ensuring a 19-nt flanking region is available on both sides.

    Args:
        mrna (mRNA): The mRNA sequence object to process.
        size (int): The desired length of the siRNAs (e.g., 19, 21).

    Returns:
        List[siRNA]: A list of siRNA objects, each representing a potential target.
    """
    sirnas = []
    if mrna.id:
        mrna_label = f"{mrna.id}_"
    else:
        mrna_label = ""

    # We stop 19 nt before the end to ensure a full 19 nt 3' flank is available for the last siRNA
    for i in range(len(mrna.seq) - size - 38):
        # The target site on the mRNA starts after the 19-nt 5' flank.
        sirna_antisense_seq = mrna.seq[i + 19 : i + 19 + size].reverse_complement_rna()
        sirna_id = f"{mrna_label}siRNA_{i+1}"
        sirna_antisense_seq_record = SeqRecord(sirna_antisense_seq, id=sirna_id)
        sirnas.append(siRNA(sirna_antisense_seq_record))
    return sirnas


def main() -> None:
    if ARGS.mrna_fasta:
        try:
            fasta_file = Path(ARGS.mrna_fasta)
            mrna_records = load_and_validate_fasta(fasta_file)
            mrnas = [mRNA(seq_record) for seq_record in mrna_records]
            output_path = fasta_file.with_suffix(".sirna_prediction.csv") 
        except (FileNotFoundError, ValueError) as e:
            print(f"Error: {e}")
            exit(1)
    elif ARGS.gene:
        db_path = Path(ARGS.db_dir)
        db_path.mkdir(exist_ok=True)
        genes.get_genome_file("human", db_path, fasta_type=genes.FastaType.CDS)
        conversion = genes.parse_conversion_file(db_path)
        gene = ARGS.gene.upper()
        target_seq = genes.get_target_seqs(conversion, gene, db_path)
        mrnas = [mRNA(target_seq)]
        output_path = Path(f"{ARGS.gene.upper()}.sirna_prediction.csv")
    results = {
        "siRNA_ID": [],
        "siRNA_Sense": [],
        "siRNA_Antisense": [],
        "Efficacy_Prediction": [],
        "Rank_Order": [],
    }
    for mrna in mrnas:
        new_sirnas = split_sirnas(mrna, ARGS.sirna_size)
        x = get_nucleotide_features(mrna, new_sirnas)
        efficacies = get_predictions(x, ARGS.model)
        rank_order = np.argsort(np.argsort(efficacies)[::-1])
        results["siRNA_ID"].extend([sirna.id for sirna in new_sirnas])
        results["siRNA_Sense"].extend(
            [sirna.seq.reverse_complement_rna() for sirna in new_sirnas]
        )
        results["siRNA_Antisense"].extend([sirna.seq for sirna in new_sirnas])
        results["Efficacy_Prediction"].extend(efficacies)
        results["Rank_Order"].extend(rank_order + 1)
    DataFrame(results).to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")


if __name__ == "__main__":
    ARGS = arguments()
    main()
