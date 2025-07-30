import gzip
import os
import requests
import shutil
import tempfile

from enum import StrEnum
from Bio import SeqIO, SeqRecord
from simple_term_menu import TerminalMenu
from typing import List, Tuple, Dict, Optional
from pathlib import Path
from pandas import DataFrame, read_csv

ENSEMBL_SPECIES_CONVERSION = {
    "human": {
        "dataset": "hsapiens_gene_ensembl",
        "cdna_genome_file_name": "Homo_sapiens.GRCh38.cdna.all.fa",
        "cdna_genome_file_link": "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "cds_genome_file_name": "Homo_sapiens.GRCh38.cds.all.fa",
        "cds_genome_file_link": "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz",
    },
    "mouse": {
        "dataset": "mmusculus_gene_ensembl",
        "cdna_genome_file_name": "Mus_musculus.GRCm39.cdna.all.fa",
        "cdna_genome_file_link": "https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
    },
    "rat": {
        "dataset": "rnorvegicus_gene_ensembl",
        "cdna_genome_file_name": "Rattus_norvegicus.mRatBN7.2.cdna.all.fa",
        "cdna_genome_file_link": "https://ftp.ensembl.org/pub/release-113/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz",
    },
    "monkey": {
        "dataset": "mfascicularis_gene_ensembl",
        "cdna_genome_file_name": "Macaca_fascicularis.Macaca_fascicularis_6.0.cdna.all.fa",
        "cdna_genome_file_link": "https://ftp.ensembl.org/pub/release-113/fasta/macaca_fascicularis/cdna/Macaca_fascicularis.Macaca_fascicularis_6.0.cdna.all.fa.gz",
    },
    "rabbit": {
        "dataset": "ocuniculus_gene_ensembl",
        "cdna_genome_file_name": "Oryctolagus_cuniculus.OryCun2.0.cdna.all.fa",
        "cdna_genome_file_link": "https://ftp.ensembl.org/pub/release-113/fasta/oryctolagus_cuniculus/cdna/Oryctolagus_cuniculus.OryCun2.0.cdna.all.fa.gz",
    },
}

class FastaType(StrEnum):
    CDS = "cds"
    CDNA = "cdna"

def download_file(url: str, out_path: Path) -> None:
        print(f"Downloading {url} to {out_path}...")
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(out_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)

def parse_conversion_file(db_dir: Path) -> Dict[str, Dict[str, str]]:
    """
    Parses the human conversion file and if it does not exist, downloads it first
    from: https://www.ensembl.org/biomart/martview and places it within the db_dir

    Args:
        db_dir (Path): Local path to save all downloaded and parsed databases called by ARGS.db_dir

    Returns:
        Dict[str, Dict[str, str]]: Human Ensembl transcript ID coversion to other values.
            Dict[ensemble transcript id, Dict[GeneName/EntrezID/GeneID/Canonical/GeneID/RefSeq: Value]]
    """
    conversion_data = human_conversion_df(db_dir)
    conversion_dict = {}
    for (
        ensembl_transcript,
        gene_name,
        ensembl_gene,
        canonical,
        entrez_id,
        refseq,
    ) in zip(
        conversion_data["Transcript stable ID version"],
        conversion_data["Gene name"],
        conversion_data["Gene stable ID"],
        conversion_data["Ensembl Canonical"],
        conversion_data["NCBI gene (formerly Entrezgene) ID"],
        conversion_data["RefSeq match transcript (MANE Select)"],
    ):
        conversion_dict[ensembl_transcript] = {}
        if str(gene_name) != "nan":
            conversion_dict[ensembl_transcript]["GeneName"] = gene_name
        elif str(ensembl_gene) != "nan":
            conversion_dict[ensembl_transcript]["GeneName"] = ensembl_gene
        else:
            conversion_dict[ensembl_transcript]["GeneName"] = ensembl_transcript

        if str(entrez_id) != "nan":
            conversion_dict[ensembl_transcript]["EntrezID"] = entrez_id
        else:
            conversion_dict[ensembl_transcript]["EntrezID"] = ""

        conversion_dict[ensembl_transcript]["GeneID"] = ensembl_gene
        conversion_dict[ensembl_transcript]["Canonical"] = bool(canonical)
        if str(refseq) == "nan":
            conversion_dict[ensembl_transcript]["RefSeq"] = False
        else:
            conversion_dict[ensembl_transcript]["RefSeq"] = True
    return conversion_dict

def get_target_seqs(
    conversion: Dict[str, Dict[str, str]], gene: str, db_dir: Path
) -> SeqRecord.SeqRecord:
    """
    Returns the sequence of the desired target transcript. 
    Args:
        conversion Dict[str, Dict[str, str]]: Human Ensembl transcript ID coversion to other values.
                Dict[ensemble transcript id, Dict[GeneName/EntrezID/GeneID/Canonical/GeneID/RefSeq: Value]]
        gene (gene): Desired target gene called by ARGS.gene
        db_dir (Path): Local path to save all downloaded and parsed databases called by ARGS.db_dir

    Returns:
        Bio.SeqRecord.SeqRecord: sequence CDS

    Raises:
        ValueError: Gene not within the conversion file
        ValueError: No transcripts selected
        ValueError: Transcript not found within the fasta file
    """
    transcripts = []
    for transcript in conversion.keys():
        if conversion[transcript]["GeneName"] == gene:
            refseq_text = (
                "and RefSeq" if conversion[transcript]["RefSeq"] else "but not RefSeq"
            )
            canonical_text = (
                "Canonical" if conversion[transcript]["Canonical"] else "Non-canonical"
            )
            transcripts.append(f"{canonical_text} {refseq_text}: {transcript}")
    if len(transcripts) == 0:
        raise ValueError("Gene not found in conversion file")
    transcripts.sort()
    terminal_menu = TerminalMenu(
        transcripts,
        title=f"Select main human variant of interest for {gene}",
    )

    chosen_transcript = transcripts[terminal_menu.show()].split(": ")[1]

    cds = None

    for seq_record in SeqIO.parse(
        get_genome_file("human", db_dir, fasta_type=FastaType.CDS),
        format="fasta",
    ):
        if seq_record.id == chosen_transcript:
            cds = seq_record
            break
    if cds is None:
        raise ValueError("Transcript not found in species cDNA FASTA")

    return cds

def download_ensembl(
    out_file: Path,
    dataset_name: str,
    attribute_list: List[str],
    filter_name_value: Optional[Tuple[str, str]] = None,
) -> None:
    """
    Using wget and the biomart API, downloads the desired biomart dataset to out_file.
    These are typically ortholog conversion, transcript ID conversion, etc. databases.
    Values to enter here can be retrieved from creating a query at http://www.ensembl.org/biomart/martview
    and clicking on the XML button.

    Args:
        out_file (Path): The location to which the database is downloaded.
        dataset_name (str): Biomart dataset, typically '<species>_gene_ensembl'.
        attribute_list (List[str]): The column data types for the columns within the database.
        filter_name_value (Optional[Tuple[str, str]]): The filter [Filter type, filter value] to select for the
                database. Within this algorithm, this is used to filter for Ensembl gene ID for SNPs as the file
                is too large to download for all Ensembl gene IDs and is prevented by Biomart.
    """
    attributes = ""
    if filter_name_value is not None:
        filter = f'<Filter name = "{filter_name_value[0]}" value = "{filter_name_value[1]}"/>'
    else:
        filter = ""
    for attribute in attribute_list:
        attributes = f'{attributes}<Attribute name = "{attribute}" />'
    fd, path = tempfile.mkstemp()
    path = Path(path)
    # download_template = f'wget -O {path} \'http://useast.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "{dataset_name}" interface = "default" >{filter}{attributes}</Dataset></Query>\''
    download_template = f'wget -O {path} \'http://ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "{dataset_name}" interface = "default" >{filter}{attributes}</Dataset></Query>\''
    os.system(download_template)
    with open(path, "r") as outputfile:
        for line_number, line in enumerate(outputfile):
            if "<title>Service unavailable</title>" in line:
                raise SystemError(
                    "Ensembl download currenlty unavailable.  Try again later."
                )
            if line_number > 6:
                break
    # Make sure file is fully downloaded
    path.rename(out_file)
    os.close(fd)

def human_conversion_df(db_dir: Path) -> DataFrame:
    """
    Downloads and returns the human conversion database from Biomart.  For each Ensembl transcript
    ID, contains Ensembl gene ID, Gene name, whether it is an Ensembl canonical transcript, and Entrez ID

    Args:
        db_dir (Path): Local path to save all downloaded and parsed databases called by ARGS.db_dir

    Returns:
        DataFrame: Human Ensembl transcript ID coversion to other values.
                ensemble transcript id to GeneName/EntrezID/GeneID/Canonical/GeneID/RefSeq
    """
    human_conversion_file = db_dir / "human_transcript_conversion.tsv"
    if not human_conversion_file.exists():
        print(f"Downloading human transcript conversion file")
        download_ensembl(
            out_file=human_conversion_file,
            dataset_name=ENSEMBL_SPECIES_CONVERSION["human"]["dataset"],
            attribute_list=[
                "ensembl_gene_id",
                "ensembl_transcript_id_version",
                "external_gene_name",
                "transcript_is_canonical",
                "entrezgene_id",
                "transcript_mane_select",
            ],
        )
    return read_csv(human_conversion_file, sep="\t")

def get_genome_file(
    species: str, db_dir: Path, fasta_type: FastaType = FastaType.CDNA
) -> Path:
    """
    Downloads the species fasta file from Ensembl.

    Args:
        species (str): Any species key found in ENSEMBL_SPECIES_CONVERSION
        db_dir (Path): Local path to save all downloaded and parsed databases called by ARGS.db_dir
        fasta_type (FastaType): Either 'CDNA' or 'CDS' to use for the creation of siRNAs

    Returns:
        Path: the path location of the genome fasta file use with Bio.SeqIO elsewhere in this program.
    """
    genome_file = (
        db_dir
        / ENSEMBL_SPECIES_CONVERSION[species][f"{fasta_type.value}_genome_file_name"]
    )
    if not genome_file.exists():
        genome_file_gz = Path(f"{genome_file}.gz")
        download_link = ENSEMBL_SPECIES_CONVERSION[species][
            f"{fasta_type.value}_genome_file_link"
        ]
        download_file(download_link, genome_file_gz)
        print(f"Unzipping {genome_file_gz.name}...")
        with gzip.open(genome_file_gz, "rb") as f_in:
            with open(genome_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        genome_file_gz.unlink()
    return genome_file
