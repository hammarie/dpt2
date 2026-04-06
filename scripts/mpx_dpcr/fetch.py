# src/mpx_dpcr/fetch.py
from Bio import Entrez, SeqIO

# Required by NCBI: identify yourself
Entrez.email = "your.email@example.com"


def _fetch_fasta_by_id(accession_or_gi: str) -> str:
    """Fetch a FASTA record directly by accession or GI."""
    with Entrez.efetch(db="nucleotide", id=accession_or_gi, rettype="fasta", retmode="text") as handle:
        # Use fasta-pearson to tolerate leading comment lines without warnings.
        seq_record = SeqIO.read(handle, "fasta-pearson")
    return str(seq_record.seq).upper()


def fetch_gene_sequence(identifier: str, organism="Homo sapiens") -> str:
    """
    Fetch the nucleotide sequence (RefSeq mRNA) for a given gene symbol or accession.

    Returns a plain DNA string (ATCG...) by first attempting to read the identifier
    directly (treating it as an accession), then falling back to a gene-symbol search.
    """
    # Try treating the identifier as an accession/GI first
    try:
        return _fetch_fasta_by_id(identifier)
    except Exception:
        pass  # Fall back to gene-symbol search

    search_query = f"{identifier}[Gene] AND {organism}[Organism] AND mRNA[Filter]"
    with Entrez.esearch(db="nucleotide", term=search_query, retmax=1) as handle:
        record = Entrez.read(handle)

    if not record["IdList"]:
        raise ValueError(f"No mRNA found for {identifier}")

    seq_id = record["IdList"][0]
    return _fetch_fasta_by_id(seq_id)
