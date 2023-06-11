# Gene Promoter Sequence Retrieval

## Description
This Python script retrieves the promoter sequence (1,000 base pairs upstream of the transcription start site) for a specified list of genes. The retrieved sequences can be from human, mouse, fly or zebrafish genomes. It outputs the promoter sequences in FASTA format, a CSV file with detailed information about each gene, and a SAM file with the mapping of each promoter sequence to the respective genome.

## Usage
```bash
python get_promoter.py -s {species} (-g {gene_names} | -f {filename})
```

## Parameters
- `-s, --species`: Species genome to use. Options are "human", "mouse", "fly", "zebrafish".
- `-g, --genes`: Comma-separated list of gene names.
- `-f, --file`: Filename of a text file containing gene names, one per line.

You should use either `-g` or `-f` option to input gene names, but not both.

## Outputs
- `output.fasta`: A FASTA file containing promoter sequences for each input gene.
- `gene_info.csv`: A CSV file containing detailed information about each gene, including the promoter and gene start and end positions, and the strand.
- `output.sam`: A SAM file containing the mapping of each promoter sequence to the genome.

## Installation
This script requires Python 3 and the following Python libraries:
- `requests`
- `csv`
- `argparse`

You can install these libraries using pip:

```bash
pip install requests csv argparse
```

## Example
```bash
python get_promoter.py -s human -g BRCA2,P53
```

This command retrieves the promoter sequences for human BRCA2 and P53 genes.

---

Remember to replace "get_promoter.py" with the actual filename of your script.
