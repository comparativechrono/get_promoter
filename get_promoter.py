#!/usr/bin/env python
# Requires: requests, csv

import argparse
import csv
import requests

SPECIES_DICT = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
    "fly": "drosophila_melanogaster",
    "zebrafish": "danio_rerio"
}

def reverse_complement(seq):
    complement_table = str.maketrans('ACGT', 'TGCA')
    return seq.translate(complement_table)[::-1]

def get_promoter_sequence(gene_names, species_common_name):
    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/{species}/{gene_name}?content-type=application/json"
    species = SPECIES_DICT[species_common_name]

    sequences = []

    with open('gene_info.csv', 'w', newline='') as csvfile:
        fieldnames = ['gene_name', 'promoter_start', 'promoter_end', 'gene_start', 'gene_end', 'strand']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        with open('output.sam', 'w', newline='') as samfile:
            fieldnames = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
            sam_writer = csv.writer(samfile, delimiter='\t')
            sam_writer.writerow(fieldnames)

            for gene_name in gene_names:
                r = requests.get(server + ext.format(gene_name=gene_name, species=species), headers={"Content-Type": "application/json"})
                if not r.ok:
                    r.raise_for_status()

                decoded = r.json()
                if decoded['strand'] == 1:
                    promoter_start = decoded['start'] - 1000
                    promoter_end = decoded['start'] - 1
                else:
                    promoter_start = decoded['end'] + 1
                    promoter_end = decoded['end'] + 1000

                seq_ext = "/sequence/region/{species}/{seq_region}:{start}:{end}:1?content-type=text/x-fasta".format(
                    species=decoded['species'],
                    seq_region=decoded['seq_region_name'],
                    start=promoter_start,
                    end=promoter_end)

                seq_r = requests.get(server + seq_ext, headers={"Content-Type": "text/x-fasta"})
                if not seq_r.ok:
                    seq_r.raise_for_status()

                sequence_header, sequence_forward = seq_r.text.split('\n', 1)
                sequence_forward = sequence_forward.replace('\n', '')
                if decoded['strand'] == -1:
                    sequence = reverse_complement(sequence_forward)
                else:
                    sequence = sequence_forward

                merged_header = '>{} | {}'.format(gene_name, sequence_header[1:])
                sequences.append('{}\n{}'.format(merged_header, sequence))

                writer.writerow({'gene_name': gene_name, 'promoter_start': promoter_start, 'promoter_end': promoter_end,
                                 'gene_start': decoded['start'], 'gene_end': decoded['end'], 'strand': decoded['strand']})

                if decoded['strand'] == 1:
                    start = decoded['start'] - 1000
                else:
                    start = decoded['end'] + 1

                sam_writer.writerow([merged_header, '0', decoded['seq_region_name'], start, '255', '1000M', '*', '0', '0', sequence_forward, '*'])

    with open('output.fasta', 'w') as f:
        f.write('\n'.join(sequences))

def main():
    parser = argparse.ArgumentParser(description='Get promoter sequences for a list of genes.')
    parser.add_argument('-s', '--species', required=True, choices=SPECIES_DICT.keys(), help='Species genome to use.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-g', '--genes', help='comma-separated list of gene names')
    group.add_argument('-f', '--file', help='file containing a list of gene names, one per line')

    args = parser.parse_args()

    if args.genes:
        gene_names = args.genes.split(',')
    elif args.file:
        with open(args.file, 'r') as f:
            gene_names = [line.strip() for line in f]

    get_promoter_sequence(gene_names, args.species)

if __name__ == "__main__":
    main()
