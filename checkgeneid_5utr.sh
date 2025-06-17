#!/bin/bash

# Step 1: Extract all gene IDs from the GFF3 file
awk '$3 == "mRNA" {
    split($9, info, ";");
    split(info[1], id, "=");
    split(id[2], gene, ".");
    print gene[1];
}' Osativa_323_v7.0.gene.primary.transcripts.gff3 | sort | uniq > all_gene_ids.txt

# Step 2: Extract all gene IDs that have five_prime_UTR
awk '$3 == "five_prime_UTR" {
    split($9, info, ";");
    split(info[1], id, "=");
    split(id[2], gene, ".");
    print gene[1];
}' Osativa_323_v7.0.gene.primary.transcripts.gff3 | sort | uniq > genes_with_5prime_utr.txt

# Step 3: Compare the two lists to find genes without five_prime_UTR
comm -23 all_gene_ids.txt genes_with_5prime_utr.txt > genes_without_5prime_utr.txt

echo "Process completed. The list of genes without five_prime_UTR can be found in genes_without_5prime_utr.txt"

