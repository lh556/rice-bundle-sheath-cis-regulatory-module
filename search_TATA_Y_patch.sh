#!/bin/bash

# Step 1: Clean up the Windows carriage returns from gene_list.txt
tr -d '\r' < gene_list.txt > cleaned_gene_list.txt

# Step 2: Extract all five_prime_UTR and mRNA information to separate BED files
awk '$3 == "five_prime_UTR" { 
    split($9, info, ";");
    split(info[1], id, "=");
    split(id[2], gene, ".");
    gene_id = gene[1];
    print $1 "\t" $4-1 "\t" $5 "\t" gene_id "\t.\t" $7; 
}' Osativa_323_v7.0.gene.primary.transcripts.gff3 > five_prime_UTR.bed

awk '$3 == "mRNA" { 
    split($9, info, ";");
    split(info[1], id, "=");
    split(id[2], gene, ".");
    gene_id = gene[1];
    print $1 "\t" $4-1 "\t" $5 "\t" gene_id "\t.\t" $7; 
}' Osativa_323_v7.0.gene.primary.transcripts.gff3 > mrna.bed

# Step 3: Sort the BED files
sort -k1,1 -k4,4 -k2,2n five_prime_UTR.bed > sorted_five_prime_UTR.bed
sort -k1,1 -k4,4 -k2,2n mrna.bed > sorted_mrna.bed

# Step 4: Merge five_prime_UTR coordinates for the same gene (take start of first occurrence, end of last, and track strand)
awk '{
    if ($4 == prev_gene && $1 == prev_chr) {
        prev_end = $3;  # Update the end position for the gene
    } else {
        if (NR > 1) {
            print prev_chr "\t" prev_start "\t" prev_end "\t" prev_gene "\t.\t" prev_strand;
        }
        prev_chr = $1; prev_start = $2; prev_end = $3; prev_gene = $4; prev_strand = $6;
    }
} END { print prev_chr "\t" prev_start "\t" prev_end "\t" prev_gene "\t.\t" prev_strand; }' sorted_five_prime_UTR.bed > merged_five_prime_UTR.bed

# Step 5: Filter gene list to those with five_prime_UTR
grep -Ff cleaned_gene_list.txt merged_five_prime_UTR.bed > filtered_five_prime_UTR.bed

# Step 6: Extract only gene IDs from filtered_five_prime_UTR.bed and find genes from the cleaned gene list that do not have a five_prime_UTR
cut -f4 filtered_five_prime_UTR.bed > genes_with_5prime_UTR.txt

# Now use grep to find genes in the cleaned_gene_list.txt that are NOT in the genes_with_5prime_UTR.txt
grep -Fwvf genes_with_5prime_UTR.txt cleaned_gene_list.txt > genes_without_5prime_UTR.txt

# Step 7: Get mRNA coordinates only for genes that do not have five_prime_UTR
grep -Ff genes_without_5prime_UTR.txt sorted_mrna.bed > genes_without_5prime_UTR.bed

# Step 8: Adjust start/end for genes with five_prime_UTR
awk '{
    if ($6 == "+") {
        new_start = ($2 - 200 >= 0) ? $2 - 200 : 0;  # Shift upstream for +
        new_end = $3;  # Keep end unchanged for +
    } else {
        new_start = $2;  # Keep start unchanged for -
        new_end = $3 + 200;  # Shift downstream for -
    }
    print $1 "\t" new_start "\t" new_end "\t" $4 "\t" $5 "\t" $6;
}' filtered_five_prime_UTR.bed > adjusted_filtered_five_prime_UTR.bed

# Step 9: Adjust start/end for genes without five_prime_UTR using mRNA coordinates
awk '{
    if ($6 == "+") {
        new_start = ($2 - 200 >= 0) ? $2 - 200 : 0;  # Shift upstream for +
        new_end = $2;  # Keep end unchanged for +
    } else {
        new_start = $3;  # Keep start unchanged for -
        new_end = $3 + 200;  # Shift downstream for -
    }
    print $1 "\t" new_start "\t" new_end "\t" $4 "\t" $5 "\t" $6;
}' genes_without_5prime_UTR.bed > adjusted_mrna.bed

# Step 10: Combine the adjusted five_prime_UTR and mRNA coordinates
cat adjusted_filtered_five_prime_UTR.bed adjusted_mrna.bed > combined_adjusted_bed.bed

# Step 11: Fetch the FASTA sequence based on the adjusted BED file with gene ID as the name
bedtools getfasta -fi Osativa_323_v7.0.fa -bed combined_adjusted_bed.bed -fo combined_adjusted_sequences.fa -s -name

echo "Process completed. The adjusted FASTA sequences are in combined_adjusted_sequences.fa"

# Step 12: Run FIMO on the FASTA sequences
fimo --oc fimo_output --thresh 1E-4 JASPR_PolII_with_Y_patch.meme combined_adjusted_sequences.fa

# Step 13: Extract gene IDs for only the Y-patch motif
# Assuming 'Y-patch' is the name of the motif in the FIMO output
awk 'NR > 1 && $2 == "Y-patch" { 
    gsub(/\(\+\)/, "", $3);  # Remove (+)
    gsub(/\(\-\)/, "", $3);  # Remove (-)
    print $3; 
}' fimo_output/fimo.tsv | sort | uniq > genes_with_Y_patch_motif.txt

echo "Gene list containing Y-patch motif saved to genes_with_Y_patch_motif.txt"
