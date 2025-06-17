#!/bin/bash

# Directory containing gene lists
genelist_dir="./genelist"

# Loop through each gene list file in the directory
for gene_list in "$genelist_dir"/*; do
    # Extract the base name of the gene list file (without path and extension)
    gene_list_name=$(basename "$gene_list" .txt)

    echo "Processing gene list: $gene_list_name"
    echo "-----------------------------------"

    # Step 1: Remove carriage return characters from the gene list
    echo "Step 1: Cleaning $gene_list_name.txt by removing carriage return characters..."
    tr -d '\r' < "$gene_list" > "cleaned_${gene_list_name}.txt"
    echo "Cleaning complete. First few lines of cleaned_${gene_list_name}.txt:"
    head "cleaned_${gene_list_name}.txt"
    echo "-----------------------------------"

    # Step 2: Extract from 2000bp upstream to 1000bp downstream of mRNA from GFF3 to a BED file
    echo "Step 2: Extracting 2000bp upstream promoter regions and 1000bp downstream 3'UTR regions from 'mRNA' features..."
    awk '$3 == "mRNA"' Osativa_323_v7.0.gene.primary.transcripts.gff3 | \
    awk 'BEGIN {OFS = "\t"}
    $3 == "mRNA" {
      split($9, attributes, ";");
      for (i in attributes) {
          if (attributes[i] ~ /^ID=/) {
              gene_id = gensub(/^ID=([^\.]+).*/, "\\1", 1, attributes[i]);
              break;
          }
      }
        if ($7 == "+") {
            start_promoter=$4-2000; if (start_promoter < 0) start_promoter=0;
            end_3utr=$5+1000;
            print $1, start_promoter, $4, gene_id, ".", $7;
        } else {
            end_promoter=$5+2000;
            start_3utr=$4-1000; if (start_3utr < 0) start_3utr=0;
            print $1, $5, end_promoter, gene_id, ".", $7;
        }
    }' > gene_regions.bed
    echo "Promoter and 3'UTR region extraction complete. First few lines of gene_regions.bed:"
    head gene_regions.bed
    echo "-----------------------------------"

    # Step 3: Filter the BED file based on gene locus from the cleaned gene list
    echo "Step 3: Filtering promoter regions based on gene list..."
    grep -Ff "cleaned_${gene_list_name}.txt" gene_regions.bed > "filtered_gene_regions_${gene_list_name}.bed"
    echo "Filtering complete. First few lines of filtered_gene_regions_${gene_list_name}.bed:"
    head "filtered_gene_regions_${gene_list_name}.bed"
    echo "-----------------------------------"

    # Step 4: Intersecting ATAC peaks with 2000bp upstream promoter and 1000bp downstream 3'UTR regions
    echo "Step 4: Intersecting ATAC peaks with gene regions..."
    bedtools intersect -a "filtered_gene_regions_${gene_list_name}.bed" -b TIGR7_DHSs.bed -wa -wb | \
    awk 'BEGIN {OFS = "\t"} {
        gene_id = $4;
        peak_id = $10;
        peak_start = $8;
        peak_end = $9;
        print gene_id, $1, peak_start, peak_end, peak_id;
    }' > "intersected_gene_peaks_${gene_list_name}.txt"

    echo "Intersection complete. First few lines of intersected_gene_peaks_${gene_list_name}.txt:"
    head "intersected_gene_peaks_${gene_list_name}.txt"
    echo "-----------------------------------"

    # Step 4.1: Sliding window of 300bp with step size 100bp for peaks larger than 300bp
    echo "Step 4.1: Applying sliding window to peaks larger than 300bp..."
    while IFS=$'\t' read -r gene_id chr start end peak_name; do
        # Calculate the peak size
        peak_size=$((end - start))

        # If the peak size is less than or equal to 300, write it directly to the output file
        if [ "$peak_size" -le 300 ]; then
            echo -e "$gene_id\t$chr\t$start\t$end\t$peak_name" >> "peaks_slide_window_${gene_list_name}.bed"
        else
            # If the peak size is greater than 300, perform a sliding window
            window_size=300
            step_size=100
            window_start=$start
            window_count=1

            # Slide the window across the peak
            while [ "$window_start" -lt "$end" ]; do
                # Calculate the window end
                window_end=$((window_start + window_size))
                if [ "$window_end" -gt "$end" ]; then
                    window_end=$end
                fi

                # Create the new peak name with _1, _2, _3, etc.
                new_peak_name="${peak_name}_${window_count}"

                # Write the window to the output file
                echo -e "$gene_id\t$chr\t$window_start\t$window_end\t$new_peak_name" >> "peaks_slide_window_${gene_list_name}.bed"

                # Move the window start by the step size and increment the window count
                window_start=$((window_start + step_size))
                window_count=$((window_count + 1))

                # Stop if the next window would end beyond the peak's end
                if [ "$window_end" -ge "$end" ]; then
                    break
                fi
            done
        fi
    done < "intersected_gene_peaks_${gene_list_name}.txt"

    echo "Processing complete. Output saved to peaks_slide_window_${gene_list_name}.bed"
    echo "-----------------------------------"

    # Step 5: Fetching FASTA sequences for each peak and saving them to a single file
    echo "Step 5: Fetching FASTA sequences for each peak..."
    output_fasta="peak_sequences_${gene_list_name}.fasta"
    > "$output_fasta"  # Clear output file if it exists

    while read -r gene_info; do
        peak_name=$(echo "$gene_info" | awk '{print $5}')
        chr=$(echo "$gene_info" | awk '{print $2}')
        start=$(echo "$gene_info" | awk '{print $3}')
        end=$(echo "$gene_info" | awk '{print $4}')

        peak_file=$(mktemp)
        echo -e "$chr\t$start\t$end\t$peak_name" > "$peak_file"

        fasta_sequence=$(bedtools getfasta -fi Osativa_323_v7.0.fa -bed "$peak_file" | grep -v "^>")
        
        if [[ -n "$fasta_sequence" ]]; then
            echo ">$peak_name" >> "$output_fasta"
            echo "$fasta_sequence" >> "$output_fasta"
            echo "FASTA sequence saved for peak $peak_name"
        else
            echo "Warning: No FASTA sequence found for $peak_name ($chr:$start-$end)"
        fi

        rm "$peak_file"
    done < "peaks_slide_window_${gene_list_name}.bed"
    echo "FASTA extraction complete. Sequences saved in $output_fasta."
    echo "-----------------------------------"

    # Step 6: Running FIMO on the combined FASTA file
    echo "Step 6: Running FIMO on all peak sequences..."
    fimo --oc "fimo_output_${gene_list_name}" --thresh 1e-3 JASPR_plant_nonredundant_deduplicate.txt "$output_fasta"
    cp "fimo_output_${gene_list_name}/fimo.tsv" "fimo_results_${gene_list_name}.txt"

    echo "FIMO motif scanning completed. Results saved in fimo_results_${gene_list_name}.txt."
    echo "-----------------------------------"

    # Step 7: Append motif family information to FIMO output based on motif cluster CSV
    echo "Step 7: Adding motif family information to FIMO results..."
    awk 'BEGIN {FS=OFS="\t"} 
         FNR==NR {motif_family[$2]=$3; next} 
         {if ($1 in motif_family) print $0, motif_family[$1]; else print $0, "Unknown"}' \
    2023JASPAR_656_motif_id_name_cluster_unnest.tsv "fimo_results_${gene_list_name}.txt" \
    > "fimo_results_with_family_${gene_list_name}.txt"
    echo "Family annotation complete. First few lines of fimo_results_with_family_${gene_list_name}.txt:"
    head "fimo_results_with_family_${gene_list_name}.txt"
    echo "-----------------------------------"

    # Step 8: Filter FIMO results based on motif family and p-value thresholds
    echo "Step 8: Filtering FIMO results based on motif family and p-value thresholds..."
    awk 'BEGIN {FS=OFS="\t"} 
         FNR==1 {print $0; next} 
         {
             motif_family = $NF; 
             p_value = $8;
             if (motif_family == "WRKY" || motif_family == "G2-like") {
                 if (p_value <= 0.001) print $0
             } else if (motif_family == "MYBR_B" || motif_family == "DOF" ||  motif_family == "IDD" || motif_family == "DPBF") {
                 if (p_value <= 0.0001) print $0
             }
         }' "fimo_results_with_family_${gene_list_name}.txt" > "fimo_filtered_results_${gene_list_name}.txt"

    if [ $? -eq 0 ]; then
        echo "Filtering complete. First few lines of fimo_filtered_results_${gene_list_name}.txt:"
        head "fimo_filtered_results_${gene_list_name}.txt"
    else
        echo "Error: Filtering step encountered a problem."
        exit 1
    fi
    echo "-----------------------------------"

    # Step 9: Check if each peak contains all specified motif families
    echo "Step 9: Checking if each peak contains all specified motif families..."

    # Sort the FIMO results by sequence_name
    sorted_fimo="fimo_sorted_results_${gene_list_name}.txt"
    sort -k3,3 -k4,4n "fimo_filtered_results_${gene_list_name}.txt" > "$sorted_fimo"

    # Initialize an output file for peaks containing all motif families
    > "peaks_with_all_motifs_${gene_list_name}.txt"
    unique_promoter_ids=()

    # Initialize variables
    current_peak=""
    current_families=()

    # Read through sorted FIMO results to collect motif families per peak
    while IFS=$'\t' read -r motif_id motif_alt_id sequence_name start stop strand score p_value q_value matched_sequence family; do
        # If processing a new peak, check the families of the previous peak
        if [[ "$sequence_name" != "$current_peak" ]]; then
            # Check if the previous peak contained all required motif families
            if [[ -n "$current_peak" && "${#current_families[@]}" -gt 0 ]]; then
                if [[ " ${current_families[@]} " == *" WRKY "* && " ${current_families[@]} " == *" G2-like "* && \
                      " ${current_families[@]} " == *" MYBR_B "* && " ${current_families[@]} " == *" IDD "* && \
                      " ${current_families[@]} " == *" DOF "* && " ${current_families[@]} " == *" DPBF "* ]]; then
                    echo "$current_peak" >> "peaks_with_all_motifs_${gene_list_name}.txt"
                    unique_promoter_ids+=("$current_peak")
                fi
            fi
            # Reset for the new peak
            current_peak="$sequence_name"
            current_families=("$family")
        else
            # Add the family to the current list if not already included
            if [[ ! " ${current_families[@]} " =~ " $family " ]]; then
                current_families+=("$family")
            fi
        fi
    done < "$sorted_fimo"

    # Handle the last peak in the file
    if [[ -n "$current_peak" && "${#current_families[@]}" -gt 0 ]]; then
        if [[ " ${current_families[@]} " == *" WRKY "* && " ${current_families[@]} " == *" G2-like "* && \
              " ${current_families[@]} " == *" MYBR_B "* && " ${current_families[@]} " == *" IDD "* && \
              " ${current_families[@]} " == *" DOF "* && " ${current_families[@]} " == *" DPBF "* ]]; then
            echo "$current_peak" >> "peaks_with_all_motifs_${gene_list_name}.txt"
            unique_promoter_ids+=("$current_peak")
        fi
    fi

    # Save unique promoter IDs to a text file with gene IDs appended
    # Sort intersected_gene_peaks.txt by the 5th column (peak name)
    sort -k5,5 "intersected_gene_peaks_${gene_list_name}.txt" > "intersected_gene_peaks_sorted_${gene_list_name}.txt"

    # Sort peaks_with_all_motifs.txt by the first column (peak name)
    sort "peaks_with_all_motifs_${gene_list_name}.txt" > "peaks_with_all_motifs_sorted_${gene_list_name}.txt"
    sort -k5,5 "peaks_slide_window_${gene_list_name}.bed" > "peaks_slide_window_sorted_${gene_list_name}.bed"

    # Join the sorted files and extract gene IDs
    join -1 5 -2 1 "peaks_slide_window_sorted_${gene_list_name}.bed" "peaks_with_all_motifs_sorted_${gene_list_name}.txt" | awk '{print $0}' | sort -u > "unique_promoter_ids_${gene_list_name}.txt"

    echo "Analysis complete. Results saved to peaks_with_all_motifs_${gene_list_name}.txt and unique IDs to unique_promoter_ids_${gene_list_name}.txt."
    head "unique_promoter_ids_${gene_list_name}.txt"
    echo "-----------------------------------"
done