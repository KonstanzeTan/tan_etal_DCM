#!/bin/bash

# Define file paths
test_set="sentinels_194_sorted.bed"
background_dir="/scratch/users/ntu/konstan7/20240820_enr_genreg/backg_bedfiles"
features_file="combined_roadmap39celltypes_dhs_histonemarks_sorted.bed"
output_file="enriched_genefeatures.txt"

# Initialize the output file with headers
echo -e "test_set_overlap_count\taverage_background_overlap_count\tcell_type\tfeature\texceeding_background_count\ttest_type" > "$output_file"

# Create a temporary file to store background overlap counts
temp_bg_counts=$(mktemp)

# Function to count exceeding background sets
count_exceeding_backgrounds() {
    local test_count=$1
    local bg_counts_file=$2
    local test_type=$3
    
    if [ "$test_type" = "enrichment" ]; then
        # Count how many background values are greater than the test value
        local exceeding_count=$(awk -v test="$test_count" '$1 > test' "$bg_counts_file" | wc -l)
    else
        # Count how many background values are less than the test value
        local exceeding_count=$(awk -v test="$test_count" '$1 < test' "$bg_counts_file" | wc -l)
    fi
    
    echo "$exceeding_count"
}

# Iterate over each unique cell type and feature in the features file
awk '{print $4 "\t" $5}' "$features_file" | sort | uniq | while read -r cell_type feature; do
    # Filter the features file for the current cell type and feature
    awk -v ct="$cell_type" -v ft="$feature" '$4 == ct && $5 == ft' "$features_file" > temp_feature.bed

    # Count overlaps with the test set
    test_overlap_count=$(bedtools intersect -a "$test_set" -b temp_feature.bed -c -sorted | awk '{sum+=$NF} END {print sum}')

    # Reset temporary background counts file
    > "$temp_bg_counts"

    # Iterate over each background file
    for bg_file in "$background_dir"/*.bed; do
        # Count overlaps with the current background file
        bg_overlap_count=$(bedtools intersect -a "$bg_file" -b temp_feature.bed -c -sorted | awk '{sum+=$NF} END {print sum}')
        echo "$bg_overlap_count" >> "$temp_bg_counts"
    done

    # Calculate the average overlap count for the background files
    avg_bg_overlap_count=$(awk '{sum+=$1} END {print sum/NR}' "$temp_bg_counts")

    # Determine test type (enrichment or depletion)
    if (( $(echo "$test_overlap_count > $avg_bg_overlap_count" | bc -l) )); then
        test_type="enrichment"
    else
        test_type="depletion"
    fi

    # Count exceeding background sets
    exceeding_count=$(count_exceeding_backgrounds "$test_overlap_count" "$temp_bg_counts" "$test_type")

    # Append the results to the output file
    echo -e "$test_overlap_count\t$avg_bg_overlap_count\t$cell_type\t$feature\t$exceeding_count\t$test_type" >> "$output_file"
done

# Clean up temporary files
rm temp_feature.bed "$temp_bg_counts"

echo "Overlap counts and exceeding background counts have been calculated and saved to $output_file"
