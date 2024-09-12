#!/bin/bash

# Define file paths
test_set="sentinels_194_sorted.bed"  # Replace with your actual test set file name
background_dir="<YOUR_DIRECTORY>/backg_bedfiles"  # Directory containing the 1000 background files
features_file="combined_roadmap39celltypes_dhs_histonemarks_sorted.bed"  # File with chr, start, end, cell type, and feature
output_file="enriched_features.txt"

# Initialize the output file with headers
echo -e "test_set_overlap_count\taverage_background_overlap_count\tcell_type\tfeature" > "$output_file"

# Create a temporary file to store background overlap counts
temp_bg_counts=$(mktemp)

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

  # Append the results to the output file
  echo -e "$test_overlap_count\t$avg_bg_overlap_count\t$cell_type\t$feature" >> "$output_file"
done

# Clean up temporary files
rm temp_feature.bed "$temp_bg_counts"

echo "Overlap counts have been calculated and saved to $output_file"
