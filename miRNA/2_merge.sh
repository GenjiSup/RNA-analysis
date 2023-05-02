#!/bin/bash

# set the directories to search for miR counts files
counts_dirs="/share/analysis/Carlo/miRNA/AZA:/share/analysis/Carlo/miRNA/DMSO:/share/analysis/Carlo/miRNA/CYC"

# set the output file name
output_file="merged_counts.csv"

# create an empty file for the merged counts
> $output_file

# loop through each directory and find all counts files
for counts_dir in ${counts_dirs//:/ }; do
    counts_files=($(find "$counts_dir" -name "miR.Counts.csv"))
    
    # loop through each counts file and append its contents to the merged counts file
    for counts_file in "${counts_files[@]}"; do
        echo "Merging counts from $counts_file"
        if [[ ! -s $output_file ]]; then
            # if the merged counts file is empty, copy the first file's contents to the merged counts file
            cut -d , -f 1,2 "$counts_file" > "$output_file"
        else
            # otherwise, join the new file with the existing merged counts file on the first column
            join -t , -a 1 -e '0' "$output_file" <(cut -d , -f 1,2 "$counts_file") > tmp_file
            mv tmp_file "$output_file"
        fi
    done
done

# delete the second row of the merged counts file
sed -i '2d' $output_file

echo "Merged counts file created: $output_file"
