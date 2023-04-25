#!/bin/bash

# set the directories to search for miR counts files
counts_dirs="/share/analysis/Carlo/miRNA/AZA:/share/analysis/Carlo/miRNA/DMSO:/share/analysis/Carlo/miRNA/CYC"

# set the output file name
output_file="merged_counts.csv"

# create an empty file for the merged counts
touch $output_file

# loop through each directory and find all counts files
    counts_files=$(find $counts_dir -name "miR.Counts.csv")
    
    # initialize a flag variable to check if it's the first file being merged
    first_file=true

    # loop through each counts file and append its contents to the merged counts file
    for counts_file in $counts_files; do
        echo "Merging counts from $counts_file"
        if [ "$first_file" = true ]; then
            # remove the suffix "_R1.fastq" from the headers
            sed -i '1s/_R1.fastq//' $counts_file
            # copy the first file's contents to the merged counts file, keeping only the first column
            cut -d , -f 1 $counts_file > $output_file
            # append the second column to the merged counts file
            paste $output_file <(cut -d , -f 2 $counts_file) > tmp_file
            mv tmp_file $output_file
            # set the flag variable to false, indicating that the first file has been merged
            first_file=false
        else
            # remove the suffix "_R1.fastq" from the headers
            sed -i '1s/_R1.fastq//' $counts_file
            # remove the suffix "_R1.fastq" from the headers
            sed -i '1s/.fastq//' $counts_file
            # append the second column to the merged counts file
            paste $output_file <(cut -d , -f 2 $counts_file) > tmp_file
            mv tmp_file $output_file
        fi
    done

# delete the second row of the merged counts file
sed -i '2d' $output_file

echo "Merged counts file created: $output_file"
