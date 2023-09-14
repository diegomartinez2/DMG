#!/bin/bash

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="\$1"
    case $key in
        -in|--input)
        input_dir="\$2"
        shift
        shift
        ;;
        -out|--output)
        output_file="\$2"
        shift
        shift
        ;;
        *)
        echo "Unknown option: \$1"
        exit 1
        ;;
    esac
done

# Check if input directory and output file are provided
if [[ -z $input_dir || -z $output_file ]]; then
    echo "Usage: program -in <input directory name> -out <output compressed filename>"
    exit 1
fi

# Compress the directory using tar and gzip
tar -czvf "$output_file" "$input_dir" --exclude="$input_dir/Downloads" --exclude="$input_dir/.cache"
gzip -9 "$output_file"

echo "Compression complete!"
