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
#tar -cvf "$output_file" "$input_dir" --exclude="$input_dir/Downloads" --exclude="$input_dir/.cache"
#gzip -9 "$output_file"
#tar -c "$input_dir" | gzip --best -c > "$output_file"
#tar cf - "$input_dir" -P | pv -s $(du -sb "$input_dir" | awk '{print $1}') | gzip --best > "$output_file"
#tar cf - /folder-with-big-files -P | pv -s $(du -sb /folder-with-big-files | awk '{print $1}') | gzip > big-files.tar.gz
(tar cf - "$input_dir" | pv -n -s 'du -sb "$input_dir" | awk '{print $1}'' | gzip -9 > "$output_file") 2>&1 | dialog --gauge 'Progress' 7 70
#(tar cf - . | pv -n -s 'du -sb . | awk '{print $1}'' | gzip -9 > out.tgz) 2>&1 | dialog --gauge 'Progress' 7 70
echo "Compression complete!"
