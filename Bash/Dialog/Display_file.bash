#!/bin/bash

# Function to monitor file changes
monitor_file() {
  tail -f "\$1" | dialog --title "File Changes" --textbox - 0 0
}

# Check if file argument is provided
if [ "\$1" ]; then
  # Monitor file and display changes
  monitor_file "\$1"
else
  # Display file selection dialog if no argument is provided
  files=("file1.txt" "file2.txt" "file3.txt")
  selected_file=$(dialog --title "Select File" --menu "Choose a file to monitor:" 0 0 0 "${files[@]}" 2>&1 >/dev/tty)
  monitor_file "$selected_file"
fi
