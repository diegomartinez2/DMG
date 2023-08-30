#!/bin/bash

# Function to monitor file changes
monitor_file() {
  tail -f "\$1" | dialog --title "File Changes" --textbox - 0 0
}

# Function to display file selection dialog
select_file() {
  files=("file1.txt" "file2.txt" "file3.txt")
  selected_file=$(dialog --title "Select File" --menu "Choose a file to monitor:" 0 0 0 "${files[@]}" 2>&1 >/dev/tty)

  # Prompt for pattern to match
  pattern=$(dialog --title "Enter Pattern" --inputbox "Enter a pattern to match in the file:" 0 0 2>&1 >/dev/tty)

  # Monitor file and display only matching lines
  monitor_file "$selected_file" | grep "$pattern" | dialog --title "Matching Lines" --textbox - 0 0
}

# Call the file selection function
select_file
