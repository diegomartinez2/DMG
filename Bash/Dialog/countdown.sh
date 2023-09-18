#!/bin/bash

# Parse command-line arguments
# while [[ $# -gt 0 ]]; do
#     key="\$1"
#     case $key in
#         --time)
#             time="\$2"
#             shift
#             shift
#             ;;
#         *)
#             echo "Unknown option: \$1"
#             exit 1
#             ;;
#     esac
# done

# Validate input
# if [[ -z $time ]]; then
#     echo "Please provide the number of seconds using --time option."
#     exit 1
# fi

# Countdown function


# Function to display countdown using 'Dialog'
countdown_dialog() {
  local seconds=$1
  dialog --infobox "Countdown: $seconds seconds" 0 0
  sleep 1
}

# Function to perform the countdown
perform_countdown() {
  local seconds=$1

  while [ $seconds -ge 0 ]; do
    countdown_dialog $seconds
    seconds=$((seconds - 1))
  done

#  dialog --msgbox "Countdown complete!" 0 0
}

# Main script
clear
echo "Enter the number of seconds to countdown:"
read seconds

perform_countdown "$seconds"
echo $seconds
