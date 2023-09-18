#!/bin/bash

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="\$1"
    case $key in
        -t|--time)
            time="\$2"
            shift
            shift
            ;;
        *)
            echo "Unknown option: \$1"
            exit 1
            ;;
    esac
done

# Validate input
if [[ -z $time ]]; then
    echo "Please provide the number of seconds using --time option."
    exit 1
fi

# Countdown function
countdown() {
    secs=\$1
    while [[ $secs -gt 0 ]]; do
        echo "Waiting $secs seconds..."
        sleep 1
        ((secs--))
    done
    echo "Countdown completed!"
}

# Display countdown using Dialog
dialog --title "Countdown" --infobox "Countdown started for $time seconds." 5 40
countdown "$time" | dialog --title "Countdown" --gauge "Countdown in progress..." 7 50
