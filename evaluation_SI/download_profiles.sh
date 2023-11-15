#!/bin/bash

# List of numbers (xxx) for file names
numbers=(83 90 91 92 93 95 102 111 116)  # Add your desired numbers

# Base URL of the HTTPS server
base_url="https://scienceweb.whoi.edu/itp/data/"

# Directory where you want to save the downloaded files
download_directory="/p/project/chhb19/mueller29/mat_files/ITP_WHOI/"

# Create the download directory if it doesn't exist
mkdir -p "$download_directory"

# Loop through the list of numbers and download the files
for number in "${numbers[@]}"; do
    filename="itp${number}grddata.zip"
    file_url="${base_url}${filename}"
    
    # Use curl to download the file
    curl -o "${download_directory}${filename}" "$file_url"
    
    # Check the exit status of curl and print a message accordingly
    if [ $? -eq 0 ]; then
        echo "Downloaded: $filename"
    else
        echo "Failed to download: $filename"
    fi
done



