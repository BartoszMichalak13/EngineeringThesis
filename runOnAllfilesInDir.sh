#!/bin/bash

# Ensure the script is passed a directory and a number as arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: bash $0 <ALGORITHMS(number)> <directory to read from>"
    echo "Notes: Optionally one can use sudo; Use this script from its directory"
    exit 1
fi

# number representing algorithms to run
ALGORITHMS=$1

# Directory containing files and the number to pass
DIR=$2

# Check if the provided directory is valid
if [ ! -d "$DIR" ]; then
    echo "Error: $DIR is not a valid directory."
    exit 1
fi

# Path to the GraphGenerator program
GRAPH_GENERATOR="./build/GraphGenerator"

# Check if GraphGenerator exists and is executable
if [ ! -x "$GRAPH_GENERATOR" ]; then
    echo "Error: $GRAPH_GENERATOR does not exist or is not executable."
    exit 1
fi

# Iterate through each file in the directory
for FILE in "$DIR"/*; do
    if [ -f "$FILE" ]; then  # Ensure it's a file (not a subdirectory, etc.)
        echo "Processing file: $FILE with ALGORITHMS: $ALGORITHMS"
        echo "$GRAPH_GENERATOR" 0 "$ALGORITHMS" "$FILE" "$RESULTS" 0
        "$GRAPH_GENERATOR" 0 "$ALGORITHMS" "$FILE" "/results/$FILE.results" 0
    fi
done