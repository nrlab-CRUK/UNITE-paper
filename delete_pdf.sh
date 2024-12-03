#!/bin/bash

# Define the directory to start searching (use current directory by default)
DIR=${1:-$(pwd)}

echo "Searching for .pdf files in $DIR and its subdirectories..."

# Find and delete .pdf files
find "$DIR" -type f -name "*.pdf" -print -exec rm -f {} \;

echo "All .pdf files have been deleted."

