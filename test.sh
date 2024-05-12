#!/bin/bash

# Navigate to the directory where the compare executable and output files are located
cd /path/to/your/src

# Define the base path for reference files
reference_path="../references"

# Array of output file names
outputs=("out1" "out2" "out3")

# Loop through each output file and compare it with the reference
for out in "${outputs[@]}"; do
    echo "[andrei.dragomir1401@fep8 src]$ ./compare $reference_path/$out $out 0.001"
    ./compare $reference_path/$out $out 0.001
done