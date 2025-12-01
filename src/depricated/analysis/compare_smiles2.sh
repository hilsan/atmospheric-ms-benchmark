#!/bin/bash
# compare_smiles_sanity.sh

file1="$1"
file2="$2"

if [ $# -ne 2 ]; then
  echo "Usage: $0 file1.csv file2.csv"
  exit 1
fi

# Extract SMILES columns (assuming column name contains 'SMILES')
awk -F, 'NR==1{for(i=1;i<=NF;i++) if(tolower($i) ~ /smiles/) col=i} NR>1{print $col}' "$file1" > smiles1.txt
awk -F, 'NR==1{for(i=1;i<=NF;i++) if(tolower($i) ~ /smiles/) col=i} NR>1{print $col}' "$file2" > smiles2.txt

# Paste side by side
echo "SMILES from $file1 vs $file2:"
paste smiles1.txt smiles2.txt | column -t

