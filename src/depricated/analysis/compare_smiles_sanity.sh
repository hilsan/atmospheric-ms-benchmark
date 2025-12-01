#!/bin/bash
# compare_smiles.sh
# Usage: ./compare_smiles.sh fun_groups.csv SMILES_arom.csv

file1="$1"
file2="$2"

if [[ -z "$file1" || -z "$file2" ]]; then
  echo "Usage: $0 fun_groups.csv SMILES_arom.csv"
  exit 1
fi

# Extract the SMILES column from each file
# Assumes the column is literally named "SMILES"
col1=$(head -1 "$file1" | tr ',' '\n' | nl -v0 | grep -w "SMILES" | awk '{print $1}')
col2=$(head -1 "$file2" | tr ',' '\n' | nl -v0 | grep -w "SMILES" | awk '{print $1}')

if [[ -z "$col1" || -z "$col2" ]]; then
  echo "Error: Could not find 'SMILES' column in one of the files."
  exit 1
fi

# Extract just the SMILES columns (skip header)
cut -d',' -f$((col1+1)) "$file1" | tail -n +2 > smiles1.txt
cut -d',' -f$((col2+1)) "$file2" | tail -n +2 > smiles2.txt

# Compare line by line

# Compare line by line
if diff -q smiles1.txt smiles2.txt > /dev/null; then
  echo "✅ SMILES columns are identical."
else
  echo "❌ SMILES columns differ. Showing differing lines:"
  diff smiles1.txt smiles2.txt | head -40
fi

# Compare line by line
if diff -q smiles1.txt smiles2.txt > /dev/null; then
  echo "✅ SMILES columns are identical (line by line)."
else
  echo "❌ SMILES columns differ. Showing first few differing lines:"
  diff --side-by-side --suppress-common-lines smiles1.txt smiles2.txt | head -40
fi


# Clean up
rm smiles1.txt smiles2.txt

