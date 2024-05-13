#!/usr/bin/env bash
set -euo pipefail

# Get text files for ProteinPaint

FILE_ARGS="data/manual/protein_paint_args.txt"

echo "KANSL1 chr17 46029916 46193429" > $FILE_ARGS
echo "MNX1 chr7 157004854 157010663" >> $FILE_ARGS