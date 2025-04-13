#!/bin/bash

# Usage: ./run_tests.sh n1 n2 n3 ...

# Check if at least one `n` value is provided
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 n1 n2 n3 ..."
  exit 1
fi

# Loop over all provided values for `n`
for n in "$@"; do
  seed=$((n)) # Calculate seed as n (+ 1)
  echo "n = $n (seed = $seed)"
  Rscript properness_test.R 10000 1000 "$n" TRUE "$seed"
done

