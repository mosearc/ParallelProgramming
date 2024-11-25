#!/bin/bash

for (( k = 4; k <= 4096; k=k*2 )); do

  gcc-9 -o SEQ SEQ.c
  for i in {1..5}
  do
    ./SEQ "$k"
  done

  gcc-9 -o OMP OMP.c -fopenmp
  for m in {1..96}
  do
    for i in {1..5}
    do
      export OMP_NUM_THREADS="$m"; ./OMP "$k" "$m"
    done
  done

done