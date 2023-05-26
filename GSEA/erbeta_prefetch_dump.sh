#!/bin/bash

while read p; do
  prefetch "$p"
  fasterq-dump -e 24 --outdir ./"$p"/ "$p"
  pigz -p 48 ./"$p"/*.fastq
  rm ./"$p"/*.sra
done <erbeta_list.txt
