#!/bin/bash

# arg 1 should be the full path to the pipeline database
# arg 2 = full path to the bacterial genome DB (made with makeblastdb)
# arg 3 = full path to the Orthologs folder
# arg 4 = type of xml blast file (genome results, transcriptome results, or ORFs present in both dbs)

pipe_path="$echo$1"
db_path="$echo$2"
outdir="$echo$3"
xml="$echo$4"

../dependencies/blast_for_uproteins/bin/tblastn -query "$pipe_path" -db "$db_path" -db_gencode 11 -evalue 0.01 -outfmt 5 -out "$outdir"/homologues"$xml".xml
