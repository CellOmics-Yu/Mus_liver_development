#!/bin/bash

# Run motif enrichment analysis with HOMER.

findMotifs.pl gene_list.txt mouse output_dir -end 100 -p 2 -fdr -nogo
