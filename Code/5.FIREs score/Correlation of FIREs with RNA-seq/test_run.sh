#!/bin/bash
#Author:yzy

samtools depth 1.R5_merge.bam > 3.DK_cover.txt
python modify_border.py -f _test_chrom.len -r 2000 -d test_boundary_down2kb_BY4741.bed
python count_depth.py -f test_chrom.len -r 2000 -d 3.DK_cover.txt -m test_HiCNormCis.bedgraph.checked
python correlation_transcription_3C_use.py chr1-2000.fire chr1-2000.rna.depth chr1-1000.border.bed DK_chr1_correlation_transcription_coverage_2K.pdf 2