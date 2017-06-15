#!/bin/bash
cd $(pwd)
grep "" sample_*/*short*|grep -v "AL" > summary.xls
grep "" sample_*/*svaba_chr4.xls > svaba_chr4_ALL.xls
grep "" sample_*/*svaba_chr8.xls > svaba_chr8_ALL.xls
grep "" sample_*/*svaba_chr11.xls > svaba_chr11_ALL.xls
grep "" sample_*/*pindel_chr4.xls > pindel_chr4_ALL.xls
grep "" sample_*/*pindel_chr8.xls > pindel_chr8_ALL.xls
grep "" sample_*/*pindel_chr11.xls > pindel_chr11_ALL.xls
grep "" sample_*/*lumpy_chr4.xls > lumpy_chr4_ALL.xls
grep "" sample_*/*lumpy_chr8.xls > lumpy_chr8_ALL.xls
grep "" sample_*/*lumpy_chr11.xls > lumpy_chr11_ALL.xls
