#!/usr/bin/env bash

cd /global/projectb/scratch/qpzhang/ViCA/100k_segments/
vica_cli.py split -i ../sorted.fa --out ./ --split_depth 6



cd /global/projectb/scratch/qpzhang/ViCA/300k_segments/
vica_cli.py split -i ../sorted.fa --out ./ --split_depth 6 --classes '{2: 300000, 2157: 300000, 2759: 300000, 10239: 300000}'


cd /global/projectb/scratch/qpzhang/ViCA/1M_segments//
vica_cli.py split -i ../sorted.fa --out ./ --split_depth 6 --classes'{2: 1000000, 2157: 1000000, 2759: 1000000, 10239: 1000000}'
