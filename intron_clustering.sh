#!/bin/bash -l




conda activate regtools

python leafcutter/clustering/leafcutter_cluster_regtools.py -j test_juncfiles.txt -m 50 -o intron_clustering -l 500000 --nochromcheck=NOCHROMCHECK

conda deactivate
