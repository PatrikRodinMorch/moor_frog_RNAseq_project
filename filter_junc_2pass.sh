#!/bin/bash -l




cat *_SJ.out.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > 2pass_SJ_out_filtered.tab
