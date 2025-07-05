#!/bin/bash -l


FILES="10_Sample_SL-2404-P1sample38_TE.cntTable
11_Sample_SL-2404-P1sample41_TE.cntTable
12_Sample_SL-2404-P1sample86_TE.cntTable
13_Sample_SL-2404-P1sample45_TE.cntTable
14_Sample_SL-2404-P1sample48_TE.cntTable
15_Sample_SL-2404-P1sample58_TE.cntTable
16_Sample_SL-2404-P1sample59_TE.cntTable
17_Sample_SL-2404-P1sample61_TE.cntTable
18_Sample_SL-2404-P1sample62_TE.cntTable
19_Sample_SL-2404-P1sample71_TE.cntTable
1_Sample_SL-2404-P1sample1_TE.cntTable
20_Sample_SL-2404-P1sample77_TE.cntTable
21_Sample_SL-2404-P1sample78_TE.cntTable
22_Sample_SL-2404-P1sample79_TE.cntTable
23_Sample_SL-2404-P2sample125_TE.cntTable
24_Sample_SL-2404-P1sample87_TE.cntTable
25_Sample_SL-2404-P1sample89_TE.cntTable
26_Sample_SL-2404-P1sample96_TE.cntTable
27_Sample_SL-2404-P2sample98_TE.cntTable
28_Sample_SL-2404-P2sample99_TE.cntTable
29_Sample_SL-2404-P2sample113_TE.cntTable
2_Sample_SL-2404-P1sample9_TE.cntTable
30_Sample_SL-2404-P2sample114_TE.cntTable
31_Sample_SL-2404-P2sample120_TE.cntTable
32_Sample_SL-2404-P2sample121_TE.cntTable
33_Sample_SL-2404-P2sample122_TE.cntTable
34_Sample_SL-2404-P2sample174_TE.cntTable
35_Sample_SL-2404-P2sample134_TE.cntTable
36_Sample_SL-2404-P2sample135_TE.cntTable
37_Sample_SL-2404-P2sample136_TE.cntTable
38_Sample_SL-2404-P2sample137_TE.cntTable
39_Sample_SL-2404-P2sample138_TE.cntTable
3_Sample_SL-2404-P1sample13_TE.cntTable
40_Sample_SL-2404-P2sample155_TE.cntTable
41_Sample_SL-2404-P2sample156_TE.cntTable
42_Sample_SL-2404-P2sample161_TE.cntTable
43_Sample_SL-2404-P2sample167_TE.cntTable
44_Sample_SL-2404-P2sample170_TE.cntTable
45_Sample_SL-2404-P1sample22_TE.cntTable
46_Sample_SL-2404-P2sample175_TE.cntTable
47_Sample_SL-2404-P2sample180_TE.cntTable
48_Sample_SL-2404-P2sample187_TE.cntTable
49_Sample_SL-2404-P2sample191_TE.cntTable
4_Sample_SL-2404-P1sample20_TE.cntTable
50_Sample_SL-2404-P1sample5_TE.cntTable
51_Sample_SL-2404-P1sample6_TE.cntTable
52_Sample_SL-2404-P1sample14_TE.cntTable
53_Sample_SL-2404-P1sample17_TE.cntTable
54_Sample_SL-2404-P1sample18_TE.cntTable
55_Sample_SL-2404-P1sample21_TE.cntTable
56_Sample_SL-2404-P1sample73_TE.cntTable
57_Sample_SL-2404-P1sample30_TE.cntTable
58_Sample_SL-2404-P1sample31_TE.cntTable
59_Sample_SL-2404-P1sample42_TE.cntTable
5_Sample_SL-2404-P1sample26_TE.cntTable
60_Sample_SL-2404-P1sample49_TE.cntTable
61_Sample_SL-2404-P1sample50_TE.cntTable
62_Sample_SL-2404-P1sample54_TE.cntTable
63_Sample_SL-2404-P1sample56_TE.cntTable
64_Sample_SL-2404-P1sample63_TE.cntTable
65_Sample_SL-2404-P1sample64_TE.cntTable
66_Sample_SL-2404-P1sample69_TE.cntTable
67_Sample_SL-2404-P2sample117_TE.cntTable
68_Sample_SL-2404-P1sample75_TE.cntTable
69_Sample_SL-2404-P1sample80_TE.cntTable
6_Sample_SL-2404-P1sample27_TE.cntTable
70_Sample_SL-2404-P1sample84_TE.cntTable
71_Sample_SL-2404-P1sample93_TE.cntTable
72_Sample_SL-2404-P2sample100_TE.cntTable
73_Sample_SL-2404-P2sample104_TE.cntTable
74_Sample_SL-2404-P2sample105_TE.cntTable
75_Sample_SL-2404-P2sample106_TE.cntTable
76_Sample_SL-2404-P2sample110_TE.cntTable
77_Sample_SL-2404-P2sample116_TE.cntTable
78_Sample_SL-2404-P2sample162_TE.cntTable
79_Sample_SL-2404-P2sample128_TE.cntTable
7_Sample_SL-2404-P1sample29_TE.cntTable
80_Sample_SL-2404-P2sample129_TE.cntTable
81_Sample_SL-2404-P2sample132_TE.cntTable
82_Sample_SL-2404-P2sample139_TE.cntTable
83_Sample_SL-2404-P2sample145_TE.cntTable
84_Sample_SL-2404-P2sample146_TE.cntTable
85_Sample_SL-2404-P2sample149_TE.cntTable
86_Sample_SL-2404-P2sample151_TE.cntTable
87_Sample_SL-2404-P2sample152_TE.cntTable
88_Sample_SL-2404-P2sample157_TE.cntTable
89_Sample_SL-2404-P2sample192_TE.cntTable
8_Sample_SL-2404-P1sample36_TE.cntTable
90_Sample_SL-2404-P2sample163_TE.cntTable
91_Sample_SL-2404-P2sample171_TE.cntTable
92_Sample_SL-2404-P2sample176_TE.cntTable
93_Sample_SL-2404-P2sample181_TE.cntTable
94_Sample_SL-2404-P2sample185_TE.cntTable
95_Sample_SL-2404-P2sample188_TE.cntTable
9_Sample_SL-2404-P1sample37_TE.cntTable"


tail -n+2 10_Sample_SL-2404-P1sample38_TE.cntTable > 10_Sample_SL-2404-P1sample38_TE.cntTable_temp
awk '{print $1}' 10_Sample_SL-2404-P1sample38_TE.cntTable_temp > TE_count_matrix_genes.txt



for f in $FILES;do
	cat ${f} | tail -n+2 | awk '{print $2}' > ${f}_cleaned
done

rm temp.txt
touch temp.txt

paste -d "\t" TE_count_matrix_genes.txt {1..95}_Sample_*_TE.cntTable_cleaned >> temp.txt

cp NewHeader_liver.txt TE_count_matrix_final.txt

cat temp.txt >> TE_count_matrix_final.txt
