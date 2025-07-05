#!/bin/bash -l

module load bioinfo-tools
module load star/2.7.9a

index=/proj/naiss2023-23-494/patrik/reference_guided/

FILES="Sample_SL-2404-P1sample1
Sample_SL-2404-P1sample9
Sample_SL-2404-P1sample13
Sample_SL-2404-P1sample20
Sample_SL-2404-P1sample26
Sample_SL-2404-P1sample27
Sample_SL-2404-P1sample29
Sample_SL-2404-P1sample36
Sample_SL-2404-P1sample37
Sample_SL-2404-P1sample38
Sample_SL-2404-P1sample41
Sample_SL-2404-P1sample45
Sample_SL-2404-P1sample48
Sample_SL-2404-P1sample58
Sample_SL-2404-P1sample59
Sample_SL-2404-P1sample61
Sample_SL-2404-P1sample62
Sample_SL-2404-P1sample71
Sample_SL-2404-P1sample77
Sample_SL-2404-P1sample78
Sample_SL-2404-P1sample79
Sample_SL-2404-P1sample86
Sample_SL-2404-P1sample87
Sample_SL-2404-P1sample89
Sample_SL-2404-P1sample96
Sample_SL-2404-P2sample98
Sample_SL-2404-P2sample99
Sample_SL-2404-P2sample113
Sample_SL-2404-P2sample114
Sample_SL-2404-P2sample120
Sample_SL-2404-P2sample121
Sample_SL-2404-P2sample122
Sample_SL-2404-P2sample125
Sample_SL-2404-P2sample134
Sample_SL-2404-P2sample135
Sample_SL-2404-P2sample136
Sample_SL-2404-P2sample137
Sample_SL-2404-P2sample138
Sample_SL-2404-P2sample155
Sample_SL-2404-P2sample156
Sample_SL-2404-P2sample161
Sample_SL-2404-P2sample167
Sample_SL-2404-P2sample170
Sample_SL-2404-P2sample174
Sample_SL-2404-P2sample175
Sample_SL-2404-P2sample180
Sample_SL-2404-P2sample187
Sample_SL-2404-P2sample191
Sample_SL-2404-P1sample5
Sample_SL-2404-P1sample6
Sample_SL-2404-P1sample14
Sample_SL-2404-P1sample17
Sample_SL-2404-P1sample18
Sample_SL-2404-P1sample21
Sample_SL-2404-P1sample22
Sample_SL-2404-P1sample30
Sample_SL-2404-P1sample31
Sample_SL-2404-P1sample42
Sample_SL-2404-P1sample49
Sample_SL-2404-P1sample50
Sample_SL-2404-P1sample54
Sample_SL-2404-P1sample56
Sample_SL-2404-P1sample63
Sample_SL-2404-P1sample64
Sample_SL-2404-P1sample69
Sample_SL-2404-P1sample73
Sample_SL-2404-P1sample75
Sample_SL-2404-P1sample80
Sample_SL-2404-P1sample84
Sample_SL-2404-P1sample93
Sample_SL-2404-P2sample100
Sample_SL-2404-P2sample104
Sample_SL-2404-P2sample105
Sample_SL-2404-P2sample106
Sample_SL-2404-P2sample110
Sample_SL-2404-P2sample116
Sample_SL-2404-P2sample117
Sample_SL-2404-P2sample128
Sample_SL-2404-P2sample129
Sample_SL-2404-P2sample132
Sample_SL-2404-P2sample139
Sample_SL-2404-P2sample145
Sample_SL-2404-P2sample146
Sample_SL-2404-P2sample149
Sample_SL-2404-P2sample151
Sample_SL-2404-P2sample152
Sample_SL-2404-P2sample157
Sample_SL-2404-P2sample162
Sample_SL-2404-P2sample163
Sample_SL-2404-P2sample171
Sample_SL-2404-P2sample176
Sample_SL-2404-P2sample181
Sample_SL-2404-P2sample185
Sample_SL-2404-P2sample188
Sample_SL-2404-P2sample192"

for f in $FILES;do
    echo $f
    base=$(basename $f)
    echo $base

    STAR --runThreadN 20 \
    --genomeDir $index \
    --readFilesIn /crex/proj/naiss2023-23-494/patrik/transcriptome_populations/200117_A00605_0098_BHV2CFDSXX/${f}/*_R1_*_paired.fastq.gz /crex/proj/naiss2023-23-494/patrik/transcriptome_populations/200117_A00605_0098_BHV2CFDSXX/${f}/*_R2_*_paired.fastq.gz --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --readFilesCommand zcat \
    --outFilterMismatchNmax 10 \
    --outFileNamePrefix $base"_"
done

echo "done!"


