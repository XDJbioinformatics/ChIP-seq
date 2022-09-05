INDEX=/stor9000/apps/users/NWSUAF/2021051320/workspeace/data/Chip-seq/random
REF=/stor9000/apps/users/NWSUAF/2021051320/workspeace/ref/bowtie2_ref/genome.fa
script=/stor9000/apps/users/NWSUAF/2021051320/workspeace/script
soft=/stor9000/apps/users/NWSUAF/2021051320/workspeace/software
r=1






if [ ! -d ${INDEX}/2.mapping ]
	then
		mkdir ${INDEX}/2.mapping
fi

cd ${INDEX}/2.mapping

bowtie2 -x ${REF} -p 4 -U ${INDEX}/raw/ssrp_random.fq.gz 2> ssrp_random.mapping.log | samtools view -o ssrp_random.bam
samtools index ${INDEX}/2.mapping/ssrp_random.bam

bowtie2 -x ${REF} -p 4 -U ${INDEX}/raw/ssrp_input_random.fq.gz 2> ssrp_input_random.mapping.log | samtools view -o ssrp_input_random.bam
samtools index ${INDEX}/2.mapping/ssrp_input_random.bam



# 5. Generate bw file

if [ ! -d ${INDEX}/5.bw ]
	then
		mkdir ${INDEX}/5.bw
fi

cd ${INDEX}/5.bw


bamCoverage -p 4 -b ${INDEX}/2.mapping/ssrp_random.bam -o ${INDEX}/5.bw/ssrp_random.bw --normalizeUsing CPM --binSize 10 --extendReads 200

bamCoverage -p 4 -b ${INDEX}/2.mapping/ssrp_input_random.bam -o ${INDEX}/5.bw/ssrp_input_random.bw --normalizeUsing CPM --binSize 10 --extendReads 200

bamCompare -b1 ${INDEX}/5.bw/ssrp_random.bw -b2 ${INDEX}/5.bw/ssrp_input_random.bw -o ${INDEX}/5.bw/ssrp_random_input_CPM.bw --operation log2 -p 4 --scaleFactorsMethod None --normalizeUsing CPM --binSize 10 --extendReads 200
