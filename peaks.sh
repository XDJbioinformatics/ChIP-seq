



INDEX=/stor9000/apps/users/NWSUAF/2021051320/workspeace/data/Chip-seq/fas2_hira_H31_H33
genomef=/stor9000/apps/users/NWSUAF/2021051320/workspeace/ref/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

# 6. calling peaks

if [ ! -d ${INDEX}/6.callpeak ]
	then
		mkdir ${INDEX}/6.callpeak
fi

cd ${INDEX}/6.callpeak

# 6.1 Callpeak for each repeat of each sample first

macs2 callpeak -t ${INDEX}/3.filter_again/Col-H31_R1.ff.bam -c ${INDEX}/3.filter_again/Col-H31_Input_R1.ff.bam -f BAM -n Col-H31_R1 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/Col-H33_R1.ff.bam -c ${INDEX}/3.filter_again/Col-H33_Input_R1.ff.bam -f BAM -n Col-H33_R1 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/fas2-4-H31_R1.ff.bam -c ${INDEX}/3.filter_again/fas2-4-H31_Input_R1.ff.bam -f BAM -n fas2-4-H31_R1 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/fas2-4-H33_R1.ff.bam -c ${INDEX}/3.filter_again/fas2-4-H33_Input_R1.ff.bam -f BAM -n fas2-4-H33_R1 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/hira-1-H31_R1.ff.bam -c ${INDEX}/3.filter_again/hira-1-H31_Input_R1.ff.bam -f BAM -n hira-1-H31_R1 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/hira-1-H33_R1.ff.bam -c ${INDEX}/3.filter_again/hira-1-H33_Input_R1.ff.bam -f BAM -n hira-1-H33_R1 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/ssrp1-1_R2.ff.bam -c ${INDEX}/3.filter_again/ssrp1-1_input_R2.ff.bam -f BAM -n ssrp1-1_R2 -g 135000000 --keep-dup all -p 0.01 --nomodel


# 6.2 Combine repeated samples and then callpeak

macs2 callpeak -t ${INDEX}/3.filter_again/ssrp1-1_R1.ff.bam ${INDEX}/3.filter_again/ssrp1-1_R2.ff.bam -c ${INDEX}/3.filter_again/ssrp1-1_input_R1.ff.bam ${INDEX}/3.filter_again/ssrp1-1_input_R2.ff.bam -f BAM -n ssrp1-1 -g 135000000 --keep-dup all -p 0.01 --nomodel


# 6.3 idr analysis

idr --samples ${INDEX}/6.callpeak/ssrp1-1_R1_peaks.narrowPeak ${INDEX}/6.callpeak/ssrp1-1_R2_peaks.narrowPeak --peak-list ${INDEX}/6.callpeak/ssrp1-1_peaks.narrowPeak --input-file-type narrowPeak --output-file ${INDEX}/6.callpeak/ssrp1-1.idr --rank p.value --soft-idr-threshold 0.05 --plot


# 6.4 Save peaks with high consistency as narrowPeak files

if [ ! -d ${INDEX}/6.callpeak/filter ]
	then
		mkdir ${INDEX}/6.callpeak/filter
fi

cd ${INDEX}/6.callpeak/filter


cut -f 1-10 ${INDEX}/6.callpeak/ssrp1-1.idr | sort -k1,1 -k2,2n -k3,3n > ${INDEX}/6.callpeak/filter/ssrp1-1.idr.f.narrowPeak


# get the coordination of summit
# only use the top 500 peaks having the highest pscore
sort -k8,8nr ${INDEX}/6.callpeak/filter/ANP32E-3Flag_R1_peaks.narrowPeak | head -n 500 | awk -v OFS="\t" '{print $1,$2+$10,$2+$10+1}' > ${INDEX}/6.callpeak/filter/summit.bed

# extend summmit to upstream 50 bp and downstream 50 bp, respectively
awk -v OFS="\t" '{print $1,$2-50,$2+50}' ${INDEX}/6.callpeak/filter/summit.bed | awk '$2>=0' > ${INDEX}/6.callpeak/filter/summit.l50.r50.bed

# get sequences
bedtools getfasta -fi ${genomef} -bed ${INDEX}/6.callpeak/filter/summit.l50.r50.bed > ${INDEX}/6.callpeak/filter/summit.l50.r50.fa







