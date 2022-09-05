#！/usr/bin/env bash
#Author: Bio_XDJ
#Created Time: 2021/09/23


# samples datas name is samplename_rep.fastq.gz    eg: H2AZ_arp6-1_mbd9-3_R1.fastq.gz,
# H2AZ_arp6-1_R1.fastq.gz

# First, you must create a sample.txt for tell the script your sample name  
# eg: H2AZ_arp6-1_mbd9-3
#	  H2AZ_arp6-1
#	  H2AZ_mbd9-3



INDEX=/stor9000/apps/users/NWSUAF/2021051320/workspeace/data/Chip-seq/fas2_hira_H31_H33
REF=/stor9000/apps/users/NWSUAF/2021051320/workspeace/ref/bowtie2_ref/genome.fa
script=/stor9000/apps/users/NWSUAF/2021051320/workspeace/script
soft=/stor9000/apps/users/NWSUAF/2021051320/workspeace/software
thread=10
r=1


# 1.filter data

if [ ! -d ${INDEX}/1.filter ]
	then
		mkdir ${INDEX}/1.filter
fi

cd ${INDEX}/1.filter

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		fastp -w ${thread} -q 20 -l 35 -i ${INDEX}/raw/${sample}_R${rep}.fq.gz -o ${INDEX}/1.filter/${sample}_R${rep}.fq.gz -h ${INDEX}/1.filter/${sample}_R${rep}.html -j ${INDEX}/1.filter/${sample}_R${rep}.json
	done
done


# 2. mapping

if [ ! -d ${INDEX}/2.mapping ]
	then
		mkdir ${INDEX}/2.mapping
fi

cd ${INDEX}/2.mapping

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		bowtie2 -x ${REF} -p ${thread} -U ${INDEX}/1.filter/${sample}_R${rep}.fq.gz 2>${sample}_R${rep}.mapping.log | samtools view -o ${sample}_R${rep}.bam
	done
done

echo -e "Sample\toverall_alignment_rate" >> ${INDEX}/2.mapping/mapping_rate.txt
for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
        echo -e "${sample}_R${rep}\t$(tail -1 ${INDEX}/2.mapping/${sample}_R${rep}.mapping.log | awk '{print $1}')" >> ${INDEX}/2.mapping/mapping_rate.txt
	done
done

# 3. filter again

if [ ! -d ${INDEX}/3.filter_again ]
	then
		mkdir ${INDEX}/3.filter_again
fi

cd ${INDEX}/3.filter_again

# 3.1 sort

if [ ! -d ${INDEX}/3.filter_again/sort ]
	then
		mkdir ${INDEX}/3.filter_again/sort
fi

cd ${INDEX}/3.filter_again/sort

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		samtools sort -@ ${thread} -o ${INDEX}/3.filter_again/sort/${sample}_R${rep}.s.bam ${INDEX}/2.mapping/${sample}_R${rep}.bam
		
		samtools index ${INDEX}/3.filter_again/sort/${sample}_R${rep}.s.bam
		
		samtools flagstat ${INDEX}/3.filter_again/sort/${sample}_R${rep}.s.bam > ${INDEX}/3.filter_again/${sample}_R${rep}.s.metrics.txt
	done
done

# 3.2 basic filter

if [ ! -d ${INDEX}/3.filter_again/basic_filter ]
	then
		mkdir ${INDEX}/3.filter_again/basic_filter
fi

cd ${INDEX}/3.filter_again/basic_filter

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		samtools view -@ ${thread}  -F 3844 -o ${INDEX}/3.filter_again/basic_filter/${sample}_R${rep}.flagF.bam ${INDEX}/3.filter_again/sort/${sample}_R${rep}.s.bam
		
		samtools index ${INDEX}/3.filter_again/basic_filter/${sample}_R${rep}.flagF.bam
		
		samtools flagstat ${INDEX}/3.filter_again/basic_filter/${sample}_R${rep}.flagF.bam > ${INDEX}/3.filter_again/${sample}_R${rep}.flagF.metrics.txt
		
	done
done


# 3.4 Remove PCR duplicates

if [ ! -d ${INDEX}/3.filter_again/pcr_filter ]
	then
		mkdir ${INDEX}/3.filter_again/pcr_filter
fi

cd ${INDEX}/3.filter_again/pcr_filter

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		java -jar /stor9000/apps/users/NWSUAF/2021051320/workspeace/software/picard.jar MarkDuplicates INPUT=${INDEX}/3.filter_again/basic_filter/${sample}_R${rep}.flagF.bam OUTPUT=${INDEX}/3.filter_again/pcr_filter/${sample}_R${rep}.pcrF.bam METRICS_FILE=${INDEX}/3.filter_again/pcr_filter/${sample}_R${rep}.pcrF.picard.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
		
		samtools index ${INDEX}/3.filter_again/pcr_filter/${sample}_R${rep}.pcrF.bam
		
		samtools flagstat ${INDEX}/3.filter_again/pcr_filter/${sample}_R${rep}.pcrF.bam > ${INDEX}/3.filter_again/${sample}_R${rep}.pcrF.metrics.txt
	done
done

# 3.5 Remove readers from organelles

if [ ! -d ${INDEX}/3.filter_again/organelles_filter ]
	then
		mkdir ${INDEX}/3.filter_again/organelles_filter
fi

cd ${INDEX}/3.filter_again/organelles_filter

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
        bam=${INDEX}/2.mapping/${sample}_R${rep}.bam
		
        samtools view -H $bam | sed -n 's/@SQ\tSN:\(.*\)\tLN:\([0-9]*\)/\1\t\2/p' >${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.txt
		
        awk '{print $1"\t0\t"$2}' ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.txt >${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.bed
		
        grep -i -v "chrM\|mt\|mito\|Pt\|mitochondria\|chrC\|cp\|chloroplast" ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.txt >${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.nog.txt
		
        grep -i -v "chrM\|mt\|mito\|Pt\|mitochondria\|chrC\|cp\|chloroplast" ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.bed >${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.nog.bed
		
		samtools view -L ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.nog.bed -o ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.organelleF.bam ${INDEX}/3.filter_again/pcr_filter/${sample}_R${rep}.pcrF.bam
		
		samtools index ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.organelleF.bam
		
		samtools flagstat ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.organelleF.bam >${INDEX}/3.filter_again/${sample}_R${rep}.organelleF.metrics.txt
	done
done

# 3.6 link

cd ${INDEX}/3.filter_again

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		ln -s ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.organelleF.bam ${INDEX}/3.filter_again/${sample}_R${rep}.ff.bam
		ln -s ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.organelleF.bam.bai ${INDEX}/3.filter_again/${sample}_R${rep}.ff.bam.bai
	done
done

# 3.7 summary

#create Rscript

>${INDEX}/3.filter_again/mapping.summary.barplot.R

echo "library(tidyverse)" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "library(cowplot)" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "args <- commandArgs(T)" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "inputf <- args[1]" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "data <- read.table(inputf, header = F, as.is = T)" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "colnames(data) <- c(\"process\", \"number\")" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "labels <- c(s=\"Raw\", flagF=\"Basic filtering\", pf=\"Impreporly paired reads filtering\", qualityF=\"Low mapping quality reads filtering\"," >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo  "      pcrF=\"PCR duplicates filteing\", organelleF=\"Organellic reads filtering\", zf=\"Mate missed reads filtering\")" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "data\$process <- labels[data\$process]" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "data\$rate <- data\$number / data[1,2]" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "data %>% mutate(process = factor(process, levels = process)) %>%" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "ggplot(aes(x = rate, y = process)) + geom_bar(stat = 'identity') +" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "  ylab(\"\") + theme_cowplot()" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "outputf <- sub(\"txt$\", \"png\", args[1])" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo "ggsave(outputf, height = 3, width = 6)" >>${INDEX}/3.filter_again/mapping.summary.barplot.R

chmod 755 ${INDEX}/3.filter_again/mapping.summary.barplot.R

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
	
		cd ${INDEX}/3.filter_again
		mkdir ${INDEX}/3.filter_again/${sample}_R${rep}
		>${INDEX}/3.filter_again/${sample}_R${rep}/${sample}_R${rep}.summary.metrics.txt
		for tp in s flagF pcrF organelleF; 
		do
			num=$(head -n 1 ${INDEX}/3.filter_again/${sample}_R${rep}.$tp.metrics.txt | awk '{print $1}')
			echo -e "$tp\t$num" >>${INDEX}/3.filter_again/${sample}_R${rep}/${sample}_R${rep}.summary.metrics.txt
			done
		cd ${INDEX}/3.filter_again/${sample}_R${rep}
		Rscript ${INDEX}/3.filter_again/mapping.summary.barplot.R ${INDEX}/3.filter_again/${sample}_R${rep}/${sample}_R${rep}.summary.metrics.txt
	done
done

# 4. Quality assessment

if [ ! -d ${INDEX}/4.quality_assessment ]
	then
		mkdir ${INDEX}/4.quality_assessment
fi

cd ${INDEX}/4.quality_assessment

# 4.1 Detect peaks

if [ ! -d ${INDEX}/4.quality_assessment/detect_peaks ]
	then
		mkdir ${INDEX}/4.quality_assessment/detect_peaks
fi

cd ${INDEX}/4.quality_assessment/detect_peaks

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		macs2 callpeak -t ${INDEX}/3.filter_again/${sample}_R${rep}.ff.bam -f BAM -n ${sample}_R${rep} -g 135000000 --keep-dup all --nomodel
	done
done

# 4.2 Library complexity assessment

if [ ! -d ${INDEX}/4.quality_assessment/library_complexity ]
	then
		mkdir ${INDEX}/4.quality_assessment/library_complexity
fi

cd ${INDEX}/4.quality_assessment/library_complexity

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		echo -e "#TotalReads\tDistinctReads\tOneRead\tTwoReads\tNRF=Distinct/Total\tPBC1=OneRead/Distinct\tPBC2=OneRead/TwoReads" >${sample}_R${rep}.pbc.metrics.txt
		bedtools bamtobed -i ${INDEX}/3.filter_again/quality_filter/${sample}_R${rep}.qualityF.bam | awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$6} else{print $1,$3,$6}}' | sort --parallel=${thread} | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' >> ${sample}_R${rep}.pbc.metrics.txt
	done
done

# 4.3 Cross correlatin analys

if [ ! -d ${INDEX}/4.quality_assessment/cross_correlation ]
	then
		mkdir ${INDEX}/4.quality_assessment/cross_correlation
fi

cd ${INDEX}/4.quality_assessment/cross_correlation


for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do 
		Rscript ${script}/run_spp.R -c=${INDEX}/3.filter_again/${sample}_R${rep}.ff.bam -p=2 -s=0:5:500 -savp=${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.pdf -out=${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.metrics.txt
	done
done

pdfunite *.pdf all.pdf




# 5. Generate bw file

if [ ! -d ${INDEX}/5.bw ]
	then
		mkdir ${INDEX}/5.bw
fi

cd ${INDEX}/5.bw


for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do 
			dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
			
			bamCoverage -p ${thread} -b ${INDEX}/3.filter_again/${sample}_R${rep}.ff.bam -o ${INDEX}/5.bw/${sample}_R${rep}.ff.bw --normalizeUsing CPM --binSize 10 --extendReads $dsize
	done
done


multiBigwigSummary bins -b ${INDEX}/5.bw/*.bw -p ${thread} -o ${INDEX}/5.bw/results.npz --outRawCounts ${INDEX}/5.bw/peak_counts.txt 

plotCorrelation --corData ${INDEX}/5.bw/results.npz --corMethod spearman --whatToPlot heatmap --plotFile ${INDEX}/5.bw/co.pdf



# 6. calling peaks

if [ ! -d ${INDEX}/6.callpeak ]
	then
		mkdir ${INDEX}/6.callpeak
fi

cd ${INDEX}/6.callpeak

# 6.1 Callpeak for each repeat of each sample first

for sample in $(cat ${INDEX}/sample2.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		macs2 callpeak -t ${INDEX}/3.filter_again/${sample}_H2A_R${rep}.ff.bam -c ${INDEX}/3.filter_again/${sample}_input_R${rep}.ff.bam -f BAM -n ${sample}_H2A_R${rep} -g 135000000 --keep-dup all -p 0.01 --nomodel
	done
done

# 6.2 Combine repeated samples and then callpeak

for sample in $(cat ${INDEX}/sample2.txt)
do
		macs2 callpeak -t ${INDEX}/3.filter_again/${sample}_H2A_R1.ff.bam ${INDEX}/3.filter_again/${sample}_H2A_R2.ff.bam -c ${INDEX}/3.filter_again/${sample}_input_R1.ff.bam ${INDEX}/3.filter_again/${sample}_input_R2.ff.bam -f BAM -n ${sample}_H2A -g 135000000 --keep-dup all -p 0.01 --nomodel
done

# 6.3 idr analysis

for sample in $(cat ${INDEX}/sample2.txt)
do
		idr --samples ${INDEX}/6.callpeak/${sample}_H2A_R1_peaks.narrowPeak ${INDEX}/6.callpeak/${sample}_H2A_R2_peaks.narrowPeak --peak-list ${INDEX}/6.callpeak/${sample}_H2A_peaks.narrowPeak --input-file-type narrowPeak --output-file ${INDEX}/6.callpeak/${sample}_H2A.idr --rank p.value --soft-idr-threshold 0.05 --plot
done

# 6.4 Save peaks with high consistency as narrowPeak files

if [ ! -d ${INDEX}/6.callpeak/filter ]
	then
		mkdir ${INDEX}/6.callpeak/filter
fi

cd ${INDEX}/6.callpeak/filter

for sample in $(cat ${INDEX}/sample2.txt)
do
		cut -f 1-10 ${INDEX}/6.callpeak/${sample}_H2A.idr | sort -k1,1 -k2,2n -k3,3n > ${INDEX}/6.callpeak/filter/${sample}_H2A.idr.f.narrowPeak
done




macs2 callpeak -t ${INDEX}/3.filter_again/HDA6_R1.ff.bam -c ${INDEX}/3.filter_again/HDA6_input_R1.ff.bam -f BAM -n HDA6_R1 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/HDA6_R2.ff.bam -c ${INDEX}/3.filter_again/HDA6_input_R1.ff.bam -f BAM -n HDA6_R2 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/HDA6_R1.ff.bam ${INDEX}/3.filter_again/HDA6_R2.ff.bam -c ${INDEX}/3.filter_again/HDA6_input_R1.ff.bam -f BAM -n HDA6 -g 135000000 --keep-dup all -p 0.01 --nomodel

idr --samples ${INDEX}/6.callpeak/HDA6_R1_peaks.narrowPeak ${INDEX}/6.callpeak/HDA6_R2_peaks.narrowPeak --peak-list ${INDEX}/6.callpeak/HDA6_peaks.narrowPeak --input-file-type narrowPeak --output-file ${INDEX}/6.callpeak/HDA6.idr --rank p.value --soft-idr-threshold 0.05 --plot

cut -f 1-10 ${INDEX}/6.callpeak/HDA6.idr | sort -k1,1 -k2,2n -k3,3n > ${INDEX}/6.callpeak/HDA6.idr.f.narrowPeak











