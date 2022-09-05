#！/usr/bin/env bash
#Author: Bio_XDJ
#Created Time: 2021/09/23


# samples datas name is samplename_rep.fastq.gz    eg: H2AZ_arp6-1_mbd9-3_R1.fq.gz,
# H2AZ_arp6-1_R1.fq.gz

# First, you must create a sample.txt for tell the script your sample name  
# eg: H2AZ_arp6-1_mbd9-3
#	  H2AZ_arp6-1
#	  H2AZ_mbd9-3

# If have sample2.txt file,must have name file



INDEX=/stor9000/apps/users/NWSUAF/2021051320/workspeace/data/Chip-seq/
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
		fastp -w ${thread} -q 20 -l 35 -i ${INDEX}/raw/${sample}_R${rep}_1.fq.gz -o ${INDEX}/1.filter/${sample}_R${rep}_1.fq.gz -I ${INDEX}/raw/${sample}_R${rep}_2.fq.gz -O ${INDEX}/1.filter/${sample}_R${rep}_2.fq.gz -h ${INDEX}/1.filter/${sample}_R${rep}.html -j ${INDEX}/1.filter/${sample}_R${rep}.json
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
		bowtie2 -x ${REF} -p ${thread} -1 ${INDEX}/1.filter/${sample}_R${rep}_1.fq.gz -2 ${INDEX}/1.filter/${sample}_R${rep}_2.fq.gz 2>${sample}_R${rep}.mapping.log | samtools view -o ${sample}_R${rep}.bam
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
		samtools view -@ ${thread} -F 3844 -o ${INDEX}/3.filter_again/basic_filter/${sample}_R${rep}.flagF.bam ${INDEX}/3.filter_again/sort/${sample}_R${rep}.s.bam
		
		samtools index ${INDEX}/3.filter_again/basic_filter/${sample}_R${rep}.flagF.bam
		
		samtools flagstat ${INDEX}/3.filter_again/basic_filter/${sample}_R${rep}.flagF.bam > ${INDEX}/3.filter_again/${sample}_R${rep}.flagF.metrics.txt
		
	done
done

#3.3 Filter improperly paired read

if [ ! -d ${INDEX}/3.filter_again/improperly_paired ]
	then
		mkdir ${INDEX}/3.filter_again/improperly_paired
fi

cd ${INDEX}/3.filter_again/improperly_paired

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		samtools view -@ ${thread}  -F 3844 -o ${INDEX}/3.filter_again/improperly_paired/${sample}_R${rep}.pF.bam ${INDEX}/3.filter_again/basic_filter/${sample}_R${rep}.flagF.bam
		
		samtools index ${INDEX}/3.filter_again/improperly_paired/${sample}_R${rep}.pF.bam
		
		samtools flagstat ${INDEX}/3.filter_again/improperly_paired/${sample}_R${rep}.pF.bam > ${INDEX}/3.filter_again/${sample}_R${rep}.pF.metrics.txt
		
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
		java -jar /stor9000/apps/users/NWSUAF/2021051320/workspeace/software/picard.jar MarkDuplicates INPUT=${INDEX}/3.filter_again/improperly_paired/${sample}_R${rep}.pF.bam OUTPUT=${INDEX}/3.filter_again/pcr_filter/${sample}_R${rep}.pcrF.bam METRICS_FILE=${INDEX}/3.filter_again/pcr_filter/${sample}_R${rep}.pcrF.picard.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
		
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
		
        grep -i -v "chrM\|mt\|mito\|mitochondria\|Pt\|chrC\|cp\|chloroplast" ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.txt >${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.nog.txt
		
        grep -i -v "chrM\|mt\|mito\|mitochondria\|Pt\|chrC\|cp\|chloroplast" ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.bed >${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.nog.bed
		
		samtools view -L ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.nog.bed -o ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.organelleF.bam ${INDEX}/3.filter_again/pcr_filter/${sample}_R${rep}.pcrF.bam
		
		samtools index ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.organelleF.bam
		
		samtools flagstat ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.organelleF.bam >${INDEX}/3.filter_again/${sample}_R${rep}.organelleF.metrics.txt
	done
done

#3.6 remove unpaired read

if [ ! -d ${INDEX}/3.filter_again/unpaired_read ]
	then
		mkdir ${INDEX}/3.filter_again/unpaired_read
fi

cd ${INDEX}/3.filter_again/unpaired_read

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		samtools sort -n ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.organelleF.bam | samtools fixmate -r - -| samtools view -b -f 2 | samtools sort -o ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam
		
		samtools index ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam
		
		samtools flagstat ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam > ${INDEX}/3.filter_again/${sample}_R${rep}.zF.metrics.txt
		
	done
done

# cheak insert size
for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		java -jar /stor9000/apps/users/NWSUAF/2021051320/workspeace/software/picard.jar CollectInsertSizeMetrics I=${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam O=${sample}.insertsize.metrics.txt H=${sample}.insertsize.pdf HISTOGRAM_WIDTH=500
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
echo "labels <- c(s=\"Raw\", flagF=\"Basic filtering\", pF=\"Impreporly paired reads filtering\", qualityF=\"Low mapping quality reads filtering\"," >>${INDEX}/3.filter_again/mapping.summary.barplot.R
echo  "      pcrF=\"PCR duplicates filteing\", organelleF=\"Organellic reads filtering\", zF=\"Mate missed reads filtering\")" >>${INDEX}/3.filter_again/mapping.summary.barplot.R
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
	
		cd ${INDEX}/3.filter_againll
		
		mkdir ${INDEX}/3.filter_again/${sample}_R${rep}
		>${INDEX}/3.filter_again/${sample}_R${rep}/${sample}_R${rep}.summary.metrics.txt
		for tp in s flagF pcrF organelleF pF zF; 
		do
			num=$(head -n 1 ${INDEX}/3.filter_again/${sample}_R${rep}.$tp.metrics.txt | awk '{print $1}')
			echo -e "$tp\t$num" >>${INDEX}/3.filter_again/${sample}_R${rep}/${sample}_R${rep}.summary.metrics.txt
			done
		cd ${INDEX}/3.filter_again/${sample}_R${rep}
		Rscript ${INDEX}/3.filter_again/mapping.summary.barplot.R ${INDEX}/3.filter_again/${sample}_R${rep}/${sample}_R${rep}.summary.metrics.txt
	done
done

# 3.8 Extract read1 mapping information

cd ${INDEX}/3.filter_again

for sample in $(cat ${INDEX}/sample.txt)
do
	for ((rep=1;rep<${r}+1;rep++))
	do
		samtools view -b -F 128 ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam | bedtools bamtobed -i - | bedtools bedtobam -i - -g ${INDEX}/3.filter_again/organelles_filter/${sample}_R${rep}.chrSize.txt > ${INDEX}/3.filter_again/${sample}_R${rep}.ff.bam
		
		samtools index ${INDEX}/3.filter_again/${sample}_R${rep}.ff.bam
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
		macs2 callpeak -t ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam -f BAM -n ${sample}_R${rep} -g 135000000 --keep-dup all --nomodel
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
		bedtools bamtobed -i ${INDEX}/3.filter_again/improperly_paired/${sample}_R${rep}.pF.bam | awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$6} else{print $1,$3,$6}}' | sort --parallel=${thread} | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' >> ${sample}_R${rep}.pbc.metrics.txt
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
		Rscript ${script}/run_spp.R -c=${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam -p=4 -s=0:5:500 -savp=${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.pdf -out=${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.metrics.txt
	done
done

pdfunite *.pdf all.pdf




# 5. Generate bw file

if [ ! -d ${INDEX}/5.bw ]
	then
		mkdir ${INDEX}/5.bw
fi

for sample in $(cat ${INDEX}/sample.txt)
	do
	for ((rep=1;rep<${r}+1;rep++))
		do
			dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
			
			bamCoverage -b ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam -o ${INDEX}/5.bw/${sample}_R${rep}.bw -p ${thread} --normalizeUsing CPM --binSize 10 --extendReads $dsize
		done
	done

cp /stor9000/apps/users/NWSUAF/2021051320/workspeace/data/Chip-seq/BigWig/trans.bed ${INDEX}/5.bw/

for sample in $(cat ${INDEX}/sample.txt)
	do
		for ((rep=1;rep<${r}+1;rep++))
		do
			computeMatrix scale-regions -R ${INDEX}/5.bw/trans.bed -S ${INDEX}/5.bw/${sample}_R${rep}.bw --startLabel TSS --endLabel TTS -m 5000 -b 2000 -a 2000 -p ${thread} -o ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz
			
			plotProfile -m ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz -out ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.profile.pdf --dpi 600 --yAxisLabel CPM --outFileSortedRegions ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.regions.txt  --samplesLabel ${sample}_R${rep}
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

macs2 callpeak -t ${INDEX}/3.filter_again/unpaired_read/ANP32E-3Flag_R1.zF.bam -c ${INDEX}/3.filter_again/unpaired_read/ANP32E-3Flag_Input_R1.zF.bam -f BAMPE -n ANP32E-3Flag_R1 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/unpaired_read/ANP32E-3Flag_R2.zF.bam -c ${INDEX}/3.filter_again/unpaired_read/ANP32E-3Flag_Input_R2.zF.bam -f BAMPE -n ANP32E-3Flag_R2 -g 135000000 --keep-dup all -p 0.01 --nomodel

macs2 callpeak -t ${INDEX}/3.filter_again/unpaired_read/ANP32E-3Flag_R1.zF.bam ${INDEX}/3.filter_again/unpaired_read/ANP32E-3Flag_R2.zF.bam -c ${INDEX}/3.filter_again/unpaired_read/ANP32E-3Flag_Input_R1.zF.bam  ${INDEX}/3.filter_again/unpaired_read/ANP32E-3Flag_Input_R2.zF.bam -f BAMPE -n ANP32E-3Flag -g 135000000 --keep-dup all -p 0.01 --nomodel

idr --samples ${INDEX}/6.callpeak/ANP32E-3Flag_R1_peaks.narrowPeak ${INDEX}/6.callpeak/ANP32E-3Flag_R2_peaks.narrowPeak --peak-list ${INDEX}/6.callpeak/ANP32E-3Flag_peaks.narrowPeak --input-file-type narrowPeak --output-file ${INDEX}/6.callpeak/ANP32E-3Flag.idr --rank p.value --soft-idr-threshold 0.05 --plot

cut -f 1-10 ${INDEX}/6.callpeak/ANP32E-3Flag.idr | sort -k1,1 -k2,2n -k3,3n > ${INDEX}/6.callpeak/ANP32E-3Flag.idr.f.narrowPeak



# 6. calling peaks

if [ ! -d ${INDEX}/6.peaks ]
	then
		mkdir ${INDEX}/6.peaks
fi

cd ${INDEX}/6.peaks

if [ -a ${INDEX}/sample2.txt ]
	then
		touch ${INDEX}/6.peaks/Have_Input_file
		
		if ((${r} > 1))
			then
			# Have biological replicate
		
			touch ${INDEX}/6.peaks/Have_replicate
		
			# 6.1 Callpeak for each repeat of each sample first
			
			
			for sample in $(cat ${INDEX}/sample2.txt)
			do
				for ((rep=1;rep<${r}+1;rep++))
				do
					macs2 callpeak -t ${INDEX}/3.filter_again/unpaired_read/${sample}_${NAME}_R${rep}.zF.bam -c ${INDEX}/3.filter_again/unpaired_read/${sample}_In_R${rep}.zF.bam -f BAMPE -n ${sample}_${NAME}_R${rep} -g 135000000 --keep-dup all -p 0.01
				done
			done



			# 6.2 Combine repeated samples and then callpeak
			
			for sample in $(cat ${INDEX}/sample2.txt)
			do
				macs2 callpeak -t ${INDEX}/3.filter_again/unpaired_read/${sample}_${NAME}_R1.zF.bam ${INDEX}/3.filter_again/unpaired_read/${sample}_${NAME}_R2.zF.bam -c ${INDEX}/3.filter_again/unpaired_read/${sample}_In_R1.zF.bam ${INDEX}/3.filter_again/unpaired_read/${sample}_In_R2.zF.bam -f BAMPE -n ${sample}_${NAME} -g 135000000 --keep-dup all -p 0.01
			done
			
			
			# 6.3 idr analysis
			
			for sample in $(cat ${INDEX}/sample2.txt)
			do
				idr --samples ${INDEX}/6.peaks/${sample}_${NAME}_R1_peaks.narrowPeak ${INDEX}/6.peaks/${sample}_${NAME}_R2_peaks.narrowPeak --peak-list ${INDEX}/6.peaks/${sample}_${NAME}_peaks.narrowPeak --input-file-type narrowPeak --output-file ${INDEX}/6.peaks/${sample}_${NAME}.idr --rank p.value --soft-idr-threshold 0.05 --plot
			done
			
			
			# 6.4 Save peaks with high consistency as narrowPeak files
			
			for sample in $(cat ${INDEX}/sample2.txt)
			do
				cut -f 1-10 ${INDEX}/6.peaks/${sample}_${NAME}.idr | sort -k1,1 -k2,2n -k3,3n > ${INDEX}/6.peaks/${sample}_${NAME}.idr.f.narrowPeak
			done
			
		else
			touch ${INDEX}/6.peaks/Have_NOT_replicate
		
			# 6.1 Callpeak for each repeat of each sample first
			
			
			for sample in $(cat ${INDEX}/sample2.txt)
			do
				for ((rep=1;rep<${r}+1;rep++))
				do
					macs2 callpeak -t ${INDEX}/3.filter_again/unpaired_read/${sample}_${NAME}_R${rep}.zF.bam -c ${INDEX}/3.filter_again/unpaired_read/${sample}_In_R${rep}.zF.bam -f BAMPE -n ${sample}_${NAME}_R${rep} -g 135000000 --keep-dup all -p 0.01
				done
			done
		fi
		
	else
	
		touch ${INDEX}/6.peaks/Have_NOT_Input_file
		
		if ((${r} > 1))
			then
			# Have biological replicate
		
			touch ${INDEX}/6.peaks/Have_replicate
		
			# 6.1 Callpeak for each repeat of each sample first
			
			
			for sample in $(cat ${INDEX}/sample2.txt)
			do
				for ((rep=1;rep<${r}+1;rep++))
				do
					macs2 callpeak -t ${INDEX}/3.filter_again/unpaired_read/${sample}_${NAME}_R${rep}.zF.bam -c ${INDEX}/3.filter_again/unpaired_read/${sample}_In_R${rep}.zF.bam -f BAMPE -n ${sample}_${NAME}_R${rep} -g 135000000 --keep-dup all -p 0.01
				done
			done
		fi
fi



			# 6.2 Combine repeated samples and then callpeak
			
			for sample in $(cat ${INDEX}/sample2.txt)
			do
				macs2 callpeak -t ${INDEX}/3.filter_again/unpaired_read/${sample}_${NAME}_R1.zF.bam ${INDEX}/3.filter_again/unpaired_read/${sample}_${NAME}_R2.zF.bam -c ${INDEX}/3.filter_again/unpaired_read/${sample}_In_R1.zF.bam ${INDEX}/3.filter_again/unpaired_read/${sample}_In_R2.zF.bam -f BAMPE -n ${sample}_${NAME} -g 135000000 --keep-dup all -p 0.01
			done
			
			
			# 6.3 idr analysis
			
			for sample in $(cat ${INDEX}/sample2.txt)
			do
				idr --samples ${INDEX}/6.peaks/${sample}_${NAME}_R1_peaks.narrowPeak ${INDEX}/6.peaks/${sample}_${NAME}_R2_peaks.narrowPeak --peak-list ${INDEX}/6.peaks/${sample}_${NAME}_peaks.narrowPeak --input-file-type narrowPeak --output-file ${INDEX}/6.peaks/${sample}_${NAME}.idr --rank p.value --soft-idr-threshold 0.05 --plot
			done
			
			
			# 6.4 Save peaks with high consistency as narrowPeak files
			
			for sample in $(cat ${INDEX}/sample2.txt)
			do
				cut -f 1-10 ${INDEX}/6.peaks/${sample}_${NAME}.idr | sort -k1,1 -k2,2n -k3,3n > ${INDEX}/6.peaks/${sample}_${NAME}.idr.f.narrowPeak
			done
fi







computeMatrix scale-regions -R ${INDEX}/5.bw/trans.bed -S ${INDEX}/5.bw/anp_H3_merge.bw ${INDEX}/5.bw/Col_H3_merge.bw ${INDEX}/5.bw/nrp_H3_merge.bw ${INDEX}/5.bw/tri_H3_merge.bw --startLabel TSS --endLabel TTS -m 5000 -b 2000 -a 2000 -p 4 -o ${INDEX}/5.bw/all_H3_merge.trans.ff.mtr.gz

plotProfile -m ${INDEX}/5.bw/all.trans.ff.mtr.gz -out ${INDEX}/5.bw/all.trans.ff.profile.pdf --dpi 600 --yAxisLabel CPM --outFileSortedRegions ${INDEX}/5.bw/all.trans.ff.regions.txt --samplesLabel anp32e-1 columbia "nrp1-1,nrp2-2" "anp32e-1,nrp1-1,nrp2-2" --colors \#E41A1C \#FF7F00 \#984EA3 \#377EB8 --plotHeight 15 --plotWidth 20 --perGroup


multiBigwigSummary bins -b ${INDEX}/bw/*.bw -p 4 -o ${INDEX}/results.npz --outRawCounts ${INDEX}/peak_counts.txt 

plotCorrelation --corData results.npz --corMethod spearman --whatToPlot heatmap --plotFile co.pdf

























# Determine whether there is biological duplication

if ((${r} > 1))
	then
		# Have biological replicate
		if [ ! -d ${INDEX}/5.bw/merge ]
			then
			mkdir ${INDEX}/5.bw/merge
		fi
		
		cd ${INDEX}/5.bw/merge
		
		touch ${INDEX}/5.bw/Have_replicate
		
		for sample in $(cat ${INDEX}/sample.txt)
		do

			samtools merge -o ${INDEX}/5.bw/merge/${sample}_merge.ff.bam ${INDEX}/3.filter_again/unpaired_read/${sample}_R1.zF.bam ${INDEX}/3.filter_again/unpaired_read/${sample}_R2.zF.bam

			samtools index ${INDEX}/5.bw/merge/${sample}_merge.ff.bam
		done
		
		cd ${INDEX}/5.bw
		
		# Determine whether there is input file
		
		if [ -a ${INDEX}/sample2.txt ]
			then
				# Have biological replicate and input
				touch ${INDEX}/5.bw/Have_Input
				
				for sample in $(cat ${INDEX}/sample.txt)
					do
					for ((rep=1;rep<${r}+1;rep++))
						do
							dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
							
							bamCoverage -b ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam -o ${INDEX}/5.bw/${sample}_R${rep}.coverage.bw -p 4 --normalizeUsing CPM --binSize 10 --extendReads $dsize
						done
					done
				
				for sample in $(cat ${INDEX}/sample2.txt)
				do
					for ((rep=1;rep<${r}+1;rep++))
					do
						dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_${NAME}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
						
						bamCompare -b1 ${INDEX}/3.filter_again/unpaired_read/${sample}_${NAME}_R${rep}.zF.bam -b2 ${INDEX}/3.filter_again/unpaired_read/${sample}_In_R${rep}.zF.bam -o ${INDEX}/5.bw/${sample}_${NAME}_R${rep}.bw --operation log2 -p 4 --scaleFactorsMethod None --normalizeUsing CPM --binSize 10 --extendReads $dsize
					done
				done
				
				
				for sample in $(cat ${INDEX}/sample2.txt)
				do
					dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_${NAME}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
					
					bamCompare -b1 ${INDEX}/5.bw/merge/${sample}_${NAME}_merge.ff.bam -b2 ${INDEX}/5.bw/merge/${sample}_${NAME}_merge.ff.bam -o ${INDEX}/5.bw/${sample}_${NAME}_merge.bw --operation log2 -p 4 --scaleFactorsMethod None --normalizeUsing CPM --binSize 10 --extendReads $dsize
				done
				
				
				for sample in $(cat ${INDEX}/sample.txt)
				do
					for ((rep=1;rep<${r}+1;rep++))
					do
						computeMatrix scale-regions -R ${INDEX}/5.bw/trans.bed -S ${INDEX}/5.bw/${sample}_R${rep}.bw --startLabel TSS --endLabel TTS -m 5000 -b 2000 -a 2000 -p 4 -o ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz
					
						plotProfile -m ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz -out ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.profile.pdf --dpi 600 --yAxisLabel CPM --outFileSortedRegions ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.regions.txt --regionsLabel ALL --samplesLabel ${sample}_R${rep}
					done
				done
				
				for sample in $(cat ${INDEX}/sample.txt)
				do
					computeMatrix scale-regions -R ${INDEX}/5.bw/trans.bed -S ${INDEX}/5.bw/${sample}_merge.bw --startLabel TSS --endLabel TTS -m 5000 -b 2000 -a 2000 -p 4 -o ${INDEX}/5.bw/${sample}_merge.trans.ff.mtr.gz
				
					plotProfile -m ${INDEX}/5.bw/${sample}_merge.trans.ff.mtr.gz -out ${INDEX}/5.bw/${sample}_merge.trans.ff.profile.pdf --dpi 600 --yAxisLabel CPM --outFileSortedRegions ${INDEX}/5.bw/${sample}_merge.trans.ff.regions.txt --samplesLabel ${sample}_merge
				done
			
			else
			# Have biological replicate and NOT input
				touch ${INDEX}/5.bw/Have_NOT_Input
			
			for sample in $(cat ${INDEX}/sample.txt)
				do
					for ((rep=1;rep<${r}+1;rep++))
					do
						dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
						
						bamCoverage -b ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam -o ${INDEX}/5.bw/${sample}_R${rep}.bw -p 4 --normalizeUsing CPM --binSize 10 --extendReads $dsize
					done
				done
				
				
				for sample in $(cat ${INDEX}/sample.txt)
				do
					dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
					
					bamCoverage -b ${INDEX}/5.bw/merge/${sample}_merge.ff.bam -o ${INDEX}/5.bw/${sample}_merge.bw -p 4 --scaleFactorsMethod None --normalizeUsing CPM --binSize 10 --extendReads $dsize
				done
				
				
				for sample in $(cat ${INDEX}/sample.txt)
				do
					for ((rep=1;rep<${r}+1;rep++))
					do
						computeMatrix scale-regions -R ${INDEX}/5.bw/trans.bed -S ${INDEX}/5.bw/${sample}_R${rep}.bw --startLabel TSS --endLabel TTS -m 5000 -b 2000 -a 2000 -p 4 -o ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz
					
						plotProfile -m ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz -out ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.profile.pdf --dpi 600 --yAxisLabel CPM --outFileSortedRegions ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.regions.txt  --samplesLabel ${sample}_R${rep}
					done
				done
				
				for sample in $(cat ${INDEX}/sample.txt)
				do
					computeMatrix scale-regions -R ${INDEX}/5.bw/trans.bed -S ${INDEX}/5.bw/${sample}_merge.bw --startLabel TSS --endLabel TTS -m 5000 -b 2000 -a 2000 -p 4 -o ${INDEX}/5.bw/${sample}_merge.trans.ff.mtr.gz
				
					plotProfile -m ${INDEX}/5.bw/${sample}_merge.trans.ff.mtr.gz -out ${INDEX}/5.bw/${sample}_merge.trans.ff.profile.pdf --dpi 600 --yAxisLabel CPM --outFileSortedRegions ${INDEX}/5.bw/${sample}_merge.trans.ff.regions.txt --samplesLabel ${sample}_merge
				done
		fi
	else
		# NOT have biological replicate
		touch ${INDEX}/5.bw/Have_NOT_replicate
		
		# Determine whether there is input file
		
		if [ -a ${INDEX}/sample2.txt ]
			then
				# NOT have biological replicate and have input
				touch ${INDEX}/5.bw/Have_Input
				
				for sample in $(cat ${INDEX}/sample.txt)
					do
					for ((rep=1;rep<${r}+1;rep++))
						do
							dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
							
							bamCoverage -b ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam -o ${INDEX}/5.bw/${sample}_R${rep}.coverage.bw -p 4 --normalizeUsing CPM --binSize 10 --extendReads $dsize
						done
					done
				
				for sample in $(cat ${INDEX}/sample2.txt)
				do
					for ((rep=1;rep<${r}+1;rep++))
					do
						dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_${NAME}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
						
						bamCompare -b1 ${INDEX}/3.filter_again/unpaired_read/${sample}_${NAME}_R${rep}.zF.bam -b2 ${INDEX}/3.filter_again/unpaired_read/${sample}_In_R${rep}.zF.bam -o ${INDEX}/5.bw/${sample}_${NAME}_R${rep}.bw --operation log2 -p 4 --scaleFactorsMethod None --normalizeUsing CPM --binSize 10 --extendReads $dsize
					done
				done
				
				
				for sample in $(cat ${INDEX}/sample.txt)
				do
					for ((rep=1;rep<${r}+1;rep++))
					do
						computeMatrix scale-regions -R ${INDEX}/5.bw/trans.bed -S ${INDEX}/5.bw/${sample}_R${rep}.bw --startLabel TSS --endLabel TTS -m 5000 -b 2000 -a 2000 -p 4 -o ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz
					
						plotProfile -m ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz -out ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.profile.pdf --dpi 600 --yAxisLabel CPM --outFileSortedRegions ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.regions.txt --regionsLabel ALL --samplesLabel ${sample}_R${rep}
					done
				done
				
			
			else
			# NOT have biological replicate and input
				touch ${INDEX}/5.bw/Have_NOT_Input
			
			for sample in $(cat ${INDEX}/sample.txt)
				do
					for ((rep=1;rep<${r}+1;rep++))
					do
						dsize=$(awk '{print $3}' ${INDEX}/4.quality_assessment/cross_correlation/${sample}_R${rep}.cc.metrics.txt | awk -F ',' '{print $1}')
						
						bamCoverage -b ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam -o ${INDEX}/5.bw/${sample}_R${rep}.bw -p 4 --normalizeUsing CPM --binSize 10 --extendReads $dsize
					done
				done
				
				
				for sample in $(cat ${INDEX}/sample.txt)
				do
					for ((rep=1;rep<${r}+1;rep++))
					do
						computeMatrix scale-regions -R ${INDEX}/5.bw/trans.bed -S ${INDEX}/5.bw/${sample}_R${rep}.bw --startLabel TSS --endLabel TTS -m 5000 -b 2000 -a 2000 -p 4 -o ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz
					
						plotProfile -m ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.mtr.gz -out ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.profile.pdf --dpi 600 --yAxisLabel CPM --outFileSortedRegions ${INDEX}/5.bw/${sample}_R${rep}.trans.ff.regions.txt  --samplesLabel ${sample}_R${rep}
					done
				done
				
		fi
fi






for sample in $(cat ${INDEX}/sample2.txt)
	do
	for ((rep=1;rep<${r}+1;rep++))
	do
		bamCompare -b1 ${INDEX}/3.filter_again/unpaired_read/${sample}_R${rep}.zF.bam -b2 ${INDEX}/3.filter_again/unpaired_read/${sample}_Input_R${rep}.zF.bam -o ${INDEX}/5.bw/${sample}_R${rep}.bw --operation log2 -p 4 --scaleFactorsMethod None --normalizeUsing CPM --binSize 10 --extendReads 200
	done
done





seqtk sample -s100 ANP32E-3Flag_R1_1.fq.gz 500000 | gzip > ../../ANP32E_random/raw/random_R1_1.fq.gz