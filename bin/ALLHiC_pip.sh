#!/bin/bash

usage()
{
	echo "    Usage: `basename $0` -r reference -1 R1.fq -2 R2.fq -k group_count [-e enzyme] [-t threads] [-b bin_size]"
	echo "          -r: reference genome"
	echo "          -1: Lib_R1.fq.gz"
	echo "          -2: Lib_R2.fq.gz"
	echo "          -k: group_count"
	echo "          -e: enzyme_sites (HindIII: AAGCTT; MboI: GATC), default: HindIII"
	echo "          -t: threads, default: 10"
	echo "          -b: bin_size for hic heatmap, can be divided with comma, default: 500k"
	exit 0
}

### get options
while getopts ':r:1:2:k:e:t:b:' OPT; do
	case $OPT in
		r)
			ref="$OPTARG";;
		1)
			R1="$OPTARG";;
		2)
			R2="$OPTARG";;
		e)
			enzyme="$OPTARG";;
		k)
			group_count="$OPTARG";;
		t)
			threads="$OPTARG";;
		b)
			bin_size="$OPTARG";;
		?)
			usage;;
	esac
done

### check required variants
if [ -z $ref ] || [ -z $R1 ] || [ -z $R2 ] || [ -z $group_count ]; then
	usage
fi

### set default values while optional variants were not set
if [ -z $threads ]; then
	threads=10
fi

if [ -z $bin_size ]; then
	bin_size=500k
fi

if [ -z $enzyme ]; then
	enzyme=AAGCTT
fi

enzyme=`echo $enzyme | tr '[a-z]' '[A-Z]'`

if [ $enzyme = HINDIII ]; then
	enzyme=AAGCTT
fi

if [ $enzyme = MBOI ]; then
	enzyme=GATC
fi

### link required files
ln -s ${ref} ./seq.fasta
ln -s ${R1} ./Lib_R1.fastq.gz
ln -s ${R2} ./Lib_R2.fastq.gz

### index reference genome
/public1/user_program/sentieon-genomics-201711/bin/bwa index seq.fasta
samtools faidx seq.fasta


### 1st round of mapping
/public1/user_program/sentieon-genomics-201711/bin/bwa mem -M -t $threads -K 10000000 seq.fasta Lib_R1.fastq.gz Lib_R2.fastq.gz | /public1/user_program/sentieon-genomics-201711/bin/sentieon util sort -r seq.fasta -o sorted.bam -t $threads --sam2bam -i - 
samtools index sorted.bam

### correct contig
ALLHiC_corrector -m sorted.bam -r seq.fasta -o seq.HiCcorrected.fasta -t $threads

### 2nd round of mapping
/public1/user_program/sentieon-genomics-201711/bin/bwa index seq.HiCcorrected.fasta
samtools faidx seq.HiCcorrected.fasta
/public1/user_program/sentieon-genomics-201711/bin/bwa mem -M -t $threads -K 10000000 seq.HiCcorrected.fasta Lib_R1.fastq.gz Lib_R2.fastq.gz | /public1/user_program/sentieon-genomics-201711/bin/sentieon util sort -r seq.HiCcorrected.fasta -o sample.bwa_mem.bam -t $threads --sam2bam -i -

### partition
ALLHiC_partition -r seq.HiCcorrected.fasta -e $enzyme -k $group_count -b sample.bwa_mem.bam

### optimize
rm cmd.list
for((K=1;K<=$group_count;K++));do echo "allhic optimize sample.bwa_mem.counts_$enzyme.${group_count}g${K}.txt sample.bwa_mem.clm" >> cmd.list;done
ParaFly -c cmd.list -CPU $group_count

### check ParaFly success
i=0
while [ -e "FailedCommands" ]
do
    i=$((i+1))
    if [[ $i -gt 100 ]]
    then
        echo "FATAL ERROR, check all files"
        exit -1
    fi  
    mv FailedCommands cmd.list
    ParaFly -c cmd.list -CPU $threads
done

### build
unset PERL5LIB
ALLHiC_build seq.HiCcorrected.fasta

### plot
#perl /public1/user_program/script/getFaLen.pl -i groups.asm.fasta -o len.txt
#grep sample len.txt > chrn.list
#ALLHiC_plot sample.bwa_mem.bam groups.agp chrn.list $bin_size pdf

