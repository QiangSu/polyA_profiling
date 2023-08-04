##1. processing preparation. This section creates several directories (folders) using the mkdir command. These directories are used to organize the output files generated during the data processing.

mkdir raw_cutadapt raw_trimmomatic mapped_reads assembled_reads qc_result rseqc_result
path_raw_read=~/raw_read
path_raw_cutadapt=~/raw_cutadapt
path_raw_trimmomatic=~/raw_trimmomatic
path_qc_result=~/qc_result
path_rseqc_result=~/rseqc_result
bam_dir=~/bam_reads

# extract the file name
ls $path_raw_read/*_001.fastq.gz | cut -d '/' -f 6 | cut -d '_' -f 1-2 > $path_raw_read/samplelist.txt
cat $path_raw_read/samplelist.txt

##2. Raw read quality control. FastQC is a quality control tool for high-throughput sequence data. It assesses the quality of the raw reads.
fastqc -t 5 -o $path_qc_result $path_raw_read/*.fastq.gz

#The fastqc command is executed with the -t option to specify the number of threads to use, and the -o option to set the output directory.
multiqc $path_qc_result

##3. cutadapt. Cutadapt is a tool used for adapter trimming in next-generation sequencing data.
for i in $(cat $path_raw_read/samplelist.txt)
#A loop is used to process each sample listed in samplelist.txt
do
        cutadapt -j 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${path_raw_cutadapt}/${i}_cutadapt_R1.fastq.gz -p ${path_raw_cutadapt}/${i}_cutadapt_R2.fastq.gz $path_raw_read/${i}_raw_1.fq.gz $path_raw_read/${i}_raw_2.fq.gz
done


##4.trimmomatic. Trimmomatic is a tool used for quality trimming of next-generation sequencing data.
for i in $(cat $path_raw_read/samplelist.txt)

#Similar to the Cutadapt step, a loop processes each sample from samplelist.txt
do
trimmomatic PE -threads 15 $path_raw_cutadapt/${i}_cutadapt_R1.fastq.gz $path_raw_cutadapt/${i}_cutadapt_R2.fastq.gz -baseout $path_raw_trimmomatic/${i}_cutadapt_trim.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:20
done


##5. hisat2. Hisat2 is a fast and sensitive alignment tool for mapping RNA-seq reads to a reference genome. 
for i in $(cat $path_raw_read/samplelist.txt)
do
        name=$(echo ${i%R1_001.fastq})
        hisat2 -p 60 --dta -x /home/Ref/Hg19/hisat/hg19_tran -1 $path_raw_trimmomatic/${i}_cutadapt_trim_1P.fastq -2 $path_raw_trimmomatic/${i}_cutadapt_trim_2P.fastq -S ${sam_dir}/${i}_cutadapt_trim.sam
done
echo hisat2 done

##6.samtools. Samtools is a suite of tools for manipulating SAM/BAM files.

for i in $(cat $path_raw_read/samplelist.txt)
do
        samtools sort -@ 10 -o ${bam_dir}/${i}_cutadapt_trim.bam ${sam_dir}/${i}_cutadapt_trim.sam
done
echo samtools done

##7. samtools index *.bam. This step creates an index file for each BAM file using the samtools index command.

for i in $(cat $path_raw_read/samplelist.txt)
do
        samtools index -@ 10 ${bam_dir}/${i}_cutadapt_trim.bam
done

##8.gene body coverage(using RSeQC)
for i in $(cat $path_raw_read/samplelist.txt)
do
        geneBody_coverage.py -i ${bam_dir}/${i}_cutadapt_trim.bam -r /home/suqiang/ref/hg19_RefSeq.bed -o ${bam_dir}/${i}_cutadapt_trim_bam_geneBody_coverage
done

##9.1 Retrieve read IDs based on bam for a specific gene
for i in $(cat ./gene_list.txt)
do
        samtools view -F 4 sample.bam ${i} | cut -f 1 >> read_ids.txt
done

##9.2 Retrieve read IDs based on specific sequence of gene for a specific gene
grep -B1 'gene_specific_sequence' cutadapt_trim_1P.fastq | sed 's/1:N:0:ATCTCGTA+ATAGAGAG//g' | awk NR%3==1 > read_ID_list.txt

##10. Extracting the poly(T)-contained R2 reads from read_ids list
for i in $(cat read_ID_list.txt)
do
        grep -A1 $i R2_cutadapt_trim_1P.fastq >> R2_cutadapt_trim_1P.txt
done

##11. Measuring the poly(T) length 
for ((i=10; i<144; i++))
do
        egrep -c A{$i} R2_cutadapt_trim_1P.txt  >> polyT_length-10-144.txt
done




