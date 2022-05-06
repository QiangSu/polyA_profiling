###this code for transcriptome-wide analysis###
#The .csv including all file names of sequencing data should be in the same directory of this script.

cat RNA-seq*.csv | sed "1d" | awk -F"," '{print $1}' > tmp_Sequenced_name.txt
cat RNA-seq*.csv | sed "1d" | awk -F"," '{print $2}' > tmp_Sequenced_path.txt
cat RNA-seq*.csv | sed "1d" | awk -F"," '{print $3}' > tmp_Sample_name.txt
cat RNA-seq*.csv | sed "1d" | awk -F"," '{print $4}' > tmp_Source.txt

##1.directory path preparing 
mkdir raw_cutadapt raw_trimmomatic mapped_reads assembled_reads qc_result rseqc_result

path_raw_read=/home/liangyongjie/practise/raw_read
##path_clean_read=/home/liangyongjie/practise/clean_read
path_qc_result=/home/liangyongjie/practise/qc_result
path_raw_cutadapt=/home/liangyongjie/practise/raw_cutadapt
path_raw_trimmomatic=/home/liangyongjie/practise/raw_trimmomatic
##path_software=/home/liangyongjie/practise/biosoft
path_rseqc_result=/home/liangyongjie/practise/rseqc_result
##index=/home/Ref/mm10/hisat/genome
bam_dir=/home/liangyongjie/practise/assembled_reads
sam_dir=/home/liangyongjie/practise/mapped_reads
##bam_sort=/home/liangyongjie/practise/bam_sort

ls $path_raw_read/*_1.fq.gz | cut -d '_' -f 2 | cut -d '/' -f 2 > $path_raw_read/samplelist.txt
cat $path_raw_read/samplelist.txt

##for i in $(cat $path_raw_read/samplelist.txt)
##do
##	filename=$(basename $i)
####	ln -s $path_raw_read/$filename $path_targetfile/$filename
##done


##2.fastqc- quality control checks on raw_read
fastqc -t 10 -o $path_qc_result $path_raw_read/*.fq.gz
##for i in $(cat samplelist.txt)
##do
##	fastqc -t 10 -o $path_qc_result $path_raw_read/${i}_1.fq.gz $path_raw_read/${i}_2.fq.gz
##done
echo fastqc done 

##3.multiqc-aggregating multiple fastqc reports
multiqc $path_qc_result
echo multiqc done

##4.cutadapt-removing adapter
for i in $(cat $path_raw_read/samplelist.txt)
do
	cutadapt -j 10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${path_raw_cutadapt}/${i}_cutadapt_R1.fastq -p ${path_raw_cutadapt}/${i}_cutadapt_R2.fastq $path_raw_read/${i}_1.fq.gz $path_raw_read/${i}_2.fq.gz
done
echo cut adapter done

##5.trimmomatic-removing low-quality reads
for i in $(cat $path_raw_read/samplelist.txt)
do
	trimmomatic PE -threads 10 $path_raw_cutadapt/${i}_cutadapt_R1.fastq $path_raw_cutadapt/${i}_cutadapt_R2.fastq -baseout $path_raw_trimmomatic/${i}_cutadapt_trim.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:20
done
echo trim done

##6.hisat2-trimmed reads mapping
for i in $(cat $path_raw_read/samplelist.txt)
do
	name=$(echo ${i%R1_001.fastq})
	hisat2 -p 10 --dta -x /home/Ref/rno/hisat/Rat_tran -1 $path_raw_trimmomatic/${i}_cutadapt_trim_1P.fastq -2 $path_raw_trimmomatic/${i}_cutadapt_trim_2P.fastq -S ${sam_dir}/${i}_cutadapt_trim.sam
done
echo hisat2 done

##7.samtools-converting file.sam to file.bam
for i in $(cat $path_raw_read/samplelist.txt)
do
        samtools sort -@ 10 -o ${bam_dir}/${i}_cutadapt_trim.bam ${sam_dir}/${i}_cutadapt_trim.sam
done
echo samtools done

##8.RSeQC-python environment
#cp ${path_software}/RSeQC-4.0.0/scripts/bam_stat.py ${bam_dir}/
#cp ${path_software}/RSeQC-4.0.0/scripts/read_distribution.py ${bam_dir}/
for i in $(cat $path_raw_read/samplelist.txt)
do
	bam_stat.py -i ${bam_dir}/${i}_cutadapt_trim.bam > ${path_rseqc_result}/${i}_bam_stat.log
	read_distribution.py -i ${bam_dir}/${i}_cutadapt_trim.bam -r /home/liangyongjie/practise/rno_ncbiRefSeq.bed > ${path_rseqc_result}/${i}_distribution.log
done
echo RSeQC done

##9.featureCounts
featureCounts -p -t exon -g gene_id -a /home/liangyongjie/practise/Rattus_norvegicus.Rnor_6.0.101.gtf -o sample_counts_gene_id.txt ${bam_dir}/*.bam
echo output gene_id done
#featureCounts -p -t exon -g gene_name -a /home/Ref/Hg19/hisat/Homo_sapiens.GRCh37.75.gtf -o sample_counts_gene_name.txt ${bam_dir}/*.bam
#echo output gene_name done

###This part is for poly(A)-tail profiling of individual gene
##10. specific-gene-identifying-grouping-from-R1
grep -B1 'specific_sequence_adjacent_polyA-tail' file.R1_cutadapt_trim_1P.fastq | sed 's/1:N:0:ATCTCGTA+ATAGAGAG//g' | awk file._R1_cutadapt_trim_1P.modif.txt 

##11. Paired-read-grouping-in-R2
for i in $(cat file.R1_cutadapt_trim_1P.modif.txt)
do
        grep -A1 $i file.R2_cutadapt_trim_1P.fastq >> file.R2_cutadapt_trim_1P.txt
done

##12. Poly(T)-length-determination-in-cummulating-profile

for ((i=10; i<144; i++))
do
egrep -c T{$i} file.R2_cutadapt_trim_1P.txt >> file.R2_cutadapt_trim_1P_polyT_accum.txt
done




