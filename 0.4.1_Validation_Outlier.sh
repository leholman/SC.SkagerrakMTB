#######Validation Outlier (Loxodonta) reads ##########
##1.Extraction of reads from lca file
#Print sample list based on name files
ls *lca.txt > sample.list

#Taxa list file
nano Loxodonta
> save it as > 'taxa.list'

#Create files based on SAMPLE.list and Taxa.list
while read -r line
do
arr=($line)
#if [ "${arr[1]}" = "$(basename $folder)" ]
lib=${arr[0]}
echo $lib
cat taxa.list | parallel -j20 "grep {} $lib | cut -f1,2,3,4,5,6,7 -d: > $lib.{}.readID.txt"
done < sample.list

#Remove lines in file were not found name_genus
wc -l *.readID.txt| awk '$1 == 0' | awk '{print $2}' > rm.list
cat rm.list | parallel -j20 "rm {}"

#Create sum-up file total sequences found per sample
wc -l *.readID.txt| paste > tot_genus_sequences.txt

#Create fastq from readIDs
seqtk subseq LV7001884464-LV7001305844-NM-MD9-2250_S3_L004.merge.filt2.sga4.fq LV7001884464-LV7001305844-NM-MD9-2250_S3_L004.lca.txt.Loxodonta.readID.txt > LV7001884464-LV7001305844-NM-MD9-2250_S3_L004.lca.txt.Loxodonta.fq &
seqtk subseq LV7001884490-LV7001305740-NM-MD9-2450_S2_L004.merge.filt2.sga4.fq LV7001884490-LV7001305740-NM-MD9-2450_S2_L004.lca.txt.Loxodonta.readID.txt > LV7001884490-LV7001305740-NM-MD9-2450_S2_L004.lca.txt.Loxodonta.fq &

#build indexes on genomes
bowtie2-build -f GCF_000001905.1_Loxafr3.0_genomic.fna Elephante_africanus 
bowtie2-build -f Mammuthus_primigenius_assisted_Loxafr3.0_HiC.fasta Mammuthus_primigenius

##2. Part B. Mapping to the references

for file in *.Loxodonta.fq
do
bname=$(basename $file)
echo $bname

for DB in Elephante_africanus    
do
echo Mapping $file against $DB
bowtie2 --threads 5 -x $DB -U $file --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in Mammuthus_primigenius
do
echo Mapping $file against $DB
bowtie2 --threads 5 -x $DB -U $file --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done
done

#sorting bam files and generating bamcov histograms and tables
for file in *.bam
do 
samtools sort -O BAM -o sort_$file $file
./bamcov/bamcov -m -w0 sort_$file | paste > bamcov_hist_$file.txt 
./bamcov/bamcov sort_$file | paste > bamcov_table_$file.txt 
done

##extract only primary alignments from bam file
samtools view -h sort_LV7001884453-LV7000630780-NM-MD9-2650_S1_L004.lca.txt.Loxodonta.fq.Mammuthus_primigenius.bam | grep -v 'XS:i:'  > Loxodonta_Mammuthus_2650_unique.bam
samtools view -h sort_LV7001884453-LV7000630780-NM-MD9-2650_S1_L004.lca.txt.Loxodonta.fq.Elephante_africanus.bam | grep -v 'XS:i:'  > Loxodonta_Elephant_2650_unique.bam 

#extract NM tag and read length (I am CHECKING this command)

