#1.reformat UMIs
for i in `find . -name "*R1.fastq"`; do bsub "python reformat_umi_fastq.py -l $i -r ${i/R1/R2} -L $i.reformated  -R ${i/R1/R2}.reformated"; done

#2a.remove rRNA reads and align with STAR from pipipes and deduplicate
for i in `find . -name "*R1.fastq.reformated"`; do bsub "piPipes rna -l $i -r ${i/R1/R2} -g mm10 -o $i.pipipes -c 8"; done
#2b.mark duplicates
module load python/2.7.9
module load gcc/5.1.0
module load htslib/1.3
module load libcurl/7.37.0
module load samtools/0.1.19
pip install --user pysam
for i in $(find . -name "*.mm10.sorted.bam"); do bsub  "python umi_mark_duplicates.py -f $i -p 4"; done
#2c.remove duplicates
for i in $(find . -name "*.deumi.sorted.bam"); do bsub "samtools view -b -F 0x400 $i > $i.dedup"; done
#2d.sort by read name
for i in $(find . -name "*.deumi.sorted.bam.dedup"); do bsub "samtools sort $i $i.temp"; done

#3a.coverage of uniquely mapping reads
module load bedtools/2.26.0
for i in $(find . -name "*.deumi.sorted.bam.dedup.temp.bam"); do bsub "bedtools bamtobed -bed12 -tag NH -i $i | awk 'BEGIN{FS=OFS=\"\\t\"}(\$5==1){if(substr(\$4,length(\$4))==1){\$6=(\$6==\"+\"?\"-\":\"+\")}; print \$0}' > $i.unique.bed12; normscale=\$(cut -f4 $i.unique.bed12 | sort | uniq | wc -l | awk '{print 2000000/\$1}'); bedtools genomecov -scale \$normscale -split -bg -strand + -i $i.unique.bed12 -g mm10.pipipes.genome > $i.plus; bedtools genomecov -scale \$normscale -split -bg -strand - -i $i.unique.bed12 -g mm10.pipipes.genome > $i.temp; awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t-\"\$4}' $i.temp > $i.minus ;cat $i.plus $i.minus | sort -k1,1 -k2,2n > $i.bedgraph; rm $i.minus $i.plus $i.temp $i.unique.bed12"; done
#3b.convert bedgraphs to bigWigs
for i in $(find . -name "*.bedgraph"); do bsub "awk '(\$4>0)' $i > $i.plus; awk '(\$4<0)' $i > $i.minus; bedGraphToBigWig $i.plus mm10.pipipes.genome $i.plus.bigWig; bedGraphToBigWig $i.minus mm10.pipipes.genome $i.minus.bigWig; rm $i.plus $i.minus"; done

#4.per gene coverage
#mm10.38.92.gtf.gene_boundaries.bed file contains gene boundaries of all ENSMUSG entries in mm10
for i in $(find . -name "*.bedgraph"); do bsub "awk '
(FNR==NR)&&(\$4>0){for(i=\$2;i<\$3;i++){covplus[\$1,i]+=\$4}}
(FNR==NR)&&(\$4<0){for(i=\$3;i>\$2;i--){covminus[\$1,i]+=(-1*\$4)}}
(FNR<NR)&&(\$6==\"+\"){gene=0;gene5=0;gene3=0;
for(i=\$2;i<\$3;i++){gene+=covplus[\$1,i]};
for(i=\$2;i<(\$2+1000);i++){gene5+=covplus[\$1,i]};
for(i=(\$2+1000);i<\$3;i++){gene3+=covplus[\$1,i]};
if ((\$3-\$2-1000)!=0){
print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$6\"\\t\"(\$3-\$2)\"\\t\"\$4\"\\t\"gene\"\\t\"gene*1000/(\$3-\$2)\"\\t\"gene5\"\\t\"gene3\"\\t\"gene3*1000/(\$3-\$2-1000)}else{
print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$6\"\\t\"(\$3-\$2)\"\\t\"\$4\"\\t\"gene\"\\t\"gene*1000/(\$3-\$2)\"\\t\"gene5\"\\t\"gene3\"\\t0\"}}
(FNR<NR)&&(\$6==\"-\"){gene=0;gene5=0;gene3=0;
for(i=\$3;i>\$2;i--){gene+=covminus[\$1,i]};
for(i=\$3;i>(\$3-1000);i--){gene5+=covminus[\$1,i]};
for(i=(\$3-1000);i>\$2;i--){gene3+=covminus[\$1,i]};
if ((\$3-\$2-1000)!=0){
print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$6\"\\t\"(\$3-\$2)\"\\t\"\$4\"\\t\"gene\"\\t\"gene*1000/(\$3-\$2)\"\\t\"gene5\"\\t\"gene3\"\\t\"gene3*1000/(\$3-\$2-1000)}else{
print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$6\"\\t\"(\$3-\$2)\"\\t\"\$4\"\\t\"gene\"\\t\"gene*1000/(\$3-\$2)\"\\t\"gene5\"\\t\"gene3\"\\t0\"}}' $i mm10.38.92.gtf.gene_boundaries.bed > $i.gro_mm10.38.92.pergene"; done
