#1.reformat UMI
module load python2/2.7.9
for i in `find . -maxdepth 1 -name "*R1.fastq"`; do bsub "python reformat_umi_fastq.py -q 10 -l $i -r ${i/R1/R2} -L $i.reformated  -R ${i/R1/R2}.reformated 2> $i.umi"; done

#2.trim Illumina adapters
module load cutadapt/4.1
for i in $(ls *R1.fastq.reformated); do bsub "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 30 --pair-filter=any -j 40 --overlap 15 -o ${i/R1.fastq.reformated/R1.rfmtd.trmd.fq.gz} -p ${i/R1.fastq.reformated/R2.rfmtd.trmd.fq.gz} $i ${i/R1.fastq.reformated/R2.fastq.reformated} 1> $i.txt"; done

#3.count ERCC reads
module load python/2.7.9
module load gcc/12.2.0
pip install --user pysam
module load samtools/0.1.19
module load bowtie/1.0.0
bowtie-build ercc.fa ercc
for i in `find . -maxdepth 1 -name "*R1.fastq.reformated"`; do bsub "bowtie -S -I 100 -X 600 --fr ercc -1 $i -2 ${i/R1/R2} > $i.ercc.sam; samtools view -Sbt ercc.fa.fai $i.ercc.sam >  $i.ercc.bam; samtools sort $i.ercc.bam $i.ercc.sorted; python umi_mark_duplicates.py -f $i.ercc.sorted.bam -p 4; samtools view -b -F 0x400 $i.ercc.deumi.sorted.bam > $i.ercc.deumi.sorted.bam.dedup; samtools view $i.ercc.deumi.sorted.bam.dedup > $i.ercc.deumi.sorted.bam.dedup.sam; rm $i.ercc.deumi.sorted.bam.dedup $i.ercc.deumi.sorted.bam.bai $i.ercc.deumi.sorted.bam $i.ercc.sorted.bam.bai $i.ercc.sorted.bam $i.ercc.bam $i.ercc.sam"; done
for i in `find . -maxdepth 1 -name "*.ercc.deumi.sorted.bam.dedup.sam"`;  do echo $i; cut -f3 $i | sort | uniq -c |awk '{print $2"\t"$1}' > temp.tmp; awk -v f="$i" '{a+=$2/2}END{print f"\t"a}' temp.tmp >> ercc.reads; done; rm temp.tmp

#4a.run piPpipes to remove rRNA and align with STAR
#mm10 genome was downloaded from https://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/
#mm10 annotation was downloaded from https://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/
for i in `find . -maxdepth 1 -name "*R1.rfmtd.trmd.fq.gz"`; do bsub "piPipes rna -l $i -r ${i/R1/R2} -g mm10 -o $i.pipipes -c 8"; done

#4b.move bam files from piPipes to current directory for deduplicating
for i in $(find . -maxdepth 3 -name "*.mm10.sorted.bam"); do mv $i .; done

#5a.mark duplicates in bam files
module load python2/2.7.9
module load gcc/12.2.0
pip install --user pysam
for i in $(find . -maxdepth 1 -name "*.mm10.sorted.bam"); do bsub  "python2 umi_mark_duplicates.py -f $i -p 8"; done
#5b.deduplicate
module load samtools/0.1.19
for i in $(find . -maxdepth 1 -name "*.deumi.sorted.bam"); do bsub "samtools view -b -F 0x400 $i > $i.dedup"; done
#5c.sort deduplicated files by chrom pos
for i in $(find . -maxdepth 1 -name "*.deumi.sorted.bam.dedup"); do bsub "samtools sort $i $i.sorted"; done


#6.calculate the number of all mapped reads from bam
module load bedtools/2.26.0
for i in $(find . -maxdepth 1 -name "*.bam"); do bsub "bedtools bamtobed -bed12 -tag NH -i $i | awk -v f=$i '{a[\$4]++}END{for (j in a){total++};print f\"\\t\"total/2}' > $i.all.mapped.reads.txt"; done

#7.stringtie quantify transcript TPMs on 3 annotations
#mm10.38.92.gtf was downloaded from https://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/
#mm10rmsk.gtf was downloaded from https://genome.ucsc.edu/cgi-bin/hgTables?&clade=mammal&org=Mouse&db=mm10&hgta_group=varRep&hgta_track=rmsk
module load stringtie/2.2.1
for i in $(ls *.deumi.sorted.bam.dedup.*.bam); do bsub "stringtie $i  --rf -o $i.stringtie_mm10 -b $i.ctab_mm10 -p 20 -G mm10.38.92.gtf -A $i.gene_abund.tab_mm10 -e; stringtie $i  --rf -o $i.stringtie_mm10_rmsk -b $i.ctab_mm10_rmsk -p 20 -G mm10.rmsk.gtf -A $i.gene_abund.tab_mm10_rmsk -e; stringtie $i  --rf -o $i.stringtie_cluster -b $i.ctab_cluster -p 20 -G mm10.cluster.corrected.gtf -A $i.gene_abund.tab_cluster -e"; done

#8.HTSeq quantify for DESeq2
module load samtools/0.1.19
module load python2/2.7.9
module load htseq/0.6.1
pip install --user pysam
for i in $(find . -maxdepth 1 -name "*.bam"); do bsub "htseq-count -f bam -r pos -s reverse -t exon -i gene_id -m union $i mm10.38.92.gtf > $i.sortedbynameGENE.htseq"; done
for i in $(find . -maxdepth 1 -name "*.bam"); do bsub "htseq-count -f bam -r pos -s reverse -t exon -i gene_id -m union $i mm10.rmsk.gtf > $i.RMSK.htseq"; done

#9a.uniquely mapping reads normalized coverage
module load bedtools/2.26.0
for i in $(find . -name "*.bam"); do bsub "bedtools bamtobed -bed12 -tag NH -i $i | awk 'BEGIN{FS=OFS=\"\\t\"}(\$5==1){if(substr(\$4,length(\$4))==1){\$6=(\$6==\"+\"?\"-\":\"+\")}; print \$0}' > $i.unique.bed12; normscale=\$(awk '{a[\$4]++}END{for (j in a){total++};print 2000000/total}' $i.unique.bed12); bedtools genomecov -scale \$normscale -split -bg -strand + -i $i.unique.bed12 -g mm10.pipipes.genome > $i.plus; bedtools genomecov -scale \$normscale -split -bg -strand - -i $i.unique.bed12 -g mm10.pipipes.genome > $i.temp; awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t-\"\$4}' $i.temp > $i.minus ;cat $i.plus $i.minus | sort -k1,1 -k2,2n > $i.bedgraph; rm $i.minus $i.plus $i.temp $i.unique.bed12"; done

#9b.convert bedgraphs to bigWigs
for i in $(find . -maxdepth 1 -name "*.bedgraph"); do bsub "awk '(\$4>0)' $i > $i.plus; awk '(\$4<0)' $i > $i.minus; bedGraphToBigWig $i.plus mm10.pipipes.genome $i.plus.bigWig; bedGraphToBigWig $i.minus mm10.pipipes.genome $i.minus.bigWig; rm $i.plus $i.minus"; done
