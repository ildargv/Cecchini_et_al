#1.trim the adapter
module load fastx_toolkit/0.0.14
for i in $(ls *.fastq); do bsub "fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 15 -c -v -i $i > $i.trimmed"; done

#2.remove duplicates and anything <18nt from all READs
for i in $(find . -maxdepth 1 -name "*.trimmed"); do bsub "awk '(NR%4==1){name=\$1}(NR%4==2){total++;if((substr(\$1,length(\$1)-11,3)==\"GTC\")&&(substr(\$1,length(\$1)-5,3)==\"TAG\")&&(((substr(\$1,4,3)==\"CGA\")&&(substr(\$1,10,3)==\"TAC\"))||((substr(\$1,4,3)==\"ATC\")&&(substr(\$1,10,3)==\"AGT\")))){umis++;if(a[\$1]!=1){nondup++;a[\$1]=1;if(length(\$1)>47){longer18++; print name; print substr (\$1,16,length (\$1)-30);getline; print; getline; print substr (\$1,16,length (\$1)-30)}}}}END{print FILENAME\"\\t\"total\"\\t\"umis\"\\t\"nondup\"\\t\"longer18 > FILENAME\".dup\"}' $i > $i.deUMI.dedup.fq" ; done 

#3.remove rRNA reads with 1 mismatch allowed
module load bowtie/1.0.0
bowtie-build rRNA.fa rRNA
for i in $(find . -name "*.fq"); do bsub "bowtie --un $i.x_filter_1mm -k 1 -v 1 rRNA $i > $i.rRNA" ; done 

#4.keep 28-32nt reads
for i in $(find . -name "*.x_filter_1mm"); do bsub "awk '(FNR%4==1){name=\$0; getline; if((length(\$1)>=28)&&(length(\$1)<=32)){print name; print; getline; print; getline; print}}' $i > $i.2832"; done

#5.align with STAR with 1 mismatch allowed
#mm10 genome was downloaded from https://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/
#mm10 annotation was downloaded from https://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/
module load star/2.5.3a
for i in $(find . -maxdepth 1 -name "*.x_filter_1mm.2832"); do bsub "STAR \
 --runMode alignReads \
 --limitOutSAMoneReadBytes 1000000 \
 --genomeDir ./mm10/STARIndex/ \
 --readFilesCommand cat \
 --readFilesIn $i \
 --runThreadN 8 \
 --outFilterScoreMin 0 \
 --outFilterScoreMinOverLread 0.72 \
 --outFilterMatchNmin 0 \
 --outFilterMatchNminOverLread 0.72 \
 --outFilterMultimapScoreRange 1 \
 --outFilterMultimapNmax -1 \
 --outFilterMismatchNmax 1 \
 --alignIntronMax 0 \
 --alignIntronMin 21 \
 --outFilterIntronMotifs None \
 --genomeLoad NoSharedMemory \
 --outFileNamePrefix $i.star. \
 --outSAMunmapped None \
 --outReadsUnmapped Fastx \
 --outSJfilterReads Unique \
 --seedSearchStartLmax 10 \
 --seedSearchStartLmaxOverLread 1.0 \
 --chimSegmentMin 0 2>&1 1> $i.STAR.log;
 InputReads=\`grep 'Number of input reads' $i.star.Log.final.out | awk '{print \$NF}'\`;
 UniquReads=\`grep 'Uniquely mapped reads number' $i.star.Log.final.out | awk '{print \$NF}'\`;
 MultiReads=\`grep 'Number of reads mapped to multiple loci' $i.star.Log.final.out | awk '{print \$NF}'\`;
 AllMapReads=\$((UniquReads+MultiReads));
 UnMapReads=\$((InputReads-UniquReads-MultiReads));
 echo -e \"genomic_mapper_reads:\t\${AllMapReads}\" >> $i.log;
 echo -e \"genomic_unique_mapper_reads:\t\${UniquReads}\" >> $i.log;
 echo -e \"genomic_multiple_mapper_reads:\t\${MultiReads}\" >> $i.log;
 echo -e \"genomic_unmappable_reads:\t\${UnMapReads}\" >> $i.log; rm $i.star.SJ.out.tab $i.star.Log.progress.out $i.star.Log.out $i.STAR.log $i.star.Unmapped.out.mate1; mv $i.star.Aligned.out.sam $i.STAR.sam; rm -r $i.star._STARtmp"; done

#6.processing sam file
module load samtools/0.1.19
for i in $(find . -maxdepth 1 -name "*.STAR.sam"); do bsub "samtools view -bS $i > $i.bam; samtools sort $i.bam $i.bam.sorted; samtools index $i.bam.sorted.bam ; rm -rf  $i.bam"; done

#7.stringtie quantitate transcript TPMs
#mm10.38.92.gtf was downloaded from https://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/
module load stringtie/1.3.4
for i in $(find . -maxdepth 1 -name "*.bam.sorted.bam"); do bsub "stringtie $i  --fr -o $i.stringtie_mm10 -b $i.ctab_mm10 -p 20 -G mm10.38.92.gtf -A $i.gene_abund.tab_mm10 -e"; done

#8a.coverage for uniquely mapping reads
module load bedtools/2.26.0
for i in $(find . -name "*.bam.sorted.bam"); do bsub "bedtools bamtobed -bed12 -tag NH -i $i | awk 'BEGIN{FS=OFS=\"\\t\"}(\$5==1)' > $i.unique.bed12; normscale=\$(awk '{a[\$4]++}END{for (j in a){total++};print 1000000/total}' $i.unique.bed12); bedtools genomecov -scale \$normscale -split -bg -strand + -i $i.unique.bed12 -g mm10.pipipes.genome > $i.plus; bedtools genomecov -scale \$normscale -split -bg -strand - -i $i.unique.bed12 -g mm10.pipipes.genome > $i.temp; awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t-\"\$4}' $i.temp > $i.minus ;cat $i.plus $i.minus | sort -k1,1 -k2,2n > $i.bedgraph; rm $i.minus $i.plus $i.temp $i.unique.bed12"; done
#8b.convert bedgraphs to bigWigs
for i in $(find . -maxdepth 1 -name "*.bam.bedgraph"); do bsub "awk '(\$4>0)' $i > $i.plus; awk '(\$4<0)' $i > $i.minus; bedGraphToBigWig $i.plus mm10.pipipes.genome $i.plus.bigWig; bedGraphToBigWig $i.minus mm10.pipipes.genome $i.minus.bigWig; rm $i.plus $i.minus"; done

#9a.coverage 5ends only
module load bedtools/2.26.0
for file in $(find . -maxdepth 1 -name "*.sorted.bam"); do bsub "bedtools bamtobed -bed12 -tag NH -i $file | awk 'BEGIN{FS=OFS=\"\\t\"}(\$5==1){\$4=1;print}' > $file.bed12"; done
#9b.5-to-5 distance
for i in $(find . -name "*.bed12"); do bsub "awk '{if(\$6==\"+\"){aplus[\$1,\$2]++}else{aminus[\$1,\$3]++}}
END{for (ij in aplus) {split(ij,indices,SUBSEP); i=indices[1]; j=indices[2]; for (k=1;k<=100;k++){res[k]+=aplus[i,j+k]}}; for (ij in aminus) {split(ij,indices,SUBSEP); i=indices[1]; j=indices[2]; for (k=1;k<=100;k++){res[k]+=aminus[i,j-k]}}; for (i=1;i<=100;i++){print i\"\\t\"res[i]}}' $i > $i.5to5_same.uniq_reads.txt"; done;
#9c.collapse 5ends
for i in $(find . -maxdepth 1  -name "*.bed12"); do bsub "awk '{if (\$6==\"+\"){aplus[\$1,\$2]++}else{aminus[\$1,\$3]++}}END{for (ij in aplus) {split(ij,indices,SUBSEP); i=indices[1]; j=indices[2]; print i\"\\t\"j\"\\t\"j+1\"\\t\"aplus[i,j]\"\\tna\\t+\"};for (ij in aminus) {split(ij,indices,SUBSEP); i=indices[1]; j=indices[2]; print i\"\\t\"j-1\"\\t\"j\"\\t\"aminus[i,j]\"\\tna\\t-\"}}' $i > $i.1"; done
#9d.normalize by sequencing depth
for i in $(find . -name "*.1"); do bsub "awk '{a[FNR]=\$0; total+=\$4; num=FNR}END{for (i=1;i<=num;i++){split (a[i],b,\"\\t\");print b[1]\"\\t\"b[2]\"\\t\"b[3]\"\\t\"1000000*b[4]/total\"\\t\"b[4]\"\\t\"b[6]}}' $i > $i.rpm"; done
#9e.make bedgraph for 5'ends
for i in $(find . -maxdepth 1 -name "*.rpm"); do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(\$6==\"+\"){print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4}(\$6==\"-\"){print \$1\"\\t\"\$2\"\\t\"\$3\"\\t-\"\$4}' $i > $i.bedgraph"; done
#9f.convert bedgraphs to bigWig
for i in $(find . -maxdepth 1 -name "*.bedgraph"); do bsub "sort -k1,1 -k2,2n $i > $i.srt; awk '(\$4>0)' $i.srt > $i.plus; awk '(\$4<0)' $i.srt > $i.minus; bedGraphToBigWig $i.plus mm10.pipipes.genome $i.plus.bigWig; bedGraphToBigWig $i.minus mm10.pipipes.genome $i.minus.bigWig; rm $i.plus $i.minus; rm $i.srt"; done

#10.HTSeq to quantify reads in CDS
module load python/2.7.9
module load python/2.7.9_packages/HTSeq/0.6.1/
for i in $(find . -maxdepth 1 -name "*.bam.sorted.bam"); do i=${i:2};  bsub "htseq-count -f bam -r pos -s yes -t CDS -i gene_id -m union $i mm10.38.92.gtf > $i.sortedbynameCDS.htseq"; done