#1.trim the adapter
module load cutadapt/4.1
for i in $(ls *.fastq); do bsub "cutadapt -a TGGAATTCTCGGGTGCCAAGG  -m 1 --overlap 15 -o ${i/.fastq/.trimmed.fq} $i  1> $i.txt"; done

#2.filter at q20p100 (Phred >=20 at all positions) 
for i in $(ls *.trimmed.fq); do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){phred[\$1]=\$2}(FNR<NR)&&(FNR%4==1){line1=\$0;getline;line2=\$0; getline;line3=\$0;getline; lowqual=0;for(i=1;i<=length(\$1);i++){if((phred[substr(\$1,i,1)]==\"\")||(phred[substr(\$1,i,1)]<30)){lowqual++}};if(lowqual==0){print line1;print line2;print line3;print}}' phred.table $i > $i.q20p100"; done

#3.UMI:remove duplicates and anything <18nt
for i in $(find . -maxdepth 1 -name "*.trimmed.fq"); do bsub "awk '(NR%4==1){name=\$1}(NR%4==2){total++;if((substr(\$1,length(\$1)-11,3)==\"GTC\")&&(substr(\$1,length(\$1)-5,3)==\"TAG\")&&(((substr(\$1,4,3)==\"CGA\")&&(substr(\$1,10,3)==\"TAC\"))||((substr(\$1,4,3)==\"ATC\")&&(substr(\$1,10,3)==\"AGT\")))){umis++;if(a[\$1]!=1){nondup++;a[\$1]=1;if(length(\$1)>47){longer18++; print name; print substr (\$1,16,length (\$1)-30);getline; print; getline; print substr (\$1,16,length (\$1)-30)}}}}END{print FILENAME\"\\t\"total\"\\t\"umis\"\\t\"nondup\"\\t\"longer18 > FILENAME\".dup\"}' $i > $i.deUMI.dedup.fq" ; done 

#4.remove rRNA with 1 mismatch allowed
module load bowtie/1.3.1
bowtie-build rRNA.fa rRNA
for i in $(find . -name "*.deUMI.dedup.fq"); do bsub "bowtie --un $i.x_filter_1mm -k 1 -v 1 rRNA $i > /dev/null" ; done 

#5a.extract spikein counts
module load bowtie/1.3.1
bowtie-build 9oligo.set.fa 9oligo.set
for i in $(find . -name "*.x_filter_1mm"); do bsub "bowtie --norc --un $i.x_spikein.fq -v 0 9oligo.set $i | awk 'BEGIN{FS=OFS=\"\\t\"}(substr(\$3,index(\$3,\"-\")+1)==length(\$5)){spike[\$3]++}END{for (i in spike){print i\"\\t\"spike[i]}}' > $i.spikein" ; done 
#5b.merge spikein count data
gawk '{a[$1][FILENAME]=$2}END{printf "name";for (i in a["114-26"]){printf "\t"i};printf("\n");for (r in a){printf r; for (f in a["114-26"]){if (a[r][f]==""){printf "\t0"}else{printf "\t"a[r][f]}};printf ("\n")}}' *.spikein > spikein.txt

#6.analyze data with Tailor
for i in `find . -name "*.x_spikein.fq"`; do bsub "Tailor/run_tailing_pipeline.sh -q 20 -T 10 -t Tailor/annotation/mm10.genomic_features -i $i -c 8 -g Tailor/indexes/mm10.fa -o $i.Tailor" ; done

#7.convert bed2 files from Tailor output to rpm files; columns mean the following:
#1-coordinates, for multimappers: first occurence only
#2-strand
#3-templated sequence
#4-rpm (total counts normalized by seq depth)
#5-rpm for reads with nontemplated nucleotides
for i in $(find . -name "*.p20.bed2"); do  bsub "awk 'BEGIN{FS=OFS=\"\\t\"}((\$9==0)&&(dealt[\$7]!=1)){nontailed[substr (\$7, 1, length (\$7)-\$9)]+=\$4; all[substr (\$7, 1, length (\$7)-\$9)]=1;coordinates[substr (\$7, 1, length (\$7)-\$9)]=\$1\":\"\$2\"-\"\$3\"\\t\"\$6; alluniq+=\$4; dealt[\$7]=1}((\$9>0)&&(dealt[\$7]!=1)){tailed[substr (\$7, 1, length (\$7)-\$9)]+=\$4;all[substr (\$7, 1, length (\$7)-\$9)]=1 ;coordinates[substr (\$7, 1, length (\$7)-\$9)]=\$1\":\"\$2\"-\"\$3\"\\t\"\$6; alluniq+=\$4; dealt[\$7]=1}END{print \"coordinates\\tstrand\\tsequenceASis\\ttotalRPM\\ttailedRPM\"; for (r in all){print coordinates[r]\"\\t\"r\"\\t\"(tailed[r]+nontailed[r])*1000000/alluniq\"\\t\"tailed[r]*1000000/alluniq}}' $i |sort -k5 -g -r > $i.rpm"; done

#8.merge reads with the same prefix 25nt-55nt (output is longest sequence followed by nt18-nt55,t18-t55,untailed,tailed,total)
for i in $(find . -maxdepth 1 -name "*.p20.bed2.rpm"); do bsub "sort -k3,3 $i > $i.sorted; gawk 'BEGIN{ch=1;started=0}(length (\$3)>24)&&(length (\$3)<56)&&(started==0){for (i=18;i<=55;i++){tempnt[i]=0;
tempt[i]=0};prefix=\$3;started=1;tempnt[length(prefix)]+=\$4;
tempt[length(prefix)]+=\$5;}
(length (\$3)>24)&&(length (\$3)<56)&&(started==1){while ((length (\$3)<56)&&(length (\$3)>24)&&(substr(\$3,1,length(prefix))==prefix)&&(ch!=0)){prefix=\$3;tempnt[length(prefix)]+=\$4;tempt[length(prefix)]+=\$5;ch=getline};
if((length(prefix)>24)&&(length(prefix)<56)){for (i=18;i<=55;i++){nt[prefix,i]=tempnt[i];
t[prefix,i]=tempt[i]};
all[prefix]=1};
for (i=18;
i<=55;
i++){tempnt[i]=0;
tempt[i]=0};
tempnt[length(\$3)]+=\$4;
tempt[length(\$3)]+=\$5;
prefix=\$3}END{for (r in all){totalnt=0;totalt=0;printf (r);
for (i=18;
i<=55;
i++){printf (\"\\t\");
printf ((nt[r,i]-t[r,i])); totalnt+=(nt[r,i]-t[r,i])};
for (i=18;
i<=55;
i++){printf (\"\\t\");
printf (t[r,i]); totalt+=t[r,i]};
printf (\"\\t\");
printf (totalnt);
printf (\"\\t\");
printf (totalt);
printf (\"\\t\");
printf (totalt+totalnt);
printf (\"\\n\")}}' $i.sorted > $i.2555.new; rm $i.sorted"; done

#9.identify the modal length of each piRNA and align to mouse genome
module load bowtie/1.3.1
for i in $(find . -name "*.bed2.rpm.2555.new"); do bsub "awk '(length(\$1)>25)&&(\$80>1){max1=0; max1n=0; for (i=8;i<=39;i++){if ((\$i+\$(i+38))>max1){max1=(\$i+\$(i+38)); max1n=i+16}};len[\$1]=max1n; ab[\$1]=\$80}END{for (i in ab){print \">\"i\";\"ab[i]\";\"len[i]; print i}}' $i > $i.fa; bowtie -f -v 0 -k 5 /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/mm10_set/mm10.pipipes $i.fa | awk '{print \$3\"\\t\"\$4\"\\t\"(\$4+length(\$5))\"\\t\"\$1\"\\t\"\$7+1\"\\t\"\$2}' > $i.bowtie; awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){places[\$4]++}(FNR<NR){\$5=places[\$4];if(printed[\$4]==\"\"){print;printed[\$4]=1}}' $i.bowtie $i.bowtie > $i.25.all1.realbed2;  rm $i.fa $i.bowtie"; done

#10a.merge data by genotype 
awk 'BEGIN{FS=OFS="\t"}($6=="+"){split($4,fields,";");abplus[$1"\t"$2"\t"$6]=abplus[$1"\t"$2"\t"$6]";"fields[2];totalabplus[$1"\t"$2"\t"$6]+=fields[2]; seqplus[$1"\t"$2"\t"$6]=seqplus[$1"\t"$2"\t"$6]";"fields[1]":"fields[3];totalnumplus[$1"\t"$2"\t"$6]++;infoplus[$1"\t"$2"\t"$6]=$1"\t"$2"\t"$2+1}($6=="-"){split($4,fields,";");abminus[$1"\t"$3"\t"$6]=abminus[$1"\t"$3"\t"$6]";"fields[2];totalabminus[$1"\t"$3"\t"$6]+=fields[2]; totalnumminus[$1"\t"$3"\t"$6]++;infominus[$1"\t"$3"\t"$6]=$1"\t"$3-1"\t"$3;seqminus[$1"\t"$3"\t"$6]=seqminus[$1"\t"$3"\t"$6]";"fields[1]":"fields[3]}END{for(plus in totalnumplus){if(totalnumplus[plus]==6){print infoplus[plus]"\t"totalabplus[plus]/6"\t"abplus[plus]"\t+\t"seqplus[plus]}};for(minus in totalnumminus){if(totalnumminus[minus]==6){print infominus[minus]"\t"totalabminus[minus]/6"\t"abminus[minus]"\t-\t"seqminus[minus]}}}' *_p17mut_*.all1.realbed2 > p17.all1.merged
awk 'BEGIN{FS=OFS="\t"}($6=="+"){split($4,fields,";");abplus[$1"\t"$2"\t"$6]=abplus[$1"\t"$2"\t"$6]";"fields[2];totalabplus[$1"\t"$2"\t"$6]+=fields[2]; seqplus[$1"\t"$2"\t"$6]=seqplus[$1"\t"$2"\t"$6]";"fields[1]":"fields[3];totalnumplus[$1"\t"$2"\t"$6]++;infoplus[$1"\t"$2"\t"$6]=$1"\t"$2"\t"$2+1}($6=="-"){split($4,fields,";");abminus[$1"\t"$3"\t"$6]=abminus[$1"\t"$3"\t"$6]";"fields[2];totalabminus[$1"\t"$3"\t"$6]+=fields[2]; totalnumminus[$1"\t"$3"\t"$6]++;infominus[$1"\t"$3"\t"$6]=$1"\t"$3-1"\t"$3;seqminus[$1"\t"$3"\t"$6]=seqminus[$1"\t"$3"\t"$6]";"fields[1]":"fields[3]}END{for(plus in totalnumplus){if(totalnumplus[plus]==7){print infoplus[plus]"\t"totalabplus[plus]/7"\t"abplus[plus]"\t+\t"seqplus[plus]}};for(minus in totalnumminus){if(totalnumminus[minus]==7){print infominus[minus]"\t"totalabminus[minus]/7"\t"abminus[minus]"\t-\t"seqminus[minus]}}}' *_p9mut_*.all1.realbed2 > p9.all1.merged
awk 'BEGIN{FS=OFS="\t"}($6=="+"){split($4,fields,";");abplus[$1"\t"$2"\t"$6]=abplus[$1"\t"$2"\t"$6]";"fields[2];totalabplus[$1"\t"$2"\t"$6]+=fields[2]; seqplus[$1"\t"$2"\t"$6]=seqplus[$1"\t"$2"\t"$6]";"fields[1]":"fields[3];totalnumplus[$1"\t"$2"\t"$6]++;infoplus[$1"\t"$2"\t"$6]=$1"\t"$2"\t"$2+1}($6=="-"){split($4,fields,";");abminus[$1"\t"$3"\t"$6]=abminus[$1"\t"$3"\t"$6]";"fields[2];totalabminus[$1"\t"$3"\t"$6]+=fields[2]; totalnumminus[$1"\t"$3"\t"$6]++;infominus[$1"\t"$3"\t"$6]=$1"\t"$3-1"\t"$3;seqminus[$1"\t"$3"\t"$6]=seqminus[$1"\t"$3"\t"$6]";"fields[1]":"fields[3]}END{for(plus in totalnumplus){if(totalnumplus[plus]==3){print infoplus[plus]"\t"totalabplus[plus]/3"\t"abplus[plus]"\t+\t"seqplus[plus]}};for(minus in totalnumminus){if(totalnumminus[minus]==3){print infominus[minus]"\t"totalabminus[minus]/3"\t"abminus[minus]"\t-\t"seqminus[minus]}}}' *_p9p17mut_*.all1.realbed2 > p9p17.all1.merged
awk 'BEGIN{FS=OFS="\t"}($6=="+"){split($4,fields,";");abplus[$1"\t"$2"\t"$6]=abplus[$1"\t"$2"\t"$6]";"fields[2];totalabplus[$1"\t"$2"\t"$6]+=fields[2]; seqplus[$1"\t"$2"\t"$6]=seqplus[$1"\t"$2"\t"$6]";"fields[1]":"fields[3];totalnumplus[$1"\t"$2"\t"$6]++;infoplus[$1"\t"$2"\t"$6]=$1"\t"$2"\t"$2+1}($6=="-"){split($4,fields,";");abminus[$1"\t"$3"\t"$6]=abminus[$1"\t"$3"\t"$6]";"fields[2];totalabminus[$1"\t"$3"\t"$6]+=fields[2]; totalnumminus[$1"\t"$3"\t"$6]++;infominus[$1"\t"$3"\t"$6]=$1"\t"$3-1"\t"$3;seqminus[$1"\t"$3"\t"$6]=seqminus[$1"\t"$3"\t"$6]";"fields[1]":"fields[3]}END{for(plus in totalnumplus){if(totalnumplus[plus]==12){print infoplus[plus]"\t"totalabplus[plus]/12"\t"abplus[plus]"\t+\t"seqplus[plus]}};for(minus in totalnumminus){if(totalnumminus[minus]==12){print infominus[minus]"\t"totalabminus[minus]/12"\t"abminus[minus]"\t-\t"seqminus[minus]}}}' *_WT*.all1.realbed2 > wt.all1.merged
awk 'BEGIN{FS=OFS="\t"}($6=="+"){split($4,fields,";");abplus[$1"\t"$2"\t"$6]=abplus[$1"\t"$2"\t"$6]";"fields[2];totalabplus[$1"\t"$2"\t"$6]+=fields[2]; seqplus[$1"\t"$2"\t"$6]=seqplus[$1"\t"$2"\t"$6]";"fields[1]":"fields[3];totalnumplus[$1"\t"$2"\t"$6]++;infoplus[$1"\t"$2"\t"$6]=$1"\t"$2"\t"$2+1}($6=="-"){split($4,fields,";");abminus[$1"\t"$3"\t"$6]=abminus[$1"\t"$3"\t"$6]";"fields[2];totalabminus[$1"\t"$3"\t"$6]+=fields[2]; totalnumminus[$1"\t"$3"\t"$6]++;infominus[$1"\t"$3"\t"$6]=$1"\t"$3-1"\t"$3;seqminus[$1"\t"$3"\t"$6]=seqminus[$1"\t"$3"\t"$6]";"fields[1]":"fields[3]}END{for(plus in totalnumplus){if(totalnumplus[plus]==1){print infoplus[plus]"\t"totalabplus[plus]/1"\t"abplus[plus]"\t+\t"seqplus[plus]}};for(minus in totalnumminus){if(totalnumminus[minus]==1){print infominus[minus]"\t"totalabminus[minus]/1"\t"abminus[minus]"\t-\t"seqminus[minus]}}}' *_p6mut_B*.all1.realbed2 > p6B.all1.merged

#10a.merge all genotypes
awk 'BEGIN{FS=OFS="\t"}(FILENAME=="wt.all1.merged"){wt[$1"\t"$2"\t"$3"\t"$6]=$0}(FILENAME=="p9.all1.merged"){p9[$1"\t"$2"\t"$3"\t"$6]=$4}(FILENAME=="p17.all1.merged"){p17[$1"\t"$2"\t"$3"\t"$6]=$4}(FILENAME=="p9p17.all1.merged"){p9p17[$1"\t"$2"\t"$3"\t"$6]=$4}(FILENAME=="p6B.all1.merged"){p6[$1"\t"$2"\t"$3"\t"$6]=$4}END{for(pirna in wt){if(p9[pirna]==""){p9[pirna]=0};if(p17[pirna]==""){p17[pirna]=0};if(p9p17[pirna]==""){p9p17[pirna]=0};if(p6[pirna]==""){p6[pirna]=0};print wt[pirna]"\t"p9[pirna]"\t"p17[pirna]"\t"p9p17[pirna]"\t"p6[pirna]}}' wt.all1.merged p9.all1.merged p17.all1.merged p9p17.all1.merged p6B.all1.merged > merged.all1.piRNAs

#11a.add cluster info and remove piRNAs a/s to clusters 
module load bedtools/2.26.0
for i in $(ls merged.all1.piRNAs); do bedtools intersect -s -wao -a $i -b piRNA.cluster.pachytene.bed6 | awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11]=substr($15,1,length($15)-2)}END{for (i in a){if(a[i]!=""){print i"\t"a[i]}else{print i"\t."}}}' > $i.cl; done

#11b.add gene info
module load bedtools/2.26.0
for i in $(ls merged.all1.piRNAs.cl); do bedtools intersect -s -wao -a $i -b mm10.38.92.gtf.exons.geneNnameNtype.bed | awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12]=$16}END{for (i in a){if(a[i]!=""){print i"\t"a[i]}else{print i"\t."}}}' > $i.g; done

#11c.add transposon info
module load bedtools/2.26.0
for i in $(ls merged.all1.piRNAs.cl.g); do bedtools intersect -wao -a $i -b mm10.rmsk.bed | awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13]=$17"\t"$19}END{for (i in a){if(a[i]!=""){print i"\t"a[i]}else{print i"\t.\t."}}}' > $i.r; done
1-chr
2-start
3-end
4-mean ab in wildtype
5-wt abundance for each replicate
6-strand
7-all seq and max length
8-p9 mean abundance
9-p17 mean abundance
10-p9p17 mean abundance
11-p6 mean abundance
12-cluster
13-GENE
14-TE
15-TE strand

#get piRNAs that are un in mutants to identify pi6, pi9, pi17-dependent
for i in $(ls merged.all1.piRNAs.cl.g.r); do awk 'BEGIN{FS=OFS="\t"}($8<0.1){print > FILENAME".p9"}($9<0.1){print > FILENAME".p17"}($10<0.1){print > FILENAME".p9p17"}($11<0.1){print > FILENAME".p6"}' $i; done

#convert merged.all1.piRNAs.cl.g.r and into the 12.inter.2555.1ppm.realbed.piRNAs format
for i in $(ls merged.all1.piRNAs.cl.g.r); do awk 'BEGIN{FS=OFS="\t"}{print $1"\t"$2"\t"$3"\t"substr($7,2,30)"\t"$4"\t"$6"\t"$12"\t"$13"\t"$14"\t"$15"\t"$8"\t"$9"\t"$10"\t"$11}' $i > $i.rfmtd; done

#OUTPUT by COLUMN
1-chr
2-start
3-end
4-seq (25nt prefix)
5-mean wt abundance
6-strand
7-cluster
8-gene
9-TE
10-TE strand
11-p9 mean abundance
12-p17 mean abundance
13-p9p17 mean abundance
14-p6 mean abundance
