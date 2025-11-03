#1.reformat UMIs
for i in `find . -maxdepth 1 -name "*.R1.fastq"`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR%4==1)&&(FNR==NR){name=\$1;getline;if((index(\$1,\"N\")==0)&&(((substr(\$1,4,3)==\"ATC\")&&(substr(\$1,10,3)==\"AGT\"))||((substr(\$1,4,3)==\"CGA\")&&(substr(\$1,10,3)==\"TAC\")))){a[FNR-1]=1; split(name,b,\" \"); c[FNR-1]=b[1]\"_\"substr(\$1,1,15);print b[1]\"_\"substr(\$1,1,15)\" \"b[2]> FILENAME\".rfmtd.fq\"; print substr(\$1,16)> FILENAME\".rfmtd.fq\"; getline;print > FILENAME\".rfmtd.fq\"; getlineprint substr(\$1,16)> FILENAME\".rfmtd.fq\"}}(FNR<NR)&&(a[FNR]==1)&&(FNR%4==1){split(\$1,b,\" \"); print c[FNR]\" \"b[2] > FILENAME\".rfmtd.fq\" ; getline;print > FILENAME\".rfmtd.fq\";getline;print > FILENAME\".rfmtd.fq\";getline;print > FILENAME\".rfmtd.fq\"}' $i ${i/.R1./.R2.}"; done

#2.run piPpipes to remove rRNA and align with STAR
#mm10 genome was downloaded from https://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/
#mm10 annotation was downloaded from https://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/
for i in `find . -maxdepth 1 -name "*R1.fastq.rfmtd.fq"`; do bsub "piPipes deg -l $i -r ${i/.R1./.R2.} -g mm10 -o $i.pipipes.deg -c 8"; done

#3.collapse 5ends for uniquely mapping reads
for i in $(find . -maxdepth 1  -name "*.unique.bed12"); do bsub "awk '{if (\$6==\"+\"){aplus[\$1,\$2]+=1/\$5}else{aminus[\$1,\$3]+=1/\$5}}END{for (ij in aplus) {split(ij,indices,SUBSEP); i=indices[1]; j=indices[2]; print i\"\\t\"j\"\\t\"j+1\"\\t\"aplus[i,j]\"\\tna\\t+\"};for (ij in aminus) {split(ij,indices,SUBSEP); i=indices[1]; j=indices[2]; print i\"\\t\"j-1\"\\t\"j\"\\t\"aminus[i,j]\"\\tna\\t-\"}}' $i > $i.1"; done

#4.normalize to sequencing depth 
for i in $(find . -name "*.1"); do bsub "awk '{a[FNR]=\$0; total+=\$4; num=FNR}END{for (i=1;i<=num;i++){split (a[i],b,\"\\t\");print b[1]\"\\t\"b[2]\"\\t\"b[3]\"\\t\"1000000*b[4]/total\"\\t\"b[4]\"\\t\"b[6]}}' $i > $i.rpm"; done

#5.select reads outside pachytene clusters
module load bedtools/2.26.0
for i in $(find . -name "*.rpm"); do bsub "bedtools intersect -wa -v -a $i -b piRNA.cluster.pachytene.merged.bed6 > $i.out.pch"; done

#6.merge combinations WT->pi mutants
for i in $(ls ../*WT*.out.pch); do si=${i:3}; si=${si%%.x_rRNA*}; si=${si%%.unique*}; for j in $(ls p9?.x_rRNA.mm10.sorted.f0x40.noS.unique.bed12.1.rpm.out.pch); do sj=${j%%.x_rRNA*}; bsub "awk '(FNR==NR){mutrpm[\$1,\$6,\$2]=\$4;mutcount[\$1,\$6,\$2]=\$5;mutcoord[\$1,\$6,\$2]=\$1\"\\t\"\$2\"\\t\"\$3;mutstrand[\$1,\$6,\$2]=\$6}(FNR<NR){if(mutrpm[\$1,\$6,\$2]==\"\"){print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\",0\\t\"\$5\",0\\t\"\$6}else{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\",\"mutrpm[\$1,\$6,\$2]\"\\t\"\$5\",\"mutcount[\$1,\$6,\$2]\"\\t\"\$6}}' $j $i > $si.$sj.unique.out.m"; done; done

#7.filter at 5+ reads in WT to speed up processing
for i in $(find . -name "*.m"); do si=${i//unique.bed12.1.rpm.out.pch./}; bsub "awk '(FS=OFS=\"\\t\"){split (\$5,a,\",\");if((a[1])>=5){split (\$4,a,\",\");if(a[1]==0){\$4=15}else{if(a[2]==0){\$4=-15}else {\$4=log(a[2]/a[1])/log(2)}}; \$5=\$5\",\"a[1]\",\"a[2]; print}}' $i > $si.5.d"; done

#8.intersect with denovo genes, rmsk, and mm10.38.92.genesNtype
module load bedtools/2.26.0
for i in $(find . -name "*.d"); do bsub "bedtools intersect -s -wao -a $i -b 4N.2samples.merged.stringtie.gtf.gene_merged_exons.bed | awk '{a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6]=\$10}END{for (i in a){print i\"\\t\"a[i]}}' > $i.gdn"; done
#mm10.rmsk.bed was downloaded from https://genome.ucsc.edu/cgi-bin/hgTables?&clade=mammal&org=Mouse&db=mm10&hgta_group=varRep&hgta_track=rmsk
for i in $(find . -name "*.gdn"); do bsub "bedtools intersect -wao -a $i -b mm10.rmsk.bed | awk '{if(\$6==\$13){a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7]=\$11\"\\t+\"}else{a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7]=\$11\"\\t-\"}}END{for (i in a){print i\"\\t\"a[i]}}' > $i.r"; done
for i in $(find . -name "*.r"); do bsub "bedtools intersect -s -wao -a $i -b mm10.38.92.gtf.exons.geneNnameNtype.bed | awk '{if (\$14==-1){b[1]=\".\";b[2]=\".\"}else{split (\$14,b,\":\")}; a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9]=\$13\";\"b[1]\";\"b[2]}END{for (i in a){print i\"\\t\"a[i]}}' > $i.g"; done

#9.get -25+15 sequence for REVERSED and -60+60 sequence for not REVERSED
for i in $(find . -name "*.d.gdn.r.g"); do bsub "awk '(\$6==\"+\"){if((\$2-60)>-1){
print \$1\"\\t\"\$2-25\"\\t\"\$2+15\"\\t\"\$1\";\"\$2\";\"\$3\";\"\$4\";\"\$5\";\"\$6\";\"\$7\";\"\$8\";\"\$9\";\"\$10\"\\t\"\$5\"\\t-\";
print \$1\"\\t\"\$2-60\"\\t\"\$2+60\"\\t\"\$1\";\"\$2\";\"\$3\";\"\$4\";\"\$5\";\"\$6\";\"\$7\";\"\$8\";\"\$9\";\"\$10\"\\t\"\$5\"\\t+\"}}(\$6==\"-\"){if((\$3-60)>-1){
print \$1\"\\t\"\$3-15\"\\t\"\$3+25\"\\t\"\$1\";\"\$2\";\"\$3\";\"\$4\";\"\$5\";\"\$6\";\"\$7\";\"\$8\";\"\$9\";\"\$10\"\\t\"\$5\"\\t+\";
print \$1\"\\t\"\$3-60\"\\t\"\$3+60\"\\t\"\$1\";\"\$2\";\"\$3\";\"\$4\";\"\$5\";\"\$6\";\"\$7\";\"\$8\";\"\$9\";\"\$10\"\\t\"\$5\"\\t-\"}}' $i > $i.200" ; done
module load bedtools/2.26.0
for i in $(find . -name "*.d.gdn.r.g.200"); do si=${i//.d.gdn.r.g.200/}; bsub "bedtools getfasta -s -name -fi mm10.pipipes.fa -bed $i -fo $si.f"; done
for i in $(find . -name "*.f"); do bsub "awk '(FNR%4==1){print; getline; print;getline; getline; print}' $i > $i.both"; done
#OUTPUT by COLUMN
#line1->deg_info::ccordinates
#line2->revcomp -25+15
#line3->sense -60+60

#10.get context sequence from transcripts 
#10a.get mm10.38.92.gtf ref for 0+fpkm transcripts
#ctab files are per transcript FPKM estimates from stringtie for 6 reps of wild-type control
#mm10.38.92.gtf.exon_coord.fasta.fa.as.form file contain antisense sequences of all mm10 transcripts accompanied by exon coordinates
awk 'BEGIN{FS=OFS="\t"}
(FILENAME=="WT1_t_data.ctab")&&(FNR>1){ab37[$6]=$12}
(FILENAME=="WT2_t_data.ctab")&&(FNR>1){ab55[$6]=$12}
(FILENAME=="WT3_t_data.ctab")&&(FNR>1){ab56[$6]=$12}
(FILENAME=="WT4_reprep_id1_t_data.ctab")&&(FNR>1){ab61[$6]=$12}
(FILENAME=="WT5_reprep_rep2_id10_t_data.ctab")&&(FNR>1){ab62[$6]=$12}
(FILENAME=="WT6_reprep_rep2_id11_t_data.ctab")&&(FNR>1){ab63[$6]=$12}
(FILENAME=="mm10.38.92.gtf.exons.geneNnameNtype"){type[$1]=$3","$2}
(FILENAME=="mm10.38.92.gtf.exon_coord.fasta.fa.as.form")&&(FNR%2==1){split(substr($1,2),gname,";");name=$1;getline; print name";"(ab37[gname[1]]+ab55[gname[1]]+ab56[gname[1]]+ab61[gname[1]]+ab62[gname[1]]+ab63[gname[1]])/6";"type[gname[2]]";"ab37[gname[1]]":"ab55[gname[1]]":"ab56[gname[1]]":"ab61[gname[1]]":"ab62[gname[1]]":"ab63[gname[1]];print $0}' *.ctab mm10.38.92.gtf.exons.geneNnameNtype mm10.38.92.gtf.exon_coord.fasta.fa.as.form > mm10.38.92.gtf.exon_coord.fasta.fa.as.form.WT1_2_3_4_5_6_st_0fpkm
#OUTPUT by COLUMN
#line1-ENSMUST
#line2-ENSMUSG
#line3-refname
#line4-chr
#line5-strand
#line6-exons (1-2,3-4,5-6) 1-based on genomic strand
#line7-segs (1-2,3-4,5-6) 1-based on transcribed strand (reverse order for - genomic strand)
#line8-cds start 1-based on transcribed strand (reverse order for - genomic strand)
#line9-cds end 1-based on transcribed strand (reverse order for - genomic strand)
#line10-ENSMUST MEAN fpkm abundance in WT1/2/3/4/5/6 by stringtie
#line11-refname,type
#line12-ENSMUST fpkm abundance in each WT1/2/3/4/5/6 by stringtie

awk 'BEGIN{FS=OFS="\ "}(FNR==NR)&&(FNR%2==1){name=substr($1,2);getline;a[name]=$1}(FNR<NR)&&(FNR%2==1){name=substr($1,2,index($1,";")-2);print;getline;print;print a[name]}' mm10.38.92.gtf.exon_coord.fasta.fa mm10.38.92.gtf.exon_coord.fasta.fa.as.form.WT1_2_3_4_5_6_st_0fpkm > mm10.38.92.gtf.exon_coord.fasta.fa.asNs.form.WT1_2_3_4_5_6_st_0fpkm


#10b.split f into 60,000 lines per file for parallel computing
#only four *.m.5.f.both files are use here to generate data for four reps of wt controls only
for i in $(find . -name "*.m.5.f.both"); do echo $i; split -a 4 -d -l 60000 $i $i.xx. ; done 
for i in $(find . -maxdepth 1 -name "*.xx.*"); do si=${i//.xx./.xy.}; bsub "awk  'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%3==1){split(substr(\$1,2),dat,\";\");
 split(dat[6],exons,\",\");
 m=split(dat[7],segs,\",\");
 for (i=1;i<=m;i++){split(exons[i],b,\"-\");
 excoor[dat[1],i,1]=b[1]-1;
 excoor[dat[1],i,2]=b[2];
 split(segs[i],b,\"-\");
 segcoor[dat[1],i,1]=b[1]-1;
 segcoor[dat[1],i,2]=b[2]};
 chr[dat[1]]=dat[4];
 strand[dat[1]]=dat[5];
 cdsst[dat[1]]=dat[8]-1;
 cdsend[dat[1]]=dat[9];
 segnum[dat[1]]=m;
 st[dat[1]]=excoor[dat[1],1,1];
 end[dat[1]]= excoor[dat[1],m,2];
 fpkm[dat[1]]= dat[10];
 getline;seq[dat[1]]=toupper(\$1);
 getline;seqse[dat[1]]=toupper(\$1)}(FNR<NR)&&(FNR%3==1){tr=\"\";ndeg=split(substr(\$1,2,index(\$1,\"::\")-2),deg,\";\");
 
 if(deg[6]==\"+\"){for (i in seq){if ((chr[i]==deg[1])&&(deg[6]==strand[i])&&(st[i]+25<=deg[2])&&(deg[2]<=end[i]-15)){ seg=0;for(j=1;j<=segnum[i];j++){if((excoor[i,j,1]<=deg[2])&&(deg[2]<excoor[i,j,2])){seg=j}};if(seg>0){utr=\"\";loc=(segcoor[i,seg,1]+deg[2]-excoor[i,seg,1]);if(cdsend[i]>0){if(loc<cdsst[i]){utr=\"5utr\"};if(cdsend[i]<loc){utr=\"3utr\"};if((loc>=cdsst[i])&&(cdsend[i]>=loc)){utr=\"cds\"}}else{utr=\"nc\"}; if(loc<100){
 tr=tr\";\"i\";\"loc\";\"utr\";\"fpkm[i]\";\"substr(seq[i],segcoor[i,segnum[i],2]-loc-15+1,40)\";\"substr(seqse[i],1,loc+100)\";\"loc\";\"cdsst[i]\";\"cdsend[i]\";\"segcoor[i,segnum[i],2]}else{
 tr=tr\";\"i\";\"loc\";\"utr\";\"fpkm[i]\";\"substr(seq[i],segcoor[i,segnum[i],2]-loc-15+1,40)\";\"substr(seqse[i],loc-100+1,200)\";\"100\";\"cdsst[i]\";\"cdsend[i]\";\"segcoor[i,segnum[i],2]}}}}};
 
 if(deg[6]==\"-\"){for (i in seq){if ((chr[i]==deg[1])&&(deg[6]==strand[i])&&(st[i]+15<=deg[3])&&(deg[3]<=end[i]-25)){
 seg=0;for(j=1;j<=segnum[i];j++){if((excoor[i,j,1]<deg[3])&&(deg[3]<=excoor[i,j,2])){seg=j}};if(seg>0){utr=\"\";loc=(segcoor[i,segnum[i]-seg+1,1]+excoor[i,seg,2]-deg[3]);if(cdsend[i]>0){if(loc<cdsst[i]){utr=\"5utr\"};if(cdsend[i]<loc){utr=\"3utr\"};if((loc>=cdsst[i])&&(cdsend[i]>=loc)){utr=\"cds\"};}else{utr=\"nc\"}; if(loc<100){
 tr=tr\";\"i\";\"loc\";\"utr\";\"fpkm[i]\";\"substr(seq[i],segcoor[i,segnum[i],2]-loc-15+1,40)\";\"substr(seqse[i],1,loc+100)\";\"loc\";\"cdsst[i]\";\"cdsend[i]\";\"segcoor[i,segnum[i],2]}else{
 tr=tr\";\"i\";\"loc\";\"utr\";\"fpkm[i]\";\"substr(seq[i],segcoor[i,segnum[i],2]-loc-15+1,40)\";\"substr(seqse[i],loc-100+1,200)\";\"100\";\"cdsst[i]\";\"cdsend[i]\";\"segcoor[i,segnum[i],2]}}}}};
 
 printf (deg[1]);for (k=2;k<=ndeg;k++){printf (\"\\t\"deg[k])}; 
 getline; printf (\"\\t\"toupper(\$1)); 
 getline;printf (\"\\t\"toupper(\$1)); printf (\"\\t\"tr); printf (\"\\n\")}' mm10.38.92.gtf.exon_coord.fasta.fa.asNs.form.WT1_2_3_4_5_6_st_0fpkm $i > $si.tr"; done
#OUTPUT by COLUMN
 #line1-6 deg fields from .d
 #line7-denovo gene
 #line8-RMSK
 #line9-RMSK_strand
 #line10-mm10_gene
 #line11-mm10_gene_type
 #line12-mm10_gene_name 
 #line13-genome revcomp for -25+15 
 #line14-genome sense for -100+100, 
#fields of 15: 
#1-empty
#2-ENSMUST
#3-loc of cleavage site in ENSMUST
#4-utr_type
#5-mean fpkm
#6-revcomp of -25+15
#7-sense of -100+100
#8-loc of cleaveage site in 7
#9-start of CDS
#10-end of CDS
#11-length of ENSMUST

#merge data after parallel computing
for i in $(find . -maxdepth 1 -name "*.xy.0000.tr"); do echo $i; first=${i%.xy.0000.*}; last=${i#*.xy.0000.}; cat $first.*.$last >  $first.merged.$last; done

#10c.keep only most abundant (and 1+fpkm) transript context for out.pch
for i in $(find . -maxdepth 1 -name "*.tr"); do bsub "awk  'BEGIN{FS=OFS=\"\\t\"}{n=split (\$15,tr,\";\");if (n>1){max=0;maxn=0;for (i=1;i<=(n-1)/10;i++){if((tr[((i-1)*10)+5]>=1)&&(tr[((i-1)*10)+5]>max)){max=tr[((i-1)*10)+5];maxn=i}}; if(maxn>0){ 
print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12\"\\t\" tr[((maxn-1)*10)+2]\"\\t\" tr[((maxn-1)*10)+3]\"\\t\"tr[((maxn-1)*10)+4]\"\\t\"tr[((maxn-1)*10)+5]\"\\t\"tr[((maxn-1)*10)+6]\"\\t\"tr[((maxn-1)*10)+7]\"\\t\"tr[((maxn-1)*10)+8]\"\\t\"tr[((maxn-1)*10)+9]\"\\t\"tr[((maxn-1)*10)+10]\"\\t\"tr[((maxn-1)*10)+11]\"\\t\"\$13\"\\t\"\$14\"\\ttranscript\"}else{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12\"\\t\" 0\"\\t\" 0\"\\t\"0\"\\t\"0\"\\t\"\$13\"\\t\"\$14\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\tnonexpressedloc\"} }else{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12\"\\t\" 0\"\\t\" 0\"\\t\"0\"\\t\"0\"\\t\"\$13\"\\t\"\$14\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\tgenome\"}}' $i > $i.1tr"; done
#OUTPUT by COLUMN
 #line1-deg_chr
 #line2-deg_start
 #line3-deg_end
 #line4-deg_log2change
 #line5-wt_count,mut_count,wt_ppm,mut_ppm
 #line6-deg_strand
 #line7-denovo gene
 #line8-RMSK
 #line9-RMSK_strand
 #line10-mm10_gene
 #line11-mm10_gene_type
 #line12-mm10_gene_name 
 #line13-ENSMUST,
 #line14-loc_position in ENSMUST, 
 #line15-utr_type, 
 #line16-fpkm, 
 #line17-revcomp transcriptseq/or/genomeseq for -25+15, 
 #line18-sense transcriptseq/or/genomeseq for -100+100, 
 #line19-locpos of cleavage site in sense transcriptseq/or/genomeseq for -100+100
 #line20-start of CDS
 #line21-end of CDS
 #line22-length of ENSMUST
 #line23-genome revcomp for -25+15 
 #line24-genome sense for -100+100, 
 #line25-annotation of the site (genome/transcript/nonexpressedloc)