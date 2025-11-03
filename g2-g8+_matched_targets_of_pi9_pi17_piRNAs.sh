#1.get transcript abundances for 7/3/3/3 reps of Sp1, Sp2, RS, ES wild-type control RNA-seq
# *ctab files contain per transcript FPKM data from stringtie
gawk 'BEGIN{FS=OFS="\t"}(FNR>1)&&(FNR==NR){annot[$5]=$8}(FNR>1)&&(FNR<NR)&&(FILENAME~/ctab/){repname=substr(FILENAME,1,index(FILENAME,"_Sp1_")-1); fpkm[$6][repname]=$12; genemeanfpkm[$6]+=$12; genename[$6]=$10; ENSMUSGname[$6]=$9; coord[$6]=$2":"$4"-"$5" "$3; fpkm["reference"][repname]="ref"}END{for (r in fpkm){printf (r"\t"ENSMUSGname[r]"\t"genename[r]"\t"coord[r]"\t"annot[genename[r]]"\t"genemeanfpkm[r]/7); for (f in fpkm["reference"]){if (fpkm[r][f]==""){printf "\t0"}else{printf "\t"fpkm[r][f]}};printf ("\n")}}' mm10.kgXref.txt *WT*Sp1*.ctab | awk '($1!="reference")' > Sp1.7reps.TRANSCRIPTfpkm.txt

gawk 'BEGIN{FS=OFS="\t"}(FNR>1)&&(FNR==NR){annot[$5]=$8}(FNR>1)&&(FNR<NR)&&(FILENAME~/ctab/){repname=substr(FILENAME,1,index(FILENAME,"_Sp2_")-1); fpkm[$6][repname]=$12; genemeanfpkm[$6]+=$12; genename[$6]=$10; ENSMUSGname[$6]=$9; coord[$6]=$2":"$4"-"$5" "$3; fpkm["reference"][repname]="ref"}END{for (r in fpkm){printf (r"\t"ENSMUSGname[r]"\t"genename[r]"\t"coord[r]"\t"annot[genename[r]]"\t"genemeanfpkm[r]/3); for (f in fpkm["reference"]){if (fpkm[r][f]==""){printf "\t0"}else{printf "\t"fpkm[r][f]}};printf ("\n")}}' mm10.kgXref.txt *WT*Sp2*.ctab | awk '($1!="reference")' > Sp2.3reps.TRANSCRIPTfpkm.txt

gawk 'BEGIN{FS=OFS="\t"}(FNR>1)&&(FNR==NR){annot[$5]=$8}(FNR>1)&&(FNR<NR)&&(FILENAME~/ctab/){repname=substr(FILENAME,1,index(FILENAME,"_RS_")-1); fpkm[$6][repname]=$12; genemeanfpkm[$6]+=$12; genename[$6]=$10; ENSMUSGname[$6]=$9; coord[$6]=$2":"$4"-"$5" "$3; fpkm["reference"][repname]="ref"}END{for (r in fpkm){printf (r"\t"ENSMUSGname[r]"\t"genename[r]"\t"coord[r]"\t"annot[genename[r]]"\t"genemeanfpkm[r]/3); for (f in fpkm["reference"]){if (fpkm[r][f]==""){printf "\t0"}else{printf "\t"fpkm[r][f]}};printf ("\n")}}' mm10.kgXref.txt *WT*RS*.ctab | awk '($1!="reference")' > RS.3reps.TRANSCRIPTfpkm.txt

gawk 'BEGIN{FS=OFS="\t"}(FNR>1)&&(FNR==NR){annot[$5]=$8}(FNR>1)&&(FNR<NR)&&(FILENAME~/ctab/){repname=FILENAME; fpkm[$6][repname]=$12; genemeanfpkm[$6]+=$12; genename[$6]=$10; ENSMUSGname[$6]=$9; coord[$6]=$2":"$4"-"$5" "$3; fpkm["reference"][repname]="ref"}END{for (r in fpkm){printf (r"\t"ENSMUSGname[r]"\t"genename[r]"\t"coord[r]"\t"annot[genename[r]]"\t"genemeanfpkm[r]/3); for (f in fpkm["reference"]){if (fpkm[r][f]==""){printf "\t0"}else{printf "\t"fpkm[r][f]}};printf ("\n")}}' mm10.kgXref.txt *WT*ES*.ctab | awk '($1!="reference")' > ES.3reps.TRANSCRIPTfpkm.txt

#2.create transcript exclusion list to remove overlaps with pachytene clusters
#mm10.38.92.gtf.3utr.bed file contains boundaries of all 3'UTRs in mm10
module load bedtools/2.26.0
bedtools intersect -wa -a mm10.38.92.gtf.3utr.bed -b piRNA.cluster.pachytene.bed | awk 'BEGIN{FS=OFS="\t"}{split ($4,a,":");print a[1]}' | sort | uniq > ENSMUST.pachyteneclusters

#3.get spliced fasta plus genomic coordinates of exons and cds
#mm10 genome was downloaded from https://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/
#mm10.38.92.gtf was downloaded from https://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/
module load cufflinks/2.2.1
gffread mm10.38.92.gtf -g mm10.pipipes.fa -W -w mm10.38.92.gtf.exon_coord.fasta
num=$(cat mm10.38.92.gtf.exon_coord.fasta | wc -l); awk -v num="$num" '(FNR==1){print}(FNR>1){full=""; while (index($1,">")==0){full=full $1; getline; if(FNR==num){full=full $1; print toupper(full);exit}}; print toupper(full); print}' mm10.38.92.gtf.exon_coord.fasta > mm10.38.92.gtf.exon_coord.fasta.fa
for i in $(find . -maxdepth 1 -name "mm10.38.92.gtf.exon_coord.fasta.fa"); do bsub "awk '(FNR%2==1){name=\$0; getline;full=\"\";\$1=toupper(\$1);for (j=length(\$1);j>=1;j--){if(substr(\$1,j,1)==\"A\"){full=full\"T\"};if(substr(\$1,j,1)==\"T\"){full=full\"A\"};if(substr(\$1,j,1)==\"G\"){full=full\"C\"};if(substr(\$1,j,1)==\"C\"){full=full\"G\"}}; print name; print full}' $i > $i.as"; done

#4.reformat as:
#line1-ENSMUST;ENSMUSG;refname;chr;strand;exons (1-2,3-4,5-6) 1-based on genomic strand;segs (1-2,3-4,5-6) 1-based on transcribed strand (reverse order for - genomic strand);cds start on transcribed strand (reverse order for - genomic strand);cds end on transcribed strand (reverse order for - genomic strand)
#line2-sequence
awk 'BEGIN{FS=OFS="\t"}(FNR==NR){gname[$1]=$2; refname[$1]=$3}
(FNR%2==1)&&(FNR<NR){n=split(substr($0,2),fields," ");strand="";cds[1]="";cds[2]="";exons="";segs="";
for (i=2; i<=n; i++){
if(index(fields[i],"gene=")==1){total++;if (substr(fields[i],6)==refname[fields[1]]){coinc++}};
if(index(fields[i],"CDS=")==1){split(substr(fields[i],5),cds,"-")};
if(index(fields[i],"loc:")==1){strand=substr(fields[i],length(fields[i]),1);chr=substr(fields[i],5,index(fields[i],"|")-5)};
if(index(fields[i],"exons:")==1){exons=substr(fields[i],7)};
if(index(fields[i],"segs:")==1){segs=substr(fields[i],6)}};
print ">"fields[1]";"gname[fields[1]]";"refname[fields[1]]";"chr";"strand";"exons";"segs";"cds[1]";"cds[2]; getline; print}END{print coinc"\t"total"\t"coinc/total > "test.txt"}' mm10.38.92.gtf.t_g_n mm10.38.92.gtf.exon_coord.fasta.fa.as > mm10.38.92.gtf.exon_coord.fasta.fa.as.form

#5.filter for 1+fpkm transcripts
awk 'BEGIN{FS=OFS="\t"}
(FNR==NR)&&(FNR>1){ab[$1]=$6}(FNR<NR)&&(FILENAME~/geneNnameNtype/){type[$1]=$2}
(FNR<NR)&&(FNR%2==1)&&(FILENAME~/form/){split(substr($1,2),gname,";");name=$1;
if(ab[gname[1]]>=1){getline; print name";"ab[gname[1]]";"type[gname[2]];print $0}}' Sp1.7reps.TRANSCRIPTfpkm.txt mm10.38.92.gtf.exons.geneNnameNtype  mm10.38.92.gtf.exon_coord.fasta.fa.as.form >  mm10.38.92.gtf.exon_coord.fasta.fa.as.form_Sp1_7rep_1fpkm
awk 'BEGIN{FS=OFS="\t"}(FNR==NR)&&(FNR>1){ab[$1]=$6}(FNR<NR)&&(FILENAME~/geneNnameNtype/){type[$1]=$2}(FNR<NR)&&(FNR%2==1)&&(FILENAME~/form/){split(substr($1,2),gname,";");name=$1;if(ab[gname[1]]>1){getline; print name";"ab[gname[1]]";"type[gname[2]];print $0}}' Sp2.3reps.TRANSCRIPTfpkm.txt mm10.38.92.gtf.exons.geneNnameNtype mm10.38.92.gtf.exon_coord.fasta.fa.as.form >  mm10.38.92.gtf.exon_coord.fasta.fa.as.form_Sp2_3rep_1fpkm
awk 'BEGIN{FS=OFS="\t"}(FNR==NR)&&(FNR>1){ab[$1]=$6}(FNR<NR)&&(FILENAME~/geneNnameNtype/){type[$1]=$2}(FNR<NR)&&(FNR%2==1)&&(FILENAME~/form/){split(substr($1,2),gname,";");name=$1;if(ab[gname[1]]>1){getline; print name";"ab[gname[1]]";"type[gname[2]];print $0}}' RS.3reps.TRANSCRIPTfpkm.txt mm10.38.92.gtf.exons.geneNnameNtype mm10.38.92.gtf.exon_coord.fasta.fa.as.form >  mm10.38.92.gtf.exon_coord.fasta.fa.as.form_RS_3rep_1fpkm
awk 'BEGIN{FS=OFS="\t"}(FNR==NR)&&(FNR>1){ab[$1]=$6}(FNR<NR)&&(FILENAME~/geneNnameNtype/){type[$1]=$2}(FNR<NR)&&(FNR%2==1)&&(FILENAME~/form/){split(substr($1,2),gname,";");name=$1;if(ab[gname[1]]>1){getline; print name";"ab[gname[1]]";"type[gname[2]];print $0}}' ES.3reps.TRANSCRIPTfpkm.txt mm10.38.92.gtf.exons.geneNnameNtype mm10.38.92.gtf.exon_coord.fasta.fa.as.form >  mm10.38.92.gtf.exon_coord.fasta.fa.as.form_ES_3rep_1fpkm

#6.get 3'UTR sequences
awk 'BEGIN{FS=OFS="\t"}(FNR==NR)&&(FNR%2==1){split(substr($1,2),gname,";");dat[gname[1]]=substr($0,2)}(FNR<NR)&&(FNR%2==1){split(substr($1,2),gname,";");if(dat[gname[1]]!=""){print $0":"dat[gname[1]]; getline; print}}' mm10.38.92.gtf.exon_coord.fasta.fa.as.form_Sp1_7rep_1fpkm mm10.38.92.gtf.3utr.bed.fa.as > mm10.38.92.gtf.3utr.bed.fa.as_Sp1_7rep_1fpkm
awk 'BEGIN{FS=OFS="\t"}(FNR==NR)&&(FNR%2==1){split(substr($1,2),gname,";");dat[gname[1]]=substr($0,2)}(FNR<NR)&&(FNR%2==1){split(substr($1,2),gname,";");if(dat[gname[1]]!=""){print $0":"dat[gname[1]]; getline; print}}' mm10.38.92.gtf.exon_coord.fasta.fa.as.form_Sp2_3rep_1fpkm mm10.38.92.gtf.3utr.bed.fa.as > mm10.38.92.gtf.3utr.bed.fa.as_Sp2_3rep_1fpkm
awk 'BEGIN{FS=OFS="\t"}(FNR==NR)&&(FNR%2==1){split(substr($1,2),gname,";");dat[gname[1]]=substr($0,2)}(FNR<NR)&&(FNR%2==1){split(substr($1,2),gname,";");if(dat[gname[1]]!=""){print $0":"dat[gname[1]]; getline; print}}' mm10.38.92.gtf.exon_coord.fasta.fa.as.form_RS_3rep_1fpkm mm10.38.92.gtf.3utr.bed.fa.as > mm10.38.92.gtf.3utr.bed.fa.as_RS_3rep_1fpkm
awk 'BEGIN{FS=OFS="\t"}(FNR==NR)&&(FNR%2==1){split(substr($1,2),gname,";");dat[gname[1]]=substr($0,2)}(FNR<NR)&&(FNR%2==1){split(substr($1,2),gname,";");if(dat[gname[1]]!=""){print $0":"dat[gname[1]]; getline; print}}' mm10.38.92.gtf.exon_coord.fasta.fa.as.form_ES_3rep_1fpkm mm10.38.92.gtf.3utr.bed.fa.as > mm10.38.92.gtf.3utr.bed.fa.as_ES_3rep_1fpkm

#7.find matching piRNAs
#7a.split into 2000 lines for parallel computing
#merged.uniq.piRNAs.cl.g.r.rfmtd is generated in 11e of SmallRNAseq_data_processing.sh
for i in $(find . -maxdepth 1 -name "merged.uniq.piRNAs.cl.g.r.rfmtd"); do split -a 4 -d -l 2000 $i $i.xx. ; done

#7b.first, get all g2-g8 matchdes for mm10.38.92.gtf.3utr.bed.fa.as_Sp1/Sp2/RS/ES_Xrep_1fpkm
for i in $(find . -maxdepth 1 -name "*fpkm"); do si=${i:2}; si=${si//mm10.38.92.gtf.exon_coord.fasta.fa.as.form_/}; for j in $(find . -maxdepth 1 -name "*.xx.*"); do sj=${j//.xx./.xy.}; sj=${sj//merged./};  for n in 8; do bsub "awk -v n=$n 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){pirnan++;
 pirna[pirnan,1]=\$4;
 pirna[pirnan,2]=\$5;
 pirna[pirnan,3]=\$1\":\"\$2\":\"\$3\":\"\$6;
 pirna[pirnan,4]=\$7;
 pirna[pirnan,5]=\$13}(FNR<NR)&&(FNR%2==1){split(substr(\$1,2),dat,\";\");
 split(dat[6],exons,\",\");
 m=split(dat[7],segs,\",\");
 for (i=1;i<=m;i++){split(exons[i],b,\"-\");
 excoor[i,1]=b[1];
 excoor[i,2]=b[2];
 split(segs[i],b,\"-\");
 segcoor[i,1]=b[1];
 segcoor[i,2]=b[2]};
 getline;for (i=1;i<=pirnan;i++){seeds=split(\$1,sprout,substr(pirna[i,1],2,n-1));
 g10[1]=length(sprout[1])+9; ;g10sense[1]=(segcoor[m,2]-g10[1]+1); for (seed=2;seed<seeds;seed++) {g10[seed]=(g10[seed-1]+length(sprout[seed])+n-1);g10sense[seed]=(segcoor[m,2]-g10[seed]+1)};
 for(seed=1;seed<seeds;seed++) {if(dat[8]!=\"\"){if(g10sense[seed]<dat[8]){annot=\"utr5\"}else{if(g10sense[seed]>dat[9]){annot=\"utr3\"}else{annot=\"cds\"}}}else{annot=\"nc\"};
 for (j=1;j<=m;j++){if((g10sense[seed]>=segcoor[j,1])&&(g10sense[seed]<=segcoor[j,2])){str=\"\";if(dat[5]==\"+\"){if ((g10[seed]-39)<1){for(l=-1;l>=(g10[seed]-39-1);l--){str=str\"N\"};str=str substr(\$1,1,90+g10[seed]-39-1)}else{str=substr(\$1,g10[seed]-39,90)};printf (dat[4]\"\\t\"(excoor[j,1]+(g10sense[seed]-segcoor[j,1])-1)\"\\t\"(excoor[j,1]+(g10sense[seed]-segcoor[j,1]))\"\\t\"dat[1]\",\"dat[2]\",\"dat[11]\",\"annot\"\\t\"dat[10]\"\\t\"dat[5]\"\\t\"str);
 for (k=1;k<=5;k++){printf (\"\\t\"pirna[i,k])}; printf (\"\\n\")}else{if ((g10[seed]-39)<1){for(l=-1;l>=(g10[seed]-39-1);l--){str=str\"N\"};str=str substr(\$1,1,90+g10[seed]-39-1)}else{str=substr(\$1,g10[seed]-39,90)};printf (dat[4]\"\\t\"(excoor[m-j+1,2]-(g10sense[seed]-segcoor[j,1])-1)\"\\t\"(excoor[m-j+1,2]-(g10sense[seed]-segcoor[j,1]))\"\\t\"dat[1]\",\"dat[2]\",\"dat[11]\",\"annot\"\\t\"dat[10]\"\\t\"dat[5]\"\\t\"str); for (k=1;k<=9;k++){printf (\"\\t\"pirna[i,k])}; printf (\"\\n\")}}}}}}' $j $i > $sj.$si.$n.ms"; done; done; done

#7c.add 4 reps of wild type control degradome to transcript matches (*.d files are generated in 7 of 5P-seq_data_initial_processing.sh)
#keep only matches with at least n=[0..19] matches outside of seed and produce table with total targeting piRNAs (pi9pi17 or non-pi9pi17) vs (non-degradome,degradome,degradome_down) vs (utr5,utr3,cds,nc): 2x3x4=24 permutations total
for j in $(find . -maxdepth 1 -name "*.d"); do sj=${j:2}; sj=${sj/.p2p9p17_E.unique.out.m.5.d/}; for i in $(find . -maxdepth 1 -name "*.xy.*_1fpkm.8.ms"); do for n in {0..19}; do 
/ms_short_32G "awk -v sj=$sj -v n=$n 'BEGIN{FS=OFS=\"\\t\"}
(FILENAME==\"merged.uniq.piRNAs.cl.g.r.rfmtd\"){pir30[substr(\$4,1,25)]=\$4}
(FILENAME~/clusters/){clust[\$1]=1}
(FILENAME~/.unique.out.m.5.d/){split(\$5,a,\",\");if(a[3]>=0.04){if (\$4>-2){degseq[\$1\$2\$3\$6]=\"0\"}else{degseq[\$1\$2\$3\$6]=\"3\"}}}
(FILENAME~/ms/){ split(\$4,g,\",\");
 if(clust[g[1]]!=1){\$8=pir30[\$8]; nf=split(FILENAME,f,\".\");
if((substr(\$7,32,f[nf-1]-1)==substr(\$8,2,f[nf-1]-1))){mtchd=0;for (i=1;i<=22;i++){if(substr(\$7,38+i,1)==substr(\$8,8+i,1)){mtchd++}};if(mtchd>=n){genes[g[2]]=1;
 if(deg[g[2]\$1\$2\$3\$6\$8]==\"\"){deg[g[2]\$1\$2\$3\$6\$8]++;
degab[g[2]\$1\$2\$3\$6\$8]=\$9;
deggene[g[2]\$1\$2\$3\$6\$8]=g[2];
degloc[g[2]\$1\$2\$3\$6\$8]=g[4];
if(degseq[\$1\$2\$3\$6]!=\"\"){degev[g[2]\$1\$2\$3\$6\$8]=degseq[\$1\$2\$3\$6]}else{degev[g[2]\$1\$2\$3\$6\$8]=\"na\"};
if(\$12<0.2){deg917[g[2]\$1\$2\$3\$6\$8]=\"r\"}else{deg917[g[2]\$1\$2\$3\$6\$8]=\"nr\"}}else{if(degloc[g[2]\$1\$2\$3\$6\$8]!=g[4]){if(g[4]==\"cds\"){degloc[g[2]\$1\$2\$3\$6\$8]=g[4]}else{if (g[4]==\"utr3\"){degloc[g[2]\$1\$2\$3\$6\$8]=g[4]}else{if (g[4]==\"utr5\"){degloc[g[2]\$1\$2\$3\$6\$8]=g[4]}}}}}}}}}END{
for (d in deg){
genedata917ab[deggene[d],deg917[d],degev[d],degloc[d]]+=degab[d]};
for (gene in genes){
if (genedata917ab[gene,\"nr\",\"na\",\"utr5\"]==\"\"){genedata917ab[gene,\"nr\",\"na\",\"utr5\"]=0};
if (genedata917ab[gene,\"nr\",\"na\",\"cds\"]==\"\"){genedata917ab[gene,\"nr\",\"na\",\"cds\"]=0};
if (genedata917ab[gene,\"nr\",\"na\",\"utr3\"]==\"\"){genedata917ab[gene,\"nr\",\"na\",\"utr3\"]=0};
if (genedata917ab[gene,\"nr\",\"na\",\"nc\"]==\"\"){genedata917ab[gene,\"nr\",\"na\",\"nc\"]=0};
if (genedata917ab[gene,\"nr\",\"0\",\"utr5\"]==\"\"){genedata917ab[gene,\"nr\",\"0\",\"utr5\"]=0};
if (genedata917ab[gene,\"nr\",\"0\",\"cds\"]==\"\"){genedata917ab[gene,\"nr\",\"0\",\"cds\"]=0};
if (genedata917ab[gene,\"nr\",\"0\",\"utr3\"]==\"\"){genedata917ab[gene,\"nr\",\"0\",\"utr3\"]=0};
if (genedata917ab[gene,\"nr\",\"0\",\"nc\"]==\"\"){genedata917ab[gene,\"nr\",\"0\",\"nc\"]=0};
if (genedata917ab[gene,\"nr\",\"3\",\"utr5\"]==\"\"){genedata917ab[gene,\"nr\",\"3\",\"utr5\"]=0};
if (genedata917ab[gene,\"nr\",\"3\",\"cds\"]==\"\"){genedata917ab[gene,\"nr\",\"3\",\"cds\"]=0};
if (genedata917ab[gene,\"nr\",\"3\",\"utr3\"]==\"\"){genedata917ab[gene,\"nr\",\"3\",\"utr3\"]=0};
if (genedata917ab[gene,\"nr\",\"3\",\"nc\"]==\"\"){genedata917ab[gene,\"nr\",\"3\",\"nc\"]=0};
if (genedata917ab[gene,\"r\",\"na\",\"utr5\"]==\"\"){genedata917ab[gene,\"r\",\"na\",\"utr5\"]=0};
if (genedata917ab[gene,\"r\",\"na\",\"cds\"]==\"\"){genedata917ab[gene,\"r\",\"na\",\"cds\"]=0};
if (genedata917ab[gene,\"r\",\"na\",\"utr3\"]==\"\"){genedata917ab[gene,\"r\",\"na\",\"utr3\"]=0};
if (genedata917ab[gene,\"r\",\"na\",\"nc\"]==\"\"){genedata917ab[gene,\"r\",\"na\",\"nc\"]=0};
if (genedata917ab[gene,\"r\",\"0\",\"utr5\"]==\"\"){genedata917ab[gene,\"r\",\"0\",\"utr5\"]=0};
if (genedata917ab[gene,\"r\",\"0\",\"cds\"]==\"\"){genedata917ab[gene,\"r\",\"0\",\"cds\"]=0};
if (genedata917ab[gene,\"r\",\"0\",\"utr3\"]==\"\"){genedata917ab[gene,\"r\",\"0\",\"utr3\"]=0};
if (genedata917ab[gene,\"r\",\"0\",\"nc\"]==\"\"){genedata917ab[gene,\"r\",\"0\",\"nc\"]=0};
if (genedata917ab[gene,\"r\",\"3\",\"utr5\"]==\"\"){genedata917ab[gene,\"r\",\"3\",\"utr5\"]=0};
if (genedata917ab[gene,\"r\",\"3\",\"cds\"]==\"\"){genedata917ab[gene,\"r\",\"3\",\"cds\"]=0};
if (genedata917ab[gene,\"r\",\"3\",\"utr3\"]==\"\"){genedata917ab[gene,\"r\",\"3\",\"utr3\"]=0};
if (genedata917ab[gene,\"r\",\"3\",\"nc\"]==\"\"){genedata917ab[gene,\"r\",\"3\",\"nc\"]=0};
print gene\"\\t\"genedata917ab[gene,\"nr\",\"na\",\"utr5\"]\"\\t\"genedata917ab[gene,\"nr\",\"na\",\"cds\"]\"\\t\"genedata917ab[gene,\"nr\",\"na\",\"utr3\"]\"\\t\"genedata917ab[gene,\"nr\",\"na\",\"nc\"]\"\\t\"genedata917ab[gene,\"nr\",\"0\",\"utr5\"]\"\\t\"genedata917ab[gene,\"nr\",\"0\",\"cds\"]\"\\t\"genedata917ab[gene,\"nr\",\"0\",\"utr3\"]\"\\t\"genedata917ab[gene,\"nr\",\"0\",\"nc\"]\"\\t\"genedata917ab[gene,\"nr\",\"3\",\"utr5\"]\"\\t\"genedata917ab[gene,\"nr\",\"3\",\"cds\"]\"\\t\"genedata917ab[gene,\"nr\",\"3\",\"utr3\"]\"\\t\"genedata917ab[gene,\"nr\",\"3\",\"nc\"]\"\\t\"genedata917ab[gene,\"r\",\"na\",\"utr5\"]\"\\t\"genedata917ab[gene,\"r\",\"na\",\"cds\"]\"\\t\"genedata917ab[gene,\"r\",\"na\",\"utr3\"]\"\\t\"genedata917ab[gene,\"r\",\"na\",\"nc\"]\"\\t\"genedata917ab[gene,\"r\",\"0\",\"utr5\"]\"\\t\"genedata917ab[gene,\"r\",\"0\",\"cds\"]\"\\t\"genedata917ab[gene,\"r\",\"0\",\"utr3\"]\"\\t\"genedata917ab[gene,\"r\",\"0\",\"nc\"]\"\\t\"genedata917ab[gene,\"r\",\"3\",\"utr5\"]\"\\t\"genedata917ab[gene,\"r\",\"3\",\"cds\"]\"\\t\"genedata917ab[gene,\"r\",\"3\",\"utr3\"]\"\\t\"genedata917ab[gene,\"r\",\"3\",\"nc\"] > FILENAME\".\"sj\".917ab.\"n\".seed30\"}}' merged.uniq.piRNAs.cl.g.r.rfmtd ENSMUST.pachyteneclusters $j $i  "; done; done; done
#merge after parallel computing
for i in $(find . -maxdepth 1 -name "*.xy.0000.*seed30"); do echo $i; first=${i%.xy.0000.*}; last=${i#*.xy.0000.}; awk 'BEGIN{FS=OFS="\t"}{fields=NF;genes[$1]=1;for (i=2;i<=NF;i++){dat[$1,i]+=$i}}END{for (g in genes){printf (g);for(i=2;i<=fields;i++){printf("\t"dat[g,i])};printf("\n")}}' $first.*.$last > $first.merged.$last; done

#8a.count ELAVL1 in 3UTRs of each transcript expressed in 7/3/3 reps of Sp1, Sp2, RS wild-type control RNA-seq
awk 'BEGIN{FS=OFS="\t"}{gene[$9]+=$12;tr[$6]+=$12;gname[$6]=$9}END{for(t in tr){if(gene[gname[t]]>0){print t"\t"tr[t]/7"\t"tr[t]/gene[gname[t]]"\t"gname[t]"\t"gene[gname[t]]/7}}}' *WT*Sp1*ctab > all.transcript_fractions_Sp1_7reps.txt
awk 'BEGIN{FS=OFS="\t"}{gene[$9]+=$12;tr[$6]+=$12;gname[$6]=$9}END{for(t in tr){if(gene[gname[t]]>0){print t"\t"tr[t]/3"\t"tr[t]/gene[gname[t]]"\t"gname[t]"\t"gene[gname[t]]/3}}}' *WT*Sp2*ctab > all.transcript_fractions_Sp2_3reps.txt
awk 'BEGIN{FS=OFS="\t"}{gene[$9]+=$12;tr[$6]+=$12;gname[$6]=$9}END{for(t in tr){if(gene[gname[t]]>0){print t"\t"tr[t]/3"\t"tr[t]/gene[gname[t]]"\t"gname[t]"\t"gene[gname[t]]/3}}}' *WT*RS*ctab > all.transcript_fractions_RS_3reps.txt

awk 'BEGIN{FS=OFS="\t";kmer[1]="TTTTTTT";kmer[2]="TTTATTT";kmer[3]="TTTGTTT";kmer[4]="TATTTAT";kmer[5]="ATTTTTA";kmer[6]="ATTTATT";kmer[7]="AATTTTA";kmer[8]="AATATTT"}(FNR==NR)&&(FNR%2==1){name=substr($1,2,18);getline;if(length($1)>=7){num[name]=0;for (i=1;i<=7;i++){num[name]+=(split($1,a,kmer[i])-1)}}}(FNR<NR)&&(num[$1]!=""){print $0"\t"num[$1] > FILENAME".aremfl"}' mm10.38.92.gtf.3utr.fasta.fa all.transcript_fractions_*reps.txt

awk 'BEGIN{FS=OFS="\t"}(FILENAME~/Sp1_7reps/){datSp1[$1]=$6;all[$1]=1}(FILENAME~/Sp2_3reps/){datsp2[$1]=$6;all[$1]=1}(FILENAME~/RS_3reps/){datrs[$1]=$6;all[$1]=1}END{for (i in all){if(datSp1[i]==""){datSp1[i]="0"};if(datsp2[i]==""){datsp2[i]="0"};if(datrs[i]==""){datrs[i]="0"};print i"\t"datSp1[i]"\t"datsp2[i]"\t"datrs[i]}}' all.transcript_fractions_Sp1_7reps.txt.aremfl all.transcript_fractions_Sp2_3reps.txt.aremfl all.transcript_fractions_RS_3reps.txt.aremfl > all.aremfl.Sp1x7_Sp2x3_RSx3.txt

#8b.add 4 reps of wild type control degradome to transcript matches (*.d files are generated in 7 of 5P-seq_data_initial_processing.sh)
#keep only matches with at least n=[0..19] matches outside of seed and 3+ ELAVL1 sites
#produce table with total targeting piRNAs (pi9pi17 or non-pi9pi17) vs (non-degradome,degradome,degradome_down) vs (utr5,utr3,cds,nc): 2x3x4=24 permutations total
for j in $(find . -maxdepth 1 -name "*.d"); do sj=${j:2}; sj=${sj/.p2p9p17_E.unique.out.m.5.d/}; for i in $(find . -maxdepth 1 -name "*.xy.*_1fpkm.8.ms" | grep -v ES_3rep ); do for n in {0..19}; do bsub "awk -v sj=$sj -v file=$i -v n=$n 'BEGIN{FS=OFS=\"\\t\"; nf=split(file,fa,\".\");if(substr(fa[fn-2],1,2)==\"Sp1\"){hurind=2};if(substr(fa[fn-2],1,2)==\"Sp\"){hurind=3};if(substr(fa[fn-2],1,2)==\"RS\"){hurind=4}}
(FILENAME==\"merged.uniq.piRNAs.cl.g.r.rfmtd\"){pir30[substr(\$4,1,25)]=\$4}
(FILENAME==\"all.aremfl.Sp1x7_Sp2x3_RSx3.txt\"){hur[\$1]=\$(hurind)}
(FILENAME==\"ENSMUST.pachyteneclusters\"){clust[\$1]=1}
(FILENAME~/.unique.out.m.5.d/){split(\$5,a,\",\");if(a[3]>=0.04){if (\$4>-2){degseq[\$1\$2\$3\$6]=\"0\"}else{degseq[\$1\$2\$3\$6]=\"3\"}}}
(FILENAME~/ms/){ split(\$4,g,\",\");
 if((clust[g[1]]!=1)&&(hur[g[1]]>2)){\$8=pir30[\$8]; nf=split(FILENAME,f,\".\");
if((substr(\$7,32,f[nf-1]-1)==substr(\$8,2,f[nf-1]-1))){mtchd=0;for (i=1;i<=22;i++){if(substr(\$7,38+i,1)==substr(\$8,8+i,1)){mtchd++}};if(mtchd>=n){genes[g[2]]=1;
 if(deg[g[2]\$1\$2\$3\$6\$8]==\"\"){deg[g[2]\$1\$2\$3\$6\$8]++;
degab[g[2]\$1\$2\$3\$6\$8]=\$9;
deggene[g[2]\$1\$2\$3\$6\$8]=g[2];
degloc[g[2]\$1\$2\$3\$6\$8]=g[4];
if(degseq[\$1\$2\$3\$6]!=\"\"){degev[g[2]\$1\$2\$3\$6\$8]=degseq[\$1\$2\$3\$6]}else{degev[g[2]\$1\$2\$3\$6\$8]=\"na\"};
if(\$12<0.2){deg917[g[2]\$1\$2\$3\$6\$8]=\"r\"}else{deg917[g[2]\$1\$2\$3\$6\$8]=\"nr\"}}else{if(degloc[g[2]\$1\$2\$3\$6\$8]!=g[4]){if(g[4]==\"cds\"){degloc[g[2]\$1\$2\$3\$6\$8]=g[4]}else{if (g[4]==\"utr3\"){degloc[g[2]\$1\$2\$3\$6\$8]=g[4]}else{if (g[4]==\"utr5\"){degloc[g[2]\$1\$2\$3\$6\$8]=g[4]}}}}}}}}}END{
for (d in deg){
genedata917ab[deggene[d],deg917[d],degev[d],degloc[d]]+=degab[d]};

for (gene in genes){
if (genedata917ab[gene,\"nr\",\"na\",\"utr5\"]==\"\"){genedata917ab[gene,\"nr\",\"na\",\"utr5\"]=0};
if (genedata917ab[gene,\"nr\",\"na\",\"cds\"]==\"\"){genedata917ab[gene,\"nr\",\"na\",\"cds\"]=0};
if (genedata917ab[gene,\"nr\",\"na\",\"utr3\"]==\"\"){genedata917ab[gene,\"nr\",\"na\",\"utr3\"]=0};
if (genedata917ab[gene,\"nr\",\"na\",\"nc\"]==\"\"){genedata917ab[gene,\"nr\",\"na\",\"nc\"]=0};
if (genedata917ab[gene,\"nr\",\"0\",\"utr5\"]==\"\"){genedata917ab[gene,\"nr\",\"0\",\"utr5\"]=0};
if (genedata917ab[gene,\"nr\",\"0\",\"cds\"]==\"\"){genedata917ab[gene,\"nr\",\"0\",\"cds\"]=0};
if (genedata917ab[gene,\"nr\",\"0\",\"utr3\"]==\"\"){genedata917ab[gene,\"nr\",\"0\",\"utr3\"]=0};
if (genedata917ab[gene,\"nr\",\"0\",\"nc\"]==\"\"){genedata917ab[gene,\"nr\",\"0\",\"nc\"]=0};
if (genedata917ab[gene,\"nr\",\"3\",\"utr5\"]==\"\"){genedata917ab[gene,\"nr\",\"3\",\"utr5\"]=0};
if (genedata917ab[gene,\"nr\",\"3\",\"cds\"]==\"\"){genedata917ab[gene,\"nr\",\"3\",\"cds\"]=0};
if (genedata917ab[gene,\"nr\",\"3\",\"utr3\"]==\"\"){genedata917ab[gene,\"nr\",\"3\",\"utr3\"]=0};
if (genedata917ab[gene,\"nr\",\"3\",\"nc\"]==\"\"){genedata917ab[gene,\"nr\",\"3\",\"nc\"]=0};
if (genedata917ab[gene,\"r\",\"na\",\"utr5\"]==\"\"){genedata917ab[gene,\"r\",\"na\",\"utr5\"]=0};
if (genedata917ab[gene,\"r\",\"na\",\"cds\"]==\"\"){genedata917ab[gene,\"r\",\"na\",\"cds\"]=0};
if (genedata917ab[gene,\"r\",\"na\",\"utr3\"]==\"\"){genedata917ab[gene,\"r\",\"na\",\"utr3\"]=0};
if (genedata917ab[gene,\"r\",\"na\",\"nc\"]==\"\"){genedata917ab[gene,\"r\",\"na\",\"nc\"]=0};
if (genedata917ab[gene,\"r\",\"0\",\"utr5\"]==\"\"){genedata917ab[gene,\"r\",\"0\",\"utr5\"]=0};
if (genedata917ab[gene,\"r\",\"0\",\"cds\"]==\"\"){genedata917ab[gene,\"r\",\"0\",\"cds\"]=0};
if (genedata917ab[gene,\"r\",\"0\",\"utr3\"]==\"\"){genedata917ab[gene,\"r\",\"0\",\"utr3\"]=0};
if (genedata917ab[gene,\"r\",\"0\",\"nc\"]==\"\"){genedata917ab[gene,\"r\",\"0\",\"nc\"]=0};
if (genedata917ab[gene,\"r\",\"3\",\"utr5\"]==\"\"){genedata917ab[gene,\"r\",\"3\",\"utr5\"]=0};
if (genedata917ab[gene,\"r\",\"3\",\"cds\"]==\"\"){genedata917ab[gene,\"r\",\"3\",\"cds\"]=0};
if (genedata917ab[gene,\"r\",\"3\",\"utr3\"]==\"\"){genedata917ab[gene,\"r\",\"3\",\"utr3\"]=0};
if (genedata917ab[gene,\"r\",\"3\",\"nc\"]==\"\"){genedata917ab[gene,\"r\",\"3\",\"nc\"]=0};
print gene\"\\t\"genedata917ab[gene,\"nr\",\"na\",\"utr5\"]\"\\t\"genedata917ab[gene,\"nr\",\"na\",\"cds\"]\"\\t\"genedata917ab[gene,\"nr\",\"na\",\"utr3\"]\"\\t\"genedata917ab[gene,\"nr\",\"na\",\"nc\"]\"\\t\"genedata917ab[gene,\"nr\",\"0\",\"utr5\"]\"\\t\"genedata917ab[gene,\"nr\",\"0\",\"cds\"]\"\\t\"genedata917ab[gene,\"nr\",\"0\",\"utr3\"]\"\\t\"genedata917ab[gene,\"nr\",\"0\",\"nc\"]\"\\t\"genedata917ab[gene,\"nr\",\"3\",\"utr5\"]\"\\t\"genedata917ab[gene,\"nr\",\"3\",\"cds\"]\"\\t\"genedata917ab[gene,\"nr\",\"3\",\"utr3\"]\"\\t\"genedata917ab[gene,\"nr\",\"3\",\"nc\"]\"\\t\"genedata917ab[gene,\"r\",\"na\",\"utr5\"]\"\\t\"genedata917ab[gene,\"r\",\"na\",\"cds\"]\"\\t\"genedata917ab[gene,\"r\",\"na\",\"utr3\"]\"\\t\"genedata917ab[gene,\"r\",\"na\",\"nc\"]\"\\t\"genedata917ab[gene,\"r\",\"0\",\"utr5\"]\"\\t\"genedata917ab[gene,\"r\",\"0\",\"cds\"]\"\\t\"genedata917ab[gene,\"r\",\"0\",\"utr3\"]\"\\t\"genedata917ab[gene,\"r\",\"0\",\"nc\"]\"\\t\"genedata917ab[gene,\"r\",\"3\",\"utr5\"]\"\\t\"genedata917ab[gene,\"r\",\"3\",\"cds\"]\"\\t\"genedata917ab[gene,\"r\",\"3\",\"utr3\"]\"\\t\"genedata917ab[gene,\"r\",\"3\",\"nc\"] > FILENAME\".\"sj\".917ab.\"n\".seedaremfl\"}}' merged.uniq.piRNAs.cl.g.r.rfmtd all.aremfl.Sp1x7_Sp2x3_RSx3.txt ENSMUST.pachyteneclusters $j $i  "; done; done; done
#merge after parallel computing
for i in $(find . -maxdepth 1 -name "*.xy.0000.*seedaremfl"); do echo $i; first=${i%.xy.0000.*}; last=${i#*.xy.0000.}; awk 'BEGIN{FS=OFS="\t"}{fields=NF;genes[$1]=1;for (i=2;i<=NF;i++){dat[$1,i]+=$i}}END{for (g in genes){printf (g);for(i=2;i<=fields;i++){printf("\t"dat[g,i])};printf("\n")}}' $first.*.$last > $first.merged.$last; done

