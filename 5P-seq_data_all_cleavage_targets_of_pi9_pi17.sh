#all p9 and p17 cleavage targets
#1a.select deg sites at >=0.1ppm and going down by >8-fold in pi9pi17 mutant
for file in `ls *.d`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(\$4<=-3){split(\$5,a,\",\");if(a[3]>=0.1){\$5=a[3];print}}' $file > $file.01" ; done

#1b.get corresponding transcript data for each putative cleavage product
# *.1tr files are generated in 10c of 5P-seq_data_initial_processing.sh
for wt in WT1 WT2 WT3 WT4; do for file in $wt.*.01; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){deg[\$1\":\"\$2\":\"\$6]=\$5}(FNR<NR){if((deg[\$1\":\"\$2\":\"\$6]!=\"\")&&(\$10!=\".\")&&(\$16>1)){print \$1\"\\t\"\$2\"\\t\"(\$2+1)\"\\t\"deg[\$1\":\"\$2\":\"\$6]\"\\t\"\$10\"\\t\"\$6\"\\t\"\$13\"\\t\"\$14\"\\t\"\$16\"\\t\"\$17\"\\t\"\$12}}' $file $wt.*.1tr > $file.t" ; done; done

#1c.find piRNAs with 19+ matches between g2..g25, as well as g3..g15 for 5ppm piRNAs, g4..g16 for 10ppm piRNAs, g5..g17 for 50ppm piRNAs, to 1+ ppm pi9/pi17 piRNAs
for file in `ls *.01.t`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FILENAME==\"merged.uniq.piRNAs.cl.g.r.p9p17\"){p9p17[\$1\"\\t\"\$2\"\\t\"\$6]=1}(FILENAME==\"merged.uniq.piRNAs.cl.g.r\"){if((p9p17[\$1\"\\t\"\$2\"\\t\"\$6]!=\"\")&&(\$4>1)){pirnas[substr(\$7,2,25)]=\$4;ab[substr(\$7,2,25)]=\$4;pos[substr(\$7,2,25)]=\$1\",\"\$2\",\"\$3\",\"\$6}}(FILENAME!=\"merged.uniq.piRNAs.cl.g.r\")&&(FILENAME!=\"merged.uniq.piRNAs.cl.g.r.p9p17\"){for(pirna in pirnas){if((substr(\$10,5+2,14)==substr(pirna,2,14))&&(ab[pirna]>5)){print \"215\\t\"pirnas[pirna]\"\\t\"\$0\"\\t\"pos[pirna]}else{if((substr(\$10,5+3,14)==substr(pirna,3,14))&&(ab[pirna]>10)){print \"316\\t\"pirnas[pirna]\"\\t\"\$0\"\\t\"pos[pirna]}else{if((substr(\$10,5+4,14)==substr(pirna,4,14))&&(ab[pirna]>50)){print \"417\\t\"pirnas[pirna]\"\\t\"\$0\"\\t\"pos[pirna]}else{matches=0;for (g=2;g<=25;g++){if(substr(\$10,5+g,1)==substr(pirna,g,1)){matches++}};if(matches>=19){print matches\"\\t\"pirnas[pirna]\"\\t\"\$0\"\\t\"pos[pirna]}}}}}}' merged.uniq.piRNAs.cl.g.r.p9p17 merged.uniq.piRNAs.cl.g.r $file > $file.01cnd" ; done

#all non-p9 and non-p17 cleavage targets
#2a.select deg sites at >=0.1ppm and changing by <2-fold in pi9pi17 mutant
for file in `ls *.d`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(\$4>=-1){split(\$5,a,\",\");if(a[3]>=0.1){\$5=a[3];print}}' $file > $file.01.nt" ; done

#2b.get corresponding transcript data for each putative cleavage product
# *.1tr files are generated in 10c of 5P-seq_data_initial_processing.sh
for wt in WT36 WT54 WT58 WT61; do for file in $wt.*.01.nt; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){deg[\$1\":\"\$2\":\"\$6]=\$5}(FNR<NR){if((deg[\$1\":\"\$2\":\"\$6]!=\"\")&&(\$10!=\".\")&&(\$16>1)){print \$1\"\\t\"\$2\"\\t\"(\$2+1)\"\\t\"deg[\$1\":\"\$2\":\"\$6]\"\\t\"\$10\"\\t\"\$6\"\\t\"\$13\"\\t\"\$14\"\\t\"\$16\"\\t\"\$17\"\\t\"\$12}}' $file $wt.*.1tr > $file.t" ; done; done

#2c.find the ones with 19+ matches between g2..g25 to a 1+ppm non-pi9/pi17 piRNA
#split for parallel computing
for i in $(find . -maxdepth 1 -name "*.01.nt.t"); do split -a 4 -d -l 20000 $i $i.xx. ; done
for file in `ls *.xx.*`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FILENAME==\"merged.uniq.piRNAs.cl.g.r.p9p17\"){p9p17[\$1\"\\t\"\$2\"\\t\"\$6]=1}(FILENAME==\"merged.uniq.piRNAs.cl.g.r\"){if((p9p17[\$1\"\\t\"\$2\"\\t\"\$6]==\"\")&&(\$4>1)){pirnas[substr(\$7,2,25)]=\$4;ab[substr(\$7,2,25)]=\$4;pos[substr(\$7,2,25)]=\$1\",\"\$2\",\"\$3\",\"\$6}}(FILENAME!=\"merged.uniq.piRNAs.cl.g.r\")&&(FILENAME!=\"merged.uniq.piRNAs.cl.g.r.p9p17\"){for(pirna in pirnas){if((substr(\$10,5+2,14)==substr(pirna,2,14))&&(ab[pirna]>5)){print \"215\\t\"pirnas[pirna]\"\\t\"\$0\"\\t\"pos[pirna]}else{if((substr(\$10,5+3,14)==substr(pirna,3,14))&&(ab[pirna]>10)){print \"316\\t\"pirnas[pirna]\"\\t\"\$0\"\\t\"pos[pirna]}else{if((substr(\$10,5+4,14)==substr(pirna,4,14))&&(ab[pirna]>50)){print \"417\\t\"pirnas[pirna]\"\\t\"\$0\"\\t\"pos[pirna]}else{matches=0;for (g=2;g<=25;g++){if(substr(\$10,5+g,1)==substr(pirna,g,1)){matches++}};if(matches>=19){print matches\"\\t\"pirnas[pirna]\"\\t\"\$0\"\\t\"pos[pirna]}}}}}}' merged.uniq.piRNAs.cl.g.r.p9p17 merged.uniq.piRNAs.cl.g.r $file > ${file//.xx./.xy.}.nt_01cnd" ; done
#2d.merge data after parallel computing
for i in $(find . -maxdepth 1 -name "*.xy.0000.*"); do echo $i; first=${i%.xy.0000.*}; last=${i#*.xy.0000.}; cat $first.*.$last > $first.merged.$last; done

#3.merge pi9/pi17 and non-pi9/non-pi17 targets 
#all.4N_7repWT_p9p17_em1_AEFGHI_em2_CDMN.deseq.txt contains TPM and log2FC data for each transcript
for i in `ls *.t.01cnd`; do echo ${i/.t.01cnd/.01both}.20; awk -v n=20 'BEGIN{FS=OFS="\t"}(FNR==NR)&&($11>1){if(($3<0.01)&&($6<0.01)){ss=1}else{ss=0};if(($2 < 0 ? -$2 : $2)<($5 < 0 ? -$5 : $5)){rsq[$7]=$2"\t"ss}else{rsq[$7]=$5"\t"ss}}(FNR<NR)&&(FILENAME~/merged/)&&($1>=n){if(nt[$7]<$2){nt[$7]=$2}}(FNR<NR)&&(FILENAME!~/merged/)&&($1>=n)&&(FILENAME~/.t.01cnd/){if(t[$7]<$2){t[$7]=$2}}(FNR<NR)&&(FILENAME!~/merged/)&&(FILENAME~/.nt.t/){deg[$5]=1}END{for (g in rsq){if((nt[g]!="")&&(t[g]!="")){print g"\t"rsq[g]"\tboth\t"t[g]"\t"nt[g]};if((nt[g]!="")&&(t[g]=="")){print g"\t"rsq[g]"\tnt\t"0"\t"nt[g]};if((nt[g]=="")&&(t[g]!="")){print g"\t"rsq[g]"\tt\t"t[g]"\t"0};if((nt[g]=="")&&(t[g]=="")&&(deg[g]==1)){print g"\t"rsq[g]"\trest\t"0"\t"0};if((nt[g]=="")&&(t[g]=="")&&(deg[g]=="")){print g"\t"rsq[g]"\trestlow\t"0"\t"0}}}' all.4N_7repWT_p9p17_em1_AEFGHI_em2_CDMN.deseq.txt $i ${i/.t.01cnd/.nt.t} ${i/.t.01cnd/.nt.t.merged.nt_01cnd} > ${i/.t.01cnd/.01both}.20; done; done

#4.get detaied deg data
for i in `ls *.01both.20`; do echo $i; awk 'BEGIN{FS=OFS="\t"}(FNR==NR)&&($4=="t"){g[$1]=$2}(FNR<NR)&&(g[$7]!=""){print $7"\t"g[$7]"\t"$2"\t"$6"\t"$11}' $i ${i/.01both.20/.t.01cnd} | sort | uniq > ${i/.01both.20/.deg_details}; done

#5.get number of statistically-significantly increased in each targeted set
for i in `ls *.??both.20`; do awk 'BEGIN{FS=OFS="\t"}($2>0)&&($3==1)&&($4=="t"){ss++}($2>0)&&($3==1)&&($4=="both"){ss-both++}($4=="t"){total++}($4=="both"){total_both++}END{print FILENAME"\t"ss"\t"total"\t"ss/total"\t"ss_both+ss"\t"total+total_both"\t"(ss+ss_both)/(total+total_both)}' $i; done > ss.txt

#6.add GROseq data
#GRO.pergene.mm10.38.92.xls files contains GRO seq coverage per gene
awk 'BEGIN{FS=OFS="\t"}(FNR==NR)&&(FNR>1){gro[$1]=($28+$37+$46)/3"\t"$28"\t"$37"\t"$46}(FNR<NR){print $0"\t"gro[$1] > FILENAME".gro"}'  GRO.pergene.mm10.38.92.xls *.??both.??