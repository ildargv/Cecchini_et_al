#1.get each target transcript sequence
#mm10.38.92.gtf.fasta.fa.t_g_n.as contains antisense sequence of all mm10 transcripts accompanied by transcript,gene,common names
# p6targets.txt, p9targets.txt, p17targets.txt are lists of ENSMUSG gene names whose abundance changes significantly (FDR<0.01) in piRNA mutants
for file in `ls *targets.txt`; do echo $file; awk 'BEGIN{FS=OFS="\t"}(FNR==NR){genes[$1]=FILENAME}(FNR<NR)&&(FNR%2==1){split(substr($1,2),a,":");if(genes[a[2]]!=""){getline; print > genes[a[2]]"."a[2]"."a[1]"."a[3]".trgt"}}' $file mm10.38.92.gtf.fasta.fa.t_g_n.as ; done

#2a.match with cleaving piRNAs
#merged.all1.piRNAs.cl.g.r.pX files are generated in 11d of SmallRNAseq_data_processing.sh
for file in `ls p6targets.txt.*.trgt`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){split(\$7,a,\";\");pirnas[substr(a[2],1,30)]=\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$12\"\\t\"\$6\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$14\"\\t\"\$15}(FNR<NR){for(pirna in pirnas){for (pos=1;pos<=length(\$1)-28;pos++){matches=0;mmpos=\"\";for (g=2;g<=30;g++){if(substr(\$1,pos+(g-2),1)!=substr(pirna,g,1)){mmpos=mmpos\";0\"}else{matches++;mmpos=mmpos\";1\"}};if(matches>=15){print matches\"\\t\"mmpos\"\\t\"pos\"\\t\"pos+8\"\\t\"length(\$1)-(pos+8)\"\\t\"pirna\"\\t\"pirnas[pirna]}}}}' merged.all1.piRNAs.cl.g.r.p6 $file > $file.all1.gds" ; done

for file in `ls p9targets.txt.*.trgt`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){split(\$7,a,\";\");pirnas[substr(a[2],1,30)]=\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$12\"\\t\"\$6\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$14\"\\t\"\$15}(FNR<NR){for(pirna in pirnas){for (pos=1;pos<=length(\$1)-28;pos++){matches=0;mmpos=\"\";for (g=2;g<=30;g++){if(substr(\$1,pos+(g-2),1)!=substr(pirna,g,1)){mmpos=mmpos\";0\"}else{matches++;mmpos=mmpos\";1\"}};if(matches>=15){print matches\"\\t\"mmpos\"\\t\"pos\"\\t\"pos+8\"\\t\"length(\$1)-(pos+8)\"\\t\"pirna\"\\t\"pirnas[pirna]}}}}' merged.all1.piRNAs.cl.g.r.p9 $file > $file.all1.gds" ; done

for file in `ls p17targets.txt.*.trgt`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){split(\$7,a,\";\");pirnas[substr(a[2],1,30)]=\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$12\"\\t\"\$6\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$14\"\\t\"\$15}(FNR<NR){for(pirna in pirnas){for (pos=1;pos<=length(\$1)-28;pos++){matches=0;mmpos=\"\";for (g=2;g<=30;g++){if(substr(\$1,pos+(g-2),1)!=substr(pirna,g,1)){mmpos=mmpos\";0\"}else{matches++;mmpos=mmpos\";1\"}};if(matches>=15){print matches\"\\t\"mmpos\"\\t\"pos\"\\t\"pos+8\"\\t\"length(\$1)-(pos+8)\"\\t\"pirna\"\\t\"pirnas[pirna]}}}}' merged.all1.piRNAs.cl.g.r.p17 $file > $file.all1.gds" ; done

#2b.keep the ones with deg evidence
#four *.1tr files are generated in 10c of 5P-seq_data_initial_processing.sh
for file in `ls *.gds`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(NR==1){split(FILENAME,f,\".\")}(FNR==NR){sitematches[f[4]\"\\t\"\$5]++;siteinfo[f[4]\"\\t\"\$5\"\\t\"sitematches[f[4]\"\\t\"\$5]]=\$0}(FNR<NR)&&(sitematches[\$13\"\\t\"\$14]!=\"\"){for(i=1;i<=sitematches[\$13\"\\t\"\$14];i++){if(siteprinted[\$13\"\\t\"\$14\"\\t\"i]==\"\"){print siteinfo[\$13\"\\t\"\$14\"\\t\"i]\"\\t\"FILENAME\"\\t\"\$0;siteprinted[\$13\"\\t\"\$14\"\\t\"i]=1}}}' $file *.1tr > $file.deg" ; done

#3a.filter out non-expressed transcripts all gene
#transcript_7reps.txt contains mean fpkm abundance of each transcript in mm10 from 7 reps of wt controls
for file in `ls *.deg`; do echo $file; awk 'BEGIN{FS=OFS="\t"}(FNR==NR){tr[$6]+=$12;g[$9]+=$12}(FNR<NR)&&(tr[$32]/g[$29]>0){print > FILENAME".0"}' transcript_7reps.txt $file; done

#3b.add all 16 permutations of deg data
#16 permutations (*.d files) are generated in 7 of 5P-seq_data_initial_processing.sh
for file in `ls *.deg.0`; do bsub "awk 'BEGIN{FS=OFS=\"\\t\"}(FILENAME~/unique/){split (\$5,ab,\",\"); degab[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$6]=degab[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$6]\";\"ab[3]; degchng[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$6]=degchng[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$6]\";\"\$4}(FILENAME~/gds.deg/){print \$1\"\\t\"\$2\"\\t\"\$20\"\\t\"\$21\"\\t\"\$22\"\\t\"degchng[\$20\"\\t\"\$21\"\\t\"\$22\"\\t\"\$25]\"\\t\"degab[\$20\"\\t\"\$21\"\\t\"\$22\"\\t\"\$25]\"\\t\"\$25\"\\t\"\$27\"\\t\"\$28\"\\t\"\$29\"\\t\"\$30\"\\t\"\$31\"\\t\"\$32\"\\t\"\$33\"\\t\"\$34\"\\t\"\$35\"\\t\"\$36\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12\"\\t\"\$13\"\\t\"\$14\"\\t\"\$15\"\\t\"\$16\"\\t\"\$17\"\\t\"\$18}' *.d  $file > $file.info" ; done

#OUTPUT by COLUMN
#line1	number of matching nucleotides
#line2	matching nucleotides starting form g2: e.g., "";1;1;1;0;1;0;0;1;0;0;1;1;0;1;1;1;1;0;0;1;0;1;0;1;1;0;0;0;0"
#degradome site data
#line3	chromosome
#line4	start
#line5	end
#line6	log2FC mut/wt
#line7	ppm abundance
#line8	strand
#line9	TE
#line10	TE strand
#line11	mm10 ENSMUSG name
#line12	gene_type (e.g., protein_coding)
#line13	gene common name (e.g., Brca2)
#line14	mm10 ENSMUST name
#line15	location in transcript
#line16	transcript annotation
#line17	transcript fpkm
#line18	target site sequence
#line19	piRNA sequence
#piRNA data
#line20	chromosome
#line21	start
#line22	end
#line23	ppm abundance in wt
#line24	piRNA locus of origin
#line25 strand
#line26	ppm abundance in p9
#line27	ppm abundance in p17
#line28	ppm abundance in p9p17
#line29	ppm abundance in p6
#line30	TE
#line31	TE strand

#3c.filter for all 4 reps of wt degradome (= all 16 permutations) to have detectable cleavage product and for pairing to have longest uninterrupted stretches
for file in `ls *.info`; do awk 'BEGIN{FS=OFS="\t"}{n=split($6,chng,";"); if(n>16){chng_down=0;for(i=2;i<=n;i++){if(chng[i]<-2){chng_down++}};if(chng_down==(n-1)){over3=0;current=0;m=split($2,pair,";");for (j=3;j<=m;j++){if((pair[j-1]==1)&&(pair[j]==1)){current++}else{current=0};if(current>3){over3++}};print over3+$1"\t"over3"\t"n-1"\t"$0}}}' $file | sort -k1,1nr > $file.flt; done
