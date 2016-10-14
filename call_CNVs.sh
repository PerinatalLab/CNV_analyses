#!/bin/bash

## 1. exporting BAFs and LRRs
# Genome Studio was used
# clustered using Illumina's .egt
# exported GType, Score, LRR and BAF for each sample x SNP
INFILE=/Julius/Documents/ARRAYS/MoBa_BAR/Data/ALL_BAF_logRR.txt
DIR=~/Desktop/MoBa_BAF/
PENNCNV=~/soft/penncnv/

## 2. lifting over to normal coordinates
# cut allele positions from the file
cut -f 1,3,4 ${INFILE} > ${DIR}hg18positions.txt
# reformat
awk -v OFS='\t' 'NR>1{print "chr" $2, $3-1, $3, $1}' ${DIR}hg18positions.txt > ${DIR}hg18positions.pos
# use UCSC's liftOver to remap hg18->hg38
~/soft/liftOver ${DIR}hg18positions.pos ~/soft/inrich/refs/hg18ToHg38.over.chain ${DIR}hg38positions.pos ${DIR}hg38_unmapped.pos -positions
# but DO NOT use the generated .pos files, because they have no SNP names
# use the latest .bedmapped file
awk -v OFS='\t' 'FNR==NR{m[$4]=$3; c[$4]=$1; next}
		$4 in m{print $4, $1, $3, c[$4], m[$4]}' liftOver_lab_workstation_22f4_5039e0.bedmapped ${DIR}hg18positions.pos > ${DIR}hg18to38positions.pos
rm liftOver_lab_workstation_*

# lift the supplied pennCNV files to hg38
awk -v OFS='\t' 'BEGIN{print "Name", "Chr", "Position", "PFB"}
		FNR==NR{c[$1]=$4; p[$1]=$5; next} $1 in p{print $1, substr(c[$1], 4), p[$1], $4}' ${DIR}hg18to38positions.pos ${PENNCNV}lib/hhall.hg18.pfb > ${PENNCNV}lib/hh660.hg38.pfb
awk -v OFS='\t' 'BEGIN{print "Name", "Chr", "Position", "GC"}
		FNR==NR{c[$1]=$4; p[$1]=$5; next} $1 in p{print $1, substr(c[$1], 4), p[$1], $4}' ${DIR}hg18to38positions.pos ${PENNCNV}lib/hhall.hg18.gcmodel > ${PENNCNV}lib/hh660.hg38.gcmodel

## 3. start splitting the superfile and analyzing
# cut a slice of 100 people, keeping the initial columns
for c in {1..34}
do
  SUBDIR=${DIR}calls${c}/
  mkdir ${SUBDIR}
  CUTSTOP=$((c*400+4))
  CUTSTART=$((CUTSTOP-399))
  cut -f 1-4,${CUTSTART}-${CUTSTOP} ${INFILE} > ${SUBDIR}SPLIT_BAF_logRR.txt

  # cut each person into a separate file
  perl ${PENNCNV}kcolumn.pl ${SUBDIR}SPLIT_BAF_logRR.txt \
	split 4 -heading 4 -tab \
	-name_by_header -out ${SUBDIR}sample

  # analyze
  perl ${PENNCNV}detect_cnv.pl -test \
	-hmm ${PENNCNV}lib/hhall.hmm \
	-pfb ${PENNCNV}lib/hh660.hg38.pfb \
	-gcmodel ${PENNCNV}lib/hh660.hg38.gcmodel \
	${SUBDIR}sample* \
	-log ${DIR}SPLIT_CALLS_${c}.log \
	-out ${DIR}SPLIT_CALLS_${c}.rawcn

  rm ${SUBDIR}sample*
  rm ${SUBDIR}SPLIT_BAF_logRR.txt
done

# 4. QC (note that no QC of samples or SNPs was done before here)
cat ${DIR}SPLIT_CALLS*.rawcn > ${DIR}ALL_CALLS.rawcn
cat ${DIR}SPLIT_CALLS*.log > ${DIR}ALL_CALLS.log
## default params are SD(LRR)<0.3, drift(BAF)<0.01, WF<0.05 
perl ${PENNCNV}filter_cnv.pl ${DIR}ALL_CALLS.rawcn \
	-qclogfile ${DIR}ALL_CALLS.log \
	-qcpassout ${DIR}ALL_CALLS_defaultqc.txt \
	-qcsumout ${DIR}ALL_CALLS_defaultqc.log \
	-output ${DIR}ALL_CALLS_defaultqc.rawcn

## make a random subset of SNPs to use for PCA
awk -v FS="\t" 'NR==1{o=""; for(i=8; i<=NF; i+=4){o=o FS $i}; print o; next}
	NR % 100==1{for(i=8; i<=NF; i+=4){
		printf("%.3f\t", $i)
	}; printf("\n")}' ${INFILE} > ${DIR}SUBSET_logRR.txt
## delete first space manually in vim and load to R.
## without the SD(LRR) filter, the qc command should look like this:
perl ${PENNCNV}filter_cnv.pl ${DIR}ALL_CALLS.rawcn \
	-qclogfile ${DIR}ALL_CALLS.log \
	-qclrrsd 2 \
	-qcpassout ${DIR}ALL_CALLS_nosdfilter.txt \
	-qcsumout ${DIR}ALL_CALLS_nosdfilter.log \
	-output ${DIR}ALL_CALLS_nosdfilter.rawcn

## restructure for R
awk 'BEGIN{print "CHR", "START", "STOP", "NUMSNP", "LENGTH", "CN", "ID", "STARTSNP", "ENDSNP"}
	{gsub("=", FS);
	split($1, a, "-");
	gsub(":", FS, a[1]);
	gsub(",", "", $7);
	print substr(a[1], 4), a[2], $3, $5, $7, gensub(/.*sample\./, "", 1, $8), $10, $12}' ${DIR}ALL_CALLS_nosdfilter.rawcn > ${DIR}ALL_CALLS_nosdfilter_nice.rawcn

## find longest possible CNP regions
sort -k1 -k2 -k3 -n ${DIR}ALL_CALLS_nosdfilter_nice.rawcn | awk -v f=${DIR}ALL_nosdfilter_cnpreg_long.txt 'NR==1{print $0, "REG";
	print "REG", "CHR", "START", "STOP">>f; next}
	NR==2{c=$1; b=$2; e=$3; n=1}
	$1==c && $2<=e{print $0, n; if($3>e)e=$3; next}
	{print n++, c, b, e>>f; c=$1; b=$2; e=$3; next}
	END{print n, c, b, e>>f}
' > ${DIR}ALL_CALLS_nosdfilter_wreg_long.txt

## OR: assign closest gene within 5Mbp 
perl ${PENNCNV}scan_region.pl ${DIR}ALL_CALLS_nosdfilter.rawcn \
	~/soft/inrich/refs/ucsc_ref_hg38.map \
	-refgene \
	-expandmax 5m > ${DIR}ALL_CALLS_nosdfilter_genes.rawcn

## reformat and load to R for regression
awk '{split($4, c, "="); split($5, a, "sample."); split($8, g, ","); for(i in g){print a[2], g[i], c[2]}}' ${DIR}ALL_CALLS_nosdfilter_genes.rawcn > ${DIR}ALL_nosdfilter_genes.txt
