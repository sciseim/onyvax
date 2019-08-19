# GFOLD.sh
#!/usr/bin/env bash
# ############################## 
#  GFOLD TIME
# ############################## 
# note: standard gsl not enough on Ubuntu to compile gfold: sudo apt-get install libgsl0-dev
#
# mv 20_2004.bam AQ0420.bam
# mv 15_2004.bam AQ0415.bam
# mv 11_2004.bam AQ0411.bam
# mv 96_2003.bam AQ0396.bam

# create counts
for BAMFILE in *.bam; do
SAMPLENAME=${BAMFILE/.bam/}
samtools view $BAMFILE | /media/sciseim/working/onyvax/GFOLDoutput/gfold/gfold.V1.1.2/gfold count -ann genes.gtf -tag stdin -o $SAMPLENAME.read_cnt; done  # here hg19

# compare all to RWPE1
for GFOLDCOUNTFILE in *.read_cnt; do
SAMPLENAME=${GFOLDCOUNTFILE/.read_cnt/}
/media/sciseim/working/onyvax/GFOLDoutput/gfold/gfold.V1.1.2/gfold diff -s1 RWPE1 -s2 $SAMPLENAME -suf .read_cnt -o RWPE1-VS-$SAMPLENAME.diff;done

