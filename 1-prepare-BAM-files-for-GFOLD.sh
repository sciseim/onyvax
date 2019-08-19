# 1-prepare-BAM-files-for-GFOLD.sh
#!/usr/bin/env bash
# ############################## 
#  PREPARE FOR GFOLD
# ############################## 
 # do this locally on ewha computer hdd
OUTPUTFOLDERNAME=GFOLDoutput
mkdir $OUTPUTFOLDERNAME ;
cd $OUTPUTFOLDERNAME ;
OUTPUTFOLDER=$(pwd)
echo $OUTPUTFOLDER ;

# on IHBIserver
cd /media/sciseim/IHBIgremlin/Inge/onyvax/bam_files/hg19 ;

for d in *; do
cd $d 
cp accepted_hits.bam $OUTPUTFOLDER/$d.bam ; 
cd ..
done
cd $OUTPUTFOLDER ;

