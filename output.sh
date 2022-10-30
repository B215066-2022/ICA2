#!/bin/bash

# Our AlignmentCounts.txt is missing header. Let's put some header so that it becomes easier when we find the means for gegenes expression levels for each group of read sequences
cat AlignmentCounts.txt | awk -vOFS="\t" '$1=$1; BEGIN { str="Chromosome Start End Name Description Tco5053_Clone1 Tco511_Clone1 Tco5132_Clone1 Tco5161_Clone2 Tco5224_WT Tco5296_Clone2 Tco5368_Clone2 Tco5375_WT Tco5220_WT Tco5572_Clone2 Tco5578_Clone2 Tco5592_Clone2 Tco5615_WT Tco5646_Clone1 Tco5774_WT Tco5797_WT Tco5814_WT Tco5864_Clone1 Tco5908_WT Tco5934_WT Tco6016_WT Tco6028_Clone2 Tco6037_WT Tco6114_Clone2 Tco6153_Clone1 Tco6195_Clone2 Tco6225_Clone2 Tco6241_WT Tco6228_Clone1 Tco6307_Clone1 Tco6320_Clone1 Tco6395_Clone1 Tco6442_Clone1 Tco6475_Clone2 Tco6449_WT Tco6520_Clone1 Tco6580_Clone1 Tco660_WT Tco6619_Clone1 Tco6624_Clone1 Tco6658_Clone1 Tco6778_Clone1 Tco6821_Clone2 Tco6890_WT Tco6894_WT Tco6950_Clone1"; split(str,arr," "); for(i in arr) printf("%s\t", arr[i]);print}' > /home/$USER/ICA1/Mapping/AlignmentCountsWH.txt
ls *.txt

# This should calculate the mean for every sample but not by group yet
cat AlignmentCounts.txt | awk '{ for(i=5;i<=NF;i++) total[i]+=$i ; } 
    END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' > /home/$USER/ICA1/Mapping/MeanBySampleName.txt
ls *.txt

# Put some header to our group data
# cat MeanBySampleName.txt | awk -vOFS="\t" '$1=$1; BEGIN { str="Chromosome Start End Name Description Tco5053_Clone1 Tco511_Clone1 Tco5132_Clone1 Tco5161_Clone2 Tco5224_WT Tco5296_Clone2 Tco5368_Clone2 Tco5375_WT Tco5220_WT Tco5572_Clone2 Tco5578_Clone2 Tco5592_Clone2 Tco5615_WT Tco5646_Clone1 Tco5774_WT Tco5797_WT Tco5814_WT Tco5864_Clone1 Tco5908_WT Tco5934_WT Tco6016_WT Tco6028_Clone2 Tco6037_WT Tco6114_Clone2 Tco6153_Clone1 Tco6195_Clone2 Tco6225_Clone2 Tco6241_WT Tco6228_Clone1 Tco6307_Clone1 Tco6320_Clone1 Tco6395_Clone1 Tco6442_Clone1 Tco6475_Clone2 Tco6449_WT Tco6520_Clone1 Tco6580_Clone1 Tco660_WT Tco6619_Clone1 Tco6624_Clone1 Tco6658_Clone1 Tco6778_Clone1 Tco6821_Clone2 Tco6890_WT Tco6894_WT Tco6950_Clone1"; split(str,arr," "); for(i in arr) printf("%s\t", arr[i]);print}' > /home/$USER/ICA1/Mapping/MeanBySampleName.txt
# ls *.txt


