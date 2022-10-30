#!/bin/sh

# Our AlignmentCounts.txt is missing header. Let's put some header so that it becomes easier when we find the means for gegenes expression levels for each group of read sequences
cat AlignmentCounts.txt | awk -vOFS="\t" '$1=$1; BEGIN { str="Chr Start End Name Description Tco-5053C2 Tco-5110C1 Tco-5132 C1 Tco-5161 C2 Tco-5224WT Tco-5296C2 Tco-5368C2 Tco-5375WT Tco-5220 WT Tco-5572C2 Tco-5578C2 Tco-5592C2 Tco-5615WT Tco-5646C1 Tco-5774WT Tco-5797WT Tco-5814WT Tco-5864C1 Tco-5908WT Tco-5934 WT Tco-6016WT Tco-6028C2 Tco-6037WT Tco-6114C2 Tco-6153C1 Tco-6195 C2 Tco-6225 C2 Tco-6241WT Tco-6228C1 Tco-6307C1 Tco-6320C1 Tco-6395C1 Tco-6442C1 Tco-6475 C2 Tco-6449WT Tco-6520C1 Tco-6580C1 Tco-6600WT Tco-6619C1 Tco-6624C1 Tco-6658C1 Tco-6778C1 Tco-6821C2 Tco-6890WT Tco-6894WT Tco-6950C1"; split(str,arr," "); for(i in arr) printf("%s\t", arr[i]);print}' > AlignmentCountsWH.txt
ls *.txt

# This should calculate the mean for every sample not by group yet
cat AlignmentCountsWH.txt | awk '{ for(i=5;i<=NF;i++) total[i]+=$i ; } 
    END { for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' > MeanBySample.txt
ls *.txt


