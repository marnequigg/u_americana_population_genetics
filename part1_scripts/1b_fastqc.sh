#make a conda to download the software
conda create -n bioconda
conda activate bioconda
conda install bioconda::fastqc #version 0.12.1
conda install bioconda::multiqc #version 0.4
conda install python=3.5 #version 3.5

#lets check out the quality
nano fastqc.sh

#!bin/bash
output_dir="fastqc_results" #set the output directory
mkdir -p "$output_dir" #create the output directory if it doesnt exist
for file in *_*.fastq.gz; do fastqc -t 40 -o "$output_dir" "$file" 
done
#press ctrl+x then Y then enter to save it

#change file permissions so you can run it
chmod +x fastqc.sh

#time to run it!
screen -L bash fastqc.sh
