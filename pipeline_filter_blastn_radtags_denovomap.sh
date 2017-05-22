#!/usr/bin/sh

#PIPELINE
#BLASTn Filter
echo "Yes: 1, No: 0"
echo "Run BLASTn filter?"; read runblastn
echo "Run Stacks process_radtags?"; read runprocessradtags
echo "Run Stacks denovo_map?"; read rundenovomap

#DIRECTORIES
echo "ALL DIRECTORIES MUST BE SPECIFIED AS /path/to/dir/"
echo "Directory outputs:"; read outputsdir
echo "Directory for FASTq (*.fastq) raw reads files:"; read fqfilesdir

#create directories for outputs
mkdir "$outputsdir"fafiles    #FASTA FILES DIRECTORY
mkdir "$outputsdir"filteredfq #FASTQ FILTERED FILES DIRECTORY
mkdir "$outputsdir"individuals   #Stacks INDIVIDUALS FILES DIRECTORY
mkdir "$outputsdir"blastnfiles   #BLASTN result FILES DIRECTORY
mkdir "$outputsdir"denovooutput   #DE NOVO map output result FILES DIRECTORY

#FILES
echo "Reference FASTA file:" read reffasta
echo "Individuals population file (format: individual <tab> int - Check Stacks manual for more information):" read individualpop

#PARAMETERS
echo "Final reads length for process_radtags (int):" read finalreadslength
echo "Restriction enzyme used during GBS sequencing for process_radtags:" read restrictionenzyme

############################# START PIPELINE ################################

#BLASTn FILTER STEP

#Convert FASTq raw reads files to FASTa files as they will be used as input for BLASTn step
cat "$individualpop" | cut -f 1 | while read line
do
  echo "converting $line to fasta"
  fastq_to_fasta -Q33 -i "$fqfilesdir""$line".fastq -o "$outputsdir"fafiles/"$line".fa  #convert fastq to fasta using phred33 quality
done

#Creating reference fasta database
makeblastdb -in "$reffasta" -dbtype nucl -out "$reffasta".blastdb

#Blasting step
cat "$individualpop" | cut -f 1 | while read line
do
  echo "Blasting $line against reference"
  blastn -query "$outputsdir"fafiles/"$line".fa -db "$reffasta".blastdb -evalue 1e-20 -outfmt 6 -out "$outputsdir"blastnfiles/"$line".blastn
done

#Creating Filtered FASTq files
cat "$individualpop" | cut -f 1 | while read line
do
  echo "Generating $line filtered FASTq file"
  python blastn2fq.py "$fqfilesdir""$line".fastq "$outputsdir"blastnfiles/"$line".blastn > "$outputsdir"filteredfq/"$line".fastq
done


#### PROCESS_RADTAGS step
cat "$individualpop" | cut -f 1 | while read line
do
  echo "executing process_radtags for $line"
  process_radtags -t "$finalreadslength" -f "$outputsdir"filteredfq/"$line".fastq -i fastq -o "$outputsdir"individuals/ -b "$outputsdir"barcodes/"$line".txt -e "$restrictionenzyme" -r -q
  mv "$outputsdir"individuals/process_radtags.log "$outputsdir"individuals/process_radtags_"$line".log  #rename log file for each fastq
done


#execute denovo_map.pl pipe
denovo_map.pl -t -T 2 -m 10 -M 5 -n 20 -H -b 1 -S -O "$individualpop" -o "$outputsdir"denovooutput/ --samples "$outputsdir"individuals/
