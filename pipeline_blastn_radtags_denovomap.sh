#!/usr/bin/sh
#Execute Stacks process_radtags and denovo_map.pl pipelines, including the option of filtering raw reads based that match to a reference genome
#Dependencies: Stacks, Blast, and FastX-toolkit

#PIPELINE
#BLASTn Filter
echo "Yes: 1, No: 0"
echo "Run BLASTn filter?"; read runblastn
echo "Run Stacks process_radtags?"; read runprocessradtags
echo "Run Stacks denovo_map?"; read rundenovomap

#MASTER VARIABLES - DIRECTORIES AND FILES
echo "ALL DIRECTORIES MUST BE SPECIFIED AS /path/to/dir/"
echo "Enter the path for the outputs directory:"; read outputsdir
echo "Enter the path to the Individuals population file (format: individual <tab> int - Check Stacks manual for more information):"; read individualpop

#CONFIGURATION FOR BLASTN FILTER STEP
if [ "$runblastn" = "1" ]; then
echo "For the BLASTn Filter Step"
echo "Enter directory for FASTq (*.fastq) raw reads files:"; read fqfilesdir
echo "Reference FASTA file:"; read reffasta
echo "Creating blastnfiles, fafiles and filterdfq directories in the output directory"
mkdir "$outputsdir"fafiles    #FASTA FILES DIRECTORY
mkdir "$outputsdir"filteredfq #FASTQ FILTERED FILES DIRECTORY
mkdir "$outputsdir"blastnfiles   #BLASTN result FILES DIRECTORY
else
echo "BLASTN FILTER STEP WONT BE EXECUTED"
fi

#CONFIGURATION FOR PROCESS_RADTAGS STEP
if [ "$runprocessradtags" = "1" ]; then
echo "For the PROCESS_RADTAGS  Step"
echo "Final reads length for process_radtags (int):"; read finalreadslength
echo "Restriction enzyme used during GBS sequencing for process_radtags:"; read restrictionenzyme
echo "Path to the barcodes file (format: seq_barcode <tab> individual):"; read barcodesfile
echo "Name for the output of Process_Radtags (processed individuals):"; read processed_individuals
echo "Creating individuals and barcodes directories in the output directory"
mkdir "$outputsdir""$processed_individuals"   #Stacks individuals FILES DIRECTORY
mkdir "$outputsdir"barcodes
else
echo "PROCESS_RADTAGS STEP WONT BE EXECUTED"
fi

#CONFIGURATION FOR DE NOVO MAP STEP
if [ "$rundenovomap" = "1" ]; then
echo "For the DE NOVO MAP  Step"
echo "Enter a name for this run to create a directory to write the outputs from DE NOVO MAP step:"; read denovooutput
echo "Creating denovoouput and barcodes directories in the output directory"
mkdir "$outputsdir""$denovooutput"   #DE NOVO map output result FILES DIRECTORY
mkdir "$outputsdir"barcodes
echo "Set the value (int) for the parameter -m:"; read mvalue
echo "Set the value (int) for the parameter -M:"; read mmvalue
echo "Set the value (int) for the parameter -n:"; read nvalue
echo "Set the number (int) of threads to use -T:"; read nthreads
echo "This pipeline will not allow SNP calling from secondary reads!"
if [ "$barcodesfile" = "" ]; then
echo "Path to the barcodes file (format: seq_barcode <tab> individual):"; read barcodesfile
fi
if [ "$processed_individuals" = "" ]; then
echo "Name for the output of Process_Radtags (processed individuals):"; read processed_individuals
fi

else
echo "DE NOVO MAP STEP WONT BE EXECUTED"
fi

############################# START PIPELINE ################################
#BLASTn FILTER STEP
if [ "$runblastn" = "1" ]; then
echo "Starting BLASTn Filter Step"

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
fi #end of BLASTn FILTER STEP


########################## PROCESS_RADTAGS step ################################

if [ "$runprocessradtags" = "1" ]; then
echo "Starting PROCESS_RADTAGS  Step"

#prepare barcodes files
cat "$barcodesfile" | while read line
do
  name=`echo "$line" | cut -f 2`
  echo "$line" > "$outputsdir"barcodes/"$name".txt
done


cat "$individualpop" | cut -f 1 | while read line
do
  echo "executing process_radtags for $line"
  process_radtags -t "$finalreadslength" -f "$outputsdir"filteredfq/"$line".fastq -i fastq -o "$outputsdir""$processed_individuals"/ -b "$outputsdir"barcodes/"$line".txt -e "$restrictionenzyme" -r -q
  mv "$outputsdir""$processed_individuals"/process_radtags.log "$outputsdir""$processed_individuals"/process_radtags_"$line".log  #rename log file for each fastq
done
fi

############################## DE NOVO MAP step #################################
if [ "$rundenovomap" = "1" ]; then
echo "Starting PROCESS_RADTAGS  Step"
denovo_map.pl -t -T "$nthreads" -m "$mvalue" -M "$mmvalue" -n "$nvalue" -H -b 1 -S -O "$individualpop" -o "$outputsdir""$denovooutput"/ --samples "$outputsdir""$processed_individuals"/
fi