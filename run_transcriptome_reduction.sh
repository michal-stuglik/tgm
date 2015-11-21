#!/bin/bash -ve

#################################################################################
###########################   Transcriptome gene modeller #######################
#################################################################################

# general settings
WORKING_DIR=$(pwd)
THREADS_NUM=2
TIME_START=`date +%s`

# INPUT file
CAP_IN='sample-data/reads10k.fasta'


# CAP3 settings
CAP_PARAM_O=40  #-o  N  specify overlap length cutoff > 15 (40)
CAP_PARAM_P=99  #-p  N  specify overlap percent CDHIT_IDENTITY cutoff N > 65 (90)

# cd-hit settings
CDHIT_IDENTITY=0.95
CDHIT_WORDSIZE=8
CDHIT_MAXMEM=2000

# clustering stage settings
CLUSTER_IDENTITY=0.96 # identity: the alignment identity between sequences was at least 96%.
CLUSTER_ALIGN2SHORTEST=0.7 #  the total length(fraction) of all high-scoring pairs (HSP) for shorter sequences covered at least 70% of the total sequence length;

# alignment assesment
ASSESMENT_IDENTITY=96

# DIR & file name settings
CAP_OUT=after_cap3.fa
TGM_OUTDIR=TGM_out
BLAST_DB=selfblast_db
BLAST_OUT=blastout.txt
RENAME_OUT=after_cap_cdhit_rename.fasta
CDHIT_OUT=after_cap3_cdhit.fa

CLUSTERING_OUT=clusters_out
ALIGNMENT_OUT=alignment_out
ASSESMENT_OUT=assesment_out
ALIGNMENT_CONSENSUS_OUT=alignment_consensus_out

mkdir -p $TGM_OUTDIR
cd $TGM_OUTDIR

#################################################################################
#####################################   pipeline    #############################
#################################################################################


echo -e '\n\n------------------------- 1. cap3 Assembly  --------------------------'
# symlink to input sequences
ln -fs $WORKING_DIR'/'$CAP_IN $WORKING_DIR'/'$TGM_OUTDIR'/cap3_seq_in'
# cap3 assembly
cap3 cap3_seq_in -o $CAP_PARAM_O -p $CAP_PARAM_P > cap3.log
# joining results from cap3
cat cap3_seq_in.cap.contigs cap3_seq_in.cap.singlets > $CAP_OUT


echo -e '\n\n---------------------- 2. cd-hit-est pipline  ------------------------'
cdhit-est -i $CAP_OUT -o $CDHIT_OUT  -c $CDHIT_IDENTITY -n $CDHIT_WORDSIZE -M $CDHIT_MAXMEM


echo -e '\n\n------------------- 3. rename sequences in python-----------------------'
python $WORKING_DIR'/'renameSeqInFasta.py -f $CDHIT_OUT -o $RENAME_OUT


echo -e '\n\n---------------------- 4. self-megablast  ---------------------------'
# make database
makeblastdb -in $RENAME_OUT -dbtype nucl -title $BLAST_DB -hash_index -out $BLAST_DB
# mega-blast
blastn -task megablast -db $BLAST_DB -query $RENAME_OUT -outfmt 6 -out $BLAST_OUT -num_threads $THREADS_NUM


echo -e '\n\n--------------------------  5. clustering  ---------------------------'
# remove it in case dir exists
rm -r $CLUSTERING_OUT 2> /dev/null
mkdir -p $CLUSTERING_OUT
python $WORKING_DIR'/'cluster_contigs_after_blast.py -b $BLAST_OUT -f $RENAME_OUT -o $CLUSTERING_OUT -v $CLUSTER_IDENTITY -k $CLUSTER_ALIGN2SHORTEST


echo -e '\n\n---------------  6. clusters multiple-alignment  ---------------------'
python $WORKING_DIR'/'align_all_fasta_in_folder.py -f $CLUSTERING_OUT -p $THREADS_NUM -s "cluster" -o $ALIGNMENT_OUT


echo -e '\n\n------------------  7. clusters assesments  ------------------------'
# remove it in case dir exists
rm -r $ASSESMENT_OUT 2> /dev/null
mkdir -p $ASSESMENT_OUT
python  $WORKING_DIR'/'align_assesment.py -f $ALIGNMENT_OUT -p $THREADS_NUM -s "cluster" -i $ASSESMENT_IDENTITY -o $ASSESMENT_OUT


echo -e '\n\n--------------  8. consensus generating, seq cleaning  ------------------'

## remove it in case dir exists
rm -r $ALIGNMENT_CONSENSUS_OUT 2> /dev/null
mkdir -p $ALIGNMENT_CONSENSUS_OUT
python $WORKING_DIR'/'consensus_and_clean.py -f $ALIGNMENT_OUT -s "cluster" -o $ALIGNMENT_CONSENSUS_OUT


echo -e '\n\n--------------------  9. Collecting all sequences    ------------------------'
## remove it in case dir exists
rm reference_transcriptome.fa 2> /dev/null
cat $ALIGNMENT_CONSENSUS_OUT/* $CLUSTERING_OUT'/'all_non_cluster.fasta > reference_transcriptome.fa


# TIME
TIME_END=`date +%s`
TIME_TOTAL=$(expr $TIME_END - $TIME_START)
TIME_TOTAL_MINUTES=$(expr $TIME_TOTAL / 60)
echo 'Total time [m]: '$TIME_TOTAL_MINUTES
