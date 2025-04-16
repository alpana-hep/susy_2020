#!/bin/bash

# Base output path
OUTDIR="/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/kalpana/SkimsUL_June2023/FR"
OUTPREFIX="skimmed_Summer20UL18_ZLLGJets_MonoPhoton_PtG-15to130_batch"

# Input file pattern via xrdfs
FILES=$(xrdfs root://cmseos.fnal.gov ls -u /store/user/lpcsusyphotons/SoftPhoton/kalpana/ULSkims_June23/v1_WithMakeClass \
  | grep 'phoID_loose_runList_Summer20UL18_ZLLGJets_MonoPhoton_PtG-15to130')

# Convert to array
readarray -t FILE_ARRAY <<< "$FILES"

# Total files and batch size
TOTAL=${#FILE_ARRAY[@]}
BATCH_SIZE=100
BATCH_NUM=0

for ((i=0; i<$TOTAL; i+=$BATCH_SIZE)); do
  # Define output file name
  OUTFILE="${OUTDIR}/${OUTPREFIX}${BATCH_NUM}.root"

  # Prepare list of files for this batch
  BATCH_FILES="${FILE_ARRAY[@]:$i:$BATCH_SIZE}"

  echo "Creating: $OUTFILE"
  hadd -fk $OUTFILE $BATCH_FILES

  ((BATCH_NUM++))
done
