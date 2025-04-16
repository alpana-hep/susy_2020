#!/bin/bash

# List of UL campaigns
UL_YEARS=("Summer20UL18" "Summer20UL17" "Summer20UL16" "Summer20UL16APV")

# Common path and filename parts
OUTDIR="/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/kalpana/SkimsUL_June2023/FR"
OUTPREFIX_BASE="skimmed_{YEAR}_ZLLGJets_MonoPhoton_PtG-15to130"
INBASE="/store/user/lpcsusyphotons/SoftPhoton/kalpana/ULSkims_June23/v1_WithMakeClass"

# Batch configuration
BATCH_SIZE=100

# Loop over each year
for YEAR in "${UL_YEARS[@]}"; do
  echo "Processing for $YEAR..."

  # Adjusted output prefix and input pattern
  OUTPREFIX=$(echo $OUTPREFIX_BASE | sed "s/{YEAR}/$YEAR/")
  INPATTERN="phoID_loose_runList_${YEAR}_ZLLGJets_MonoPhoton_PtG-15to130"

  # Fetch input files
  FILES=$(xrdfs root://cmseos.fnal.gov ls -u $INBASE | grep "$INPATTERN")
  readarray -t FILE_ARRAY <<< "$FILES"

  TOTAL=${#FILE_ARRAY[@]}
  BATCH_NUM=0

  # Loop over files in batches
  for ((i=0; i<$TOTAL; i+=$BATCH_SIZE)); do
    OUTFILE="${OUTDIR}/${OUTPREFIX}_batch${BATCH_NUM}.root"
    BATCH_FILES="${FILE_ARRAY[@]:$i:$BATCH_SIZE}"

    echo "  Creating: $OUTFILE"
    hadd -fk $OUTFILE $BATCH_FILES

    ((BATCH_NUM++))
  done
done
