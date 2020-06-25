#!/bin/bash
# M.S. 14 jan 16

if [[ "x$1" == "x" ]]; then
  echo "Missing file list. Use '$0 list_file'"
  exit 1
fi

input=$1

runlist=""

while IFS='' read -r runnum
  do
     runlist+=" cosmicAnaRun_$runnum.root"
done < "$input"

echo "Merging the following runs"
cat $input
echo "to cosmicAnaRun_merged.root"
echo

hadd cosmicAnaRun_merged.root $runlist

echo
status=$?

if [[ $status -eq 0 ]]; then
  echo "Done."
else
  echo "There was an error. Check and retry"
fi
