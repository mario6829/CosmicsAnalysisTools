#!/bin/bash

if [ "x$3" == "x" ]; then
  echo "Usage: mvToOld.sh <period> <filelist> <suffix>"
  echo "(ex. mvToOld.sh LHC15a listOfRuns.txt old)"
  exit 1
fi

PERIOD=$1
RUNLIST=$2
SUFFIX=$3

exec 6>&1
exec 6<$RUNLIST
while read -u 6 nRun
  do
    alien_mv /alice/cern.ch/user/s/sitta/anacosmic/$PERIOD/$nRun /alice/cern.ch/user/s/sitta/anacosmic/$PERIOD/${nRun}_${SUFFIX}
  done
