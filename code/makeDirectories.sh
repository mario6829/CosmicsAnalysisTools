#!/bin/bash

if [ "x$2" == "x" ]; then
  echo "Usage: deleteDirectories.sh <period> <filelist>"
  echo "(ex. deleteDirectories.sh LHC15a listOfRuns.txt)"
  exit 1
fi

PERIOD=$1
RUNLIST=$2

exec 6>&1
exec 6<$RUNLIST
while read -u 6 nRun
  do
    alien_mkdir /alice/cern.ch/user/s/sitta/anacosmic/$PERIOD/$nRun
    alien_mkdir /alice/cern.ch/user/s/sitta/anacosmic/$PERIOD/$nRun/output
  done
