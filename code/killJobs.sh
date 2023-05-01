#!/bin/bash

if [ "x$1" == "x" ]; then
  echo "Usage: killJobs.sh <joblist>"
  echo "(ex. killJobs.sh listOfRuns.jobs)"
  exit 1
fi

JOBLIST=$1

sed 1d $JOBLIST | while read runJob
  do
    runNum=`echo $runJob | awk '{print $1}'`
    jobId=`echo $runJob | awk '{print $2}'`

    echo "Killing master job $runJob (Run $runNum)"
    alien.py kill $jobId
  done
