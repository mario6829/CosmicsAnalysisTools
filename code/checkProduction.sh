#!/bin/bash

firstArg=$1
if [[ "${firstArg:0:1}" == "-" ]]
  then
    option=${firstArg:1:2}
    lhcPeriod=$2
    jobList=$3
  else
    option=""
    lhcPeriod=$1
    jobList=$2
  fi

if [[ "x$jobList" == "x" || "x$option" == "xh" ]]; then
  echo "Usage: $0 [-v|-h] <period> <joblist>"
  echo "(ex. $0 LHC15a listOfRuns.jobs)"
  echo "-v : verbose - gives details on single master jobs"
  echo "-c : check production files"
  echo "-h : this help"
  exit 1
fi

sed 1d $jobList | while read runJob
  do
    runNum=`echo $runJob | awk '{print $1}'`
    jobId=`echo $runJob | awk '{print $2}'`

    masterJobStatus=`alien_ps -j $jobId | awk '{print $6}'`
    masterJobStatus=${masterJobStatus:0:1}
    case $masterJobStatus in
      D)
        jobStatus="Done"
        ;;
      I)
        jobStatus="Inserting"
        ;;
      S)
        jobStatus="Split"
        ;;
      *)
        jobStatus="Unknown"
        ;;
    esac

    totalJobs=`alien_ps -m $jobId | wc -l`

    if [[ "$option" != "c" ]]; then
      echo "Master job "$jobId" (run "$runNum") is in state "$jobStatus

      if [[ "$masterJobStatus" != "D" && "$masterJobStatus" != "S" ]]; then
        continue
      fi

      if [[ "$option" == "v" ]]; then
        jobsDone=`alien_ps -m $jobId | grep D | wc -l`
        jobsRun=`alien_ps -m $jobId | grep R | wc -l`
        jobsWait=`alien_ps -m $jobId | grep W | wc -l`
        jobsErr=`alien_ps -m $jobId | grep E | wc -l`

        echo -n "   $totalJobs jobs --"
        if [[ "$masterJobStatus" == "S" ]]; then
          echo -n " $jobsRun running -- $jobsWait waiting --"
        fi
        echo " $jobsDone done -- $jobsErr in error"
      fi

    else

      rootFiles=`alien_find /alice/cern.ch/user/s/sitta/anacosmic/$lhcPeriod/$runNum/output/000 "*.root" | head -n -2 | wc -l`
      if [ $totalJobs -eq $rootFiles ]; then
        prodMessage="OK"
      else
        prodMessage="Missing files!"
      fi

      echo -n "Master job "$jobId" (run "$runNum") : "
      echo "$totalJobs jobs -- $rootFiles Root files -- $prodMessage"
    fi

  done
