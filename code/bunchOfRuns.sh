#!/bin/bash

if [[ "x$2" == "x" ]]; then
  echo "Usage: $0 <period> <filelist>"
  echo "(ex. $0 LHC15a listOfRuns.txt)"
  exit 1
fi

lhcPeriod=$1
runList=$2

jobTable=`echo $runList | awk -F "." '{print $1}'`.jobs

exec 6>&1
exec 6<$runList

rm -f jobs.tmp
echo " Run    MasterJob" >jobs.tmp

while read -u 6 nRun
  do
    aliroot -b -q -x runGrid.C\($nRun,\"F\",\"$lhcPeriod\"\);
    rm *.so *.d *.xml anaCosmic* 
    rm AliCosmics_cxx_ACLiC_dict.cxx
    rm AliCosmics_cxx_ACLiC_dict.h
    rm AliCosmics_cxx_ACLiC_linkdef.h
    rm AliCosmics_cxx_ACLiC_map.in
  done

rm -f $jobTable
mv jobs.tmp $jobTable
