#!/bin/bash
# M.S. 14 jan 16, 22 jan 2016

if [[ "x$2" == "x" ]]; then
  echo "Usage: $0 <period> <filelist>"
  echo "(ex. $0 LHC15a listOfRuns.txt)"
  exit 1
fi

period=$1
input=$2

rm -f list.txt

while IFS='' read -r runnum
  do
     alien_find /alice/cern.ch/user/s/sitta/anacosmic/$period/$runnum interaction_events.out | grep interaction >list.txt

     echo
     echo "Getting interaction_events for run $runnum (`wc -l list.txt | awk '{print $1}'` files)"

     touch interaction_events_$runnum.out

# alien_cat seems not to work, so copy locally and then cat local copy
     while read -r evfile
       do
	 alien_cp alien://$evfile /tmp
	 cat /tmp/interaction_events.out >>interaction_events_$runnum.out
	 rm -f /tmp/interaction_events.out
     done < "list.txt"

     rm -f list.txt

###
     alien_find /alice/cern.ch/user/s/sitta/anacosmic/$period/$runnum highnumtracks.out | grep highnumtrack >list.txt

     echo
     echo "Getting highnumtracks.out for run $runnum (`wc -l list.txt | awk '{print $1}'` files)"

     touch highnumtracks_$runnum.out

# alien_cat seems not to work, so copy locally and then cat local copy
     while read -r evfile
       do
	 alien_cp alien://$evfile /tmp
	 cat /tmp/highnumtracks.out >>highnumtracks_$runnum.out
	 rm -f /tmp/highnumtracks.out
     done < "list.txt"

     rm -f list.txt

done < "$input"

echo "Done."
