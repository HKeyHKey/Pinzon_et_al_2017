#!/bin/sh

if test "$1" = ""
then echo "Please enter seed definition (e.g., '2-8')."
     read seed_def
else seed_def=$1
fi

input=seeds_to_be_analyzed_mouse-centered_$seed_def
reference=mmu

while test `grep -vc ' processed$' $input` -ne 0
do time_stamp=`date +%s`
   grep -m 1 -v 'processed' $input > todo_$time_stamp
   seed=`awk '{print $1}' todo_$time_stamp`
   clade=`awk '{print $2}' todo_$time_stamp`
   sed -i '/^'$seed' '$clade'$/ s|$| processed|' $input

   ./Module_extracts_hairpin_from_seed.pl $seed species_list_$clade $reference
   case "$reference" in "hsa") mv Hits_for_`echo $seed | tr T U`_in_additional_genomes.txt Slice_$seed_def'_human-centered/';;
                        "mmu") mv Hits_for_`echo $seed | tr T U`_in_additional_genomes.txt Slice_$seed_def'_mouse-centered/';;
   esac
done

