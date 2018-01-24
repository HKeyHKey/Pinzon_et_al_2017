#!/bin/sh

if test "$1" = ""
then echo "Please enter 3-letter code for reference species ('hsa' or 'mmu')."
     read ref
else ref=$1
fi
if test "$2" = ""
then echo "Please enter seed definition ('2-7' or '2-8')."
     read seed_def
else seed_def=$2
fi

case "$ref" in "hsa") clade_list='Hominidae Catarrhini Boreoeutheria Euteleostomi';directory=$seed_def'_human-centered';regular_name='human';if test "$seed_def" = "2-7";then output_Blaise='Conservation_miRNA_targets_seed_2-7';else output_Blaise='Conservation_miRNA_targets';fi;;
    "mmu") clade_list='Murinae Boreoeutheria Euteleostomi';directory=$seed_def'_mouse-centered';regular_name='mouse';if test "$seed_def" = "2-7";then output_Blaise='Conservation_miRNA_targets_mouse_centered_seed_2-7';else output_Blaise='Conservation_miRNA_targets_mouse_centered';fi;; # It's not worth analyzing the "Muroidea" clade: the only additional species (compared to Murinae) is the Chinese hamster (genome assembly: criGri1), which is not in the whole genome alignment ".maf" files anyway
esac


echo "Clade_of_seed Analyzed_clade Seed Number_of_seed_matches_in_UTR_alignment Number_of_seed_matches_in_UTR_alignment_outside_of_analyzed_clade" > Control_$directory'.dat'
for clade1 in `echo $clade_list` # this is the clade where the seed belongs
do for clade2 in `echo $clade_list` # these are the clades where we will look for conserved seed matches
   do for seed in `grep ' name: '$clade1' ' /mnt/data/home/herve.seitz/Blaise/miRNA_absence_verification_by_Herve/Seed_selection/out_$seed_def'_'$ref | awk '{print $1}'`
      do withU=`echo $seed | tr T U`
          if test `cat /mnt/data/home/herve.seitz/Blaise/miRNA_absence_verification_by_Herve/Moulinette_$directory'.dat' | grep 'Hits_for_'$withU'_' | awk '{print $1}'` -eq 0 # Check that my (HMMer or Blast) + RNAsubopt search did not find potential orthologs for that miRNA family outside the clade of interest            
          then for species in `cat /mnt/data/home/herve.seitz/Blaise/miRNA_absence_verification_by_Herve/species_outside_$clade2`
               do assembly_list=`grep -w $species /mnt/data/home/herve.seitz/Blaise/miRNA_absence_verification_by_Herve/Seed_selection/Full_correspondance_table | awk -F ',' '{print $2}'`
                  for assembly in `echo $assembly_list`
                  do if test `grep -c ','$assembly',' ../Species_in_100species_alignment.txt` -eq 1
                     then grep -v '^[ACGT][ACGT]* *$' /mnt/data/home/herve.seitz/Blaise/$output_Blaise/sites_db/sites_in_UTR_all/$seed'.txt' | grep -w $assembly
                     fi
                 done
               done | sort | uniq -c > $clade1'-specific_'$seed_def'_seed_families_'$regular_name'-centered/Control_'$clade2'_'$seed'_matches.dat'
               echo $clade1 $clade2 $seed `grep -v '^[ACGT][ACGT]* *$' /mnt/data/home/herve.seitz/Blaise/$output_Blaise/sites_db/sites_in_UTR_all/$seed'.txt' | wc -l` `cat $clade1'-specific_'$seed_def'_seed_families_'$regular_name'-centered/Control_'$clade2'_'$seed'_matches.dat' | wc -l` >> Control_$directory'.dat'
         fi
    done
  done
done

#Below: analysis of "non-seeds" (heptamers not found at nt 2-8 of any vertebrate miRNA in miRBase v21, and hexamers not found at nt 2-7):
for species in `grep ';Vertebrata;' /mnt/data/home/herve.seitz/Blaise/miRBase_data/organisms.txt | awk '{print $1}'`;do grep -A 1 '^> *'$species'\-' /mnt/data/home/herve.seitz/Blaise/miRNA_absence_verification_by_Herve/matureJun14.fa | grep -v '\-\-';done | grep -v '^>' | cut -c $seed_def | sort | uniq | tr U T > vertebrate_seeds_$directory
if test "$seed_def" = "2-7"
then for nt1 in A C G T;do for nt2 in A C G T;do for nt3 in A C G T;do for nt4 in A C G T;do for nt5 in A C G T;do for nt6 in A C G T;do echo $nt1''$nt2''$nt3''$nt4''$nt5''$nt6;done;done;done;done;done;done > all_possible_$directory'_seeds'
fi
if test "$seed_def" = "2-8"
then for nt1 in A C G T;do for nt2 in A C G T;do for nt3 in A C G T;do for nt4 in A C G T;do for nt5 in A C G T;do for nt6 in A C G T;do for nt7 in A C G T;do echo $nt1''$nt2''$nt3''$nt4''$nt5''$nt6''$nt7;done;done;done;done;done;done;done > all_possible_$directory'_seeds'
fi
cat all_possible_$directory'_seeds' all_possible_$directory'_seeds' vertebrate_seeds_$directory | sort | uniq -c | grep '^ *2 ' | awk '{print $2}' > non_seeds_$directory
for clade2 in `echo $clade_list` # these are the clades where we will look for conserved seed matches
do for seed in `cat non_seeds_$directory`
   do withU=`echo $seed | tr T U`
       for species in `cat /mnt/data/home/herve.seitz/Blaise/miRNA_absence_verification_by_Herve/species_outside_$clade2`;do assembly_list=`grep -w $species /mnt/data/home/herve.seitz/Blaise/miRNA_absence_verification_by_Herve/Seed_selection/Full_correspondance_table | awk -F ',' '{print $2}'`;for assembly in `echo $assembly_list`;do if test `grep -c ','$assembly',' ../Species_in_100species_alignment.txt` -eq 1;then grep -v '^[ACGT][ACGT]* *$' /mnt/data/home/herve.seitz/Blaise/$output_Blaise/sites_db/sites_in_UTR_all/$seed'.txt' | grep -w $assembly;fi;done;done | sort | uniq -c > 'non_seed_'$directory'_families/Control_'$clade2'_'$seed'_matches.dat'
           echo "non_seed" $clade2 $seed `grep -v '^[ACGT][ACGT]* *$' /mnt/data/home/herve.seitz/Blaise/$output_Blaise/sites_db/sites_in_UTR_all/$seed'.txt' | wc -l` `cat 'non_seed_'$directory'_families/Control_'$clade2'_'$seed'_matches.dat' | wc -l` >> Control_$directory'.dat'
   done
done

echo "Clade Seed Number_of_seed_matches_in_UTR_alignment Number_of_seed_matches_in_UTR_alignment_outside_of_clade" > Parsed_$directory'.dat'
if test "$ref" = "hsa"
then grep '^Hominidae Hominidae ' Control_$directory'.dat' | awk '{print $2,$3,$4,$5}' >> Parsed_$directory'.dat'
     grep '^Catarrhini Catarrhini ' Control_$directory'.dat' | awk '{print $2,$3,$4,$5}' >> Parsed_$directory'.dat'
     grep '^Boreoeutheria Boreoeutheria ' Control_$directory'.dat' | awk '{print $2,$3,$4,$5}' >> Parsed_$directory'.dat'
     grep '^Euteleostomi Euteleostomi ' Control_$directory'.dat' | awk '{print $2,$3,$4,$5}' >> Parsed_$directory'.dat'
fi
if test "$ref" = "mmu"
then grep '^Murinae Murinae ' Control_$directory'.dat' | awk '{print $2,$3,$4,$5}' >> Parsed_$directory'.dat'
     grep '^Boreoeutheria Boreoeutheria ' Control_$directory'.dat' | awk '{print $2,$3,$4,$5}' >> Parsed_$directory'.dat'
     grep '^Euteleostomi Euteleostomi ' Control_$directory'.dat' | awk '{print $2,$3,$4,$5}' >> Parsed_$directory'.dat'
fi
