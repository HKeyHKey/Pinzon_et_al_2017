#!/bin/sh

nb_Huang=`cat edited_journal.pgen.1001154.s001.txt | wc -l`
awk '{print $4}' edited_journal.pgen.1001154.s001.txt | sed 's?\(.*\)|\(.*\)|.*?\1 \2?' > tmp_Huang

#Below: for the 3442 genes in the Huang et al., 2010 list, that are not in the TargetScan list (some of them may actually be in the TargetScan list, but under another name): gene nomenclature synchronization):
awk '{print $1}' tmp_Huang | sort | uniq > list_Huang2010
cat /home/herve/Blaise/TargetScan_prediction_for_absent_miRNAs/Pan-vertebrate_alignment/Parsing_Blaise_results/For_elimination_of_overconserved_sites_in_other_analyses/vert_70_Conserved_Family_Info_overconservation_flagged_chr*.txt > TargetScan_predictions.txt
awk '{print $3}' TargetScan_predictions.txt | sort | uniq > list_TargetScan
cat list_Huang2010 list_TargetScan list_TargetScan | sort | uniq -c | grep '^ *1 ' | awk '{print $2}' > list_for_synchronization

echo "Gene_name GenBank_ID Alternative_name" > Gene_name_synonyms.txt
for gene in `cat list_for_synchronization`
do wget -O tmp_download http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene\&term=Homo+sapiens[orgn]$gene''[gene]live
   sleep 1
   for ID in `grep '^<Id>' tmp_download | sed -e 's|^<Id>||' -e 's|</Id>||'`
   do wget -O tmp_download2 http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene\&id=$ID
      sleep 0.5
      sed '/^ *gene /,/^ *prot / !d' tmp_download2 | grep '^ *locus "' | sed 's|.*"\(.*\)".*|\1|' > tmp_alternative_names
      sed '/^ *gene /,/^ *prot / !d' tmp_download2 | sed '/^ *syn {/,/^ *}/ !d' | grep -v '{' | grep -v '}' | sed 's|.*"\(.*\)".*|\1|' >> tmp_alternative_names
      sed 's|^|'$gene' '$ID' |' tmp_alternative_names >> Gene_name_synonyms.txt
   done
done



echo "Gene miRNA PCT haploinsufficiency_probability Overconservation_flag" > Comparison_Huang_Friedman.dat
for gene in `seq 1 $nb_Huang`
do line=`head -$gene tmp_Huang | tail -1`
   proba=`echo $line | awk '{print $2}'`
   name=`echo $line | awk '{print $1}'`
   grep -P '^\S+\t\S+\t'$name'\t' TargetScan_predictions.txt > tmp_grep
   if test `cat tmp_grep | wc -l` -ne 0
   then sed -e 's|\tNULL\t|\t0\t|' -e 's|^\([a-zA-Z0-9/_\.-]*\)\t.*\t\([0-9\.]*\)\t\([a-z_]*\)\t[a-z0-9XY:_+-]*$|'$name' \1 \2 '$proba' \3|' tmp_grep
   else for alternative_name in `grep '^'$name' ' Gene_name_synonyms.txt | awk '{print $3}'`
        do awk '$3=="'$alternative_name'" {print}' TargetScan_predictions.txt | sed -e 's|\tNULL\t|\t0\t|' -e 's|^\([a-zA-Z0-9/_\.-]*\)\t.*\t\([0-9\.]*\)\t\([a-z_]*\)\t[a-z0-9XY:_+-]*$|'$name' \1 \2 '$proba' \3|'
        done
   fi
done >> Comparison_Huang_Friedman.dat

awk '$5!="yes" {print}' Comparison_Huang_Friedman.dat > Comparison_Huang_Friedman_without_overconserved.dat
mv Comparison_Huang_Friedman.dat Comparison_Huang_Friedman_with_overconserved.dat

Rscript R_commands_TargetScan_vs_Huang2010_with_or_without_overconserved without_overconserved
Rscript R_commands_TargetScan_vs_Huang2010_with_or_without_overconserved with_overconserved

echo "Done. See 'Histogram_PCTmax_Conserved.ps', 'Histogram_Haploinsufficiency_Conserved.ps', 'Huang_vs_Friedman_Conserved.ps', 'Binned_Huang_vs_Friedman_Conserved.ps' and 'Correlation_test_result'."
