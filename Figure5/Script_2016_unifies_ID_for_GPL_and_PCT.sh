#!/bin/sh

tail -n +2 tmp_PCT_including_overconserved.dat | awk '{print $1}' | sed 's|\..*||' | sort | uniq > Ensembl_IDs_in_PCT


echo "Ensembl_gene_ID gene_ID_in_Genbank" > Ensembl_gene_annotation.dat
for gene in `cat Ensembl_IDs_in_PCT`
do wget http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene\&term=Ensembl:$gene
   ID=`grep -P '^\t*<Id>' esearch.fcgi?db=gene\&term=Ensembl:$gene | sed 's|^\t*<Id>\([0-9]*\)</Id>|\1|'`
   echo $gene $ID >> Ensembl_gene_annotation.dat
   sleep 0.5
done

