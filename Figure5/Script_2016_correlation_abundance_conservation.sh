#!/bin/sh

if test "$1" = ""
then echo "Please enter biological tissue name, as given in 'Organs/Organ_sample_correspondence' (e.g., 12_month_striatum)."
     read tissue
else tissue=$1
fi

for strategy in including excluding #including or excluding over-conserved seed matches
do echo "library(plotrix)" > R_cmd
   nb_samples=0
   for sample in `grep '^'$tissue' ' /home/herve/Bhowali_theory/Abundance_vs_conservation/Organs/Organ_sample_correspondence | awk '{print $2}'`
   do echo $sample"=read.table('/home/herve/Bhowali_theory/Abundance_vs_conservation/Organs/Gene_"$sample".dat',header=T)" >> R_cmd
       echo $sample"_signal="$sample"\$Intensity[!is.na("$sample"\$Gene)]" >> R_cmd
       nb_samples=`echo $nb_samples"+1" | bc`
   done
   nb_samples_1=`echo $nb_samples"+1" | bc`
   for sample in `grep '^'$tissue' ' /home/herve/Bhowali_theory/Abundance_vs_conservation/Organs/Organ_sample_correspondence | awk '{print $2}' | head -1`
   do echo "gene_list="$sample"\$Gene[!is.na("$sample"\$Gene)]" >> R_cmd
   done
   echo "sample_array=array(c(gene_list[order(gene_list)]," >> R_cmd
   for sample in `grep '^'$tissue' ' /home/herve/Bhowali_theory/Abundance_vs_conservation/Organs/Organ_sample_correspondence | awk '{print $2}'`
   do echo $sample"_signal[order(gene_list)],"
   done |sed '$ s|,$|),dim=c(length(gene_list),'$nb_samples_1'))|' >> R_cmd
   echo "" >> R_cmd

#   for SEED in `cat mmu_matureJun14.fa  | grep -v '>' | cut -c 2-8 | sort | uniq` #Rather take seed list in TargetScan files: not every miRBase seed is in there.
   for SEED in `ls *_PCT_including_overconserved.dat | sed -e 's|_PCT_including_overconserved\.dat||'`
   do grep -B 1 '^[ACGU]'$SEED mmu_matureJun14.fa | grep -v '\-\-' | grep '^>' | sed -e 's|^> *mmu-||' -e 's| .*||' | sort -n > tmp
       if test `cat tmp | wc -l` -gt 1
       then label=`sed 's|^miR-||' tmp | sort -g | head -1 | sed -e 's|^|miR-|' -e 's|$| family|'`
           MIRNA=`sed 's|^miR-||' tmp | sort -g | head -1 | sed 's|^|miR-|'`
       else label=`cat tmp`
           MIRNA=$label
       fi
       MIRNA=`echo $MIRNA | sed 's|\*|star|'`

       echo "conservation=read.table('"$SEED"_PCT_"$strategy"_overconserved.dat',header=F)" >> R_cmd
       echo "" >> R_cmd
       echo "if (length(unique(conservation\$V6))>1)"  >> R_cmd #When PCT is 0 for every recorded site, don't waste time: cor.test won't work anyway
       echo "{" >> R_cmd
       echo "x=c()" >> R_cmd
       echo "xe=c()" >> R_cmd
       echo "y=c()" >> R_cmd
       echo "for (i in 1:length(conservation\$V1[!is.na(conservation\$V1)]))" >> R_cmd
       echo "{" >> R_cmd
       echo "index=c(1:length(sample_array[,1]))[sample_array[,1]==conservation\$V1[!is.na(conservation\$V1)][i]]" >> R_cmd
       echo "if (length(index) > 0)" >> R_cmd
       echo "{" >> R_cmd
       
       echo "x=append(x,mean(c(" >> R_cmd
       
       for sample in `seq 2 $nb_samples_1`
       do echo "sample_array[index,"$sample"],"
       done | sed '$ s|,$|)))|' >> R_cmd
       echo "xe=append(xe,std.error(c(" >> R_cmd
       for sample in `seq 2 $nb_samples_1`
       do echo "sample_array[index,"$sample"],"
       done | sed '$ s|,$|)))|' >> R_cmd
       echo "y=append(y,conservation\$V6[!is.na(conservation\$V1)][i])" >> R_cmd
       echo "}" >> R_cmd
       echo "}" >> R_cmd
       echo "} else x=c()"  >> R_cmd #When PCT is 0 for every recorded site, don't waste time: cor.test won't work anyway
       echo "if (length(x[!is.na(x)])>0)" >> R_cmd
       echo "{" >> R_cmd
#       echo "postscript('predicted_targets_for_"$MIRNA"_in_"$tissue"_"$strategy"_over-conserved_sites.ps',horizontal=F,paper='special',width=6,height=5)" >> R_cmd
#       echo "plotCI(x,y,xe,err='x',ylim=c(0,1),xlab='Gene expression in "$tissue"',ylab=expression(paste('Conservation score (',P[CT],') of binding sites to "$label"',sep='')),main=paste(\"Kendall's tau=\",signif(cor.test(x,y,method='kendall')\$estimate,digits=3),\"\\\n(p-value=\",signif(cor.test(x,y,method='kendall')\$p.value,digits=3),\")\"))" | sed 's|_| |g' >> R_cmd
#       echo "dev.off()" >> R_cmd
       echo "sink('Correlation_"$MIRNA"_in_"$tissue"_"$strategy"_over-conserved_sites')" >> R_cmd
       echo "print(c('"$SEED"','"$label"'))" >> R_cmd
       echo "if (is.na(cor.test(x,y,method='kendall')\$p.value))" >> R_cmd
       echo "{" >> R_cmd
       echo "p='NA'" >> R_cmd
       echo "} else" >> R_cmd
       echo "if (cor.test(x,y,method='kendall')\$p.value < 2.2e-16)" >> R_cmd
       echo "{" >> R_cmd
       echo "p=2.2e-16" >> R_cmd
       echo "} else" >> R_cmd
       echo "{" >> R_cmd
       echo "p=signif(cor.test(x,y,method='kendall')\$p.value,digits=3)" >> R_cmd
       echo "}" >> R_cmd
       echo "print (c(signif(cor.test(x,y,method='kendall')\$estimate,digits=3),p))" >> R_cmd
       echo "sink()" >> R_cmd
       echo "} else" >> R_cmd
       echo "{" >> R_cmd
       echo "sink('Correlation_"$MIRNA"_in_"$tissue"_"$strategy"_over-conserved_sites')" >> R_cmd
       echo "print(c('"$SEED"','"$label"'))" >> R_cmd
       echo "print (c('No data available'))" >> R_cmd
       echo "sink()" >> R_cmd
       echo "}" >> R_cmd
   done
   
   R CMD BATCH R_cmd
   mv R_cmd.Rout save_R_cmd_Rout_$tissue'_'$strategy
   for RNA in `ls Correlation_*_in_$tissue"_"$strategy"_over-conserved_sites" | sed -e 's|^Correlation_||' -e 's|_in_'$tissue'_.*||'`
   do grep -A 1 -w tau Correlation_$RNA'_in_'$tissue'_'$strategy'_over-conserved_sites' | tail -1 | sed 's|^|'$RNA' |'
   done | sed '1 s|.*|RNA tau p_value\
&|' > correlation_summary_in_$tissue'_'$strategy'_over-conserved_sites.dat'
   
   echo "for_volcano=read.table('correlation_summary_in_"$tissue"_"$strategy"_over-conserved_sites.dat',header=T)" > R_volcano_plotting
   echo "postscript('Scaled_correlation_coefficient_distribution_in_"$tissue"_"$strategy"_over-conserved_sites.ps',paper='special',horizontal=F,width=5,height=5)" >> R_volcano_plotting
   echo "qval=p.adjust(for_volcano\$p_value,method='BH')" >> R_volcano_plotting
   echo "plot(for_volcano\$tau,qval,xlim=c(-.25,.25),ylim=c(1e-16,1),log='y',xlab=expression(paste(\"Kendall's \",tau)),ylab='Adjusted p-value',axes=F,">> R_volcano_plotting
   echo "main='"$tissue"')" | sed 's|_| |g' >> R_volcano_plotting
   echo "axis(1,labels=c(-0.25,-0.1,0,0.1,0.25),at=c(-0.25,-0.1,0,0.1,0.25))" >> R_volcano_plotting
   echo "axis(2,labels=c(1e-16,1e-12,1e-8,1e-4,1),at=c(1e-16,1e-12,1e-8,1e-4,1))" >> R_volcano_plotting
   echo "lines(c(-0.25,0.25),c(0.05,0.05),col='red',lty=3,lwd=2)" >> R_volcano_plotting
   echo "lines(c(0,0),c(1e-16,1),lty=3,lwd=2)" >> R_volcano_plotting
   echo "dev.off()" >> R_volcano_plotting
   R CMD BATCH R_volcano_plotting
   echo "Done. Output file is named: 'Correlation_coefficient_distribution_in_"$tissue".ps'."
   
done
