#!/usr/bin/perl

open(MIR_LIST,"miRNAs_for_predictions");
while (<MIR_LIST>)
{
    chomp;
    push(@miR_list,$_);
}
close(MIR_LIST);

open(MIR,"miR_Family_Info.txt");
while (<MIR>)
{
    chomp;
    @array=split('\t',$_);
    $name=$array[0];
### Below: to clean up inconsistencies in the TargetScan miRNA nomenclature:
    if ($name eq 'miR-140-3p.2/497b')
    {
	$name='miR-140-3p.2/497';
    }
    if ($name eq 'miR-200bc-3p/429')
    {
	$name='miR-200bc-3p/429-3p';
    }
    if ($name eq 'miR-291a-3p/294-3p/295-3p/302abd-3p')
    {
	$name='miR-291-3p/294-3p/295-3p/302-3p';
    }
### Above: to clean up inconsistencies in the TargetScan miRNA nomenclature:

    if ($array[2]==10090) #Selects murine miRNAs; problem: some miRNAs (e.g., miR-142-5p) are not in that list, while they appear in the predictions...
    {
	$seedm8{$name}=$array[1];
    }
    else
    {
	if (($seedm8{$name} eq '') && ($array[2]==10116)) #if no seed sequence is recorded for that miRNA: take at least the rat seed sequence...
	{
	    $seedm8{$name}=$array[1];
	}
	if ($seedm8{$name} eq '') #if still no seed sequence is recorded for that miRNA: take whatever seed sequence you find...
	{
		$seedm8{$name}=$array[1];
	}
    }
}
close(MIR);


### Below: for verification
foreach $miR (@miR_list)
{
    if ($seedm8{$miR} eq '')
    {
	print "$miR no recorded seed\n";
    }
}
### Above: for verification


open(OUT_EX,">tmp_PCT_excluding_overconserved.dat");
print OUT_EX "Gene_Ensembl Gene_name Transcript_Ensembl miRNA_m8seed maximal_PCT\n";
close(OUT_EX);
open(OUT_IN,">tmp_PCT_including_overconserved.dat");
print OUT_IN "Gene_Ensembl Gene_name Transcript_Ensembl miRNA_m8seed maximal_PCT\n";
close(OUT_IN);


for $chr (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,M,X,Y)
#$chr=1;
{
#    `cat "mmu_71_Conserved_Family_Info_overconservation_flagged_chr$chr"".txt" "mmu_71_Nonconserved_Family_Info_overconservation_flagged_chr$chr"".txt" > "mmu_71_Conserved_and_Nonconserved_Family_Info_overconservation_flagged_chr$chr"".txt"`; #(already done)
    $file="mmu_71_Conserved_and_Nonconserved_Family_Info_overconservation_flagged_chr$chr".".txt";
#    $file='mini.txt';
    {
	open(DATA,$file);
	while (<DATA>)
	{
	    chomp;
	    @array=split('\t',$_);

	    @uniqu=($seedm8{$array[0]});
	    if (push(@uniqu) != 1) # Because sometimes TargetScan's 'miR_Family_Info.txt' uses the concatenated names of all family members, sometimes it uses one member name
	    {
		@names=split('/',$array[0]);
		foreach $name (@names)
		{
		    if (($name!~/^miR-/) && ($name!~/^let-/))
		    {
			$name='miR-'.$name;
		    }
		    push(@seeds,$seedm8{$name});
		}
		%seen = ();
		@uniqu = grep { ! $seen{$_} ++ } @seeds;
	    }
	    
	    if (push(@uniqu) != 1)
	    {
		print "Problem with $array[0]! seeds=@uniqu\n";
	    }
	    
	    $seed=$uniqu[0];
	    
#	if ($array[10] ne 'NULL') #Exclude 'NULL' PCT values (which apparently mean that no multiple alignment was found around the seed match
#	{
#	    if ($array[11] ne 'yes')
#	    {
#		push(@{$PCT_list_excluding{$array[1]."_".$array[2]."_".$array[3]}{$seed}},$array[10]);
#	    }
#	    push(@{$PCT_list_including{$array[1]."_".$array[2]."_".$array[3]}{$seed}},$array[10]);
#	}
	    
	    if ($array[10] eq 'NULL') #Convert "NULL" PCT values into 0
	    {
		$array[10]=0;
	    }
	    if (($array[11] ne 'yes') && ($array[11] ne 'yes_unchecked'))
	    {
		push(@{$PCT_list_excluding{$array[1]."_".$array[2]."_".$array[3]}{$seed}},$array[10]);
	    }
	    push(@{$PCT_list_including{$array[1]."_".$array[2]."_".$array[3]}{$seed}},$array[10]);
	}
	close(DATA);
    }
}
#	$gene='ENSMUSG00000066595.4_Mfsd7b_ENSMUST00000085635.4';
#	$seed=$seedm8{'miR-338-3p'};
#	print "file=$file gene=$gene seed=$seed list_including=@{$PCT_list_including{$gene}{$seed}}\n";

open(OUT_EX,">>tmp_PCT_excluding_overconserved.dat");

foreach $gene (keys %PCT_list_excluding)
{
    for $seed (keys %{$PCT_list_excluding{$gene}})
    {
	$m=&max(@{$PCT_list_excluding{$gene}{$seed}});
	$gene_display=$gene;
	$gene_display=~s/_/ /g;
	print OUT_EX "$gene_display $seed $m\n";
    }
}
close(OUT_EX);
open(OUT_IN,">>tmp_PCT_including_overconserved.dat");
foreach $gene (keys %PCT_list_including)
{
#    if ($gene eq 'ENSMUSG00000066595.4_Mfsd7b_ENSMUST00000085635.4')
#    {
#	print "gene=$gene\n";
#	for $seed (keys %{$PCT_list_including{$gene}})
#	{
#	    print "seed=$seed PCT_list=@{$PCT_list_including{$gene}{$seed}}\n";
#	}
#	print "\n";
#    }
    for $seed (keys %{$PCT_list_including{$gene}})
    {
	$m=&max(@{$PCT_list_including{$gene}{$seed}});
	$gene_display=$gene;
	$gene_display=~s/_/ /g;
	print OUT_IN "$gene_display $seed $m\n";
    }
}
close(OUT_IN);


sub max
{
    $m=$_[0];
    foreach $elem (@_)
    {
	if ($elem > $m)
	{
	    $m=$elem;
	}
    }
    $m;
}
