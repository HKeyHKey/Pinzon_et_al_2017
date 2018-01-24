#!/usr/bin/perl

open(CORRESP,"Ensembl_gene_annotation.dat");
while(<CORRESP>)
{
    chomp;
    if ($_ ne 'Ensembl_gene_ID gene_ID_in_Genbank')
    {
	@array=split(' ',$_);
	$genbank{$array[0]}=$array[1];
    }
}
close(CORRESP);

foreach $file ('tmp_PCT_excluding_overconserved.dat','tmp_PCT_including_overconserved.dat')
{
    $outfile=$file;
    $outfile=~s/^tmp_//;
    open(IN,$file);
    open(OUT,">$outfile");
    while(<IN>)
    {
	chomp;
	if ($_ eq 'Gene_Ensembl Gene_name Transcript_Ensembl miRNA_m8seed maximal_PCT')
	{
	    print OUT "GenBank_ID $_\n";
	}
	else
	{
	    $memo=$_;
	    s/ .*//;
	    s/\..*//;
	    if ($genbank{$_})
	    {
		print OUT "$genbank{$_} $memo\n";
	    }
	    else #try to get a GenBank ID number using the HUGO name...
	    {
		@array=split(' ',$memo);
		$hugo=$array[1];
		if ($add{$hugo} eq '')
		{
		    $id=`./Script_2016_unifies_with_HUGO.sh $hugo`;
		    chomp($id);
		    $add{$hugo}=$id;
		}
		if ($add{$hugo})
		{
		    print OUT "$add{$hugo} $memo\n";
		}
		else
		{
		    print OUT "NA $memo\n";
		}
	    }
	}
    }
    close(IN);
    close(OUT);
}
