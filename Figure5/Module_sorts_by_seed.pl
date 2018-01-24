#!/usr/bin/perl

if ($ARGV[0] eq '')
{
    print "Please enter script input file name (e.g., ./Module_sorts_by_seed.pl PCT_excluding_overconserved.dat).\n";
}
else
{
    open(IN,$ARGV[0]);
    while(<IN>)
    {
	chomp;
	if ($_ ne 'GenBank_ID Gene_Ensembl Gene_name Transcript_Ensembl miRNA_m8seed maximal_PCT')
	{
	    @array=split(' ',$_);
	    open(OUT,">>$array[4]"."_$ARGV[0]");
	    print OUT "$_\n";
	    close(OUT);
	}
    }
    close(IN);
}
