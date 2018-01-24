#!/usr/bin/perl
use Math::Round;

$TIME_STAMP=`date +%s`;
chomp $TIME_STAMP;
$TIME_STAMP=$TIME_STAMP."_".$ARGV[0];

print "TIME_STAMP=$TIME_STAMP\n";

$LIMIT_TO_1HIT=1; #set to 1 if you just want to know whether there's at least one hit in the scanned genomes (and you don't need the exhaustive list of every hit): faster

$MACHINE='Daisy';
#$MACHINE='garonne';

#$DISTANCE_CUTOFF=2; # This is the distance threshold (calculated by muscle on hairpin alignments) to decide whether two hairpins are similar enough to be used to generate a common HMMer pattern (not used finally: muscle is too irreproducible)
$EVAL_CUTOFF=0.001; # This is the E-value threshold (calculated by blastn on hairpin alignments) to decide whether two hairpins are similar enough to be used to generate a common HMMer pattern
$HMMER_CUTOFF=1; # This is the cutoff on HMMer E-value for the identification of orthologs in the scanned genome
$BLAST_CUTOFF=1; # This is the cutoff on Blast E-value for the identification of orthologs in the scanned genome
$UNBRANCHED_CUTOFF=0.75; # This is the cutoff on predicted DeltaG for the most stable unbranched hairpin, relatively to the most stable predicted structure for the candidate (e.g., a value of 0.75 means that the most stable unbranched hairpin has to have a DeltaG < 0.75* DeltaG of the most stable predicted structure); only the structure falling within 5 kcal/mol of the most stable structure will be analyzed
$MINIMAL_LENGTH=50; # This is the minimal accepted length for pre-miRNA orthologs in the scanned genome

$MATURE_SEQ='matureJun14.fa';
$HAIRPIN_SEQ='Fused_hairpinJun14.fa';
$ANNOT='miRNAJun14.dat';
$ORGANISMS='vertebrate_list';
#$ORGANISM_CORRESPONDENCE='Correspondance_table';

print "\nThis is 'Module_extracts_hairpin_from_seed.pl'. I will scan genomes for candidate orthologs for known miRNAs from specific vertebrate seed families. Criteria:\n1. HMMer search with an E-value cutoff of $HMMER_CUTOFF, Blast search with an E-value cutoff of $BLAST_CUTOFF.\n2. Candidate pre-miRNA locus must fold in a stable unbranched hairpin, whose predicted DeltaG is at most $UNBRANCHED_CUTOFF times that of the predicted DeltaG of its most stable structure.\n3. The miRNA seed has to be perfectly conserved.\n\n";

if ($ARGV[2] eq '')
{
    print "Please enter inputs: seed sequence, then text file containing species code list then 3-letter code for reference species (e.g., ./Module_extracts_hairpin_from_seed.pl CCAGCAT species_list hsa). I will search for orthologous pre-miRNAs in these species using HMMer.\n";
}
else
{
    $CODE_3LETTERS=$ARGV[2];
    open(ANNOT,$ANNOT);
    while(<ANNOT>)
    {
	chomp;
	if (/^AC /)
	{
	    s/^AC *//;
	    s/;//;
	    $hairpin=$_;
	}
	if (/^FT  *miRNA /)
	{
	    s/.* //;
	    $position=$_;
	}
	if (/^FT  *\/accession=/)
	{
	    s/^FT  *\/accession="//;
	    s/"$//;
	    $origin{$_}=$hairpin;
	    $position_in_hairpin{$_}=$position;
#	    if ($_ eq 'MIMAT0000416')
#	    {
#		print "position_in_hairpin{MIMAT0000416}=$position\n";
#	    }
	}
    }
    close(ANNOT);

    open(HAIRPIN,$HAIRPIN_SEQ);
    while(<HAIRPIN>)
    {
	chomp;
	if (/^>/)
	{
	    s/^>[a-zA-Z\d-]* //;
	    s/ .*//;
	    $ID=$_;
	}
	else
	{
	    $hairpin_sequence{$ID}=$_;
	}
    }
    close(HAIRPIN);

    open(MATURE,$MATURE_SEQ);
    while(<MATURE>)
    {
	chomp;
	if (/^>/)
	{
	    ($name,$ID)=($_,$_);
	    $name=~s/^> *//;
	    $name=~s/ .*//;
	    $ID=~s/^>[a-zA-Z\d-]* //;
	    $ID=~s/ .*//;
	    $miRNA_name{$ID}=$name;
	    $miRNA_ID{$name}=$ID;
	}
	else
	{
	    $mature_sequence{$ID}=$_;
	}
	
    }
    close(MATURE);

    open(SPECIES,$ORGANISMS);
    while(<SPECIES>)
    {
	chomp;
	push(@species_list,$_);
    }
    close(SPECIES);

    open(SPECIES_SCAN,$ARGV[1]);
    while(<SPECIES_SCAN>)
    {
	chomp;
	push(@species_to_scan,$_);
    }
    close(SPECIES_SCAN);



	$seed=$ARGV[0];
	$seed=~tr/acgtT/ACGUU/;
        @clusters=();
	open(HITS,">Hits_for_$seed"."_in_additional_genomes.txt");
	@for_alignment=();
	open(FOR_ALIGNMENT,">tmp_for_alignment_$TIME_STAMP".".fa");
	for $mature (keys %mature_sequence)
	{
	    $species=$miRNA_name{$mature};
	    $species=~s/-.*//;
	    if (($mature_sequence{$mature}=~/^.$seed/) && (grep {$_ eq $species} @species_list))
	    {
#		print "$seed $miRNA_name{$mature} $mature $origin{$mature} $hairpin_sequence{$origin{$mature}}\n";
		print FOR_ALIGNMENT "> pre-$miRNA_name{$mature}\n$hairpin_sequence{$origin{$mature}}\n";
		push(@for_alignment,'pre-'.$miRNA_name{$mature});
	    }
	}
	close(FOR_ALIGNMENT);
#	`muscle -in tmp_for_alignment.fa -tree2 tmp.out;python parser.py > tmp_parsed;./parser2.sh > tmp_parsed2`; NOT REPRODUCIBLE!!!
	if ($MACHINE eq 'garonne')
	{
	    `formatdb -p F -i tmp_for_alignment_$TIME_STAMP.".fa";blast2 -p blastn -d tmp_for_alignment_$TIME_STAMP.".fa" -e $EVAL_CUTOFF -W 4 -i tmp_for_alignment_$TIME_STAMP.".fa" -B 1 -o tmp_aligned_$TIME_STAMP`; # This is for older blast installations (e.g., on garonne)
	}
	if ($MACHINE eq 'Daisy')
	{
	    `makeblastdb -dbtype=nucl -in=tmp_for_alignment_$TIME_STAMP".fa";blastn -db=tmp_for_alignment_$TIME_STAMP".fa" -evalue=$EVAL_CUTOFF -word_size=4 -query=tmp_for_alignment_$TIME_STAMP".fa" -outfmt=6 -out=tmp_aligned_$TIME_STAMP`; # This is for newer blast installations (e.g., on Daisy)
	}
	@selected=();
        %partner=();
	open(ALIGNED,"tmp_aligned_$TIME_STAMP");
	while(<ALIGNED>)
	{
	    chomp;
	    if ($_ !~ /^#/)
	    {
		@array=split('\t',$_);
		if ($array[0] ne $array[1])
		{
		    @to_sort=($array[0],$array[1]);
		    @sorted=sort(@to_sort);
		}

### Below: merging partner lists for $sorted[0] and $sorted[1]
		$receiving=$sorted[0];
		for $parent (keys %partner)
		{
		    if (grep {$_ eq $sorted[0]} @{$partner{$parent}})
		    {
			$receiving=$parent;
		    }
		}
		push(@{$partner{$receiving}},$sorted[1]);
		push(@{$partner{$receiving}},@{$partner{$sorted[1]}});
		delete $partner{$sorted[1]};
### Above: merging partner lists for $sorted[0] and $sorted[1]
	    }
	}
	close(ALIGNED);
### Below: removing duplicates in partner lists
	foreach $premiRNA (keys %partner)
	{
	    if ($premiRNA)
	    {
		@partner_list=@{$partner{$premiRNA}};
		push(@partner_list,$premiRNA);
		%seen = ();
		@uniqu = grep { ! $seen{$_} ++ } @partner_list;
		push(@clusters,[ @uniqu ]);
	    }
	}
### Above: removing duplicates in partner lists

### Below: adding the pre-miRNAs that were not even in 'tmp_aligned' because they didn't have any Blast hit with a correct E-value
	foreach $premiRNA (@for_alignment)
	{
	    $new=1;
	    for $aref ( @clusters )
	    {
		if (grep {$_ eq $premiRNA} @$aref)
		{
		    $new=0;
		}
	    }
	    if ($new)
	    {
		@add=($premiRNA);
		push(@clusters,[ @add ]);
	    }
	}

	for $i ( 0 .. $#clusters )
	{
	    @human_ID=();
	    $file_name="pre-miRNAs_".$i."_for_".$seed.".fa";
	    open(OUT,">$file_name");
	    for $j ( 0 .. $#{$clusters[$i]} )
	    {
		$name=$clusters[$i][$j];
		$name=~s/^pre-//;
		$this_species=$name;
		$this_species=~s/-.*//;
		if ($this_species eq $CODE_3LETTERS)
		{
		    push(@human_ID,$miRNA_ID{$name}); # Record the ID's of human orthologs for that miRNA (or: mouse orthologs, if the analysis is mouse-centered)
		}
		print OUT "> $clusters[$i][$j]\n$hairpin_sequence{$origin{$miRNA_ID{$name}}}\n";
	    }
	    close(OUT);

	    if (($LIMIT_TO_1HIT==0) || ($found_hits{$seed}==0))
	    {
		$cluster_count=$#{$clusters[$i]}+1;
		if ($cluster_count==1)
		{
		    print "   cluster #$i has only 1 member: can't build an HMMer profile with that. I will use Blast instead...\n";
		    
		    if ($MACHINE eq 'Daisy') #Anyway: I don't have enough disk space on Garonne to put all the genome blast databases...
		    {
			foreach $scanned_species (@species_to_scan)
			{
			    $genome_file='/mnt/data/home/herve.seitz/Genomes/Blast_indexes/'.$scanned_species.'.fa';
			    $output_file='blast_output_for_pre-miRNAs_'.$i.'_for_'.$seed.'_on_'.$scanned_species;
			    `blastn -db=$genome_file -evalue=$BLAST_CUTOFF -word_size=6 -query=$file_name -out=$output_file`;
			    open(CANDIDATES,$output_file);
			    $cand_sequence='';
			    while(<CANDIDATES>)
			    {
				chomp;
				if (/^>/)
				{
				    s/^> *//;
				    $cand_chr=$_;
				}
				if (/^ Score =/)
				{
				    $cand_start=0;
				    $cand_sequence='';
				}
				if (/^Sbjct /)
				{
				    @array=split(' ',$_);
				    if ($cand_start==0)
				    {
					$cand_start=$array[1];
				    }
				    $cand_sequence=$cand_sequence.$array[2];
				    $cand_end=$array[3];
				}
				if (/^$/)
				{
				    ++$empty_lines;
				    if (($empty_lines==2) && ($cand_sequence))
				    {
					$cand_sequence=~s/-//g;
					if (($LIMIT_TO_1HIT==0) || ($found_hits{$seed}==0))
					{
					    &evaluates_candidate($cand_sequence,$cand_start,$cand_end);
					}
				    }
				}
				else
				{
				    $empty_lines=0;
				}
			    }
			    close(CANDIDATES);
			    if (($LIMIT_TO_1HIT==1) && ($found_hits{$seed}>0))
			    {
				print "(stopping search after finding a potential ortholog)\n";
				last;
			    }
			}
		    }
		}
		if ($#{$clusters[$i]}>0) # do the following only if there are at least two hairpin sequences in the sequence cluster (so an HMMer profile can be built)
		{
		    print "   now scanning genomes with HMMer for orthologs to pre-miRNA sequence cluster $i...\n";
		    $file_name="pre-miRNAs_".$i."_for_".$seed.".fa";
		    $file_name1="pre-miRNAs_".$i."_for_".$seed.".aln";
		    $file_name2="pre-miRNAs_".$i."_for_".$seed.".hmm";
		    
		    if ($MACHINE eq 'garonne')
		    {
			`clustalw $file_name;hmmbuild $file_name2 $file_name1`;
		    }
		    if ($MACHINE eq 'Daisy')
		    {
			$output_clustalw=system("clustalw2 $file_name > /dev/null");
                        if ($output_clustalw)
                        {
                            print "Warning! clustalw problem for ARGV[0]=$ARGV[0] ARGV[1]=$ARGV[1] ARGV[2]=$ARGV[2]\nError message:\n$output_clustalw\n";
                        }
                        `hmmbuild $file_name2 $file_name1`;
		    }
		    foreach $scanned_species (@species_to_scan)
		    {
			print "      now scanning genome $scanned_species...\n";
			$genome_file='/mnt/data/home/herve.seitz/Genomes/'.$scanned_species.'.fa';
			$output_file='hmmer_output_for_pre-miRNAs_'.$i.'_for_'.$seed.'_on_'.$scanned_species;
			`nhmmer -o $output_file $file_name2 $genome_file`;
			open(CANDIDATES,$output_file);
			$read=0;
			$header_read=0;
			$count=0;
			%hit_ends=();
			while(<CANDIDATES>)
			{
			    chomp;
			    if ($read==0)
			    {
				if (/^ +------- ------ ----- +-------- +----- +----- +-----------$/)
				{
				    $header_read=1;
				}
				else
				{
				    if ($_ eq '')
				    {
					$header_read=0;
				    }
				    if ($header_read)
				    {
					@array=split(' ',$_);
					push(@{$hit_ends{$array[3]}},$array[5]);
				    }
				}
			    }
			    if ($_ eq 'Annotation for each hit  (and alignments):')
			    {
				$read=1;
			    }
			    else
			    {
				if ($read)
				{
				    ++$count;
				    if (/^>>/)
				    {
					$count=0;
					$cand_sequence='';
					$cand_start='';
				    }
				    if ($count==3)
				    {
					@array=split(' ',$_);
					$evalue=$array[3];
				    }
				    if (($count>=9) && (($count-9)/5==round(($count-9)/5)) && ($evalue<=$HMMER_CUTOFF)) # If this is the 9th, or 14th, or 19th, ... line of the report for that hit, and the HMMer E-value is below the cutoff
				    {
					@array=split(' ',$_);
					($cand_chr,$add_cand_sequence,$cand_end)=($array[0],$array[2],$array[3]);
					if ($cand_start eq '')
					{
					    $cand_start=$array[1];
					}
					$cand_sequence=$cand_sequence.$add_cand_sequence;
					if (grep {$_ eq $cand_end} @{$hit_ends{$cand_chr}}) # Tests if the report for that hit has been completely read (they can span over several blocks of 4 lines)
					{
					    if ((($LIMIT_TO_1HIT==0) || ($found_hits{$seed}==0)) && (length($cand_sequence)>=$MINIMAL_LENGTH))
					    {
						&evaluates_candidate($cand_sequence,$cand_start,$cand_end);
					    }
					}
				    }
				}
			    }
			}
			close(CANDIDATES);
                        if (($LIMIT_TO_1HIT==1) && ($found_hits{$seed}>0))
                        {
                            print "(stopping search after finding a potential ortholog)\n";
                            last;
                        }
		    }
                }
	    }
	}
	close(HITS);
}



sub evaluates_candidate
{
    ($cand_sequence,$cand_start,$cand_end)=($_[0],$_[1],$_[2]);
#    print "Now folding $cand_sequence (for seed $ARGV[0])...\n";
    $most_stable=`echo $cand_sequence | RNAfold | tail -1`;
    $most_stable=~s/.* \(//;
    $most_stable=~s/\)\n//;
    $unbranched=`echo $cand_sequence | RNAsubopt -s -e 5 | tail -n +2 | grep '^[\.(]*[\.)]* ' | sed 's|.* ||' | head -1`;
    $unbranched=~s/\n//;
#    print "$cand_chr: cand_sequence=$cand_sequence unbranched=$unbranched UNBRANCHED_CUTOFF=$UNBRANCHED_CUTOFF most_stable=$most_stable";
    if ($unbranched <= $UNBRANCHED_CUTOFF*$most_stable)
    {
	$OK=1;
    }
    else
    {
	$OK=0;
    }
#    print " $cand_chr cand_sequence=$cand_sequence OK=$OK\n";
    if ($OK)
    {
### Now verifies that the seed is perfectly conserved in the candidate
	`cp $file_name tmp_seed_check_$TIME_STAMP;echo "> candidate" >> tmp_seed_check_$TIME_STAMP;echo $cand_sequence >> tmp_seed_check_$TIME_STAMP;clustalw2 tmp_seed_check_$TIME_STAMP`;
#	print "Position of human miRNA $human_ID[0] in its hairpin: $position_in_hairpin{$human_ID[0]} (human_ID=@human_ID; miRNA_name=$miRNA_name{$human_ID[0]})\n";
	($seed_start,$seed_end)=($position_in_hairpin{$human_ID[0]},$position_in_hairpin{$human_ID[0]});
	$seed_start=~s/\.\..*//;
	$seed_start=$seed_start+1;
	$seed_end=$seed_start+6;
	$block=0;
	@human_sequence=();
	@candidate_sequence=();
	open(ALIGN_SEED,"tmp_seed_check_".$TIME_STAMP.".aln");
	while (<ALIGN_SEED>)
	{
	    chomp;
	    if (/^pre-$miRNA_name{$human_ID[0]} /) #for the first human ortholog in @human_ID, find the coordinates of its seed in the alignment
	    {
		s/.* //;
		@seq_array=split('',$_);
		push(@human_sequence,@seq_array);
	    }
	    if (/^candidate /)
	    {
		s/.* //;
		@seq_array=split('',$_);
		push(@candidate_sequence,@seq_array);
	    }
	}
	close(ALIGN_SEED);
#	print "human_sequence=@human_sequence\n";
	$human_length=push(@human_sequence);
	$n_align=0;
	for ($n=0;$n<$human_length;++$n)
	{
	    if ($human_sequence[$n] ne '-')
	    {
		++$n_align;
		if ($n_align==$seed_start)
		{
		    $start_in_align=$n;
		}
		if ($n_align==$seed_end)
		{
		    $end_in_align=$n;
		}
	    }
	}
#	print "seed_start=$seed_start start_in_align=$start_in_align seed_end=$seed_end end_in_align=$end_in_align\n";
	$candidate_seed='';
	for ($n=$start_in_align;$n<=$end_in_align;++$n)
	{
	    $candidate_seed=$candidate_seed.$candidate_sequence[$n];
	}
	$candidate_seed=~tr/T/U/;
        $candidate_seed=~s/-//g;
#	print "candidate_sequence=@candidate_sequence\n";
#	print "candidate_seed=$candidate_seed miRNA_seed=$seed\n";
	if ($candidate_seed eq $seed)
	{
	    $seq_display=$cand_sequence;;
	    $seq_display=~s/-//g;
	    print HITS "reference: $ref_sequence (pre-$name, mature sequence: $mature_sequence{$miRNA_ID{$name}})\ncandidate: $seq_display ($scanned_species, chromosome $cand_chr $cand_start-$cand_end)\n\n";
	    print "Found a hit!\nreference: $ref_sequence ($name, pre-miRNAs_".$i."_for_".$seed.", $mature_sequence{$miRNA_ID{$name}})\ncandidate: $seq_display ($scanned_species chromosome $cand_chr $cand_start-$cand_end). Most stable structure DeltaG: $most_stable; most stable unbranched hairpin DeltaG: $unbranched\n\n";
	    ++$found_hits{$seed};
	}
    }
}
