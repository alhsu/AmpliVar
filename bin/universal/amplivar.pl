#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use File::Copy;
use File::Basename;

=head1 NAME

amplivar.pl

=head1 VERSION

0.1

=head1 DESCRIPTION

Process amplicon read pairs from fastq.gz to genotypes, variant calls, quality and coverage

=head1 INPUT

FLAGS:
    B<-i> input sequence directory
    B<-o> output directory
    B<-j> usual suspects file
    B<-k> amplicon flanks
    B<-h> (optional) help


=head1 OUTPUT

uncompressed fastq files (.fastq)
filtered fasta files (.fna)
quality files (.csv)
matches against usual suspects (.txt)
grouped reads (_grp)

=head1 OTHER REQUIREMENTS

SeqPrep in path, bash shell

=head1 EXMAPLE USAGE

amplivar.pl -i source_directory -o output_directory -j usual_suspects  -k amplicon_flanks

=head1 AUTHOR

Graham Taylor University of Melbourne

=cut

# START OF CODE

use vars qw( $opt_i $opt_o $opt_j $opt_k $opt_h );


my $VERSION = 0.1;

our $source_dir;		#this is where the fastq.gz files are
our $output_dir;
our $usual_suspects;	#this is the list of known mutations: 4 columns, 3 tabs
our $amplicon_flanks;	#this is the flanking sequences to capture $1 (.*) 4 columns, 3 tabs
our $grouped_dir;		#grouped reads from fasta files
our $genotype_dir;		#genotypes
our $flanked_dir;		#flanked read counts
our $locus_dir ;		#sorted by locus
our $sorted_locus_dir ;	#each amplicon sorted by abundance
our $fna_dir ;			#fasta files from merged fastq
our $quals_dir ;		#qual score files from merged fastq


our $fastq ;
our $fasta ;
our $qual ;
our $grouped ;
our $flanked ; 

my $correct_usage = "USAGE: amplivar.pl [ flags ]\n\n" .
                    "  FLAGS\n" .
                    "    Required:\n" .
                    "		-i <source directory>\n" .
                    "		-o <output directory>\n" .
                    " 		-j <path to usual suspects>\n" .
                    "		-k <path to amplicon flanks>\n" .
                    "    Optional:\n" .
                    "      -h\n";

my $help =
"This program will read in all fastq.gz format sequence files\n".
"in the source directory, a file containing a list of usual suspect mutations\n".
"and file of amplicon flanking primers. It will output all matches to the usual suspect\n".
"file, matches to the flanking sequences and it will group reads within each amplicon\n". 
"for further analysis\nRequires bash shell, requires SeqPrep in path" ;

# Get command line args
getopts("i:o:j:k:h");

# Print out help
if ( $opt_h )
{
  print "\n$help\n\n$correct_usage\n\n\n";
  exit;
}

# Check that there is a source directory
if ( ! $opt_i )
{
  die "\nNo source directory given.\n\n$correct_usage\n\n\n";
}
else
{
  # Set the input directory
 $source_dir = $opt_i ;
  if (-d $source_dir)
  {
    print "getting all fastq.gz files from $source_dir\n";
  }
  else
  {
    print "\n$source_dir not found\n\n$correct_usage\n\n\n";
    die;
  }    
}

# Check that there is an output directory
if ( ! $opt_o )
{
  die "\nNo output directory given.\n\n$correct_usage\n\n\n";
}
else
{
  # Set the output directory
 $output_dir = $opt_o ;
  if (-d $output_dir)
  {
    print "writing all outputs to $output_dir\n";
  }
  else
  {
    print "\n$output_dir not found\n\n$correct_usage\n\n\n";
    die;
  }    
}

# Read usual suspects file
if ( ! $opt_j )
{
  die "\nNo usual suspects file given\n\n$correct_usage\n\n\n";
}
else
{
$usual_suspects = $opt_j ;
}
if (-e $usual_suspects)
{
print "\n$usual_suspects file will be used for genotyping\n" ;
}
else
  {
    print "\n$usual_suspects file not found\n\n$correct_usage\n\n\n";
    die;
  }    

# Read amplicon flanks file
if ( ! $opt_k )
{
  die "\nNo amplicon flanks file name given.\n\n$correct_usage\n\n\n";
}
else
{
$amplicon_flanks = $opt_k ;
}
if (-e $amplicon_flanks)
{
print "\n$amplicon_flanks file will be used for mutation scanning\n" ;
}
else
  {
    print "\n$amplicon_flanks file not found\n\n$correct_usage\n\n\n";
    die;
  }

#log timestamp
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
our $date = sprintf("%4d%02d%02d%02d%02d%02d",($year + 1900),($mon+1),$mday,$hour,$min,$sec);
our $log ="$output_dir/amplivar_$date.log";
our $LOG ;
open ($LOG, ">", "$log") ;


print $LOG "source directory\t$source_dir\noutput directory\t$output_dir\ngenotype list\t$usual_suspects\namplicon flanks\t$amplicon_flanks\n";
fastq2fasta($source_dir) ;
group_fasta_reads();
genotype_grouped_reads($grouped_dir) ;
grab_by_flanks($grouped_dir) ;
sort_by_locus($flanked_dir) ;
sort_each_locus() ;
tidy_up() ;
exit ;

# Some output to say what has gone on
# print_summary();

################################################################################
#
# Subroutines below
#
################################################################################

=item *

fastq2fasta

    Title   : fastq2fasta
    Usage   : fastq2fasta
    Function: writes filtered fasta file and csv quality file
    Returns : nothing
    Args    : $source_dir
=cut

sub fastq2fasta
{
print "\nfastq to fasta\n";
my @file_list = glob "$source_dir/*fastq.gz";
for $fastq (@file_list)
	{
	print "\n$fastq\n";
	print $LOG "$fastq\n";
	`gunzip $fastq`;
	}
@file_list = glob "$source_dir/*fastq";
for $fastq (@file_list)
	{
	print "\n$fastq\n";
	print $LOG "$fastq\n";	
	quality_filter("$fastq") ;
	}
} # END sub fastq2fasta 


=item *

quality_filter

    Title   : quality_filter
    Usage   : quality_filter
    Function: applies preset quality cutoffs and converts fastq to fasta
    Returns : nothing
    Args    : fastq
=cut

sub quality_filter
{
print "\nquality filter\n";
my $header ; # first line of fastq set of 4 lines
my $sequence ; # second line of fastq set of 4 lines
my $quality ; # fourth line of fastq set of 4 lines
my @qual ; # split $quality into array
my $qual ; # each array element
my $count ; # count number of reads output
my $output ; # file_name of reads output
my $dna ;
my $line ;
my $bquality ;
my $count_base = 0;
my $ave_qual ;
my $count_out  = 0;
my $seq_len ;
my $count_N ;
my $max_count_N ;
my $base_qual ;
my $sum_base_qual ;
my $min_mean_qual ;
my $mean_qual ;
my $low_qual ;
my $count_low_qual ;
my $max_count_low_qual ;
my $quals ;
my $pass = 0 ;
my $reject_on_low_qual = 0;
my $reject_on_mean_qual = 0;
my $reject_on_count_N = 0;
# get the fastq file
unless ( open(FASTQ,"<", $fastq) ) 
	{
    print "Cannot open file \"$fastq\"\n\n";
    exit;
	}
# set the quality filters
$max_count_N = 1 ;
$min_mean_qual = 30 ;
$low_qual = 20 ;
$max_count_low_qual = 3;

# while in <fastq> count to four (4 line file), capture header (line 1), 
# sequence (line 2), quality stats (line 4)
# as $header ; $sequence ; $quality ;
# in fastq line 3 is either just '+' or '+' then same as line one

$fastq =~ s/\_R1.fastq// ; # remove _R1.fastq
$output = "$fastq\.fna" ; # inherits name from fastq file
open FASTA, ">", $output ;
$output =~ s/\.fna// ;
$quals = "$output\.txt" ;
$quals =~ s/\.txt/_quals/ ;
open QUALS, ">", $quals ;
#loop through fastq file to get the values
$count =1 ;
while (<FASTQ>) 
	{
	chomp $_ ;
	if ($count == 1) 
		{
		$header = $_  ;
		$header =~ s/@/>/ ;
		$count ++ ;
		next ;
		}
	elsif ($count == 2) 
		{
		$sequence = $_ ;
		# get some quality stats here: read_length ($seq_len), number on Ns ($count_N)
		$seq_len = length($sequence);
		$count_N = $sequence ;
		$count_N = ($count_N = tr/N//) ;
		$count ++ ;
		next ;
		}
	elsif ($count == 3) 
	{
	$count ++ ;
	next;
	}
	elsif ($count == 4) 
	{
	$quality = $_;
	# get some quality sats here: mean_qual ($mean_qual), count_low_qual ($count_low_qual)
	# ASCII_2_phred
	@qual = split(//, $quality) ;
	$count_base = 0 ;
	$count_low_qual = 0 ;
	# convert ASCII qscore for each base to numeric, -33 for Sanger phred offset
	$sum_base_qual = 0 ;
	$bquality = "";
	while ($count_base < $seq_len) 
		{
		$base_qual = (ord($qual[$count_base]))-33 ; #convert ASCII to numeric
		$bquality = "$bquality"."$base_qual," ;		
		if ($base_qual <= $low_qual) 
			{
			$count_low_qual ++ ;
			}
		$sum_base_qual = $base_qual+$sum_base_qual ; 
		$count_base ++ ;
		}
		$mean_qual = int(10*$sum_base_qual/($count_base + 1)) ;
		$mean_qual = $mean_qual/10 ;
	$count = 1;
	$line =  "$header\n$sequence\t$bquality\t$seq_len\t$count_N\t$count_low_qual\t$mean_qual";
	#print "$line\n" ;
	if ($count_low_qual > $max_count_low_qual) 
		{
		$reject_on_low_qual ++ ;
		next ;
		}
	elsif ($mean_qual < $min_mean_qual) 
		{
		$reject_on_mean_qual ++ ;
		next;
		}
	elsif ($count_N > $max_count_N) 
		{
		$reject_on_count_N ++ ;
		next;
		}
	else 
		{
		$pass ++ ;
		print FASTA "$header\n$sequence\n" ;
		print QUALS "$bquality\n" ;
		}
	}
}
print "Rejected on low qual:\t$reject_on_low_qual\n";
print "Rejected on mean qual:\t$reject_on_mean_qual\n";
print "Rejected on count Ns:\t$reject_on_count_N\n";
print "Passed:\t$pass\n" ;
print $LOG "Rejected on low qual:\t$reject_on_low_qual\n";
print $LOG "Rejected on mean qual:\t$reject_on_mean_qual\n";
print $LOG "Rejected on count Ns:\t$reject_on_count_N\n";
print $LOG "Passed:\t$pass\n" ;
close FASTQ ;
close FASTA ;
close QUALS ;
} # END sub quality_filter 


=item *

group_fasta_reads

    Title   : group_fasta_reads
    Usage   : group_fasta_reads(directory)
    Function: groups_reads as 2 column list: read counts and read
    Returns : nothing
    Args    : none (can set min read set to include by editing $limit: default is 1)

=cut

sub group_fasta_reads 
{
print "\ngrouping\n";
my @file_list ;

# get the fasta files into an array

@file_list = glob "$source_dir/*fna" ;

for $fasta (@file_list)
	{
	print "\n$fasta\n";
	group_variants("$fasta") ; 
	}
} #END sub group_fasta_reads

#group variants makes a hash table of the variants where the key is the sequence and the value is the count

=item *

group_variants

    Title   : group_variants
    Usage   : group_variants(fasta)
    Function: groups_reads as 2 column list from fasta file
    Returns : nothing
    Args    : none

=cut

sub group_variants 
{
my $dna_filename = $fasta;
my $OUT ;
my $outname ;
my $variant ;
my %count_variant ;
my @count_variant ;
my $count_variant ;
my $limit ;
my $read_count ;
my $grouped_reads_count=0;
my $grouped_reads_sets ;
my $lower_limit ;
my @line ;
my $line ;

# Get the DNA sequence data
chomp $dna_filename;
# Does the file exist?
unless ( -e $dna_filename) 
	{
    print "File \"$dna_filename\" doesn\'t seem to exist!!\n";
    exit;
	}
# Can we open the file?
unless ( open(DNAFILE, "<", $dna_filename) ) 
	{
    print "Cannot open file \"$dna_filename\"\n\n";
    exit;
	}
$limit = 1 ; #edit this to change lowest group size to include
			 #future version make selectable?
chomp $limit ;
while (<DNAFILE>)
	{
	next if ($_ =~ m/^>/) ; 
	$_ =~ s/(.*\t?)(.*\t.*$)/$2/ ;
	$read_count ++ ;
	$count_variant{$_} ++  ;
	}

#$dna_filename =~ s/\./_/ ;
$outname = "$dna_filename\_$limit\_num_grp" ;
unless (open($OUT, '>', $outname )) 
	{
   	print "cannot open the file \"$outname\"\n\n" ;
    exit ;
	}
while (($variant, $count_variant) = each %count_variant ) 
	{
    if ($count_variant >= $limit) 
    	{
     	print $OUT "$count_variant\t$variant";
     	$grouped_reads_count = $grouped_reads_count + $count_variant ;
     	$grouped_reads_sets ++ ; 
		}
	}
close DNAFILE;
print $LOG "total fasta reads = $read_count\n" ;
print $LOG "grouped reads of $limit or more = $grouped_reads_count\n" ;
print $LOG "grouped reads sets = $grouped_reads_sets\n" ; 
print "total fasta reads = $read_count\n" ;
print "grouped reads of $limit or more = $grouped_reads_count\n" ;
print "grouped reads sets = $grouped_reads_sets\n" ;
$grouped_dir ="$output_dir/grouped";
mkdir $grouped_dir, 0755 ;
# move grouped to own directory
my @file_list = glob "$source_dir/*_grp" ;
for $line(@file_list) 
	{
	move($line, $grouped_dir) or die $!;
	}
} #END sub group_variants  

=item *

genotype_grouped_reads

    Title   : genotype_grouped_reads
    Usage   : genotype_grouped_reads(directory)
    Function: take grouped read file (read_counts; read)
    		  report genotypes from lookup table
    Returns : nothing
    Args    : none

=cut

sub genotype_grouped_reads{
# take grouped read file: read counts\tread
# report reads against each genotype target and identify source fragment from lookup table

################
# lookup_table
################
my @file_list ;
my $list ;
my $look_up ;
my $line ;
my @line ;
my @grep ;
my @list ;
print "\ngenotyping\n";
# Get the look up table: 4 column 3 tab separated, no headers, Unix line endings
# columns: gene name refseq /c.DNA/coding/sequence

$look_up = $usual_suspects;

# Does the file exist?
unless ( -e $look_up) {
    print "File \"$look_up\" doesn\'t seem to exist!!\n";
    exit;
	}
 #Can we open the file?
	unless ( open(LOOKUP, $look_up) ) 
	{
    print "Cannot open file \"$look_up\"\n\n";
    exit;
	}
print "\nlook up\t$look_up\n";
@list = <LOOKUP> ; #@list is the list of usual suspects
close LOOKUP ;
##################
#count_reads
##################

my $read = "";
my $revcomp ;
my $name = "";
my $count_regex ; 
my $count_rcregex ;
my %count_regex ;
my %count_rcregex ;
my %total ;
my $total ;
my $fasta ;
my @fasta ;
my $dna_filename ;
my @positive_f ;
my $positive_f ;
my @positive_r ;
my $positive_r ;
my @negative ;
my $negative ;
my $output_positive_f = 'PosF' ;
my $output_positive_r = 'PosR' ;
my $output_negative  = 'Neg' ;
my $genotype_report ; 
my $HANDLE_POS_F ;
my $HANDLE_POS_R ;
my $HANDLE_NEG ;
my $OUTPUT_GENOTYPE ;
my $total_reads = 0;
my @unmatchedfasta ;
my $unmatchedfasta ;
my @sequence ;
my $sequence ;
my $readcount ;
my $gene ;
my $codon ;
my $seq ;
my @seq ;
my $mapped_reads = 0 ;
my $percent_mapped ;

# Get the grouped read files as glob

@file_list = glob "$grouped_dir/*_grp" ; #@file_list is the list of grouped reads to match with @list
print "\nfile list\n@file_list\n" ;
$genotype_dir ="$output_dir/genotype";
mkdir $genotype_dir, 0755 ;  

#loop through the file list

for $dna_filename(@file_list) # loop through @file_list 
	{ 
	chomp $dna_filename;
	# Does the file exist?
	unless ( -e $dna_filename) 
		{
		print "File \"$dna_filename\" doesn\'t seem to exist!!\n";
		exit;
		}
	# Can we open the file?
	unless ( open(DNAFILE, $dna_filename) ) 
		{ 
		print "Cannot open file \"$dna_filename\"\n\n";
		exit;
		} 
	print "\nfile being processed\n$dna_filename\n";
	#name the output files
	$genotype_report = "$dna_filename"."_Genotypes.txt" ;
	$output_positive_f = "$dna_filename"."_output_positive_f\.fna" ;
	$output_positive_r = "$dna_filename"."_output_positive_r\.fna" ;
	$output_negative = "$dna_filename"."_output_negative\.txt" ;

	@grep = <DNAFILE>;

	close DNAFILE;

	@fasta = grep (/^[1-9]/,@grep) ; #select grouped sequence lines only: must begin with a number
	for $seq (@fasta) 
	   { 
        @sequence = split (/\t/, $seq ) ;
        $readcount = $sequence[0] ;
        $fasta = $sequence[1] ;
        $total_reads= $total_reads + $readcount ;
        } 

	@unmatchedfasta = @fasta ; # bin for unmatched reads

	unless (open($OUTPUT_GENOTYPE, '>', $genotype_report )) 
		{ 
		print "cannot open the genotype file \"$genotype_report\"\n\n" ;
		exit ;
		} 

	unless (open ($HANDLE_POS_F, '>',$output_positive_f )) 
		{
		print "cannot open the positive f file \"$output_positive_f\"\n\n" ;
		exit ;
		}
		
	unless (open ($HANDLE_POS_R, '>',$output_positive_r ))  
		{
		print "cannot open the positive r file \"$output_positive_r\"\n\n" ;
		exit ;
		}
	unless (open ( $HANDLE_NEG, '>',$output_negative )) 
		{
		print "cannot open the negative file \"$output_negative\"\n\n";
		exit ;
		}
		
	print "gene\ttotal\tforward\treverse\n" ; #header for output
	print $OUTPUT_GENOTYPE  "gene\ttotal\tforward\treverse\n" ;
	foreach $line (@list) #loop through usual suspects list
		{ 
		chomp $line ;
		@line = split (/\t/, $line ) ;
		$read = $line[3] ; #read from lookup table
		$name = $line[2] ;
		$codon = $line[1] ;
		$gene = $line[0] ;
		$name = "$gene "."$codon "."$name" ;
		$revcomp = reverse $read ;
		$revcomp =~ tr/ACGTacgt\[\]\.\*/TGCAtgca\]\[\*\./ ;
		@unmatchedfasta = grep (!/$read/i,@unmatchedfasta) ; #remove matches 
		@unmatchedfasta = grep (!/$revcomp/i,@unmatchedfasta) ; #remove matches 
		$count_regex{$name} = 0; 
		$count_rcregex{$name} = 0 ;
		for $seq(@fasta) 
			{
			@sequence = split (/\t/, $seq ) ;
			$readcount = $sequence[0] ;
			$fasta = $sequence[1] ; #read output from run
				if ($fasta =~ m/$read/i ) 
				{
				$count_regex{$name} = $count_regex{$name} + $readcount ;  #count matches. 
				#push @positive_f, "$name\t$readcount\>\n$fasta\n" ; # store matches
				$positive_f = "$name\t$readcount\t$fasta\n" ;
				print $HANDLE_POS_F "$positive_f" ;
				$mapped_reads = $mapped_reads + $readcount ;
				}
			   elsif ($fasta =~ m/$revcomp/i ) 
			   	{
				$count_rcregex{$name} = $count_rcregex{$name} + $readcount ; #count matches. 
				#push  @positive_r, "$name\t$readcount\>\n$fasta\n" ; #store matches
				$positive_r = "$name\t$readcount\t$fasta\n" ;
				print $HANDLE_POS_R "$positive_r" ;
				#consider revcomp and pushing into a single array
				$mapped_reads = $mapped_reads + $readcount ;
				}		
			}# end for $seq(@fasta)
		$total{$name} = $count_regex{$name} + $count_rcregex{$name}; 
		print "$name\t$total{$name}\t$count_regex{$name}\t$count_rcregex{$name}\n";
		print  $OUTPUT_GENOTYPE "$name\t$total{$name}\t$count_regex{$name}\t$count_rcregex{$name}\n";
		} #end loop through usual suspects list
	$percent_mapped = 100*$mapped_reads/$total_reads; 
	print "mapped reads\t$mapped_reads\ttotal reads\t$total_reads\tpercent mapped\t$percent_mapped\n" ;
	print $LOG "$dna_filename\nmapped reads\t$mapped_reads\ntotal reads\t$total_reads\npercent mapped\t$percent_mapped\n" ;
	close ($OUTPUT_GENOTYPE ) ;
	close ($HANDLE_POS_F) ;
	close ($HANDLE_POS_R) ;
	close ($HANDLE_NEG) ;
	} # end loop through @file_list

# move genotype files to genotype directory
@file_list = glob "$grouped_dir/*_grp_*" ;
	for $line(@file_list) 
	{
	move($line, $genotype_dir) or die $!;
	}
} # END sub genotype_grouped_reads

=item *

grab_by_flanks

    Title   : grab_by_flanks
    Usage   : grab_by_flanks(directory)
    Function: count matches to amplicon flanks
    Returns : nothing
    Args    : none
=cut

sub grab_by_flanks 
{

################
# lookup_table
################
my @file_list ;
my $list ;
my $look_up ;
my $line ;
my @line ;
my @grep ;
my @list ;

# Get the look up table: 4 column tab separated, no headers, Unix line endings
# columns: gene name /chr/positions/flanking sequence with capture brackets


$look_up = $amplicon_flanks ;

chomp $look_up;

# Does the file exist?
unless ( -e $look_up) 
	{
    print "File \"$look_up\" doesn\'t seem to exist!!\n";
    exit;
	}
 #Can we open the file?
	unless ( open(LOOKUP, $look_up) ) 
	{
    print "Cannot open file \"$look_up\"\n\n";
    exit;
	}
@list = <LOOKUP> ;
close LOOKUP ;
##################
#count_reads
##################

my $read = "";
my $revcomp ;
my $name = "";
my $count_regex ; 
my $count_rcregex ;
my %count_regex ;
my %count_rcregex ;
my %total ;
my $total ;
my $fasta ;
my @fasta ;
my $dna_filename ;
my @positive_f ;
my $positive_f ;
my @positive_r ;
my $positive_r ;
my @negative ;
my $negative ;
my $output_positive_f = 'PosF' ;
my $output_positive_r = 'PosR' ;
my $output_negative  = 'Neg' ;
my $flanked_report ; 
my $HANDLE_POS_F ;
my $HANDLE_POS_R ;
my $HANDLE_NEG ;
my $OUTPUT_GENOTYPE ;
my $total_reads = 0;
my @unmatchedfasta ;
my $unmatchedfasta ;
my @sequence ;
my $sequence ;
my $readcount ;
my $gene ;
my $codon ;
my $seq ;
my @seq ;
my $mapped_reads = 0 ;
my $percent_mapped ;

# Get the grouped read files as glob

@file_list = glob "$grouped_dir/*_grp" ; #@file_list is the list of grouped reads to match with @list
print "\nfile list\n@file_list\n" ;
$flanked_dir ="$output_dir/flanked";
mkdir $flanked_dir, 0755 ; 

# loop through file list;

for $dna_filename(@file_list) # loop through file list;
	{
	chomp $dna_filename;
	# Does the file exist?
	unless ( -e $dna_filename)
		{   
		print "File \"$dna_filename\" doesn\'t seem to exist!!\n";
		exit;
		}
	# Can we open the file?
	unless ( open(DNAFILE, $dna_filename) )
		{
		print "Cannot open file \"$dna_filename\"\n\n";
		exit;
		}
	print "\nfile being processed\t$dna_filename\n";
	$flanked_report = "$dna_filename"."_flanked.txt" ;
	$output_positive_f = "$dna_filename"."_flanked_positive_f" ;
	$output_positive_r = "$dna_filename"."_flanked_positive_r" ;
	$output_negative = "$dna_filename"."_flanked_negative" ;

	@grep = <DNAFILE>;

	close DNAFILE;

	@fasta = grep (/^[1-9]/,@grep) ; #select grouped sequence lines only: must begin with a number
	for $seq (@fasta) 
	   { 
        @sequence = split (/\t/, $seq ) ;
        $readcount = $sequence[0] ;
        $fasta = $sequence[1] ;
        $total_reads= $total_reads + $readcount ;
        } 

	@unmatchedfasta = @fasta ; # bin for unmatched reads

	unless (open($OUTPUT_GENOTYPE, '>', $flanked_report )) 
		{ 
		print "cannot open the genotype file \"$flanked_report\"\n\n" ;
		exit ;
		} 

	unless (open ($HANDLE_POS_F, '>',$output_positive_f )) 
		{
		print "cannot open the positive f file \"$output_positive_f\"\n\n" ;
		exit ;
		}
		
	unless (open ($HANDLE_POS_R, '>',$output_positive_r ))  
		{
		print "cannot open the positive r file \"$output_positive_r\"\n\n" ;
		exit ;
		}
	unless (open ( $HANDLE_NEG, '>',$output_negative )) 
		{
		print "cannot open the negative file \"$output_negative\"\n\n";
		exit ;
		}
		
	print "gene\ttotal\tforward\treverse\n" ; #header for output
	print $OUTPUT_GENOTYPE  "gene\ttotal\tforward\treverse\n" ;
	foreach $line (@list) #loop through flanks list
		{ 
		chomp $line ;
		@line = split (/\t/, $line ) ;
		$read = $line[3] ; #read from lookup table
		$name = $line[2] ;
		$codon = $line[1] ;
		$gene = $line[0] ;
		$name = "$gene "."$codon "."$name" ;
		$revcomp = reverse $read ;
		$revcomp =~ tr/ACGTacgt\[\]\.\*\(\)/TGCAtgca\]\[\*\.\)\(/ ;
		@unmatchedfasta = grep (!/$read/i,@unmatchedfasta) ; #remove matches 
		@unmatchedfasta = grep (!/$revcomp/i,@unmatchedfasta) ; #remove matches 
		#
		### Need to reset all regexes here or reset @positive_f and @positive_r ###
		$count_regex{$name} = 0; 
		$count_rcregex{$name} = 0 ;
			for $seq(@fasta) 
			{
			@sequence = split (/\t/, $seq ) ;
			$readcount = $sequence[0] ;
			$fasta = $sequence[1] ; #read output from run #should be $1
				if ($fasta =~ m/$read/i ) 
				{
				#$fasta=$1; remove
				$count_regex{$name} = $count_regex{$name} + $readcount ;  #count matches. 
				$positive_f = "$name\t$readcount\t$1\n" ;
				print $HANDLE_POS_F "$positive_f" ;
				$mapped_reads = $mapped_reads + $readcount ;
				}
			   elsif ($fasta =~ m/$revcomp/i ) 
			   	{
			   	$fasta=reverse($1);
			   	$fasta=~tr/ACGTacgt/TGCAtgca/ ;
				$count_rcregex{$name} = $count_rcregex{$name} + $readcount ; #count matches. 
				$positive_r = "$name\t$readcount\t$fasta\n" ;
				print $HANDLE_POS_F "$positive_r" ; #putting all matches in POS_F file
				$mapped_reads = $mapped_reads + $readcount ;
				}		
			}# end for $seq(@fasta)
		$total{$name} = $count_regex{$name} + $count_rcregex{$name}; 
		print "$name\t$total{$name}\t$count_regex{$name}\t$count_rcregex{$name}\n";
		print  $OUTPUT_GENOTYPE "$name\t$total{$name}\t$count_regex{$name}\t$count_rcregex{$name}\n";
		} #end loop through flanks list
	$percent_mapped = 100*$mapped_reads/$total_reads; 
	print "mapped reads\t$mapped_reads\ttotal reads\t$total_reads\tpercent mapped\t$percent_mapped\n" ;
	print $LOG "$dna_filename\nmapped reads\t$mapped_reads\ntotal reads\t$total_reads\npercent mapped\t$percent_mapped\n" ;
	close ($OUTPUT_GENOTYPE ) ;
	close ($HANDLE_NEG) ;
	} # end loop through @file_list
# move genotype files to genotype directory
@file_list = glob "$grouped_dir/*_grp_*" ;
	for $line(@file_list) 
		{
		move($line, $flanked_dir) or die $!;
		}
} # END sub grab_by_flanks



=item *

sort_by_locus

    Title   : sort_by_locus
    Usage   : take three tabbed read file output from genotype_grouped_read
    Function: Sort by locus and write new file for each locus as read count and sequence
    Returns : nothing
    Args    : nothing

=cut

sub sort_by_locus
{
# take three tabbed read file output from genotype_grouped_reads
# locus\tcount\tsequence
# sort by locus and write new file for each locus as read count and sequence


################
# lookup_table
################
my @file_list ;
my @list ;
my @look_up ;
my $look_up ;
my @line ;
my $line ;
my $read ;
my $count ;
my $name ;
my $filename ;
my $shortname ;

print "\nsorting\n";

# Get the look up table: 4 column tab separated, no headers, Unix line endings
# columns: gene name /codon/description/sequence
#print "Please type the filename of the look up table ";

$look_up = $amplicon_flanks ; #default

chomp $look_up;

print "\nsorting\n";

# Does the file exist?
unless ( -e $look_up) 
	{
    print "File \"$look_up\" doesn\'t seem to exist!!\n";
    exit;
	}
 #Can we open the file?
	unless ( open(LOOKUP, $look_up) ) 
	{
    print "Cannot open file \"$look_up\"\n\n";
    exit;
	}
@look_up = <LOOKUP> ;
close LOOKUP ;
print "\nlook up table\t$look_up\n"; #table for flanking reads

@file_list = glob "$flanked_dir/*_grp_flanked_positive_f"; #?fixed this by revcomp and merging
# positive_f and positive_r
print "\nfile list\t@file_list\n";
$locus_dir = "$output_dir/locus";
print "\nlocus dir\t$locus_dir\n";
mkdir $locus_dir, 0755 ;

for $filename(@file_list) #loop through file list
	{
	chomp $filename ;
	$shortname = $filename ;
	$shortname =~ s/_fna_//;
	$shortname =~s/merged_//;
	print "\nfilename\n";

	for $line(@look_up) 
		{
		chomp $line ;
		print "\n$line\n";
    	@line = split (/\t/, $line ) ;
    	$read = $line[2] ;
    	$count = $line[1] ;
    	$name = $line[0] ;
    	`grep $name $filename >$shortname\_$name\_locus` ;
    			}
	} #end for filename
@file_list = glob "$flanked_dir/*_locus" ;
	for $line(@file_list) 
	{
	move($line, $locus_dir) or die $!;
	}
} #end sub sort by locus

=item *
sort_each_locus

  	Title   : sort_each_locus
    Usage   : sort_each_locus(directory)
    Function: Sorts the locus by read count
    Returns : nothing
    Args    : none
=cut

sub sort_each_locus
{
my @file_list;
my $dna_filename ;

# Get the DNA sequence data
print "\ngetting list of variants\n";


@file_list = glob "$locus_dir/*_locus" ;

for $dna_filename(@file_list)
	{
	print "$dna_filename\n";
	sort_by_count($dna_filename) ;
	}
}
sub sort_by_count
{
my $dna_filename ;
my $OUT ;
my $LOG ;
my $logname ;
my $outname ;
my $variant ;
my %count_variant ;
my @count_variant ;
my $count_variant ;
my $limit ;
my $read ;
my $read_count= 0;
my $grouped_reads_count=0;
my $grouped_reads_sets = 0 ;
my $lower_limit ;
my @line ;
my $line ;
my $locus ;
my @variants ;
my %variants ;
my $variants ;
my @key ;
my $key ;
my %key ;
my %value ;
my @value;
my $value ;

$dna_filename = $_[0];
chomp $dna_filename;
	# Does the file exist?
		unless ( -e $dna_filename) 
		{
		print "File \"$dna_filename\" doesn\'t seem to exist!!\n";
		exit;
		}
		# Can we open the file?
		unless ( open(DNAFILE, "<", $dna_filename) ) 
		{
		print "Cannot open file \"$dna_filename\"\n\n";
		exit;
		}
	$limit = 1	;
	print "\nlower limit of reads included is $limit\n" ;
 	while (<DNAFILE>)
		{
		@line = split (/\t/ , $_) ;
		$locus = $line[0];
		$count_variant = $line[1] ;
		$variant = $line[2] ;
		$read_count ++ ;
		$count_variant{$variant} =  $count_variant{$variant} += $count_variant ;
		}
#		$dna_filename =~ s/\./_/ ;
		$outname = "$dna_filename\_$limit\_num_grp\_srt\.txt" ;
		$logname = "$dna_filename\_$limit\_num_grp\_log" ;
		unless (open($OUT, '>', $outname )) 
		{
		print "cannot open the file \"$outname\"\n\n" ;
		exit ;
		}
		unless (open($LOG, '>', $logname )) 
		{
		print "cannot open the file \"$logname\"\n\n" ;
		exit ;
		}
		while (($variant,$count_variant) = each (%count_variant)) 
		{
		if ($count_variant >= $limit) 
			{
			chomp $variant; chomp $count_variant ; 
			push (@variants,"$count_variant\t$variant\n") ;
			$key{$variant} = $count_variant ;
			$grouped_reads_count = $grouped_reads_count += $count_variant ;
			$grouped_reads_sets ++ ; 	
			}
		}	
	
	@key = sort { $key{$b} <=> $key{$a} or $a cmp $b} keys %key ;
	foreach $key(@key)
		{
		#print "\$locus\t$key{$key}\n$key\n";
		print $OUT "\>$locus\t$key\t$key{$key}\n"; #one line or fasta
		}
close DNAFILE;
    print $LOG "total reads = $read_count\n" ;
    print $LOG "grouped reads of $limit or more = $grouped_reads_count\n" ;
    print $LOG "grouped reads sets = $grouped_reads_sets\n" ; 
    print "total reads = $read_count\n" ;
    print "grouped reads of $limit or more = $grouped_reads_count\n" ;
    print "grouped reads sets = $grouped_reads_sets\n" ; 
close $LOG ;
close $OUT ;
}

=head1

tidy_up

    Title   : tidy_up
    Usage   : tidy_up()
    Function: Moves files into suitable folders leaving everything neat
    Returns : Nothing
    Args    : top level directory
=cut

sub tidy_up
{
my $line ;
my @glob_list ;
$grouped_dir ="$output_dir/grouped";
$locus_dir = "$output_dir/locus" ;
$sorted_locus_dir = "$output_dir/sorted" ;
$fna_dir = "$output_dir/fasta" ;
$quals_dir = "$output_dir/qual_scores" ;


mkdir $sorted_locus_dir, 0755 ;
mkdir $quals_dir, 0755 ;
mkdir $fna_dir, 0755 ;

#move fna to own directory
@glob_list = glob "$source_dir/*fna" ;
for $line(@glob_list) 
	{
	print "moving $line to $fna_dir\n";
	print $LOG "moving $line to $fna_dir\n";
	move($line, $fna_dir) or die $!;
	}

#move sorted loci to own directory
@glob_list = glob "$locus_dir/*_srt\.txt" ;
for $line(@glob_list) 
	{
	print "moving $line to $grouped_dir\n";
	print $LOG "moving $line to $sorted_locus_dir\n";
	move($line, $sorted_locus_dir) or die $!;
	}
	
#move quals to own directory
@glob_list = glob "$source_dir/*_quals" ;
for $line(@glob_list) 
	{
	print "moving $line to $quals_dir\n";
	print $LOG "moving $line to $quals_dir\n";
	move($line, $quals_dir) or die $!;
	}

} # end sub tidy files

