#!/usr/bin/perl
use strict;
use warnings;

=head1 NAME

amplisort.pl

=head1 VERSION

0.8
Last updated on 2014/02/17

=head1 DESCRIPTION

Process amplicon output from amplivar to output relevant FASTA files
After amplivar has run there should be a folder "sorted" containing all files sorted by
locus, based on the flanking file table used in amplivar.
This tool enables selectable, rapid interrogation of the sorted folder to produce fasta
files filtered by locus, sample and relative abundance compared to the most abundant
transcript for that locus

=head1 INPUT

 uses @ARGV to input six parameters as space separated variables
	1: locus directory: text; full path to directory of sorted read files from amplivar
	2: gene ID: text; any abbreviation that pulls out the gene or amplicon of interest. Or "." for all genes
	3: sample: text; any unique identifier for samples in the directory. Or "." for all samples
	4: min reads: integer; lowest number of most abundant read set to analyse.  Lower numbers are rejected
	5: minimum percentage: integer;  minimum percentage of minor allele to output
	6: suppress: integer; 1 (default) suppress 100% allele output, 0 include all 100% alleles above minimum percentage


=head1 OUTPUT

fasta files for selected reads &
fasta files for reads not passing minimum read criteria &
tsv file for read count distribution

=head1 OTHER REQUIREMENTS

none

=head1 EXMAPLE USAGE

perl amplisort.pl  <locus_directory gene_ID> <gene_ID> <sample_ID> <min_reads> <min_percent> <suppress_100%[1/0]>

=head1 AUTHOR

Graham Taylor, University of Melbourne

=cut

# START OF CODE

our @file_list ;
our $locus_list ;
our $locus_directory ;
our $gene_ID ; # gene abbreviation
our $min_reads ; # min reads in most abundant read
our $min_pc ; # min percentage reads
our $sample ; #input sample ID
our $suppress_100=1 ;
our $dir ;
our @read_set ;
our $LOW_COVERAGE ; # file handle to print low coverage loci to
our $LOG ; # file handle for log
our $OUTPUT ; #file handle for output
our $READCOUNT ; #file handle for read counts
our $sample_file ;
our $gene_ID_file ;
our $hundred_percent ="" ;
our $hundred_percent_reads="" ;
our $header_plus_read = "";

# Get the DNA sequence data
my $arg_count = scalar @ARGV ;
print "perl amplisort.pl arguments passed = $arg_count\n";
if ($arg_count != 6)
	{
	print "amplisort.pl outputs reads above a miniumum threshold selectable by gene abbreviation\n";
	print "6 arguments needed: first three as text, second three as numbers\n" ;
	print "enter perl amplisort.pl\n1\t\<locus_directory\>\n2\t\<gene abbreviation\> \(use \. for all genes\)\n3\t\<sample identifer\> \(. for all samples\)\n4\t\<lowest read count\>\n5\t\<lowest percentage reads\>\n6\t<suppress 100% reads 1 yes or 0 no\>\n";
	print "type \"perl amplisort.pl \<directory\> \<gene abbrev\> \<sample ID\> \<min read count\> \<min percent minor allele\> \<1 or 0 to skip 100% alleles>\"\n" ; 
	print "\n\n";
	exit ;
	}

$locus_directory = $ARGV[0] ;
$gene_ID = $ARGV[1] ;
$sample = $ARGV[2] ;
$min_reads = $ARGV[3] ;
$min_pc = $ARGV[4] ;
$suppress_100 = $ARGV[5] ;
($suppress_100 = 1) if $suppress_100 !=0 ;
$sample_file = $sample ;
$sample_file = "all_samples" if ($sample eq ".") ;
$gene_ID_file = $gene_ID ;
$gene_ID_file = "all_genes" if ($gene_ID eq ".") ;

#log timestamp
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
our $date = sprintf("%4d%02d%02d%02d%02d%02d",($year + 1900),($mon+1),$mday,$hour,$min,$sec);
our $log ="amplisort_$date.log";
open ($LOG, ">", "$log") ;
open ($LOW_COVERAGE, ">", "$sample_file\.max_read_count_below\_$min_reads\_$gene_ID_file\.fna" ) ;
open ($OUTPUT, ">", "$sample_file\.amplivar_$gene_ID_file\_depth_over_$min_reads\_minor_alleles_over_$min_pc\_percent\.fna") ;
open ($READCOUNT, ">", "$sample_file\.amplivar_$gene_ID_file\_depth_over_$min_reads\_minor_alleles_over_$min_pc\_percent\.tsv") ;

print "\nlocus directory is $locus_directory\n";
print "gene_ID is $gene_ID\n";
print "sample is $sample\n";
print "min depth is $min_reads\n";
print "min percent minor allele = $min_pc\n";
print "skipping 100% alleles set to $suppress_100\n" ;

print $LOG "\nlocus directory is $locus_directory\n";
print $LOG "gene_ID is $gene_ID\n";
print $LOG "sample is $sample\n";
print $LOG "min depth is $min_reads\n";
print $LOG "min percent minor allele = $min_pc\n";
print $LOG "skipping 100% alleles set to $suppress_100\n" ;

print "\ngetting list of variants\n";
print $LOG "\ngetting list of variants\n";

get_file_list() ;

for $locus_list(@file_list)
	{
	counts_of_reads($locus_list) ;
	}
	
close $LOG ;
close $OUTPUT ;
close $LOW_COVERAGE ;
close $READCOUNT ;
exit;

################
##subroutines##
###############


=head2 get_file_list

    Title   : get_file_list
    Usage   : get_file_list($source_dir)
    Function: globs file list from source directory
    Returns : @file_list
    Args    : $source_dir
=cut


sub get_file_list
	{
	print "sub get_file_list\n";
	print $LOG "sub get_file_list\n";
	@file_list = glob "$locus_directory/*txt" ;
	}#end sub get_file_list
	
=head2 counts_of_reads

    Title   : counts_of_reads
    Usage   : counts_of_reads($locus_list)
    Function: counts reads from each amplicon
    Returns : prints FASTA files
    Args    : @file_list
=cut


sub counts_of_reads
	{
	my $dnafile ;
	my $line ;
	my @line ;
	my $read ;
	my $read_count ;
	my $line_count ;
	my %line_counter ; #$line_counter{$read_count}
	my $ratio ; #ratio of read value to highest value
	my $sample_ID ;
	my @sample_description ;
	my $run_date ;
	my $library_name ;
	my $sample_name ;
	my $locus_ID ;
	my $read_length ;
	$locus_list = $_[0];
	chomp $locus_list;
	# Does the file exist?
	unless ( -e $locus_list) 
		{
		print "File \"$locus_list\" doesn\'t seem to exist!!\n";
		exit;
		}
	# Can we open the file?
	unless ( open($dnafile, "<", $locus_list) ) 
		{
		print "Cannot open file \"$locus_list\"\n\n";
		exit;
		}
#	$sample_ID = $locus_list ;
#	$sample_ID =~ s/$locus_directory\/// ;
#	@sample_description = split('_',$sample_ID) ; #New file naming edit
#	$run_date = $sample_description[0] ; #New file naming edit
#	$library_name = $sample_description[3] ; #New file naming edit
#	$sample_name = $sample_description[4] ; #New file naming edit
#	$sample_ID = "$run_date\_$library_name\_$sample_name" ; #New file naming edit
	

	my @sorted_list = <$dnafile>;
	$line_count = 0 ;
	for $line (@sorted_list)
		{
		chomp $line ;
#		next unless ($sample =~ m/($run_date)*($library_name)*($sample_name)/) ;
		next unless ($line =~ m/$gene_ID/) ;
		@line = split (/\t/, $line) ;
		$read_count = $line[2] ;
		$read = $line[1] ;
		$read_length = length($read) ;
		next if ($read_length<22) ; #Careful you don't throw the baby out with the bathwater
		$line_count ++ ;
		$locus_ID = $line[0] ;
		$line_counter{$line_count} = $read_count ; #add to hash
		$ratio = 100000*$line_counter{$line_count}/$line_counter{1} ;
		$ratio =int($ratio+0.5)/1000 ;
		if ($line_counter{1}<$min_reads)
			{
			$locus_ID =~s/\>/$line_count\_/ ;
			print "$locus_ID\n" ;
			print $LOW_COVERAGE "\>$sample_file $locus_ID $ratio\% reads $read_count\n$read\n" if ($ratio ==100) ;
			next ;
			}
		next if ($ratio<$min_pc) ;
			$locus_ID =~s/\>/$line_count\_/ ;
			if ($suppress_100 == 1 && $line_count == 1) 
				{
#				$locus_ID =~s/\>/$line_count\_/ ;
				$hundred_percent = "\>$sample_file $locus_ID reads $read_count $ratio\% length $read_length\n$read\n" ;
				$hundred_percent=~ s/ /_/g ;
				$hundred_percent_reads = "\>$sample_file\_$locus_ID\t$read_count\t$read_length\n" ;
				next ;
				}
			if ($line_count==2)
				{
				print $hundred_percent;
				print $OUTPUT $hundred_percent;
				print $READCOUNT $hundred_percent_reads ;
				}
			$locus_ID =~ s/\s+$//g;
			#my $LID = $locus_ID;
			my @tokens = split(/\s+/, $locus_ID);
			my $LID = "$tokens[0]\_$tokens[1]\_$tokens[2]";
			$header_plus_read = "\>$sample_file $LID reads $read_count $ratio\% length $read_length\n$read\n";
			$header_plus_read =~ s/ /_/g;
			print $header_plus_read ;
			print $OUTPUT $header_plus_read ;
			print $READCOUNT "\>$sample_file\_$LID\t$read_count\t$read_length\n" ;
			}
	}#end sub counts_of_reads
