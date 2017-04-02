#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 03.04.2017

###############################################################################
###############################################################################
### Remove genes with low expression.                                       ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my (
    # variable for parameters which are read in from commandline
    $help,
    $infile,
    $minExpress,
    $prefix,
    $verbose,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"           => \$help,         # Print this help
    "in=s"            => \$infile,       # Matrix with abundance data.
    "minExpression=f" => \$minExpress,   # Maximum abundance of a transcript in all samples, in order to be removed. E.g. 5 removes 95% of olfactory receptor proteins.
    "prefix:s"        => \$prefix,       # Include only genes with this prefix.
    "verbose!"        => \$verbose,      # output additional info
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################

##############################################################################
### SUBROUTINES ##############################################################
############################################################################## 

###############################################################################
sub printHelp {
###############################################################################
    # prints a help about the using and parameters of this scripts 
    # (execute if user types commandline parameter -h)
    # param:  no paramaters
    # return: no return value

    my (
	$usage,
	$sourceCode,
	@rows,
	$row,
	$option,
	$scriptInfo,
	$example,
       );

    $usage = "$0\n";


    print "\nUsage: " .  $usage . "\n";

    print "Valid options are:\n\n";
    open(MYSELF, "$0") or
      die "Cannot read source code file $0: $!\n";
    $sourceCode .= join "", <MYSELF>;
    close MYSELF;
    $sourceCode =~ s/^.+?\&GetOptions\(\n//s;
    $sourceCode =~ s/\n\).+$//s;
    @rows = split /\n/, $sourceCode;
    foreach $row (@rows){
        $option = $row;
	$option =~ s/\s+\"//g;
	$option =~ s/\"\s.+\#/\t\#/g;
	$option =~ s/=./\t<value> [required]/;
	$option =~ s/:./\t<value> [optional]/;
	$option =~ s/!/\t<non value> [optional]/;

	$row =~ s/^.*//;
	print "\t";
	printf("%-1s%-30s%-30s\n", "-",$option,$row);

    } # end of foreach $row (@rows)
    print "\n";
    print "Options may be abreviated, e.g. -h for --help\n\n";

    $example  = "$0";
}
##############################################################################
sub readCountFile {

    open(F2,"zcat $infile|") or die "\nERROR: Failed to open $infile: $!\n\n";
    my $h=<F2>;
    chomp($h);
    my @h=split(/\t/, $h);

    # print out header
    print "$h\n";

    my %data = ();
    my %relevantGenes = ();
    while(my $l = <F2>){
	next if($l !~ /^$prefix/);
	chomp($l);
	my @a=split(/\t/, $l);
	for(my $i = 1; $i < @a; $i++) {
	    $data{$a[0]}->{$i} = $a[$i];
	    $relevantGenes{$a[0]} = 1 if($a[$i] > $minExpress);
	}
    }
    close(F2);
    return (\%data, \%relevantGenes);
}
##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############
if(!defined $infile) {
    die "\nMissing infile argument.\n\tUsage: $0 -in freeze3_v2.kallisto.lib.trans.raw.wl.aliquot_id_GTEX_1.4.pcawg.transcripts.raw.tsv.gz -prefix ENS\n\n";
}
if(!defined $minExpress) {
    die "\nMissing minExpress argument.\n\tUsage: $0 -in freeze3_v2.kallisto.lib.trans.raw.wl.aliquot_id_GTEX_1.4.pcawg.transcripts.raw.tsv.gz -prefix ENS\n\n";
}
if(!defined $prefix) {
    die "\nMissing prefix argument.\n\tUsage: $0 -in freeze3_v2.kallisto.lib.trans.raw.wl.aliquot_id_GTEX_1.4.pcawg.transcripts.raw.tsv.gz -prefix ENS\n\n";
}

print STDERR "2Reading in $infile ... " if($verbose);
my ($data, $relevantGenes) = &readCountFile();
print STDERR "done\n" if($verbose);

print STDERR "Filtering ...\n" if($verbose);
foreach my $ensg (sort keys %$data) {
    my $maxAbs = 0;
    my $maxRel = 0;
    next if(!exists $relevantGenes->{$ensg});
    print $ensg;
    foreach my $i (sort {$a<=>$b} keys %{$data->{$ensg}}) {
	print "\t$data->{$ensg}->{$i}";
    }
    print "\n";
}
print STDERR "All DONE\n";
