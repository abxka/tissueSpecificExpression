#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 13.03.2017

###############################################################################
###############################################################################
### Assignes tissue specific information to set of genes.                   ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;

my (
    # variable for parameters which are read in from commandline
    $help,
    $geneFile,
    $tissueFile,
    $mapFile,
    $minEnrichment,
    $maxPvalue,
    );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"           => \$help,          # print this help
    "geneFile=s"      => \$geneFile,      # e.g. genes.txt
    "tissueFile=s"    => \$tissueFile,    # e.g. tissueSpecificGeneExpression_meta_mean_mutatedRegions.tsv.gz
    "mapFile=s"       => \$mapFile,       # e.g. ensg_ensp_enst_ense_geneName_v75.tsv.gz
    "minEnrichment=f" => \$minEnrichment, # e.g. 5
    "maxPvalue=f"     => \$maxPvalue,     # e.g. 0.01
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
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############


if(!defined $tissueFile) {
    print STDERR "\n\tPlease provide a tissue file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $geneFile) {
    print STDERR "\n\tPlease provide a gene file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $mapFile) {
    print STDERR "\n\tPlease provide a ENSG2genename file. Try $0 -help to get more information\n\n";
    exit;
}


open(F1, "zcat $mapFile|") or die "\nERROR: Failed to open $mapFile: $!\n\n";
my %m = ();
my %d = ();
while(my $l = <F1>) {
    chomp($l);
    my @a = split(/\t/, $l);
    my $ensg = $a[0];
    my $name = $a[4];
    my $db = $a[5];
    if(exists $m{$ensg}) {
	if($db eq "HGNC Symbol") {
	    $m{$ensg} = $a[4];
	    $d{$ensg} = $a[5];
	}
    } else {
	$m{$ensg} = $a[4];
	$d{$ensg} = $a[5];
    }
}
close(F1);

open(F2, "zcat $tissueFile|") or die "\nERROR: Failed to open $tissueFile: $!\n\n";
my %h = ();
my %p = ();
while(my $l = <F2>) {
    next if($l =~ /^#/);
    chomp($l);
    my @a = split(/\t/, $l);
    my $ensg = $a[2];
    $ensg =~ s/\..*//;
    $h{$m{$ensg}} = "$a[0]\t$a[1]\t$a[5]\t$a[3]\t$a[4]\t$a[6]\t".(($a[5] > $minEnrichment and $a[6] < $maxPvalue) ? 1 : 0);
    $p{$m{$ensg}} = $a[6];
}
close(F2);

open(F3, $geneFile) or die "\nERROR: Failed to open $geneFile: $!\n\n";
print "#TestedGene\tTissueWithHighestMeanExpression\tTissueWith2ndHighestMeanExpression\tMeanExpressionOfTissueWithHighestMeanExpression\tMeanExpressionOfTissueWith2ndHighestMeanExpression\tEnrichment\tWilcoxRankSumTestPvalue\tIsTissueSpecific(Enrichment>$minEnrichment,Pvalue<$maxPvalue)\n";
while(my $l = <F3>){
    next if($l =~ /^#/);
    chomp($l);
    my @a = split(/;/, $l);
    my $minS = "";
    my $minP;
    foreach my $a(@a){
	if(exists $h{$a}){
	    if($minS eq "") {
		$minS = "$a\t$h{$a}";
		$minP = $p{$a};
	    } else {
		if($p{$a} < $minP) {
		    $minS = "$a\t$h{$a}";
		    $minP = $p{$a};
		}
	    }
	} else {
	    if($minS eq "") {
		$minS = "$a\t-\t-\t-\t-\t-";
		$minP = 1;
	    }
	}
    }
    print "$minS\n";
}
close(F3);
