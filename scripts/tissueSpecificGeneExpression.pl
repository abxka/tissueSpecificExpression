#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 19.02.2017

###############################################################################
###############################################################################
### Determines tissue specific expression of genes.                         ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
use List::Util qw[min max sum];
use Statistics::R;

my (
    # variable for parameters which are read in from commandline
    $help,
    $metaFile,
    $metaCohortFile,
    $expressFile,
    $doMean,
    $outputAll,
    $exclude,
    $verbose,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"           => \$help,           # print this help
    "metaFile=s"      => \$metaFile,       # e.g. rnaseq_metadata.tsv.gz
    "cohortFile=s"    => \$metaCohortFile, # e.g. metaCohorts.tsv
    "expressFile=s"   => \$expressFile,    # e.g. freeze3_v2.kallisto.lib.trans.fpkm.wl.aliquot_id.tsv.gz
    "exclude:s"       => \$exclude,        # exclude these comma seperated tissues from analysis.
    "mean!"           => \$doMean,         # uses mean expression per tissue rather than median.
    "verbose!"        => \$verbose,        # print out additional information on calculation progress plus warning messages
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################
$| = 1;

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
    $sourceCode =~ s/.*\&GetOptions\(//s;
    $sourceCode =~ s/\).+//s;
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

###############################################################################
sub readMetaCohortFile {
    my %metaCohort = (); 
    open(F1, $metaCohortFile) or die "\nERROR: Failed to open $metaCohortFile: $!\n";
    while(my $l = <F1>) {
	next if($l =~ /^#/);
	chomp($l);
	my @a = split(/\t/, $l);
	$metaCohort{$a[0]} = $a[1];
    }
    close(F1);
    return \%metaCohort;
}
###############################################################################
sub readMetaDataFile {
    my $metaCohort = ();
    $metaCohort = &readMetaCohortFile() if(defined $metaCohortFile);

    my %sample2tissue = (); 
    open(F2, "zcat $metaFile|") or die "\nERROR: Failed to open $metaFile: $!\n";
    while(my $l = <F2>) {
	next if($l =~ /^#/);
	chomp($l);
	my @a = split(/\t/, $l);
	$sample2tissue{$a[1]} = $a[3];
	$sample2tissue{$a[1]} = $metaCohort->{$a[3]} if(exists $metaCohort->{$a[3]});
    }
    close(F2);
    return \%sample2tissue;
}
###############################################################################
sub readExpressionFile {
    my $sample2tissue = &readMetaDataFile();
    my %tissueN = ();
    my %express = ();
    open(F3, "zcat $expressFile|") or die "\nERROR: Failed to open $expressFile: $!\n";
    my $l = <F3>;
    chomp($l);
    my @samples = split(/\t/, $l);

    while(my $l = <F3>) {
	chomp($l);
	my @a = split(/\t/, $l);
	my $ensg = $a[0];

	for(my $i = 1; $i < @a; $i++) {
	    if(exists $sample2tissue->{$samples[$i]}) {
		my $tissue = $sample2tissue->{$samples[$i]};

		next if(defined $exclude and $exclude =~ /\b$tissue\b/i);

		my $expression = $a[$i];
		$expression = 0 if($expression eq "NA");

		if(exists $express{$ensg}->{$tissue}) {
		    push(@{$express{$ensg}->{$tissue}}, $expression);
		} else {
		    my @a = ($expression);
		    $express{$ensg}->{$tissue} = \@a;
		}
	    } else {
		print STDERR "WARNING: Sample $samples[$i] could not be found in $metaFile\n" if($verbose);
	    }
	}
    }
    close(F3);

    my %tissueAvgExpress = ();
    foreach my $ensg (sort keys %express) {
	foreach my $tissue (sort keys %{$express{$ensg}}) {
	    next if(defined $exclude and $exclude =~ /\b$tissue\b/i);

	    $tissueN{$tissue} = 1;

	    my $avgExpress;
	    if(defined $doMean) {
		$avgExpress = &mean(@{$express{$ensg}->{$tissue}});
	    } else {
		$avgExpress = &median(@{$express{$ensg}->{$tissue}});
	    }
	    $tissueAvgExpress{$ensg}->{$tissue} = $avgExpress;
	}
    }

    return (\%tissueAvgExpress, \%express, \%tissueN);
}
##############################################################################
sub median {
    my @values = sort{$a<=>$b} @_;

    my $med = 0;
    # if odd number of elements take middle element. Array size is given in size-1
    if(@values%2==1){
	$med = $values[int(@values/2)];
    }
    # if even number of elements take average of middle two elements
    else{
	my $n= @values;	
	$med = ($values[int(@values/2)-1] + $values[int(@values/2)]) / 2;
    }
    return $med;
}

##############################################################################
sub mean {
    my $avg = sum(@_)/@_;
    return $avg;
}
##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

if(!defined $metaFile) {
    print STDERR "\n\tPlease provide an META data file for the expression file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $expressFile) {
    print STDERR "\n\tPlease provide an expression table file. Try $0 -help to get more information\n\n";
    exit;
}

my ($tissueAvgExpress, $express, $tissueN) = &readExpressionFile();
my $R = Statistics::R->new();
my $tn = keys %$tissueN;

print "#Tissue1st\tTissue2nd\tTissueSpecificGene\tExpression1st\tExpression2nd\tEnrichment\tPvalue\n";
foreach my $ensg (sort keys %$tissueAvgExpress) {
    my $preTissueAvgExpress = -1;
    my $preTissue = "";
    # as avg tissue expression will be handled in sorted order, with the tissue that has the highest expression being first, 
    # we only need to test whether the 2nd tissue has significantly less expression, as only than we can talk about 
    # tissue specific expression.
    foreach my $tissue (sort {$tissueAvgExpress->{$ensg}->{$b} <=> $tissueAvgExpress->{$ensg}->{$a}} keys %{$tissueAvgExpress->{$ensg}}) {
	my $tissueAvgExpress = $tissueAvgExpress->{$ensg}->{$tissue};
	if($preTissueAvgExpress == -1) {
	    $preTissueAvgExpress = $tissueAvgExpress;
	    $preTissue = $tissue;
	    next;
	} else {
	    my $max = max(@{$express->{$ensg}->{$tissue}});
	    my $pvalue = -1;
	    if($max > 0) {
		$R->set('tissue', \@{$express->{$ensg}->{$tissue}});
		$R->set('preTissue', \@{$express->{$ensg}->{$preTissue}});
		$R->run(q`pvalue <- wilcox.test(tissue, preTissue)$p.value`);
		$pvalue = $R->get('pvalue');
	    }
	    my $enrichment = $preTissueAvgExpress;
	    $enrichment /= $tissueAvgExpress if($tissueAvgExpress > 0);
	    printf("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4e\n", $preTissue, $tissue, $ensg, $preTissueAvgExpress, $tissueAvgExpress, $enrichment, $pvalue);
	    last;
	}
    }
}
$R->stop();
