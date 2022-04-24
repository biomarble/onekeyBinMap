#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","F=s","M=s","od=s","ref","h" );
if(!defined($opts{i}) || !defined($opts{F}) || !defined($opts{M}) || defined($opts{h})){
                    &help();
                    exit;
}

my $infile=$opts{i};
my $Fname =$opts{F};
my $Mname =$opts{M};
my $outdir=$opts{od};
my $ref=$opts{ref};
$outdir ="outdir" if(!defined $outdir);
mkdir $outdir if (!-e $outdir);

open (INFILE, "<$infile")  || die "Cannot open input VCF file!";

my $last_chromo = '';
my @header;

my $Find;
my $Mind;
while (<INFILE>){
    chomp;

    if ($_ =~ /^##/){
        next;
    }elsif ($_ =~ /^#CHROM/){
        my @h= split (/\t/, $_);
        splice @h, 0, 9;
        for(my $i=0;$i<@h;$i++){
            if($h[$i] eq $Fname ){
                $Find=$i;
            }elsif($h[$i] eq $Mname ){
                $Mind=$i;
            }else{
                push @header,$h[$i];
            }
        }
    }else{
            my ($chromo,$pos,undef,undef,undef,undef,undef,undef,undef,@rows) = split (/\t/, $_);
            my $marker_name = $chromo . "_" . $pos;
            $marker_name="C".$marker_name if($marker_name =~/^\d/);
            
            my @parentb = split (/:/, $rows[$Find]);
            my @parenta = split (/:/, $rows[$Mind]);
            
            #disregard SNPs with missing data or heterozygous in one of the parents.
            if ($parenta[0] eq "0/1" || $parenta[0] eq "./." || $parentb[0] eq "0/1" || $parentb[0] eq "./." || $parenta[0] eq $parentb[0]){ next;}
            
            if ($chromo eq $last_chromo){
                if($ref){
                    print OUTFILE $marker_name . "\t" . $pos;
                }else{
                    print OUTFILE $marker_name;
                }
                #    foreach my $row (@row){
                for(my $i=0;$i<@rows;$i++){
                    next if($i == $Find);
                    next if($i == $Mind);
                    my $row=$rows[$i];
                    my @geno = split (/:/, $row);
                    
                    if ($geno[0] eq $parenta[0]){
                        print OUTFILE "\ta";
                    }
                    
                    elsif ($geno[0] eq "0/1"){
                        print OUTFILE "\th";
                    }
                    
                    elsif ($geno[0] eq $parentb[0]){
                        print OUTFILE "\tb";
                    }
                    
                    else {print OUTFILE "\t-";
                    }
                }
                print OUTFILE "\n";
            }else{
                if($ref){
                    open (OUTFILE, ">$outdir/geno.$chromo.tsv") || die "Cannot open output file!";
                    print OUTFILE "marker_name\tposition\t";
                    print OUTFILE join("\t", @header);
                    print OUTFILE "\n";
                    print OUTFILE $marker_name . "\t" . $pos;
                }else{
                    if(!defined(fileno(OUTFILE))){
                        open (OUTFILE, ">$outdir/geno.tsv") || die "Cannot open output file!";
                        print OUTFILE "marker_name\t";
                        print OUTFILE join("\t", @header);
                        print OUTFILE "\n";
                        print OUTFILE $marker_name;
                    }
                }
                #foreach my $row (@row){
                for(my $i=0;$i<@rows;$i++){
                    next if($i == $Find);
                    next if($i == $Mind);
                    my $row=$rows[$i];
                    my @geno = split (/:/, $row);
                    if ($geno[0] eq $parenta[0]){
                        print OUTFILE "\ta";
                    }elsif ($geno[0] eq "0/1"){
                        print OUTFILE "\th";
                    }elsif($geno[0] eq $parentb[0]){
                        print OUTFILE "\tb";
                    }else{
                        print OUTFILE "\t-";
                    }
                }
                print OUTFILE "\n";
                $last_chromo = $chromo;
        }
    }
}

sub help{
print <<"Usage.";
convert VCF to genotype csv file for inbred populations.

Usage:
-i                vcf file;
-F                F sample name in vcf
-M                M sample name in vcf
-od               Out dir, default out


-h                help document
Usage.
exit;
}

