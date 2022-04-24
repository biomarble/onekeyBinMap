use strict;
use warnings;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","F=s","M=s","od=s","r=s","minBin=s","h" );
if(!defined($opts{i}) || !defined($opts{F}) || !defined($opts{M}) || defined($opts{h})){
                    &help();
                    exit;
}

my $vcf=$opts{i};
my $Fname=$opts{F};
my $Mname=$opts{M};
my $outdir=$opts{od};
my $r=$opts{r};
my $minBin=$opts{minBin};

$r||=0.02;
$outdir||="binResult";
$minBin||=10000;

mkdir $outdir  if(!-e $outdir);

$outdir=abs_path($outdir);
$vcf=abs_path($vcf);

mkdir "$outdir/1.bin" if(!-e "$outdir/1.bin");
mkdir "$outdir/2.plot" if(!-e "$outdir/2.plot");
mkdir "$outdir/3.final" if(!-e "$outdir/3.final");

my $cmd="perl $Bin/vcf_to_tsv.pl -i $vcf -F   $Fname -M  $Mname -od  $outdir/0.input -ref ";
system $cmd;

my @chr;
foreach my $file(glob("$outdir/0.input/geno.*.tsv")){
	if($file=~/geno\.(\w+)\.tsv/){
		push @chr,$1;
	}
}

open CMD1,">$outdir/s1.binning.sh" or die $!;
open CMD2,">$outdir/s2.visualize.sh" or die $!;
open CMD3,">$outdir/s3.convert.sh" or die $!;
my ($l1,$l2);
foreach my $chr(sort @chr){
	print CMD1 "snpbinner crosspoints -i $outdir/0.input/geno.$chr.tsv -o $outdir/1.bin/crosspoint.LG.$chr.csv -r $r && ";
	print CMD1 "snpbinner bins  -i $outdir/1.bin/crosspoint.LG.$chr.csv  -o  $outdir/1.bin/cbins.LG.$chr.csv  -l $minBin\n";
	print CMD2 "snpbinner visualize --out $outdir/2.plot/LG$chr --bins $outdir/1.bin/cbins.LG.$chr.csv --crosspoints $outdir/1.bin/crosspoint.LG.$chr.csv  --snps $outdir/0.input/geno.$chr.tsv \n";
	print CMD3 "csvtk transpose $outdir/1.bin/cbins.LG.$chr.csv |sed '2,\$ {s/^/C".$chr."P/}'|sed '1 {s/bin center/BinID/}' >$outdir/3.final/LG.$chr.binGeno.csv && ";
	print CMD3 "csvtk transpose -C '\$' $outdir/1.bin/cbins.LG.$chr.csv | csvtk cut -f 3,1,2 | sed 's/^/C".$chr."P/g'|sed 's/,/,$chr,/'|csvtk add-header -n 'BinID,Chr,Start,End' >$outdir/3.final/LG.$chr.binInfo.csv \n";
	$l1.=" $outdir/3.final/LG.$chr.binGeno.csv ";
	$l2.=" $outdir/3.final/LG.$chr.binInfo.csv ";

}
close CMD1;
close CMD2;
close CMD3;

open CMD4,">$outdir/s4.finalize.sh" or die $!;
print CMD4 "csvtk concat $l1 > $outdir/3.final/all.binGeno.csv \n";
print CMD4 "csvtk concat $l2 > $outdir/3.final/all.binInfo.csv \n";
close CMD4;

my $batchcmd="bash $outdir/s1.binning.sh && bash $outdir/s2.visualize.sh && bash $outdir/s3.convert.sh && bash $outdir/s4.finalize.sh";
open RUN,">$outdir/runSNPbinner.sh" or die $!;
print RUN "$batchcmd\n";
close RUN;

system $batchcmd;

sub help{
print <<"Usage.";
snp binning using SNPBinner using VCF file.

Usage:
-i                vcf file;
-F                F sample name in vcf
-M                M sample name in vcf
-od               Out dir, default MSToutdir

-minBin           minimum bin size, default 10000;

-r                Minimum distance between crosspoints as a ratio, default 0.02
                  0.02 means 2% of the chromosome length

-h                help document
Usage.
exit;
}

