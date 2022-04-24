use strict;
use warnings;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use Getopt::Long;
my %opts;
GetOptions(\%opts,"g=s","info=s","pop=s","od=s","method=s","miss=s","pcut=s","h" );
if(!defined($opts{g}) ||!defined($opts{info}) || !defined($opts{pop}) || defined($opts{h})){
                    &help();
                    exit;
}

my $key=time();
my $infile=$opts{g};
my $info=$opts{info};
my $pop=$opts{pop};
my $outdir=$opts{od};
my $method=$opts{method};
my $maxMis=$opts{miss};
my $pcut=$opts{pcut};

$pcut=0.001 if(!defined $pcut);
$maxMis=0.5 if(!defined $maxMis);

open IN,$info or die $!;
<IN>;
my %info;
while(<IN>){
	chomp;
	my ($bid,$chr,$s,$e)=split /,/,$_;
	$info{$chr}{$bid}{'s'}=$s;
	$info{$chr}{$bid}{'e'}=$e;
}
close IN;


$pop="RIL2" if($pop eq "F2");
if($pop =~/RIL(\d*)/){
	&help if(!defined $1 or $1<2);
	$pop="RIL10" if($1>10);
}

$method||="ML";
$outdir||="MSToutdir";

mkdir $outdir if(!-e $outdir);
$outdir=abs_path($outdir);

open IN,$infile or die $!;
my $header=<IN>;
chomp $header;
my (undef,@sid)=split /,/,$header;
my $nind=scalar @sid;
my %geno;
my $allNum=0;
my $mis=0;
open OUT,">$outdir/$key.tmp.stat" or die $!;
while(<IN>){
	chomp;
	next if($_=~/^\s*$/);
	next if($_=~/^#+$/);
	my ($id,$line)=split /,/,$_,2;
	$allNum++;
	my $a= ($line =~ s/a/a/g);
	my $b= ($line =~ s/b/b/g);
	my $h= ($line =~ s/h/h/g);
	$a||=0;
	$b||=0;
	$h||=0;
	
	my $miss;
	if($pop =~/^RIL2$/){
		$miss=$nind-$a-$b-$h;
	}elsif($pop=~/RIL/){
		$miss=$nind-$a-$b;
	}else{
		die "pop not valid\n";
	}
	if($miss/$nind >$maxMis){
		$mis++;
	}
	print OUT "$id\t$a\t$b\t$h\t$miss\n";
	$line=~s/,/\t/g;
	$line=~s/\bh\b/X/g;
	$geno{$id}=$line;
}
close IN;
close OUT;

system "Rscript $Bin/filter.R  $outdir/$key.tmp.stat $outdir/filter $pcut $maxMis $pop ";

my %valid;
open IN,"$outdir/filter.valid" or die $!;
while(<IN>){
	chomp;
	$valid{$_}=1;
}
close IN;
`rm $outdir/$key.tmp.stat `;


foreach my $chr(sort keys %info){
	my @valMarker;
	foreach my $mid(sort keys %{$info{$chr}}){
		next if(!exists $valid{$mid});
		push @valMarker,$mid;
	}
	my $nloc=scalar @valMarker;
my $out=<<HEAD;
population_type $pop
population_name LG
distance_function kosambi
cut_off_p_value 2.0
no_map_dist 15.0
no_map_size 0
missing_threshold 1
estimation_before_clustering no
detect_bad_data yes
objective_function $method
number_of_loci $nloc
number_of_individual $nind
HEAD
	open OUT,">$outdir/LG.$chr.MSTmapIn.txt";
	print OUT "$out\nlocus_name\t";
	print OUT join("\t",@sid);
	print OUT "\n";
	foreach my $mid(@valMarker){
		print OUT "$mid\t$geno{$mid}\n";
	}
	close OUT;
	my $cmd="mstmap  $outdir/LG.$chr.MSTmapIn.txt    $outdir/LG.$chr.MSTmapOut.txt  >$outdir/LG.$chr.MSTmapOut.log\n";
	system $cmd;
}

open OUT,">$outdir/MST.map.txt" or die $!;
foreach my $chr(sort keys %info){
open IN,"$outdir/LG.$chr.MSTmapOut.txt" or die $!;
while(<IN>){
	chomp;
	next if($_=~/;/);
	next if($_=~/^\s*$/);
	next if($_=~/^group\s+/);
	my ($mid,$pos)=split /\s+/,$_;
	print OUT "$mid\t$chr\t$pos\n";
}
close IN;
}
close OUT;

my $nused=scalar keys %valid;
print STDERR "TotalMarker: $allNum\nMissing > ".($maxMis*100)."%: $mis\nSegregation Disorder Markers:".($allNum-$nused-$mis)." \nUsed Markers: $nused\n";

sub help{
print <<"Usage.";

MSTmap linkage map construction using SNPBinner output.

Usage:
-i                bin genotype file;
-info             bin information file, containing BinID CHR START END;

-pop              population type:  DH,BC1,F2,RIL6
                  RIL generation number should between 3~10;

-method           map construction objective function : ML(default), COUNT
                  ML:    mamaximum likelihood
                  COUNT: sum of recombination events

-od               Out dir, default MSToutdir
-h                help document
Usage.
exit;
}
