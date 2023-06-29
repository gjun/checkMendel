#!/usr/bin/perl -w

# Check Mendelian Errors in VCF genotypes using PED pedigree files

use strict;
use Getopt::Long;

my $vcf = "";
my $fam = "";
my $chr = "";
my $out = "mendelErr";
my %trios = ();
my $filter	= 0;

my $cmd = GetOptions("vcf=s"=>\$vcf, "ped=s"=>\$fam, "out=s"=>\$out, "chr=s" => \$chr);

die "Usage: perl check_mendel.pl --vcf [VCF] --ped [PED] --out [output prefix] \n" unless ($cmd && $vcf && $fam );

if ($vcf =~ m/\.gz$/)
{
	die "Error: Cannot find $vcf.tbi\n" unless ($chr eq "" || -e "$vcf.tbi");
	open(IN, "bcftools view -h $vcf|grep ^#CHROM|") || die "Cannot open $vcf\n";
}
elsif ($vcf =~ m/\.vcf$/)
{
	die "Error: --chr option is available only with tabixed vcf.gz\n" unless ($chr eq "");
	open(IN, "bcftools view -h $vcf|grep ^#CHROM|") || die "Cannot open $vcf\n";
}
else
{
	die "$vcf should be in .vcf or .gz format\n";
}

my %vcf_ids = ();
my $n_sample = 0;

while(<IN>)
{
	chomp;
	my ($c, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @ids) = split;

	$n_sample = scalar(@ids);

	for(my $i=0;$i<$n_sample;++$i)
	{
		$vcf_ids{$ids[$i]} = $i; # column number in VCF 
	}
}
close IN;

my %pats = ();
my %mats = ();

open(IN, "$fam") || die "Cannot open $fam\n";
while(<IN>)
{
	chomp;
	my ($fam, $id, $pat, $mat, $sex) = split;
	if (defined($vcf_ids{$id}) && (defined($vcf_ids{$pat}) || defined($vcf_ids{$mat})))
	{
		$trios{$id} = $vcf_ids{$id};
		if (defined($vcf_ids{$pat}))
		{
			$pats{$id} = $vcf_ids{$pat};
		}
		else
		{
			$pats{$id} = -1;
		}

		if (defined($vcf_ids{$mat}))
		{
			$mats{$id} = $vcf_ids{$mat};
		}
		else
		{
			$mats{$id} = -1;
		}
	}
}
close IN;


my %genos = ();

my @mendelErrs = (0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,1,1,0);   
my @genoCnts = ();

$genos{"./."} = 0;
$genos{"./0"} = 0;
$genos{"0/."} = 0;
$genos{"./1"} = 0;
$genos{"1/."} = 0;
$genos{".|."} = 0;
$genos{".|0"} = 0;
$genos{"0|."} = 0;
$genos{".|1"} = 0;
$genos{"1|."} = 0;

$genos{"0/0"} = 1;
$genos{"0|0"} = 1;

$genos{"0/1"} = 2;
$genos{"0|1"} = 2;
$genos{"1/0"} = 2;
$genos{"1|0"} = 2;

$genos{"1/1"} = 3;
$genos{"1|1"} = 3;

$genos{"."} = 0;
$genos{"0"} = 1;
$genos{"1"} = 2;


my %errormap = ();

foreach my $p (sort keys %genos)
{
	foreach my $m (sort keys %genos)
	{
		foreach my $s (sort keys %genos)
		{
			$errormap{"$p:$m:$s"} = $mendelErrs[$genos{$p}*16 + $genos{$m}*4 + $genos{$s}];
			$genoCnts[$genos{$p}*16 + $genos{$m}*4 + $genos{$s}] = 0;
		}
	}

}


if ($vcf =~ m/\.gz$/)
{
	if ($chr eq "")
	{
		open(IN, "bcftools view -H $vcf|") || die "Cannot open $vcf\n";
	}
	else
	{
		open(IN, "tabix $vcf $chr:0|") || die "Cannot open $vcf\n";
	}
}
elsif ($vcf =~ m/\.vcf$/)
{
	open(IN, "bcftools view -H $vcf|") || die "Cannot open $vcf\n";
}
else
{
	die "$vcf should be in .vcf or .gz format\n";
}

my %ac_cnt  = ();
my %n_var = ();

open(OUT, ">$out.var");
while(<IN>)
{
	my ($c, $pos, $id, $ref, $alt, $qual, $filt, $info, $format, @G) = split;
	my @g = ();
	my $bParseGT = 0;
	my $GTidx = -1;

	unless ($ref =~ /,/ || $alt =~ /,/)
	{
		#if ($format ne "GT")
		#{
			my (@fields) = split(/:/,$format);
			for(my $j=0;$j<@fields;++$j)
			{
				if ($fields[$j] eq "GT")
				{
					$GTidx = $j;
				}
			}
			if ($GTidx>=0)
			{
				for(my $j=0;$j<@G;++$j)
				{
					my (@gtfield) = split(/:/,$G[$j]);
					$g[$j] = $gtfield[$GTidx];
				}
			}
		#}
		unless ($GTidx <0)
		{

			die "Error, $c:$pos has ".scalar(@g)." GT fields, while header has $n_sample IDs\n" unless (scalar(@g) == $n_sample);

			my $n_error = 0;
			my $ac = 0;
			$ac = $1 if ($info =~ /AC=(\d+)/);

			my $err_trio = "";

			foreach my $k (sort keys %trios)
			{
				my ($p, $m, $s);

				$p = ($pats{$k}<0) ? "./." : $g[$pats{$k}];
				$m = ($mats{$k}<0) ? "./." : $g[$mats{$k}];
				$s = $g[$trios{$k}];
				if (!defined($errormap{"$p:$m:$s"}))
				{
					print "$p:$m:$s is not defined\n"
				}
				else
				{
					my $E = $errormap{"$p:$m:$s"};
					$n_error += $E;
					if ($E > 0 )
					{
						$err_trio .= "$G[$pats{$k}]___$G[$mats{$k}]___$G[$trios{$k}]\t";
					}
					$genoCnts[$genos{$p}*16 + $genos{$m}*4 + $genos{$s}] ++;
				}
			}
			print OUT "$c\t$pos\t$id\t$ref\t$alt\t$ac\t$n_error\t$info\t$err_trio\n";

			$ac_cnt{$ac} = (defined($ac_cnt{$ac})) ? $ac_cnt{$ac}+$n_error : $n_error;
			$n_var{$ac} = (defined($n_var{$ac})) ? $n_var{$ac}+1 : 1;
		}
	}
}
close IN;
close OUT;

open(OUT, ">$out.AC");

print OUT "#AC\tN_VAR\tMENDEL_ERR\tERR_RATE\n";
foreach my $ac (sort { $a <=> $b } keys %ac_cnt)
{
	print OUT "$ac\t$n_var{$ac}\t$ac_cnt{$ac}\t".($ac_cnt{$ac}/($n_var{$ac}*$ac))."\n";
}
close OUT;

open(OUT, ">$out.table");


print OUT "#FAT\tMAT\t./.\t0/0\t0/1\t1/1\n";
my @L = ("./.","0/0","0/1","1/1");
my $genosum = 0;
my $errsum = 0;

for(my $i=0;$i<4;++$i)
{
	for(my $j=0;$j<4;++$j)
	{
		print OUT "$L[$i]\t$L[$j]";
		for(my $k=0;$k<4;++$k)
		{
            my $idx = $i*16+$j*4+$k;
			print OUT "\t$genoCnts[$idx]";

            if (($i+$j)>0 && $k>0)
            {
                $genosum += $genoCnts[$idx];
                $errsum += $genoCnts[$idx] * $mendelErrs[$idx];
            }
		}
		print OUT "\n";
	}
}
$genosum -= $genoCnts[5] + $genoCnts[17] + $genoCnts[21];

print OUT "Non-ref error rate: ".($errsum/$genosum)."\n";
close OUT;
