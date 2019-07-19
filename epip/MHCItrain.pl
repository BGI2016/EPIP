#!/usr/bin/perl
use strict;

die "Training MHCI Allele PSSM...\t
\tperl $0 <database> <species>  <2PSSMBG> <outdir> <MinNum>\n" if(@ARGV < 2);
my ($database, $species, $pssm2bg, $outdir, $MinNum) = @ARGV;
$pssm2bg ||= "NA";
$outdir ||= "./";
$MinNum ||= 10;

open PSSMFILE,">> $outdir/pssm_file.lst" or die $!;
 
my %db;
open(FH, $database) || die $!;
while(<FH>)
{
	chomp;
	my @tmp = split /\s+/;
	next if($tmp[1] ne $species);
	push @{$db{$tmp[3]}{$tmp[6]}}, $tmp[5];
}
close FH;

for my $allele (keys %db)
{
	for my $len (keys %{$db{$allele}})
	{
		my @seq = @{$db{$allele}{$len}};
		next if(($#seq + 1) < 10);

		open(OUT, ">$outdir/train.seq") || die $!;
		print OUT ">$allele $len\n";
		foreach my $x (@seq)
		{
			print OUT "$x\n";
		}
		close OUT;

		if($pssm2bg eq "NA")
		{
            open FH2, "> $outdir/$allele\_$len.pssm" or die $!;
            my ($input, $IC50) = "$outdir/train.seq", "0";

            my @lib = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');
            
            # Define background frequency in uniprot database
            my %bg;
            $bg{'A'} = 0.0826622168195282;
            $bg{'C'} = 0.0136850332007961;
            $bg{'D'} = 0.0545883107520953;
            $bg{'E'} = 0.0675987909234391;
            $bg{'F'} = 0.0386658831073075;
            $bg{'G'} = 0.0708488499195514;
            $bg{'H'} = 0.0227226416223757;
            $bg{'I'} = 0.0597571006180889;
            $bg{'K'} = 0.0585272272907145;
            $bg{'L'} = 0.096681402591599;
            $bg{'M'} = 0.0242260080622711;
            $bg{'N'} = 0.0406364484957737;
            $bg{'P'} = 0.0469865264985876;
            $bg{'Q'} = 0.0393742846376247;
            $bg{'R'} = 0.0553599274523758;
            $bg{'S'} = 0.0654447212017059;
            $bg{'T'} = 0.0533960101143632;
            $bg{'V'} = 0.0687266172099021;
            $bg{'W'} = 0.0108407347094559;
            $bg{'Y'} = 0.0292109734641572;
            
            #Set Dirchlet mix distrubution to avoid overfit
            my ($dirchlet, $dirchlet_all, $dirchlet_ic50) = (1, 20, 1000);
            
            # Read peptides from input file
            my @seq;
            my @ic50;
            local $/ = "\n>";
            open(FH, $input) || die $!;
            while(<FH>)
            {
            	chomp;
            	s/^>//;
            	my @tmp = split /\n/, $_;
            	# Get MHC type and peptides lenngth.
            	my ($allel, $len) = split /\s+/, $tmp[0];
            	my $num = $#tmp;
            	for my $i (1..$num)
            	{
            		($seq[$i-1], $ic50[$i-1]) = split /\s+/, $tmp[$i];
            		$seq[$i-1] = uc($seq[$i-1]);
            		die "The $i peptides has a wrong length!\n" if(length($seq[$i-1]) != $len);
            	}
            	# Build PSSM
            	my %pssm;
            	for my $j (0..$len-1)
            	{
            		foreach my $char (@lib)
            		{
            			$pssm{$char}[$j] = 0;
            		}
            		my $all = 0;
            		for my $i (0..$num-1)
            		{
            			my $char = substr($seq[$i], $j, 1);
            			my $frq = 1;
            			$frq = $frq / $ic50[$i] if($IC50 == 1);
            			$all += $frq;
            			$pssm{$char}[$j] += $frq;
            		}
            
            		foreach my $char (@lib)
            		{
            			my $tmp = log(($pssm{$char}[$j]+$dirchlet)/($all+$dirchlet_all)/$bg{$char});
            			$tmp = log(($pssm{$char}[$j]+$dirchlet/$dirchlet_ic50)/($all+$dirchlet_all/$dirchlet_ic50)/$bg{$char}) if($IC50 == 1);
            			$pssm{$char}[$j] = $tmp;
            		}
            	}
            	# OUtput 
            	print FH2 ">$tmp[0]\n";
            	foreach my $char (@lib)
            	{
            		for my $j(0..$len-1)
            		{
            			print FH2 "$pssm{$char}[$j]\t";
            		}
            		print FH2 "\n";
            	}
            }
            close FH;
            local $/ = "\n";

          
		}

		`rm $outdir/train.seq`;
                print PSSMFILE "PSSM/PSSMx/$allele\_$len.pssm\t$allele\t$len\n"
	}
}
