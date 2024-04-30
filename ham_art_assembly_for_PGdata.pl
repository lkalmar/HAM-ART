#!/usr/bin/perl

#Read in the config file


$progpath = $0;
if ($progpath =~ /\//){
	($path) = $progpath =~ /(.+\/)ham_art\.pl/;
}
open(C, $path."deps.txt");
while(<C>){
	if(!/^#/){
		chomp;
		$_ =~ s/\r$//;
		($key, $opt, $val) = split(' ', $_);
		$dep{$key}{'prog'} = $val;
		$dep{$key}{'opt'} = $opt;
	}
}
close(C);

$step = shift(@ARGV);

# Assembly steps

if ($step eq "assembly"){
	@prj_ids = parseIDs(shift(@ARGV));
	for $prj_id (@prj_ids){
		%prj = GetPrjDet($prj_id);
		# Start RSTRIM
		for $i (1..2){
		    $old = "pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta.notCombined_".$i.".fastq";
		    $new = "pre_assembly/".$prj{'PRJF'}."/LIB_PEL".$i."_".$prj{'PRJID'}.".fastq";
		    `mv $old $new`;
		}
		$old = "pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta.extendedFrags.fastq";
		$new = "pre_assembly/".$prj{'PRJF'}."/LIB_SEL_".$prj{'PRJID'}.".fastq";
		`mv $old $new`;

		open(SEL, ">>pre_assembly/".$prj{'PRJF'}."/LIB_SEL_".$prj{'PRJID'}.".fastq");
		open(EXC, ">pre_assembly/".$prj{'PRJF'}."/EXCL_".$prj{'PRJID'}.".fastq");

		@files = ('_hic.notCombined_1.fastq', '_hic.notCombined_2.fastq', '_hic.extendedFrags.fastq');

		for $fend (@files){
		    open(F, "pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}.$fend);
		    $fc = 0; #counter for fastq 4 line blocks
		    $cntr = $cnt1 = $cnt2 = $cntex = $cnt0 = 0;
		    while(<F>){
		        $fc++;
		        chomp;
		        if ($fc == 1){ #ID line
		            $id = $_;
		        }
		        if ($fc == 2){ #seq line
                #print $_."\n";
                $cntr++;
                @hits = $_ =~ /(GATCGATC|GATCAATT|AATTGATC|AATTAATT)/g;
                if (scalar(@hits) == 0){
                     if ($fend =~ /notCombined/){
                        $tosel = 1;
                        @seqs = ($_);
                        $cnt0++;
                    }
                    if ($fend =~ /extendedFrags/){
                        $exclude = 1;
                        @seqs = ($_);
                        $cntex++;
                    }
                }
                if (scalar(@hits) == 1){
                    $tosel = 2;
                    ($ligsite1, $ligsite2) = $hits[0] =~ /(....)(....)/;
                    @seqs = $_ =~ /(.*$ligsite1)($ligsite2.*)/;
                    $cnt1++;
                }
                if (scalar(@hits) == 2){
                    $tosel = 3;
                    ($ligsite1, $ligsite2) = $hits[0] =~ /(....)(....)/;
                    ($ligsite3, $ligsite4) = $hits[1] =~ /(....)(....)/;
                    @seqs = $_ =~ /(.*$ligsite1)($ligsite2.*$ligsite3)($ligsite4.*)/;
                    $cnt2++;
                }
            }
		        if ($fc == 4){
		            $qual = $_;
		            if ($exclude == 0){
		                if ($tosel > 0){
		                    $fraglen = 0;
		                    for ($i = 1; $i <= $tosel; $i++){
		                        if (length($seqs[$i-1]) > 30){
		                            print SEL $id."_".$i."\n";
		                            print SEL $seqs[$i-1]."\n";
		                            print SEL "+\n";
		                            print SEL (substr($qual, $fraglen, length($seqs[$i-1]))."\n");
		                        }
		                    }
		                }
		            }
		            if ($exclude == 1){
		                print EXC $id."\n";
		                print EXC $seqs[0]."\n";
		                print EXC "+\n";
		                print EXC $qual."\n";
		            }
		            $fc = $tosel = $topel = $exclude = 0;
		        }
		    }
		    close(F);
		    print STDERR "[RT]Processed ".$prj{'PRJID'}.$fend."\n";
		    print STDERR "[RT]Number of reads processed: ".$cntr."\n";
		    if ($fend =~ /notCombined/){
		        print STDERR "[RT]Reads with no detected restriction site: ".$cnt0."\n";
		    }
		    if ($fend =~ /extendedFrags/){
		        print STDERR "[RT]Excluded reads due to the lack of restriction site: ".$cntex++."\n";
		    }
		    print STDERR "[RT]Number of reads with 1 detected restriction site: ".$cnt1."\n";
		    print STDERR "[RT]Number of reads with 2 detected restriction site: ".$cnt2."\n";
		}

		close(SEL);
		close(EXC);

		for $i (1..2){
		    $old = "pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dehost.r_".$i.".fq";
		    $new = "pre_assembly/".$prj{'PRJF'}."/ARIBA_".$prj{'PRJID'}."_meta_dehost.r_".$i.".fq";
		    `mv $old $new`;
		}

		$oldfiles = "pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."*";
		`rm -f $oldfiles`;
		# Start MetaSpades
		$folder = "assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."/";
		$runopt = "--pe1-1 pre_assembly/".$prj{'PRJF'}."/LIB_PEL1_".$prj{'PRJID'}.".fastq --pe1-2 pre_assembly/".$prj{'PRJF'}."/LIB_PEL2_".$prj{'PRJID'}.".fastq --pe1-s pre_assembly/".$prj{'PRJF'}."/LIB_SEL_".$prj{'PRJID'}.".fastq -o ".$folder;
		($cpu, $mem) = split(/\,/, $dep{'METASPADES'}{'opt'});
		$runpr = $dep{'METASPADES'}{'prog'}." -k 21,33,55,77,99,121 -t ".$cpu." -m ".$mem;
	    `$runpr $runopt`;
	}
}

sub parseIDs{
	my $str = $_[0];
	my @arr = ();
	if ($str =~ /\,/){
		@arr = split(/\,/, $str);
	}
	else{
		@arr = ($str);
	}
	return @arr;
}
sub GetPrjDet {
    my $id = shift(@_);
    my $perm = '';
    if (scalar(@_) > 0){
        $perm = shift(@_);
    }
    my %retarr = ();
    my $prj_fn = '';
    my $pass = 0;
    open(F, 'project_defs/projects_db.txt');
    while(<F>){
        chomp;
        my @fl = split(/\t/, $_);
        if ($fl[0] eq $id){
            $prj_fn = $fl[2];
            my @fl2 = split(/;/, $fl[4]);
            $retarr{'ALL_STEPS'} = $fl[4];
            $retarr{'LAST_STEP'} = shift(@fl2);
            $retarr{'PRJ_FN'} = $fl[2];
            $pass = 1;
        }
    }
    close(F);
    if ($pass == 0){
        die "There is no project with this identifier in the database!\n";
    }
    open(F, $prj_fn);
    while(<F>){
        chomp;
        if (!/^#/){
            my ($k, $v) = split(' ', $_);
            $retarr{$k} = $v;
        }
    }
    close(F);
    return %retarr;
}

sub AppendLog {
    my ($id, $str) = @_;
    open(F, 'project_defs/projects_db.txt');
    open(FT, '>project_defs/projects_db.txt.temp');
    while(<F>){
        chomp;
        my @fl = split(/\t/, $_);
        print FT $_;
        if ($fl[0] eq $id){
            print FT $str.";";
        }
        print FT "\n";
    }
    close(F);
    close(FT);
    `mv project_defs/projects_db.txt.temp project_defs/projects_db.txt`;
}
sub checkPWD {
	if (!-d 'raw_seq_data' && !-d 'project_defs'){
		die "Please run the script by being in the root folder of the pipeline!\n";
	}
}

sub RevCompl {
	my $s = $_[0];
	$s =~ tr/ACGT/TGCA/;
	$r = reverse($s);
	return $r;
}

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

sub mean {
    my @arr = @_;
    my $sum = 0;
    for (@arr){
        $sum += $_;
    }
    return $sum / scalar(@arr);
}

sub stdev {
    my @arr = @_;
    my $avg = mean(@_);
    my $sqtot = 0;
    for (@arr){
        $sqtot += ($avg - $_) ** 2;
    }
    return (($sqtot / (scalar(@arr) - 1)) ** 0.5);
}

sub median{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2){
        return $vals[int($len/2)];
    }
    else{
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub GetLastTaxa{
	my $txstr = $_[0];
	my %taxdef = ("d" => "Domain", "p" => "Phylum", "c" => "Class", "o" => "Order", "f" => "Family", "g" => "Genus", "s" => "Species");
	my @f = split(';', $txstr);
	my @retarr = ();
	for my $t (@f){
		my ($tax, $str) = split('__', $t);
		if ($str ne ''){
			@retarr = ($taxdef{$tax}, $str);
		}
	}
	return @retarr;
}
