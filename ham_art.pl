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
		($key, $opt, $val) = split(/\t/, $_);
		$dep{$key}{'prog'} = $val;
		$dep{$key}{'opt'} = $opt;
	}
}
close(C);

$step = shift(@ARGV);

# project initialisation
# uses only 1 argument, the project description file name or multiple filenames separated by comma(s)
if ($step eq 'initialise'){ 
	@fl = split(/\//, $ENV{'PWD'});
	if (pop(@fl) ne 'project_defs'){
	    die "Please run this script from the project_defs folder!\n";
	}
	
	@fns = parseIDs(shift(@ARGV));
	for $fn (@fns){
	    open(F, "projects_db.txt");
	    while(<F>){
	        chomp;
	        ($id, undef) = split(' ', $_);
	        $old_ids{$id}=1;
	    }
	    close(F);
	    open(F, $fn);
	    while(<F>){
	        chomp;
	        if(!/^#/){
	            ($k, $v) = split(' ', $_);
	            $param{$k} = $v;
	        }
	    }
	    close(F);
	    @mandatory_params = ("PRJID", "PRJF", "PRJH", "METF1", "METF2", "HICF1", "HICF2");
	    for $p (@mandatory_params){
	    	if (!exists($param{$p})){
			die "Mandatory parameter ".$p." is missing from the project description file!\n";
		}
	    }
	    #Project id check
	    if (exists($old_ids{$param{'PRJID'}})){
	        die "The project ID (".$param{'PRJID'}.") has already been used before or the project is already initiated, please use a different ID\n";
	    }
	    else{
	        print STDERR "Project ID: ".$param{'PRJID'}." --> VALID\n";
	    }
	    #Project folder check
	    @folders = ('../assembly', '../post_assembly', '../pre_assembly', '../results');
	    if (-d "../raw_seq_data/".$param{'PRJF'}){
	        print STDERR "Project folder: ".$param{'PRJF'}." --> VALID\n";
	        for $f (@folders){
	            if (-d $f."/".$param{'PRJF'}."/"){
	                print STDERR "Folder ".$f."/".$param{'PRJF'}."/ was already created\n";
	            }
	            else{
	                $fold = $f."/".$param{'PRJF'}."/";
	                `mkdir $fold`;
	                print STDERR "Folder ".$f."/".$param{'PRJF'}."/ is created\n";
	            }
	        }
	    }
	    else{
	        die "Project folder ".$param{'PRJF'}." does not exist in raw_seq_data folder, please check the PRJF parameter in the project definition file\n";
	    }

	    #Input files check
	    if (-f "../raw_seq_data/".$param{'PRJF'}."/".$param{'METF1'}){
	        print STDERR "Metagenomics file 1 (raw_seq_data/".$param{'PRJF'}."/".$param{'METF1'}.") was found --> VALID\n";
	    }
	    else{
	        die "Metagenomics file 1 (raw_seq_data/".$param{'PRJF'}."/".$param{'METF1'}.") does not exist, please check the METF1 parameter in the project definition file\n";
	    }
	    if (-f "../raw_seq_data/".$param{'PRJF'}."/".$param{'METF2'}){
	        print STDERR "Metagenomics file 2 (raw_seq_data/".$param{'PRJF'}."/".$param{'METF2'}.") was found --> VALID\n";
	    }
	    else{
	        die "Metagenomics file 2 (raw_seq_data/".$param{'PRJF'}."/".$param{'METF2'}.") does not exist, please check the METF2 parameter in the project definition file\n";
	    }
	    #if ($param{'HICPRE'} ne 'NONE'){
	        if (-f "../raw_seq_data/".$param{'PRJF'}."/".$param{'HICF1'}){
	            print STDERR "Hi-C file 1 (raw_seq_data/".$param{'PRJF'}."/".$param{'HICF1'}.") was found --> VALID\n";
	        }
	        else{
	            die "Hi-C file 1 (raw_seq_data/".$param{'PRJF'}."/".$param{'HICF1'}.") does not exist, please check the HICPRE parameter in the project definition file\n";
	        }
	        if (-f "../raw_seq_data/".$param{'PRJF'}."/".$param{'HICF2'}){
	            print STDERR "Hi-C file 2 (raw_seq_data/".$param{'PRJF'}."/".$param{'HICF2'}.") was found --> VALID\n";
	        }
	        else{
	            die "Hi-C file 2 (raw_seq_data/".$param{'PRJF'}."/".$param{'HICF2'}.") does not exist, please check the HICPRE parameter in the project definition file\n";
	        }
	    #}

	    print STDERR "\n########################################################################################################\n";
	    print STDERR "#                                                                                                      #\n";
	    print STDERR "# Everything is VALID within the project file, you are now allowed to initiate the project.            #\n";
	    print STDERR "# Once you initiate it, you will be able to go through the steps by using only the project identifier. #\n";
	    print STDERR "#                                                                                                      #\n";
	    print STDERR "########################################################################################################\n";
	    print STDERR "\nWould you like to initiate it now (y/n)?";
	    $resp = <STDIN>;
	    chomp($resp);
	    if ($resp eq 'y'){
	        print STDERR "Initiating...\n";
	        open(OUT, ">>projects_db.txt");
	        (undef,undef,undef,$mday,$mon,$year,undef,undef,undef) = localtime();
	        $datestr = $mday."_".($mon+1)."_".($year+1900);
	        $md5s = `md5sum $fn`;
	        chomp($md5s);
	        ($md5, undef) = split(' ', $md5s);
	        print OUT $param{'PRJID'}."\t".$datestr."\tproject_defs/".$fn."\t".$md5."\tinit;\n";
	        close(OUT);
	        print STDERR "\n########################################################################################################\n";
	        print STDERR "#                                                                                                      #\n";
	        print STDERR "# Project ".$param{'PRJID'}." initialised successfully, please use your project ID in all further steps #\n";
	        print STDERR "#                                                                                                      #\n";
	        print STDERR "########################################################################################################\n";
	    }
	    else{
	        print STDERR "You chose not to initiate the project yet, you can do it anytime by running the same command as you just used\n";
	    }
	}
}

# Pre assembly steps

if ($step eq 'pre_assembly'){
	@prj_ids = parseIDs(shift(@ARGV));
	for $prj_id (@prj_ids){
	    #Get the project details from the project description file
	    %prj = GetPrjDet($prj_id);
	    #Check for the appropriate working directory
	    checkPWD();

	    #Run de-duplication on metagenomics dataset
	    print STDERR "Running de-duplication on the metagenomics dataset\n";
	    if (-f "raw_seq_data/".$prj{'PRJF'}."/".$prj{'METF1'} && -f "raw_seq_data/".$prj{'PRJF'}."/".$prj{'METF2'}){
	        $runopt = "in=raw_seq_data/".$prj{'PRJF'}."/".$prj{'METF1'}." in2=raw_seq_data/".$prj{'PRJF'}."/".$prj{'METF2'}." out=pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dedup.r_1.fq out2=pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dedup.r_2.fq";
	        $runpr = $dep{'CLUMPIFY'}{'prog'}." -Xmx".$dep{'CLUMPIFY'}{'opt'}."g dedupe=t";
	        `$runpr $runopt`;
	    }
	    else{
	        die "Raw data fastq files cannot be found!\n";
	    }

	    #Run de-duplication for Hi-C dataset if there is a Hi-C dataset associated to the project
	    print STDERR "Running de-duplication on the Hi-C dataset\n";
        if (-f "raw_seq_data/".$prj{'PRJF'}."/".$prj{'HICF1'} && -f "raw_seq_data/".$prj{'PRJF'}."/".$prj{'HICF2'}){
            $runopt = "in=raw_seq_data/".$prj{'PRJF'}."/".$prj{'HICF1'}." in2=raw_seq_data/".$prj{'PRJF'}."/".$prj{'HICF2'}." out=pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_hic_dedup.r_1.fq out2=pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_hic_dedup.r_2.fq";
            $runpr = $dep{'CLUMPIFY'}{'prog'}." -Xmx".$dep{'CLUMPIFY'}{'opt'}."g dedupe=t";
	        `$runpr $runopt`;
        }
        else{
            die "Hi-C raw data fastq files cannot be found!\n";
        }

	    #Run de-hosting step for metagenomics data
	if ($prj{'PRJH'} eq 'NOHOST'){
		print STDERR "Skipping de-hosting step.\n";
		$runopt = "pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dedup.r_1.fq pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dehost.r_1.fq";
		`mv $runopt`;
		$runopt	= "pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dedup.r_2.fq pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dehost.r_2.fq";
                `mv $runopt`;
	}
	else{    
		print STDERR "Running de-hosting on the metagenomics dataset\n";
	    @hostlist = `ls host_ref_seqs/ | cut -f 1 -d "_" | sort | uniq`;
	    for $h (@hostlist){
	        chomp($h);
	        $vh{$h} = 1;
	        push(@h_arr, $h);
	    }
	    if (!exists($vh{$prj{'PRJH'}})){
	        die "No built-in database for this host, please use one of these or install the new database: ".join(',', @h_arr)."\n";
	    }
	    $runopt = "-x host_ref_seqs/".$prj{'PRJH'}."_bt2 -1 pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dedup.r_1.fq -2 pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dedup.r_2.fq --un-conc pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dehost.r_%.fq > /dev/null";
	    $runpr = $dep{'BOWTIE2'}{'prog'}." -p ".$dep{'BOWTIE2'}{'opt'}." --no-mixed --fast";
	    `$runpr $runopt`;
	}
	    #Run FLASH read merging on the metagenomics dataset
	    print STDERR "Megring overlapping reads in the metagenomics dataset\n";
	    $runopt = "-o ".$prj{'PRJID'}."_meta -d pre_assembly/".$prj{'PRJF'}."/ pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dehost.r_1.fq pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_meta_dehost.r_2.fq";
	    $runpr = $dep{'FLASH'}{'prog'}." -t ".$dep{'FLASH'}{'opt'}." -m 30";
	    `$runpr $runopt`;

	    #Run FLASH read merging for Hi-C dataset 
	    print STDERR "Megring overlapping reads in the Hi-C dataset\n";
        $runopt = "-o ".$prj{'PRJID'}."_hic -d pre_assembly/".$prj{'PRJF'}."/ pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_hic_dedup.r_1.fq pre_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_hic_dedup.r_2.fq";
        $runpr = $dep{'FLASH'}{'prog'}." -t ".$dep{'FLASH'}{'opt'}." -m 30";
	    `$runpr $runopt`;
	}
}

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
		            @hits = $_ =~ /(ACGCGT)/g;
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
		                @seqs = $_ =~ /(.*ACG)(CGT.*)/;
		                $cnt1++;
		            }
		            if (scalar(@hits) == 2){
		                $tosel = 3;
		                @seqs = $_ =~ /(.*ACG)(CGT.*ACG)(CGT.*)/;
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

# Post assembly steps

if ($step eq "post_assembly"){
	@prj_ids = parseIDs(shift(@ARGV));
	for $prj_id (@prj_ids){
		%prj = GetPrjDet($prj_id);
		#Re-alingment 1, copy contigs fasta and make bowtieDB
	    if (-f "assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."/contigs.fasta"){
	        $old = "assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."/contigs.fasta";
	        $new = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contigs.fasta";
	        `mv $old $new`;
	        $runopt = "-f post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contigs.fasta post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contigs";
	        $runpr = $dep{'BOWTIE2-B'}{'prog'}." -q";
	        `$runpr $runopt`;
	        
	    }
	    else{
	        if (-f "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contigs.fasta"){
	            $runopt = "-f post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contigs.fasta post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contigs";
	            $runpr = $dep{'BOWTIE2-B'}{'prog'}." -q";
	        	`$runpr $runopt`;
	        }
	        else{
	            die "This step has to be run after the metagenomic assembly, but the final fasta file wasn't found!\n";
	        }
	    }
	    #Re-alignment 2, doing the Bowtie2 alignment
	    $runopt = "-x post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contigs -U pre_assembly/".$prj{'PRJF'}."/LIB_SEL_".$prj{'PRJID'}.".fastq -S post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_SEL.sam";
	    $runpr = $dep{'BOWTIE2'}{'prog'}." --fast -p ".$dep{'BOWTIE2'}{'opt'};
	    `$runpr $runopt`;
	    #Extract contigs with amr genes
	    $infile = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contigs.fasta";
		$outfile = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_amr_contigs.fasta";
		$logfile = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_amr_genes.txt";
		print STDERR "Reading in contigs_fasta file for project ".$prj{'PRJF'}."...\n";
		open(F, $infile);
		while(<F>){
			chomp;
			if (/^>/){
				$id = $_;
				$id =~ s/^>//;
			}
			else{
				$seq{$id} .= $_;
			}
		}
		close(F);
		print STDERR "Performing ABRICATE search...\n";
		$runpr = $dep{'ABRICATE'}{'prog'}." --mincov 50 --threads ".$dep{'ABRICATE'}{'opt'}." --db resfinder $infile"; 
		@res = `$runpr`;
		open(OUT, ">".$logfile);
		for $line (@res){
			chomp($line);
			@fl = split(' ', $line);
			$amr{$fl[1]} = 1;
			print OUT $line."\n";
		}
		close(OUT);

		print STDERR "Writing output file...\n";
		open (OUT, ">".$outfile);
		for $id (keys %amr){
			print OUT ">".$id."\n".$seq{$id}."\n";
		}
		close(OUT);
		#Sorint SAM file
		$runopt = "-o post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_SEL_sorted.sam post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_SEL.sam";
		$runpr = $dep{'SAMTOOLS'}{'prog'}." sort -n -@ ".$dep{'SAMTOOLS'}{'opt'}." -O sam";
		`$runpr $runopt`;
		#Extracting contacts
		open(F, "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_SEL_sorted.sam");
		open(O, "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_hic_contacts");
		$preid = '';
		$pass = 0;
		$cnt = 0;
		while (<F>){
			$cnt++;
			chomp;
			if (!/^@/){
				@fl = split(/\t/, $_);
				if ($fl[1] != 4){
					$id = $fl[0];
					$seq = $fl[9];
					if ($fl[1] == 16 || $fl[1] == 272){
						$seq = RevCompl($fl[9]);
					}
					if ($id ne $preid && $pass == 1){
						$hitcnt++;
						if (scalar(keys(%hits)) <= 1){
							$uniqhit++;
						}
						else{
							@targets = ();
							for $seq (keys %hits){
								if (scalar(@{$hits{$seq}} == 1)){
									$multim = 0;
									@fl2 = split(/\t/, $hits{$seq}[0]);
									$pass2 = 0;
									if ($fl2[5] =~ /M/){
										$pass2 = 1;
										push(@targets, $fl2[2]);
									}
								}
								else{
									$pass2 = 0;
									$multim = 0;
									for ($i=0;$i<scalar(@{$hits{$seq}});$i++){
										@fl2 = split(/\t/, $hits{$seq}[$i]);
										if ($fl2[5] =~ /M/){
											$multim++;
											$tmptarget = $fl2[2];
										}
									}
		                            if ($multim == 1){
		                                push(@targets, $tmptarget);
		                                $pass2 = 1;
		                            }
								}
							}
							@utargets = uniq(@targets);
							if (scalar(@utargets) > 1){
								$valid++;
		                        for($t1=0;$t1<scalar(@utargets)-1;$t1++){
		                            for($t2=$t1+1;$t2<scalar(@utargets);$t2++){
		                                print O $utargets[$t1]."\t".$utargets[$t2]."\n";
		                            }
		                        }
							}
		                    if (scalar(@utargets) == 1){
		                        $close++;
		                    }
		                    if (scalar(@utargets) == 0){
		                        $excluded++;
		                    }
						}
						%hits = ();
					}
					$pass = 1;
					$preid = $id;
					push(@{$hits{$seq}}, $_);
				}
				else{
					$nonmapped++;
				}
			}
		}
		close(F);
		close(O);
		print STDERR "\n";
		print STDERR "Number of read groups: ".$hitcnt."\n";
		print STDERR "Non-mapped: ".$nonmapped."\n";
		print STDERR "Unique hit: ".$uniqhit."\n";
		print STDERR "Excluded reads: ".$excluded."\n";
		print STDERR "Read groups on same contig fragment: ".$close."\n";
		print STDERR "Read groups providing contact information: ".$valid."\n";
	    #Louvain network resolution
	    open(F, "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_hic_contacts");
		$cnt = 1;
		while(<F>){
		    chomp;
		    ($n1, $n2) = split(' ', $_);
	        if (!exists($lab{$n1})){
	            $lab{$n1} = $cnt;
	            $rlab{$cnt} = $n1;
	            $cnt++;
	        }
	        if (!exists($lab{$n2})){
	            $lab{$n2} = $cnt;
	            $rlab{$cnt} = $n2;
	            $cnt++;
	        }
	        $pass = 0;
	        if (exists($conn{$lab{$n1}}{$lab{$n2}})){
	            $conn{$lab{$n1}}{$lab{$n2}}++;
	            $pass = 1;
	        }
	        if (exists($conn{$lab{$n2}}{$lab{$n1}})){
	            $conn{$lab{$n2}}{$lab{$n1}}++;
	            $pass = 1;
	        }
	        if ($pass == 0){
	            $conn{$lab{$n1}}{$lab{$n2}}++;
	        }
		}
		close(F);
		$indfn = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_ind_cont.txt";
		$contfn = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contacts.bin";
		$wfn = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contacts.weights";
		$tmptfn = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_tmp.tree";
		$tmpmfn = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_tmp.memb";
		open(NF, ">".$indfn);
		for $c1 (keys %conn){
		    (undef, undef, undef, $len1, undef, undef, undef) = split(/_/, $rlab{$c1});
		    for $c2 (keys %{$conn{$c1}}){
		        (undef, undef, undef, $len2, undef, undef, undef) = split(/_/, $rlab{$c2});
		        if ($len1 >= $len2){
		            print NF ($c1."\t".$c2."\t".($conn{$c1}{$c2}/$len2)."\n");
		        }
		        else{
		            print NF ($c1."\t".$c2."\t".($conn{$c1}{$c2}/$len1)."\n");
		        }
		    }
		}
		close(NF);
		$runpr = $dep{'LOUVAIN'}{'prog'}."convert -i $indfn -o $contfn -w $wfn";
		`$runpr`;

		for $iter (1..100){
		    $runpr = $dep{'LOUVAIN'}{'prog'}."louvain $contfn -l -1 -w $wfn > $tmptfn";
		    `$runpr`;
		    $runpr = $dep{'LOUVAIN'}{'prog'}."hierarchy $tmptfn -m > $tmpmfn";
		    `$runpr`;
		    open(F, $tmpmfn);
		    $hd = <F>; #the algorithm makes a community index 0, we don't need that
		    while(<F>){
		        chomp;
		        ($id, $comm) = split(' ', $_);
		        $data{$id} .= $comm."_";
		    }
		    close(F);
		}

		for $cont (keys %data){
		    push(@{$data2{$data{$cont}}}, $rlab{$cont});
		}
		$cnt = 1;
		open(OUT, ">post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_cc_defs.txt");
		for $cc (keys %data2){
		    $sumlen = 0;
		    for $c (@{$data2{$cc}}){
		        (undef, undef, undef, $len, undef, undef, undef) = split(/_/, $c);
		        $sumlen += $len;
		    }
		        print OUT ("CC_".sprintf('%05d', $cnt)."\t".$sumlen."\t".scalar(@{$data2{$cc}})."\t".join(';', @{$data2{$cc}})."\n");
		        $cnt++;
		}
		close(OUT);
		#Cleanup
		`rm -f $indfn $contfn $wfn $tmptfn $tmpmfn`;

	    #Final core community forming

	    $contigs = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_contigs.fasta";
		$contacts = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_hic_contacts";
		$ccdef = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."_cc_defs.txt";
		$of = "results/".$prj{'PRJF'}."/".$prj{'PRJID'}."_clusters/";
		if (!-d $of){
		    `mkdir $of`;
		}
		$|++;

		print STDERR "Processing project ".$prj_id."...\n";
		print STDERR "Reading in contigs.fasta file...\n";
		open(F, $contigs);
		while(<F>){
		    chomp;
		    if(/^>/){
		        $id = $_;
		        $id =~ s/^>//;
		        (undef, undef, undef, $len, undef, $cov) = split('_', $id);
		        if ($cov > 0){
		            $lens{$id} = $len;
		            $covs{$id} = $cov;
		        }
		    }
		    else{
		        $allseq{$id}.=$_;
		    }
		}
		close(F);

		print STDERR "Reading in Hi-C contacts file...\n";
		open(F, $contacts);
		while(<F>){
		    chomp;
		    ($c1, $c2) = split(' ', $_);
		    $allc{$c1}++;
		    $allc{$c2}++;
		    $cont{$c1}{$c2}++;
		    $cont{$c2}{$c1}++;
		}
		close(F);

		%CCset = ();
		%CClen = ();
		open(F, $ccdef);
		print STDERR "Splitting CCs if necessary...\n";
		while(<F>){
		    chomp;
		    @fl = split(/\t/, $_);
		    @cs = split(/;/, $fl[3]);
		    $CClen{$fl[0]} = 0;
		    for $id (@cs){
		        if ($allc{$id} > 1 && exists($covs{$id})){
		            push(@{$CCset{$fl[0]}}, $id);
		            $CClen{$fl[0]} += $lens{$id};
		        }
		    }
		    #Split the CC if necessary
		    if ($CClen{$fl[0]} > 1000000){
		        $totlen = $CClen{$fl[0]};
		        @sorted_ids = ();
		        foreach $id (sort { $covs{$a} <=> $covs{$b} } @{$CCset{$fl[0]}}) {
		            push(@sorted_ids, $id);
		        }
		        $mincov = $covs{$sorted_ids[0]};
		        $maxcov = $covs{$sorted_ids[$#sorted_ids]};
		        $diff = log10($maxcov) - log10($mincov);
		        $bin = $diff/50;
		        %bins = ();
		        $akt = $mincov;
		        for $i (1..50){
		            if ($akt == $mincov){
		                $bins{$i}{'min'} = 0;
		            }
		            else{
		                $bins{$i}{'min'} = $akt;
		            }
		            $akt = 10**(log10($akt)+$bin);
		            if ($i == 50){
		                $bins{$i}{'max'} = 10000000;
		            }
		            else {
		                $bins{$i}{'max'} = $akt;
		            }
		        }
		        for $id (@sorted_ids){
		            for $b (keys %bins){
		                if ($covs{$id} > $bins{$b}{'min'} && $covs{$id} <= $bins{$b}{'max'}){
		                    push(@{$bins{$b}{'ids'}}, $id);
		                }
		            }
		        }
		        for $b (keys %bins){
		            $binlen = 0;
		            for $id (@{$bins{$b}{'ids'}}){
		                $binlen += $lens{$id};
		            }
		            if ($binlen > $totlen/100){
		                $bins{$b}{'incr'} = 1;
		            }
		            else{
		                $bins{$b}{'incr'} = 0;
		                if ($binlen == 0){
		                    $bins{$b}{'incr'} = 2;
		                }
		            }
		        }
		        $aktlen = 0;
		        $subnum = 1;
		        $pass = 0;
		        %sublens = ();
		        %subccid = ();
		        for $b (1..50){
		            for $id (@{$bins{$b}{'ids'}}){
		                if ($bins{$b}{'incr'} == 1){
		                    $subccid{$id} = $subnum;
		                    $pass = 1;
		                    $sublens{$subnum} += $lens{$id};
		                }
		                if ($bins{$b}{'incr'} == 0){
		                    if ($pass == 1){
		                        $subnum++;
		                        $pass = 0;
		                    }
		                }
		            }
		        }
		        if ($pass == 0){
		            $subnum--;
		        }
		        if ($subnum != 1){
		            $max = 0;
		            for $sn (1..$subnum){
		                if ($sublens{$sn}/$totlen > $max){
		                    $max = $sublens{$sn}/$totlen;
		                }
		            }
		            if ($max < 0.85){
		                print STDERR $fl[0]." needs splitting...\n";
		                $pass = 0;
		                $subsubnum = 1;
		                for $b (1..50){
		                    for $id (@{$bins{$b}{'ids'}}){
		                        if ($bins{$b}{'incr'} == 1){
		                            if ($pass == 0){
		                                $nid = $fl[0]."_s".$subsubnum;
		                                $CClen{$fl[0]} = 0;
		                            }
		                            push(@{$CCset{$nid}} ,$id);
		                            $CClen{$nid} += $lens{$id};
		                            $pass = 1;
		                        }
		                        if ($bins{$b}{'incr'} == 0){
		                            if ($pass == 1){
		                                $subsubnum++;
		                                $pass = 0;
		                            }
		                        }
		                    }
		                }
		            }
		        }
		    }
		}
		close(F);
		print STDERR "Splitting check finished, start CC extension...\n";
		for $ccid (keys %CClen){
		    if ($CClen{$ccid} > 250000 && scalar(@{$CCset{$ccid}}) > 2){
		        print STDERR "Extending ".$ccid."...\n";
		        $initsd = $initmean = 0;
		        $sdmulti = 2;
		        $totlen = $CClen{$ccid};
		        $temp_fasta = "temp/".$prj_id."_".$ccid."_temp.fasta";
		        %seq = ();
		        for $id (@{$CCset{$ccid}}){
		            $seq{$id} = 1;
		        }
		        %banned = ();
		        for $i (1..9){
		            $extended = 1;
		            while($extended == 1 && $totlen < 6000000){
		                # write out the current CC for the blast search
		                open (OUT, ">".$temp_fasta);
		                for $id (keys %seq){
		                    print OUT ">".$id."\n".$allseq{$id}."\n";
		                }
		                close(OUT);
		                # Creating blastdb from the actual CC
		                $runpr = $dep{'MAKEBLDB'}{'prog'}." -dbtype nucl -in $temp_fasta";
		                `$runpr`;
		                # Extract the core for finding contacts
		                $extended = 0;
		                $sublen = 0;
		                %subCClist = ();
		                @akt_covs = ();

		                while($sublen < $totlen*0.8 || scalar(@akt_covs) < 2){
		                    $max = 0;
		                    $maxid = '';
		                    for $id (keys %seq){
		                        if ($lens{$id} > $max && !exists($subCClist{$id})){
		                            $max = $lens{$id};
		                            $maxid = $id;
		                        }
		                    }
		                    $sublen += $max;
		                    $subCClist{$maxid} = 1;
		                    push(@akt_covs, $covs{$maxid});
		                }
		                # Excluding coverage outliers from the core
		                %tmpCClist = ();
		                $tmplen = $sublen;
		                $avcov = mean(@akt_covs);
		                $sdcov = stdev(@akt_covs);
		                $excl = 0;
		                @akt_covs = ();
		                for $id (%subCClist){
		                    if ($covs{$id} > $avcov-$sdmulti*$sdcov && $covs{$id} < $avcov+$sdmulti*$sdcov){
		                        $tmpCClist{$id} = 1;
		                        push(@akt_covs, $covs{$id});
		                    }
		                    else {
		                        $excl++;
		                        $tmplen = $tmplen - $lens{$id};
		                    }
		                }
		                $avcov = mean(@akt_covs);
		                $sdcov = stdev(@akt_covs);
		                if ($initsd == 0){
		                    $initsd = $sdcov;
		                    $initmean = $avcov;
		                }
		                print STDERR "Core mean coverage: ".$avcov." with sd: ".$sdcov.", excluded ".$excl." contigs. Length: ".$sublen."-->".$tmplen."\n";
		                $fail = 0;
		                if ($sdcov > $initsd*2){
		                    if ($sdmulti == 1){
		                        $extended = 0;
		                    }
		                    if ($sdmulti == 2){
		                        $sdmulti = 1;
		                        $extended = 1;
		                    }
		                    print STDERR "Sudden increase in coverage SD detected, reducing SD multiplier for contig involvement to 1. Recover core from previous iteration.\n";
		                    $sdmulti = 1;
		                    %tmpCClist = %preCClist;
		                    %seq = %pre_seq;
		                    $totlen = $pretotlen;
		                    $fail = 1;
		                }
		                if ($avcov < $initmean*0.8 || $avcov > $initmean*1.2){
		                    if ($sdmulti == 1){
		                        $extended = 0;
		                    }
		                    if ($sdmulti == 2){
		                        $sdmulti = 1;
		                        $extended = 1;
		                    }
		                    print STDERR "Sudden change in coverage mean detected, reducing SD multiplier for contig involvement to 1. Recover core from previous iteration.\n";
		                    
		                    %tmpCClist = %preCClist;
		                    %seq = %pre_seq;
		                    $totlen = $pretotlen;
		                    $fail = 1;
		                }
		                if ($fail == 0){
		                    $sublen = $tmplen;
		                    %subCClist = %tmpCClist;
		                    %preCClist = %subCClist;
		                    %pre_seq = %seq;
		                    # Finding connected contigs with the given connection ratio threshold
		                    %connected = ();
		                    %toact = ();
		                    $connlen = 0;
		                    for $c1 (keys %cont){
		                        if (exists($subCClist{$c1})){
		                            for $c2 (keys %{$cont{$c1}}){
		                                if (!exists($seq{$c2})){
		                                    $toact{$c2} += $cont{$c1}{$c2};
		                                }
		                            }
		                        }
		                    }
		                    for $id (keys %toact){
		                        if ($toact{$id} >= (1-$i/10)*$allc{$id} && $covs{$id} > $avcov-$sdmulti*$sdcov && $covs{$id} < $avcov+$sdmulti*$sdcov && !exists($banned{$id})){
		                            $tm_contig = "temp/".$prj_id."_".$ccid."_temp_contig.fasta";
		                            open (TMF, ">".$tm_contig);
		                            print TMF ">".$id."\n".$allseq{$id}."\n";
		                            close(TMF);
		                            $runpr = $dep{'BLASTN'}{'prog'}."-db $temp_fasta -query $tm_contig -perc_identity 85 -outfmt '6 qacc sacc length qlen qstart qend' -culling_limit 1";
		                            @res = `$runpr`;
		                            chomp(@res);
		                            $mstr = '0' x $lens{$id};
		                            for $line (@res){
		                                @fl = split(/\t/, $line);
		                                if ($fl[4] > $lens{$id}){
		                                    print STDERR "Error in blast search: ".$line."\n";
		                                }
		                                else{
		                                    substr($mstr, $fl[4]-1, $fl[2]) = '1' x $fl[2];
		                                }
		                            }
		                            $matches = $mstr =~ tr/1//;
		                            if ($matches < $lens{$id}*0.5){
		                                $connected{$id} = 1;
		                                $connlen += $lens{$id};
		                            }
		                            else{
		                                $banned{$id} = 1;
		                                # print STDERR $id." is banned..\n";
		                            }
		                        }
		                    }
		                    print STDERR "Found ".scalar(keys %connected)." connected contigs (thr=".(1-$i/10)."), additional length:".$connlen."\n";
		                    $pretotlen = $totlen;
		                    if (scalar(keys %connected) > 0){
		                        $extended = 1;
		                        for $id (keys %connected){
		                            $seq{$id} = 1;
		                            $totlen += $lens{$id};
		                        }
		                    }
		                    print STDERR ">>> Total length of the CC: ".$totlen." <<<\n";
		                    if ($i == 9){
		                        $i = 8;
		                    }
		                }
		                
		            }
		        }
		        print STDERR "Writing final fasta file for ".$ccid."...\n";
		        open (O, ">".$of.$ccid.".mfa") or die "Error writing file ".$of.$ccid."!\n";
		        for $id (keys %seq){
		            print O ">".$id."\n".$allseq{$id}."\n";
		        }
		        close(O);
		        `rm $tm_contig`;
		        `rm $temp_fasta*`;
		    }

		    else{
		        print STDERR "Consensus-cluster ".$ccid." was excluded due to small size (".$CClen{$ccid}.")\n";
		    }
		}


		#Do the cleanup in the post-processing directory

		@files = ('_contigs.1.bt2', '_contigs.2.bt2', '_contigs.3.bt2', '_contigs.4.bt2', '_contigs.rev.1.bt2', '_contigs.rev.2.bt2', '_SEL.sam', '_SEL_sorted.sam');
		for $file (@files){
		    $ff = "post_assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}.$file;
		    if (-f $ff){
		        `rm -f $ff`;
		    }
		}
		$folder = "assembly/".$prj{'PRJF'}."/".$prj{'PRJID'}."/";
		if (-d $folder){
		    `rm -rf $folder`;
		}
	}
}

# Clade refinement

if ($step eq 'clade_refinement'){
	@prj_ids = parseIDs(shift(@ARGV));
	%prj = GetPrjDet($prj_id[0]);
	$prjcat = $prj{'PRJF'};

	$bf = "results/".$prjcat."/";
	for $pid (@prj_ids){
		if (!-d $bf.$pid."_clusters"){
			die "Error: ".$bf.$pid."_clusters/ folder cannot be found!\n";
		}
	}

	@sorted_pids = sort @prj_ids;
	$of = $bf."combined_results_".$sorted_pids[0]."_to_".$sorted_pids[$#sorted_pids]."/";
	$oft = $of."temporary_fasta_files/";
	if (-d $of){
		`rm -rf $of`;
	}
	`mkdir $of`;
	`mkdir $oft`;

	for $pid (@prj_ids){
		$fromfd = $bf.$pid."_clusters/";
		@ls = `ls $fromfd`;
		chomp(@ls);
		for $fn (@ls){
			$from = $fromfd.$fn;
			$str = "_".$pid;
			$fn =~ s/\.mfa/$str/ee;
			$fn .= ".mfa";
			$to = $oft.$fn;
			`cp $from $to`;
		}
	}

	$flist = $of."CC_file_list.txt";
	`ls $oft* > $flist`;
	@fls = `cat $flist`;
	chomp(@fls);
	$sketchfile = $of."cc_sketches";
	$runpr = $dep{'MASH'}{'prog'};
	$runopt = "sketch -l $flist -o $sketchfile -s 10000";
	`$runpr $runopt`;
	$sketchfile .= ".msh";
	if (!-f $sketchfile){
		die "Sketchfile wasn't created!\n";
	}
	$dist_res = $of."cc_distance_results.txt";
	open(OUT, ">".$dist_res);
	for $fn (@fls){
		$runpr = $dep{'MASH'}{'prog'};
		$runopt = "dist $sketchfile $fn";
		@res = `$runpr $runopt`;
		chomp(@res);
		for $line (@res){
			@fl = split(/\t/, $line);
			@sfl0 = split(/\//, $fl[0]);
			@sfl1 = split(/\//, $fl[1]);
			@ssfl0 = split(/\./, $sfl0[$#sfl0]);
			@ssfl1 = split(/\./, $sfl1[$#sfl1]);
			print OUT $ssfl0[0]."\t".$ssfl1[0]."\t".$fl[2]."\n";
		}
	}
	close(OUT);

	print STDERR "Constructing distance tree...\n";
	$runpr = $path."upgma.jl";
	$runopt = $dist_res;
	`$runpr $runopt`;
	$restree = $of."full_combined_tree.txt";
	$resclann = $of."clade_text_annotation.txt";
	`mv pairs_upgma.tre $restree`;
	`mv clade_text_annotation.txt $resclann`;
	print STDERR "Performing cleanup...\n";
	`rm -rf $oft`;

	$dbfd = $of."databases/";
	if (!-d $dbfd){
		`mkdir $dbfd`;
	}

	for $p (@prj_ids){
		$db_path = $of."databases/".$p."_blastdb";
		$fasta_path = "post_assembly/".$prjcat."/".$p."_contigs.fasta";
		if (!-f $db_path.".nhr"){
			$runpr = $dep{'MAKEBLDB'}{'prog'};
			$runopt = "-dbtype nucl -in $fasta_path -out $db_path";
			`$runpr $runopt`;
			print STDERR "BLAST database ".$db_path." was created...\n";
		}
		else{
			print STDERR "BLAST database ".$db_path." was alredy created in the target folder...\n";
		}
	}

	# clade refinement 1
	#($name) = $folder =~ /\/combined_results_(.+)\//;
	$file = $of."clade_text_annotation.txt";
	# Create blast databases from contigs fasta
	%plist = ();
	open(F, $file);
	while(<F>){
		if (/^CC/){
			@fl = split(/\t/, $_);
			@fl1 = split('_', $fl[0]);
			$plist{pop(@fl1)} = 1;
		}
	}
	close(F);
	@arr = `grep \"^CC_\" $file | grep -v \"SMALL\" | cut -f 2 | sort | uniq -c | sort -gr | awk \'\$1 > 2\' | awk \'{print \$2}\'`;
	chomp(@arr);

	# clade refinement 2
	$fd = $of;
	$prj_name = $prjcat;
	#$clade_list = @prj_ids;
	@fl = split(/\//, $ENV{'PWD'});
	for $c (@arr){
		$target{$c} = 1;
	}
	$dbfd = $fd."databases/";
	$resfd = $fd."refined_CC_files/";
	$tmpd = $fd."temporary_refinement_fasta_files/";
	$gtdbtkfd = $fd."temporary_gtdbtk_folders/";
	if (!-d $dbfd){
		`mkdir $dbfd`;
	}
	if (!-d $resfd){
		`mkdir $resfd`;
	}
	if (!-d $gtdbtkfd){
		`mkdir $gtdbtkfd`;
	}
	%preg = ();
	open(F, $fd."clade_text_annotation.txt");
	while(<F>){
		if (/^CC/){
			if (/^CC_\d+_s\d_\w\d+/){
				($prjid) = $_ =~ /CC_\d+_s\d_(\w\d+)\s/;
				$preg{$prjid} = 1;
				@fl = split(/\t/, $_);
				push(@{$clades{$fl[1]}}, $fl[0]);
				$cc_data{$fl[0]}{'clade'} = $fl[1];
			}
			else{
				($prjid) = $_ =~ /CC_\d+_(\w\d+)\s/;
				$preg{$prjid} = 1;
				@fl = split(/\t/, $_);
				push(@{$clades{$fl[1]}}, $fl[0]);
				$cc_data{$fl[0]}{'clade'} = $fl[1];
			}
		}
	}
	close(F);
	@plist = keys %preg;
	#choose the best CC for the clade
	%refined_contigs = ();
	for $cl (keys %clades){
		if (exists($target{$cl})){
			@akt_arr = @{$clades{$cl}};
			print STDERR "Finding exampler in clade ".$cl." and do BLAST search in other projects...\n";
			$res = `perl_scripts/post_pipeline/find_clade_best.pl $fd $prj_name $cl`;
			
			$folder = $fd;
			$prj_name = $prjcat;
			$clade = $cl;
			open(F, $folder."clade_text_annotation.txt");
			while(<F>){
				if (/$clade/ && /^CC/ && !/SMALL/){
					($cc, undef) = split(/\t/, $_);
					push(@ccs, $cc);
				}
			}
			close(F);

			if (scalar(@ccs) <= 2 ){
				print "not_enough_CCs";
				die $clade." has less than 2 members, choosing example is not reliable!\n";
			}
			@lens = ();
			$gtdbtkroot = $folder."temporary_gtdbtk_folders/".$clade;
			if (!-d $gtdbtkroot){
				`mkdir $gtdbtkroot`;
			}
			$gtdbtkfd = $folder."temporary_gtdbtk_folders/".$clade."/genomes/";
			if (!-d $gtdbtkfd){
				`mkdir $gtdbtkfd`;
			}
			$gtdbtkfdres = $folder."temporary_gtdbtk_folders/".$clade."/res/";
			if (!-d $gtdbtkfdres){
				`mkdir $gtdbtkfdres`;
			}
			$gtdbtkfdscr = $folder."temporary_gtdbtk_folders/".$clade."/scratch/";
			if (!-d $gtdbtkfdscr){
				`mkdir $gtdbtkfdscr`;
			}
			print STDERR "Copying files into a temporary folder to perform GTDBtk analysis on them...\n";
			for $cc_fid (@ccs){
				if ($cc_fid =~ /_s\d_/){
					(undef, $cc_tmpid, $sub, $sid) = split(/_/, $cc_fid);
					$cc_id = $cc_tmpid."_".$sub;
				}
				else{
					(undef, $cc_id, $sid) = split(/_/, $cc_fid);
				}
				$from = "results/".$prj_name."/".$sid."_clusters/CC_".$cc_id.".mfa";
				$to = $gtdbtkfd.$cc_fid.".mfa";
				open(F, $from);
				%seq = ();
				%size = ();
				while(<F>){
					chomp;
					if (/^>/){
						$id = $_;
						($s) = $id =~ /length_(\d+)_cov/;
						$size{$id} = $s;
					}
					else{
						$seq{$id} = $_;
					}
				}
				close(F);
				open(OUT, ">".$to);
				foreach $id (sort { $size{$b} <=> $size{$a} } keys %size) {
				    print OUT $id."\n".$seq{$id}."\n";
				}
				close(OUT);
				print STDERR ".";
			}

			print STDERR "\n";
			$runpr = $dep{'GTDBTK'}{'prog'};
			$runopt = "identify --cpus ".$dep{'GTDBTK'}{'opt'}." --genome_dir $gtdbtkfd -x mfa --prefix $clade --out_dir $gtdbtkfdres";
			`$runpr $runopt`;
			$bacmarkerfile = $gtdbtkfdres.$clade.".bac120.markers_summary.tsv";
			$armarkerfile = $gtdbtkfdres.$clade.".ar122.markers_summary.tsv";
			@bacmiss = `cut -f 4 $bacmarkerfile`;
			chomp(@bacmiss);
			@armiss = `cut -f 4 $armarkerfile`;
			chomp(@armiss);
			if (average(@bacmiss) == 120 && average(@bacmiss) == 122){
				#No markers were found
				print "No markers were found, potentially no bacterial or archeal assembly!";
			}
			if ((120-average(@bacmiss))/120 > (122-average(@armiss))/122){
				#Mostly bacterial genomes
				@res = `cut -f 1-3 $bacmarkerfile | grep \"^CC\" | sort -grk2 | head`;
				print STDERR "Found mostly bacterial markers...\n";
			}
			else{
				#mostly archeal genomes
				@res = `cut -f 1-3 $armarkerfile | grep \"^CC\" | sort -grk2 | head`;
				print STDERR "Found mostly archeal markers...\n";
			}
			chomp(@res);

			$max = $maxu = 0;
			for $line (@res){
				print STDERR $line."\n";
				@fl = split(/\t/, $line);
				if ($fl[2] == 0){
					$fl[2] = 1;
				}
				if ($fl[1] / $fl[2] > $max && $fl[1] > $maxu * 0.8){
					$max = $fl[1] / $fl[2];
					$maxu = $fl[1];
					$best = $fl[0];
				}
			}

			print STDERR $best." was chosen as the example with ".$maxu." unique single copy genes and ".$max." unique/multiple ratio.\n";

			`rm -rf $gtdbtkroot`;
			$res = $best;
			if ($res =~ /^CC/){
				if ($res =~ /_s\d_/){ # Changed from /_s\d_/
					(undef, $tmpcid, $sub, $pid) = split(/_/, $res);
					$cid = $tmpcid."_".$sub;
				}
				else{
					(undef, $cid, $pid) = split(/_/, $res);
				}
				$cc_path = "results/".$prj_name."/".$pid."_clusters/CC_".$cid.".mfa";
				#Use the core of the example for blast (the collection of the longest contigs that take up the 90% of the whole nucleotide content)
				print STDERR "Extracting the core of the exemplar CC...\n";
				open(F, $cc_path);
				%seq = ();
				$totlen = 0;
				%lens = ();
				while(<F>){
					chomp;
					if (/^>/){
						$id = $_;
						($len) = $id =~ /length_(\d+)_cov/;
						$totlen += $len;
						$lens{$id} = $len;
					}
					else{
						$seq{$id} = $_;
					}
				}
				close(F);
				print STDERR "Total nucleotide content of the CC: ".$totlen."...\n";
				$sublen = 0;
				%subCClist = ();
				while($sublen < $totlen*0.9){
					$max = 0;
					$maxid = '';
					for $id (keys %lens){
						if ($lens{$id} > $max && !exists($subCClist{$id})){
							$max = $lens{$id};
							$maxid = $id;
						}
					}
					$sublen += $max;
					$subCClist{$maxid} = 1;
					print STDERR "Including ".scalar(keys %subCClist)." contigs with a total of ".$sublen." nucleotides...\n";
				}
				$ccc_path = $dbfd."query_".$cl.".fasta";
				open (OUT, ">".$ccc_path);
				for $id (keys %subCClist){
					print OUT $id."\n".$seq{$id}."\n";
				}
				close(OUT);
				$to_path = $resfd."CC_".$pid."_".$cl."_ex.mfa";
				`cp $cc_path $to_path`;
				print STDERR "Start the BLAST cycles...\n";
				for $tcc (sort @plist){
					if ($tcc ne $pid){
						$fasta_path = "post_assembly/".$prj_name."/".$tcc."_contigs.fasta";
						$db_path = $dbfd.$tcc."_blastdb";
						$runpr = $dep{'BLASTN'}{'prog'};
						$runopt = '-num_threads '.$dep{'BLASTN'}{'opt'}.' -query $ccc_path -db $db_path -perc_identity 98 -outfmt \'6 qacc sacc length qlen slen sstart send pident\' | awk \'\$3 > 100\'';
						@blres = `$runpr $runopt`;
						$totlen = 0;
						chomp(@blres);
						%tmp_arr = ();
						for $line (@blres){
							@fl = split(/\t/, $line);
							if (!exists($tmp_arr{$fl[1]})){
								($len) = $fl[1] =~ /length_(\d+)_cov/;
								$tmp_arr{$fl[1]} = '0' x $len;
							}
							substr($tmp_arr{$fl[1]}, $fl[5], $fl[2]) = '1' x $fl[2];
						}
						@tmp_arr2 = ();
						for $cont (keys %tmp_arr){
							$len = length($tmp_arr{$cont});
							$ones = $tmp_arr{$cont} =~ tr/1//;
							if ($ones/$len > 0.3){
								$totlen += $len;
								push(@tmp_arr2, $cont);
							}
						}
						if ($totlen > $sublen*0.2){
							@{$refined_contigs{$tcc}{$cl}} = @tmp_arr2;
							print STDERR "Total length of contigs found by BLAST for project ".$tcc." for clade ".$cl.": ".$totlen."\n";
						}
						else {
							print STDERR "Blast search resulted low amount of DNA (".$totlen.") in contigs from project ".$tcc." no CC will be refined.\n";
						}
					}
				}
			}
			print $cl." -> ".$res."\t(".scalar(@{$clades{$cl}}).")\n";
		}
	}
	print STDERR "Finished with BLAST search...\n";

	for $prj (sort @plist){
		print STDERR "Excluding contigs with no contact to the main block in project ".$prj."\n";
		open(F, "post_assembly/".$prj_name."/".$prj."_hic_contacts");
		%cc_cont = ();
		%cc_cnum = ();
		%allc = ();
		while(<F>){
			chomp;
			($c1, $c2) = split(/\t/, $_);
			$cc_cont{$c1}{$c2} = 1;
			$cc_cont{$c2}{$c1} = 1;
			$cc_cnum{$c1}{$c2}++;
			$cc_cnum{$c2}{$c1}++;
			$allc{$c1}++;
			$allc{$c2}++;
		}
		close(F);
		%seq = ();
		open(F, "post_assembly/".$prj_name."/".$prj."_contigs.fasta");
		while(<F>){
			chomp;
			if (/^>/){
				$id = $_;
				$id =~ s/^>+//;
				$seq{$id} = '';
			}
			else{
				$seq{$id} .= $_;
			}
		}
		close(F);
		%amr_seq = ();
		open(F, "post_assembly/".$prj_name."/".$prj."_amr_contigs.fasta");
		while(<F>){
			chomp;
			if (/^>/){
				$id = $_;
				$id =~ s/^>+//;
				$amr_seq{$id} = '';
			}
			else{
				$amr_seq{$id} .= $_;
			}
		}
		close(F);
		for $cl (keys %{$refined_contigs{$prj}}){
			#Filtering predicted CCs for contigs without contact to the others and the main block
			@akt_arr = @{$refined_contigs{$prj}{$cl}};
			@ref_arr = ();
			%int = ();
			for $c1 (@akt_arr){
				for $c2 (@akt_arr){
					if($c1 ne $c2 && exists($cc_cont{$c1}{$c2})){
						$int{$c1}++;
						$int{$c2}++;
					}
				}
			}
			$totlen = 0;
			for $c (@akt_arr){
				if(exists($int{$c})){
					push(@ref_arr, $c);
					($len) = $c =~ /length_(\d+)_cov/;
					$totlen += $len;
				}
			}
			print STDERR "Length after filtering out not connected contigs from refined CC for project ".$prj.": ".$totlen."\n";
			@akt_arr = @ref_arr;
			$blnum = 1;
			%blocks = ();
			while(scalar(@ref_arr) > 0){
				$c = shift(@ref_arr);
				($len) = $c =~ /length_(\d+)_cov/;
				$blocks{$blnum}{$c} = 1;
				$blockslen{$blnum} = $len;
				$found = 1;
				while($found == 1){
					$found = 0;
					@cas = keys %{$blocks{$blnum}};
					for $ca (@cas){
						@tmp = ();
						for $cb (@ref_arr){
							if ($cc_cont{$ca}{$cb} == 1){
								$found = 1;
								$blocks{$blnum}{$cb} = 1;
								($len) = $cb =~ /length_(\d+)_cov/;
								$blockslen{$blnum} += $len;
							}
							else{
								push(@tmp, $cb);
							}
						}
						@ref_arr = @tmp;
					}
				}
				$blnum++;
			}
			$max = 0;
			for $blnum (keys %blocks){
				print STDERR "Block number: ".$blnum."\tNumber of contigs: ".scalar(keys %{$blocks{$blnum}})."\tLength of the block: ".$blockslen{$blnum}."\n";
				if (scalar(keys %{$blocks{$blnum}}) > 1 && $blockslen{$blnum} > $max){
					$max = $blockslen{$blnum};
					$winner = $blnum;
				}
			}
			print STDERR "Keeping contigs from block: ".$winner."\n";
			$totlen = 0;
			for $c (keys %{$blocks{$winner}}){
				($len) = $c =~ /length_(\d+)_cov/;
				$totlen += $len;
			}
			$sublen = 0;
			%subCClist = ();
			while($sublen < $totlen*0.9){
				$max = 0;
				$maxc = '';
				for $c (keys %{$blocks{$winner}}){
					($len) = $c =~ /length_(\d+)_cov/;
					if ($len > $max && !exists($subCClist{$c})){
						$max = $len;
						$maxc = $c;
					}
				}
				$sublen += $max;
				$subCClist{$maxc} = 1;
			}
			%toact = ();
			# Core contigs are in the winner hash
			for $c1 (keys %cc_cnum){
				for $c2 (keys %{$cc_cnum{$c1}}){
					if(!exists($subCClist{$c1}) && exists($subCClist{$c2}) && exists($amr_seq{$c1})){
						$toact{$c1} += $cc_cnum{$c1}{$c2};
					}
					if(exists($subCClist{$c1}) && !exists($subCClist{$c2}) && exists($amr_seq{$c2})){
						$toact{$c2} += $cc_cnum{$c1}{$c2};
					}
				}
			}

			$totlen = 0;
			open(NF, ">".$resfd."CC_".$prj."_".$cl."_rf.mfa");
			for $c (keys %{$blocks{$winner}}){
				print NF ">".$c."\n".$seq{$c}."\n";
				($len) = $c =~ /length_(\d+)_cov/;
				$totlen += $len;
			}
			print STDERR "Refined CC size after keeping main block for project ".$prj.": ".$totlen."\n";
			for $c (keys %toact){
				print NF ">".$c."\n".$amr_seq{$c}."\n";
				($len) = $c =~ /length_(\d+)_cov/;
				$totlen += $len;
			}
			close(NF);
			print STDERR "Final size of refined CC after keeping main block and extending CC for project ".$prj.": ".$totlen."\n";
		}
	}



	$fd = $resfd;
	$thr = 250000;

	if ($fd !~ /\/$/){
		$fd = $fd."/";
	}

	$exfd = $fd."examples/";
	$allfd = $fd."allbig_by_clade/";
	$smallfd = $fd."small_CCs/";


	if (!-d $exfd){
		`mkdir $exfd`;
	}
	if (!-d $allfd){
		`mkdir $allfd`;
	}
	if (!-d $smallfd){
		`mkdir $smallfd`;
	}

	$s = $fd."*.mfa";
	@files = `ls $s`;
	chomp(@files);
	print STDERR "Calculating CC sizes...\n";
	$cnt = 0;
	%data = ();
	for $f (@files){
		$cnt++;
		if ($cnt =~ /[05]0$/){
			print STDERR "Processed ".$cnt." files...\n";
		}
		$size = `grep -v "^>" $f | wc -m`;
		chomp($size);
		$ksize = int($size/1000);
		($prj, $cl) = $f =~ /CC_(.+)_clade_(\d+)_/;
		$data{$cl}{$prj}{'fn'} = $f;
		$data{$cl}{$prj}{'size'} = $size;
		$data{$cl}{$prj}{'ksize'} = $ksize;
	}

	for $cl (keys %data){
		print STDERR "Processing clade_".$cl."...\n";
		$clfd = $allfd.$cl."/";
		if (!-d $clfd){
			`mkdir $clfd`;
		}
		%big = ();
		for $prj (keys %{$data{$cl}}){
			$fn = $data{$cl}{$prj}{'fn'};
			if ($data{$cl}{$prj}{'size'} < $thr){
				if ($data{$cl}{$prj}{'size'} == 0){
					`rm $fn`;
					print STDERR "Empty CC (".$fn.") was deleted\n";
				}
				else {
					`mv $fn $smallfd`;
				}
			}
			else{
				$big{$prj} = $data{$cl}{$prj}{'size'};
			}
		}
		$median = median(values %big);
		$min = 1000000000000;
		$medp = '';
		for $prj (keys %big){
			if (abs($big{$prj}-$median) < $min && $big{$prj} > $big{$medp}){
				$min = abs($big{$prj}-$median);
				$medp = $prj;
			}
		}
		for $prj (keys %big){
			$fn = $data{$cl}{$prj}{'fn'};
			@fnl = split(/\//, $data{$cl}{$prj}{'fn'});
			$nfn = pop(@fnl);
			$end = "_".$data{$cl}{$prj}{'ksize'}."k";
			$nfn =~ s/\.mfa/$end/ee;
			$nfn .= ".mfa";
			if ($prj eq $medp){
				$to = $exfd.$nfn;
				`cp $fn $to`;
			}
			$to = $clfd.$nfn;
			`mv $fn $to`;
		}
	}

	$fd = $resfd;

	$db = "~/rds/rds-mah1-vet_amr/mah1_scripts/resfinder_blast_db/resfinder.fas";
	@colours = ('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928');
	if ($fd !~ /\/$/){
		$fd = $fd."/";
	}
	$exfd = $fd."examples/";
	$bigfd = $fd."allbig_by_clade/";
	$tmpfd = $fd."tmpfd/";
	if (!-d $tmpfd){
		`mkdir $tmpfd`;
	}
	else{
		`rm -f $tmpfd/*`;
	}
	$resfd = $fd."analysis_results/";
	if (!-d $resfd){
		`mkdir $resfd`;
	}
	else{
		`rm -f $resfd/*`;
	}
	$restxtfd = $resfd."clade_amr_tables/";
	if (!-d $restxtfd){
		`mkdir $restxtfd`;
	}
	else{
		`rm -f $restxtfd/*`;
	}
	@ls = `ls $exfd`;
	chomp(@ls);
	for $fn (@ls){
		$fn =~ s/\.mfa//;
		($akt_cl) = $fn =~ /clade_([0-9]+)_/;
		$clade_ex{$akt_cl} = $fn;
	}

	print STDERR "Looking for AMR genes in all clades...\n";
	@clades = `ls $bigfd`;
	chomp(@clades);
	%cl_amr_data = ();
	%clnum = ();
	for $cl (@clades){
		print STDERR "Analysing clade ".$cl."...\n";
		$cldir = $bigfd.$cl."/";
		@ls = `ls $cldir`;
		$clnum{$cl} = scalar(@ls);
		chomp(@ls);
		for $fn (@ls){
			($prjid, $s) = $fn =~ /CC_(.+)_clade_\d+_[er][xf]_(\d+)k\.mfa/;
			$prj_data{$cl}{$prjid} = 1;
			$size{$cl}{$prjid} = $s;
			$query = $cldir.$fn;
			$runpr = $dep{'ABRICATE'}{'prog'}." --mincov 50 --threads ".$dep{'ABRICATE'}{'opt'}." --db resfinder $infile"; 
			@res = `$runpr`;
			for $line (@res){
				chomp($line);
				@fl = split(' ', $line);
				$cl_samr_data{$fl[5]}{$cl}{$prjid} = 1;
			}
		}
	}
	open(OUT, ">".$resfd."amr_annotations.txt");
	print OUT "DATASET_MULTIBAR\nSEPARATOR TAB\nDATASET_LABEL\tAntimicrobial resistance genes\nCOLOR\t#000000\n";
	print OUT "FIELD_COLORS";
	$cnt = $cnt1 = 0;
	$label_str = "FIELD_LABELS";
	$scale_str = "DATASET_SCALE";
	$genes_str = "CC_ID\tCLADE\tSIZE_(KB)";
	for $gene (sort(keys %cl_samr_data)){
		print OUT "\t".$colours[$cnt]."\t#ffffff";
		$genes_str .= "\t".$gene;
		$label_str .= "\t".$gene."\tNA";
		$gene =~ s/\-/_/g;
		$scale_str .= "\t".$cnt1."-0-#000000-1-0-1"."\t".($cnt1 + 0.5)."-".$gene."-".$colours[$cnt]."-2-0-2";
		$cnt1++;
		$cnt++;
		if ($cnt > $#colours){
			$cnt = 0;
		}
	}
	print OUT "\n";
	print OUT $label_str."\n";
	print OUT $scale_str."\n";
	print OUT "ALIGN_FIELDS\t0\nBORDER_WIDTH\t0.5\nMARGIN\t2\nWIDTH\t2000\nHEIGHT_FACTOR\t0.7\n";
	print OUT "DATA\n";
	$num_ccs_str = "DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\tNumber of CCs\nDATA\n";
	open(OUT3, ">".$resfd."full_amr_table.txt");
	print OUT3 $genes_str."\n";
	for $cl (keys %clade_ex){
		open(OUT2, ">".$restxtfd.$cl.".txt");
		print OUT2 $genes_str."\n";
		$num_ccs_str .= $clade_ex{$cl}."\tn=".$clnum{$cl}."\t-1\t#000000\tnormal\t1\t0\n";
		print OUT $clade_ex{$cl};
		for $gene (sort(keys %cl_samr_data)){
			if (exists($cl_samr_data{$gene}{$cl})){
				print OUT "\t".(scalar(keys %{$cl_samr_data{$gene}{$cl}}) / $clnum{$cl})."\t".(1-(scalar(keys %{$cl_samr_data{$gene}{$cl}}) / $clnum{$cl}));
			}
			else{
				print OUT "\t0\t1";
			}
		}
		for $prj (keys %{$prj_data{$cl}}){
			print OUT2 $prj."\t".$cl."\t".$size{$cl}{$prj};
			print OUT3 $prj."\t".$cl."\t".$size{$cl}{$prj};
			for $gene (sort(keys %cl_samr_data)){
				if (exists($cl_samr_data{$gene}{$cl}{$prj})){
					print OUT2 "\t1";
					print OUT3 "\t1";
				}
				else{
					print OUT2 "\t0";
					print OUT3 "\t0";
				}
			}
			print OUT2 "\n";
			print OUT3 "\n";
		}
		close(OUT2);
		print OUT "\n";
	}
	close(OUT3);
	close(OUT);
	open(OUT, ">".$resfd."ccnum_annotations.txt");
	print OUT $num_ccs_str;
	close(OUT);



	$file = $fd."analysis_results/full_amr_table.txt";
	@prjids = `cut -f 1 $file | sort | grep -v \"CC_ID\" | uniq`;
	chomp(@prjids);
	for $prjid (@prjids){

		$path = $fd;
		$prjf = $prjcat;
		$clfd = $path."allbig_by_clade/";
		$of = $path."allbig_by_clade_extended/";
		$pafolder = "post_assembly/".$prjf."/";

		if (!-d $of){
		    `mkdir $of`;
		}

		@clades = `ls $clfd`;
		chomp(@clades);
		for $sclfd (@clades){
		    $tmp = $clfd.$sclfd;
		    @fdlist = `ls $tmp`;
		    chomp (@fdlist);
		    for $fn (@fdlist){
		        ($tpid, $clid) = $fn =~ /^CC_(.+)_clade_(\d+)_[er][xf]/;
		        $ccid = $fn;
		        $ccid =~ s/\.mfa$//;
		        if ($tpid eq $pid){
		            $ccpath{$fn} = $tmp."/".$fn;
		            open (F, $ccpath{$fn});
		            while (<F>){
		                if (/^>/){
		                    chomp;
		                    $id = $_;
		                    $id =~ s/^>//;
		                    ($len) = $id =~ /length_(\d+)_cov/;
		                    $CClen{$ccid} += $len;
		                    push(@{$CCset{$ccid}}, $id);
		                    $clade_id{$ccid} = $clid;
		                }
		            }
		            close(F);
		        }
		    }
		}

		print STDERR "Processing project ".$pid."...\n";
		print STDERR "Reading in contigs.fasta file...\n";
		open(F, $pafolder.$pid."_contigs.fasta");
		while(<F>){
		    chomp;
		    if(/^>/){
		        $id = $_;
		        $id =~ s/^>//;
		        (undef, undef, undef, $len, undef, $cov) = split('_', $id);
		        if ($cov > 0){
		            $lens{$id} = $len;
		            $covs{$id} = $cov;
		        }
		    }
		    else{
		        $allseq{$id}.=$_;
		    }
		}
		close(F);
		print STDERR "Reading in amr_contigs.fasta file...\n";
		open(F, $pafolder.$pid."_amr_contigs.fasta");
		while(<F>){
		    chomp;
		    if (/^>/){
		        $id = $_;
		        $id =~ s/^>//;
		        $amr_contigs{$id} = 1;
		    }
		}
		close(F);

		print STDERR "Reading in Hi-C contacts file...\n";
		open(F, $pafolder.$pid."_hic_contacts");
		while(<F>){
		    chomp;
		    ($c1, $c2) = split(' ', $_);
		    $allc{$c1}++;
		    $allc{$c2}++;
		    $cont{$c1}{$c2}++;
		    $cont{$c2}{$c1}++;
		}
		close(F);


		for $ccid (keys %CClen){
		    if ($CClen{$ccid} > 250000 && scalar(@{$CCset{$ccid}}) > 2){
		        print STDERR "Extending ".$ccid."...\n";
		        $initsd = $initmean = 0;
		        $sdmulti = 2;
		        $totlen = $CClen{$ccid};
		        $temp_fasta = "temp/".$pid."_".$ccid."_temp.fasta";
		        %seq = ();
		        for $id (@{$CCset{$ccid}}){
		            $seq{$id} = 1;
		        }
		        %banned = ();
		        for $i (1..9){
		            $extended = 1;
		            while($extended == 1 && $totlen < 6000000){
		                open (OUT, ">".$temp_fasta);
		                for $id (keys %seq){
		                    print OUT ">".$id."\n".$allseq{$id}."\n";
		                }
		                close(OUT);
		                # Creating blastdb from the actual CC
		                $runpr = $dep{'MAKEBLDB'}{'prog'};
		                $runopt = "-dbtype nucl -in $temp_fasta";
		                `$runpr $runopt`;
		                # Extract the core for finding contacts
		                $extended = 0;
		                $sublen = 0;
		                %subCClist = ();
		                @akt_covs = ();

		                while($sublen < $totlen*0.8 || scalar(@akt_covs) < 2){
		                    $max = 0;
		                    $maxid = '';
		                    for $id (keys %seq){
		                        if ($lens{$id} > $max && !exists($subCClist{$id})){
		                            $max = $lens{$id};
		                            $maxid = $id;
		                        }
		                    }
		                    $sublen += $max;
		                    $subCClist{$maxid} = 1;
		                    push(@akt_covs, $covs{$maxid});
		                }
		                # Excluding coverage outliers from the core
		                %tmpCClist = ();
		                $tmplen = $sublen;
		                $avcov = mean(@akt_covs);
		                $sdcov = stdev(@akt_covs);
		                $excl = 0;
		                @akt_covs = ();
		                for $id (%subCClist){
		                    if ($covs{$id} > $avcov-$sdmulti*$sdcov && $covs{$id} < $avcov+$sdmulti*$sdcov){
		                        $tmpCClist{$id} = 1;
		                        push(@akt_covs, $covs{$id});
		                    }
		                    else {
		                        $excl++;
		                        $tmplen = $tmplen - $lens{$id};
		                    }
		                }
		                $avcov = mean(@akt_covs);
		                $sdcov = stdev(@akt_covs);
		                if ($initsd == 0){
		                    $initsd = $sdcov;
		                    $initmean = $avcov;
		                }
		                print STDERR "Core mean coverage: ".$avcov." with sd: ".$sdcov.", excluded ".$excl." contigs. Length: ".$sublen."-->".$tmplen."\n";
		                $fail = 0;
		                if ($sdcov > $initsd*2){
		                    if ($sdmulti == 1){
		                        $extended = 0;
		                    }
		                    if ($sdmulti == 2){
		                        $sdmulti = 1;
		                        $extended = 1;
		                    }
		                    print STDERR "Sudden increase in coverage SD detected, reducing SD multiplier for contig involvement to 1. Recover core from previous iteration.\n";
		                    $sdmulti = 1;
		                    %tmpCClist = %preCClist;
		                    %seq = %pre_seq;
		                    $totlen = $pretotlen;
		                    $fail = 1;
		                }
		                if ($avcov < $initmean*0.8 || $avcov > $initmean*1.2){
		                    if ($sdmulti == 1){
		                        $extended = 0;
		                    }
		                    if ($sdmulti == 2){
		                        $sdmulti = 1;
		                        $extended = 1;
		                    }
		                    print STDERR "Sudden change in coverage mean detected, reducing SD multiplier for contig involvement to 1. Recover core from previous iteration.\n";
		                    
		                    %tmpCClist = %preCClist;
		                    %seq = %pre_seq;
		                    $totlen = $pretotlen;
		                    $fail = 1;
		                }
		                if ($fail == 0){
		                    $sublen = $tmplen;
		                    %subCClist = %tmpCClist;
		                    %preCClist = %subCClist;
		                    %pre_seq = %seq;
		                    # Finding connected contigs with the given connection ratio threshold
		                    %connected = ();
		                    %toact = ();
		                    $connlen = 0;
		                    for $c1 (keys %cont){
		                        if (exists($subCClist{$c1})){
		                            for $c2 (keys %{$cont{$c1}}){
		                                if (!exists($seq{$c2})){
		                                    $toact{$c2} += $cont{$c1}{$c2};
		                                }
		                            }
		                        }
		                    }
		                    for $id (keys %toact){
		                        if ($toact{$id} >= (1-$i/10)*$allc{$id} && $covs{$id} > $avcov-$sdmulti*$sdcov && $covs{$id} < $avcov+$sdmulti*$sdcov && !exists($banned{$id}) && !exists($amr_contigs{$id})){
		                            $tm_contig = "temp/".$pid."_".$ccid."_temp_contig.fasta";
		                            open (TMF, ">".$tm_contig);
		                            print TMF ">".$id."\n".$allseq{$id}."\n";
		                            close(TMF);
		                            $runpr = $dep{'BLASTN'}{'prog'};
		                            $runopt = '-db $temp_fasta -query $tm_contig -perc_identity 85 -outfmt \'6 qacc sacc length qlen qstart qend\' -culling_limit 1';
		                            `$runpr $runopt`;
		                            # print Dumper(@res);
		                            chomp(@res);
		                            $mstr = '0' x $lens{$id};
		                            for $line (@res){
		                                @fl = split(/\t/, $line);
		                                if ($fl[4] > $lens{$id}){
		                                    print STDERR "Error in blast search: ".$line."\n";
		                                }
		                                else{
		                                    substr($mstr, $fl[4]-1, $fl[2]) = '1' x $fl[2];
		                                }
		                            }
		                            $matches = $mstr =~ tr/1//;
		                            if ($matches < $lens{$id}*0.5){
		                                $connected{$id} = 1;
		                                $connlen += $lens{$id};
		                            }
		                            else{
		                                $banned{$id} = 1;
		                            }
		                        }
		                    }
		                    print STDERR "Found ".scalar(keys %connected)." connected contigs (thr=".(1-$i/10)."), additional length:".$connlen."\n";
		                    $pretotlen = $totlen;
		                    if (scalar(keys %connected) > 0){
		                        $extended = 1;
		                        for $id (keys %connected){
		                            $seq{$id} = 1;
		                            $totlen += $lens{$id};
		                        }
		                    }
		                    print STDERR ">>> Total length of the CC: ".$totlen." <<<\n";
		                    if ($i == 9){
		                        $i = 8;
		                    }
		                }
		                
		            }
		        }
		        print STDERR "Writing final fasta file for ".$ccid."...\n";
		        $ofsd = $of.$clade_id{$ccid};
		        if (!-d $ofsd){
		            `mkdir $ofsd`;
		        }
		        open (O, ">".$ofsd."/".$ccid.".mfa") or die "Error writing file ".$ofsd."/".$ccid."!\n";
		        for $id (keys %seq){
		            print O ">".$id."\n".$allseq{$id}."\n";
		        }
		        close(O);
		        `rm $tm_contig`;
		        `rm $temp_fasta*`;
		    }

		    else{
		        print STDERR "Core-community ".$ccid." was not extended due to small size or too few contigs (".$CClen{$ccid}.")\n";
		        $ofsd = $of.$clade_id{$ccid};
		        if (!-d $ofsd){
		            `mkdir $ofsd`;
		        }
		        open (O, ">".$ofsd."/".$ccid.".mfa") or die "Error writing file ".$ofsd."/".$ccid."!\n";
		        for $id (@{$CCset{$ccid}}){
		            print O ">".$id."\n".$allseq{$id}."\n";
		        }
		        close(O);
		    }
		}
	}

	@sflist = `ls $of`;
	chomp(@sflist);
	print STDERR "Preparing fasta files for GTDBtk analysis...\n";
	for $sfd (@sflist){
		print STDERR "Processing clade ".$sfd."\n";
		$s = $fd.$sfd."/*.mfa";
		@flist = `ls $s`;
		chomp(@flist);
		for $fn (@flist){
			open(F, $fn);
			%seq = ();
			%size = ();
			while(<F>){
				chomp;
				if (/^>/){
					$id = $_;
					($s) = $id =~ /length_(\d+)_cov/;
					$size{$id} = $s;
				}
				else{
					$seq{$id} = $_;
				}
			}
			close(F);
			open(OUT, ">".$fn);
			foreach $id (sort { $size{$b} <=> $size{$a} } keys %size) {
			    print OUT $id."\n".$seq{$id}."\n";
			}
			close(OUT);
			print STDERR ".";
		}
		print STDERR "\n";
	}
	print STDERR "Performing GTDBtk analysis...\n";
	$tfd = $of;
	$tfd =~ s/allbig_by_clade_extended\/$//;
	$tfd .= "tmpfd/";

	@clades = `ls $fd`;
	chomp(@clades);
	for $cl (@clades){
		if ($cl =~ /\d{5}/){
			print STDERR "Creating output folder and start up process for clade: ".$cl."\n";
			$ctfd = $tfd.$cl."/";
			if (!-d $ctfd){
				`mkdir $ctfd`;
			}
			else {
				`rm -rf $ctfd*`;
			}
			$gfd = $of.$cl."/";
			$runpr = $dep{'GTDBTK'}{'prog'}." classify_wf";
			$runopt = "--genome_dir $gfd -x mfa --prefix $cl --cpus ".$dep{'GTDBTK'}{'opt'}." --out_dir $ctfd";
			`$runpr $runopt`;
		}
	}	

	$ofd = $fd."analysis_results/clade_summary_files/";
	if (!-d $ofd){
		`mkdir $ofd`;
	}
	@colours = ('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928');

	#Read in the full amr table data
	open(F, $fd."analysis_results/full_amr_table.txt");
	open(FF, ">".$fd."analysis_results/full_amr_table_filtered.txt");
	$hd = <F>;
	print FF $hd;
	chomp($hd);
	@amr_hd = split(/\t/, $hd);
	while(<F>){
		chomp;
		@fl = split(/\t/, $_);
		$pre_ext_size{$fl[1]}{$fl[0]} = $fl[2];
		$full_line{$fl[1]}{$fl[0]} = $_;
		for $i (3..$#fl){
			$amr_data{$fl[1]}{$fl[0]}{$amr_hd[$i]} = $fl[$i];
			$cl_amr_pres{$fl[1]}{$amr_hd[$i]} += $fl[$i];
			
		}
	}
	close(F);

	$datafd = $fd."tmpfd/";
	@cls = `ls $datafd`;
	chomp(@cls);

	%add_genomes = ();

	for $cl (@cls){
		print STDERR "Processing clade ".$cl."\n";
		$bac_sum = $filt = $ar_sum = 0;
		%bac_data = %ar_data = %tax_data = %filt_data = %cc_full_names = ();
		%clade_genomes = ();
		if (!-f $datafd.$cl."/gtdbtk.log"){
			print STDERR "ERROR - GTDBtk log wasn't found in folder ".$datafd.$cl."!\n";
		}
		else{
			$genomes_tmp_folder = $datafd."mash_tmp/";
			if (!-d $genomes_tmp_folder){
				`mkdir $genomes_tmp_folder`;
			}
			else {
				`rm -f $genomes_tmp_folder*`;
			}
			if (-f $datafd.$cl."/".$cl.".bac120.summary.tsv"){
				open(F, $datafd.$cl."/".$cl.".bac120.summary.tsv");
				$hd = <F>;
				while(<F>){
					if (/^CC/){
						chomp;
						@fl = split(/\t/, $_);
						($prj) = $fl[0] =~ /CC_(.+)_clade/;
						$cc_full_names{$prj} = $fl[0];
						($bac_data{$prj}{'tax'}, $bac_data{$prj}{'class'}) = GetLastTaxa($fl[1]);
						$bac_data{$prj}{'aa_perc'} = $fl[15];
						$bac_data{$prj}{'warn'} = $fl[18];
						$tax_data{$fl[1]}++;
						$bac_sum++;
						$fasta_path = $fd."allbig_by_clade_extended/".$cl."/".$fl[0].".mfa";
						$size = `grep -v \"^>\" $fasta_path | wc -m`;
						`cp $fasta_path $genomes_tmp_folder`;
						# if (exists($target{$cl})){
						# 	`cp $fasta_path $target_dir`;
						# }
						chomp($size);
						$bac_data{$prj}{'size'} = int($size/1000);
						print FF $full_line{$cl}{$prj}."\n";
						@add_arr = split('; ', $fl[14]);
						for $add_s (@add_arr){
							@sfl = split(', ', $add_s);
							if ($sfl[1] =~ /^s__/){
								$sfl[1] =~ s/^s__//;
								$clade_genomes{$sfl[0]} = $sfl[1];
								if (exists($target{$cl})){
									$add_genomes{$sfl[0]} = $sfl[1];
								}
							}
						}
					}
					
				}
				close(F);
			}
			if (-f $datafd.$cl."/".$cl.".bac120.filtered.tsv"){
				open(F, $datafd.$cl."/".$cl.".bac120.filtered.tsv");
				while(<F>){
					if (/^CC/){
						chomp;
						@fl = split(/\t/, $_);
						($prj) = $fl[0] =~ /CC_(.+)_clade/;
						$filt_data{$prj} = $fl[1];
						$filt++;
					}
				}
				close(F);
			}

			if (-f $datafd.$cl."/".$cl.".ar122.summary.tsv"){
				open(F, $datafd.$cl."/".$cl.".ar122.summary.tsv");
				$hd = <F>;
				while(<F>){
					if (/^CC/){
						chomp;
						@fl = split(/\t/, $_);
						($prj) = $fl[0] =~ /CC_(.+)_clade/;
						$cc_full_names{$prj} = $fl[0];
						($ar_data{$prj}{'tax'}, $ar_data{$prj}{'class'}) = GetLastTaxa($fl[1]);
						$ar_data{$prj}{'aa_perc'} = $fl[15];
						$ar_data{$prj}{'warn'} = $fl[18];
						$tax_data{$fl[1]}++;
						$ar_sum++;
						$fasta_path = $fd."allbig_by_clade_extended/".$cl."/".$fl[0].".mfa";
						$size = `grep -v \"^>\" $fasta_path | wc -m`;
						`cp $fasta_path $genomes_tmp_folder`;
						chomp($size);
						$ar_data{$prj}{'size'} = int($size/1000);
						print FF $full_line{$cl}{$prj}."\n";
						@add_arr = split('; ', $fl[14]);
						for $add_s (@add_arr){
							@sfl = split(', ', $add_s);
							if ($sfl[1] =~ /^s__/){
								$sfl[1] =~ s/^s__//;
								$clade_genomes{$sfl[0]} = $sfl[1];
								if (exists($target{$cl})){
									$add_genomes{$sfl[0]} = $sfl[1];
								}
							}
						}
					}
					
				}
				close(F);
			}
			if (-f $datafd.$cl."/".$cl.".ar122.filtered.tsv"){
				open(F, $datafd.$cl."/".$cl.".ar122.filtered.tsv");
				while(<F>){
					if (/^CC/){
						chomp;
						@fl = split(/\t/, $_);
						($prj) = $fl[0] =~ /CC_(.+)_clade/;
						$filt_data{$prj} = $fl[1];
						$filt++;
					}
				}
				close(F);
			}

			#making clade specific tree and annotations
			if ((scalar(keys %bac_data) + scalar(keys %ar_data) + scalar(keys %clade_genomes)) > 2){
				for $id (keys %clade_genomes){
					$from = $dep{'GTDBF'}{'prog'}.$id."_genomic.fna.gz";
					`cp $from $genomes_tmp_folder`;
					$tmp_fn = $genomes_tmp_folder.$id."_genomic.fna.gz";
					`gunzip $tmp_fn`;
				}
				$flist = $datafd.$cl."_glist.txt";
				`ls $genomes_tmp_folder* > $flist`;
				@g_filelist = `cat $flist`;
				chomp(@g_filelist);
				$sketchfile = $datafd.$cl."_sketches";
				$runpr = $dep{'MASH'}{'prog'};
				$runopt = "sketch -l $flist -o $sketchfile -s 10000";
				`$runpr $runopt`;
				$sketchfile .= ".msh";
				$dist_res = $datafd.$cl."_distance_results.txt";
				open(OD, ">".$dist_res);
				for $g_fn (@g_filelist){
					$runpr = $dep{'MASH'}{'prog'};
					$runopt = "dist $sketchfile $g_fn";
					@res = `$runpr $runopt`;
					chomp(@res);
					for $line (@res){
						@fl = split(/\t/, $line);
						@sfl0 = split(/\//, $fl[0]);
						@sfl1 = split(/\//, $fl[1]);
						@ssfl0 = split(/\./, $sfl0[$#sfl0]);
						@ssfl1 = split(/\./, $sfl1[$#sfl1]);
						print OD $ssfl0[0]."\t".$ssfl1[0]."\t".$fl[2]."\n";
					}
				}
				close(OD);
				print STDERR "Constructing distance tree...\n";
				$jpath = $path."upgma.jl";
				`julia $jpath $dist_res`;
				$resfile = $ofd.$cl."_tree.txt";
				`mv pairs_upgma.tre $resfile`;
				`rm clade_text_annotation.txt`;
				open(OD, ">".$ofd.$cl."_tree_labels.txt");
				print OD "LABELS\nSEPARATOR TAB\nDATA\n";
				for $id (keys %clade_genomes){
					($tid, undef) = split(/\./, $id);
					print OD $tid."\t".$clade_genomes{$id}."\n";
				}
				close(OD);
				`rm $sketchfile`;
				`rm $dist_res`;
				`rm $flist`;
			}
			else{
				print STDERR $cl." has less then 3 assemblies, no tree created...\n";
			}

			open (O, ">".$ofd.$cl.".txt");
			print O "---- Final report for clade ".$cl." -----\n";
			print O "Location of the final assemblies: ".$fd."allbig_by_clade_extended/".$cl."/\n";
			if ($bac_sum > 0){
				print O "\n---- Identified bacterial genomes ----\n";
				print O "PRJ_ID : PRE_EXT_SIZE : FINAL_SIZE : AA_PERC : TAX_LEVEL : CLASSIFICATION\t(!WARNINGS!)\n";
				foreach $prj (sort keys %bac_data){
					print O sprintf('%-7s', $prj).": ".sprintf('%-13s', $pre_ext_size{$cl}{$prj}."kb").": ".sprintf('%-11s', $bac_data{$prj}{'size'}."kb").": ".sprintf('%-8s', $bac_data{$prj}{'aa_perc'}).": ".sprintf('%-10s', $bac_data{$prj}{'tax'}).": ".$bac_data{$prj}{'class'}."\t(!".$bac_data{$prj}{'warn'}."!)\n";
				}
			}
			if ($ar_sum > 0){
				print O "\n---- Identified archeal genomes ----\n";
				print O "PRJ_ID : PRE_EXT_SIZE : FINAL_SIZE : AA_PERC : TAX_LEVEL : CLASSIFICATION\t(!WARNINGS!)\n";
				foreach $prj (sort keys %ar_data){
					print O sprintf('%-7s', $prj).": ".sprintf('%-13s', $pre_ext_size{$cl}{$prj}."kb").": ".sprintf('%-11s', $ar_data{$prj}{'size'}."kb").": ".sprintf('%-8s', $ar_data{$prj}{'aa_perc'}).": ".sprintf('%-10s', $ar_data{$prj}{'tax'}).": ".$ar_data{$prj}{'class'}."\t(!".$ar_data{$prj}{'warn'}."!)\n";
				}
			}


			print O "\n---- Classification strings in the clade, ordered by prevalence ----\n";
			foreach $str (sort {$tax_data{$b} <=> $tax_data{$a}} keys %tax_data){
				print O $tax_data{$str}."\t".$str."\n";
			}

			if ($filt > 0){
				print O "\n---- Filtered assemblies due to low recognisable AA content ----\n";
				foreach $prj (sort keys %filt_data){
					print O $prj."\t".$filt_data{$prj}."\n";
				}
			}

			print O "\n---- AMR profiles of samples in the clade ----\n";
			@amr_list = ();
			foreach $gene (sort keys %{$cl_amr_pres{$cl}}){
				if ($cl_amr_pres{$cl}{$gene} > 0){
					$cnt = 0;
					for $prj (keys %filt_data){
						$cnt += $amr_data{$cl}{$prj}{$gene};
					}
					if ($cl_amr_pres{$cl}{$gene} > $cnt){
						push (@amr_list, $gene);
					}
				}
			}
			if (scalar(@amr_list) > 0){
				print O "PRJID";
				open(OD, ">".$ofd.$cl."_tree_amr_annot.txt");
				print OD "DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\tresfinder AMR\nFIELD_LABELS\t".join("\t", @amr_list)."\nFIELD_COLORS";
				for $gene (@amr_list){
					print O "\t".$gene;
					@tmp_col = @colours;
					for (@amr_list){
						print OD "\t".shift(@tmp_col);
						if (scalar(@tmp_col) == 0){
							@tmp_col = @colours;
						}
					}
				}
				print O "\n";
				print OD "\nFIELD_SHAPES";
				for (@amr_list){
					print OD "\t1";
				}
				print OD "\nDATA\n";
				foreach $prj (sort keys %{$amr_data{$cl}}){
					if (!exists($filt_data{$prj})){
						print O $prj;
						print OD $cc_full_names{$prj};
						for $gene (@amr_list){
							print O "\t".$amr_data{$cl}{$prj}{$gene};
							print OD "\t".$amr_data{$cl}{$prj}{$gene};
						}
						print O "\n";
						print OD "\n";
					}
				}
				for $id (keys %clade_genomes){
					($tid, undef) = split(/\./, $id);
					print OD $tid;
					for (@amr_list){
						print OD "\t0";
					}
					print OD "\n";
				}
				close(OD);
			}
			else{
				print O "No AMR genes detected in the clade!\n";
			}
			close(O);
			
		}
	}

	close(FF);

	`rm -rf $genomes_tmp_folder`;


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
