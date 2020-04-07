use strict;
use warnings;

my$file = "complex_motif_more.fsa";
my$motiflength = 8;
#my$motiflength = 16;
if($#ARGV > -1){
	$file = $ARGV[0];
	if($#ARGV > 0){
		$motiflength = $ARGV[1];
	}
}

my$seqref = read_fsa($file);
my@seqs = @{$seqref};
print "Found ".scalar(@seqs)." sequences\n";



my$instanceref = get_random_instances(\@seqs,$motiflength);
my@instances = @{$instanceref};

print "Got random instances:\n";
print join("\n",@instances)."\n\n";
my$profileref = create_profile(\@instances);
my%profile = %{$profileref};
print "Starting profile:\n";
	my@nts = ('A','C','G','T');
	foreach my$nt (@nts){
		print $nt;
		foreach my$pos (0 .. $motiflength-1){
			print "\t".$profile{$pos}{$nt};
		}
		print "\n";
	}
	print "\n\n";

my$finalprofileref = recursive_random(\@instances);

sub recursive_random{
	my@instances = @{shift @_};

	my$oldtotal = check_solution(\@instances);

	for my$j (0 .. $#instances){

		my$leftout = int(rand($#instances));

		print "Leaving out sequence ".$leftout." : ".$instances[$leftout]."\n";

		my@traininstances;
		foreach my$i (0 .. $#instances){
			if($i == $leftout){
			} else{
				push(@traininstances,$instances[$i])
			}
		}

		my$profileref = create_profile(\@traininstances);
		my%profile = %{$profileref};

		my@leftseqs;
		push(@leftseqs,$seqs[$leftout]);
		my$newinstref = get_best_matches(\@leftseqs,\%profile,$motiflength);
		my@new_instances = @{$newinstref};

		print "New best instance:\n";
		foreach my$i (0 .. $#new_instances){
			print $new_instances[$i]."\n\n";
			$instances[$leftout] = $new_instances[$i];
		}

	}

	my$total = check_solution(\@instances);
	my$profileref = create_profile(\@instances);
	my%profile = %{$profileref};

	print "New solution: ".$total."\n";
	print "New profile:\n";
	my@nts = ('A','C','G','T');
	foreach my$nt (@nts){
		print $nt;
		foreach my$pos (0 .. $motiflength-1){
			print "\t".$profile{$pos}{$nt};
			#print "\t".int(exp(-$profile{$pos}{$nt})*(scalar(@instances)+4));
		}
		print "\n";
	}
	print "\n\n";

	if($total < $oldtotal){
		$profileref = recursive_random(\@instances);
	} else {
		print "Final profile:\n";
		my@nts = ('A','C','G','T');
		foreach my$nt (@nts){
			print $nt;
			foreach my$pos (0 .. $motiflength-1){
				#print "\t".$profile{$pos}{$nt};
				print "\t".int(exp(-$profile{$pos}{$nt})*(scalar(@instances)+4));
			}
			print "\n";
		}
		print "\n\n";
	}


	return $profileref;

}

sub check_solution{
	my@instances = @{shift @_};


	my$profileref = create_profile(\@instances);
	my%profile = %{$profileref};

	my@old_scores;
	foreach my$instance (@instances){
		push(@old_scores,score_profile($instance,\%profile));
	}
	my$oldtotal = 0;
	foreach my$score (@old_scores){
		$oldtotal += $score;
	}

	return $oldtotal;
}

sub read_fsa{
	my$file = shift @_;
	my@seqs;

	my$i = -1;
	open(FSA, $file) or die "Cannot open ".$file."\n";
	while(<FSA>){
		if(m/>/){
			$i++;
		} elsif (m/([ACTG]+)/){
			$seqs[$i] .= $1;
		}
	}
	close FSA;

	return \@seqs;

}

sub get_random_instances{
	my@seqs = @{shift @_};
	my$motiflength = shift @_;

	my@instances;
	foreach my$seq (@seqs) {
		my$pos = int(rand(length($seq)-$motiflength));
		push(@instances,substr($seq,$pos,$motiflength))
	}

	return \@instances;
}

sub create_profile{
	my@instances = @{shift @_};

	my$profilelength = length($instances[0]);

	my$pseudocount = 1;
	my@nts = ('A','C','G','T');

	my%profile;
	for my$i (0 .. $profilelength-1){
		foreach my$nt (@nts){
			$profile{$i}{$nt} = $pseudocount;
			$profile{$i}{'all'}++;
		}
		foreach my$seq (@instances){
			$profile{$i}{substr($seq,$i,1)}++;
			my$temp=substr($seq,$i,1);
			$profile{$i}{'all'}++;
		}
		foreach my$nt (@nts){
			$profile{$i}{$nt} = -log($profile{$i}{$nt}/$profile{$i}{'all'});
			$profile{$i}{$nt} = int($profile{$i}{$nt}*100)/100
		}
	}

	return \%profile;

}

sub score_profile{
	my$seq = shift @_;
	my%profile = %{shift @_};

	my$score=0;
	my@seqpos = split("",$seq);

	for my$i (0 .. $#seqpos){
		$score += $profile{$i}{$seqpos[$i]};
	}

	return $score;
}

sub get_best_matches{
	my@seqs = @{shift @_};
	my$profileref = shift @_;
	my$motiflength = shift @_;

	my@best_matches;
	foreach my$seq (@seqs){
		my$bestscore = 999;
		my$bestmatch = '';
		for my$i (0 .. length($seq)-1-$motiflength){
			my$dnaseq = substr($seq,$i,$motiflength);
			my$score = score_profile($dnaseq,$profileref);

			if($score < $bestscore){
				$bestscore = $score;
				$bestmatch = $dnaseq;
			}
		}
		push(@best_matches,$bestmatch);
	}

	return \@best_matches;
}