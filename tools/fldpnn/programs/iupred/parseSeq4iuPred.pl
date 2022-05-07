#!/usr/bin/env perl

use strict;
use warnings;

my  @sequence_file  = @ARGV; 

my $rdir = 'seqData/';
my $iupredPath = '/data0/drugtarget/iupred/';
$ENV{'IUPred_PATH'} = $iupredPath;
my %sequences;
my %seqHeads;
my($id, $info);
mkdir( $rdir) unless -e $rdir;

foreach my $sequence_file (@sequence_file){
	open ( my $input, '<', $sequence_file ) or die $!;
	{
		while ( my $stLine = <$input> ) {
			chomp( $stLine);
			if ($stLine =~ /^>+/) {
				if( $stLine =~ /\|/){
					$id = (split/\|/,$stLine)[1];
					$info = $stLine;
				}else{
					$id = substr( $stLine, 1);
					$info = $stLine;
				}
				$seqHeads{$id} .= $info;
			}else{
				$sequences{$id} .= $stLine;
			}
		}
	}

	close( $sequence_file );
}

chdir( $rdir);

foreach my $id (keys %sequences) {
	if ( $sequences{$id} ) {
		open( my $hFileSingleFasta, '>' . 'seq-' . $id . '.fasta') 
			|| die "fail open file ${id}.fasta for write $!\n";
		print $hFileSingleFasta "$seqHeads{$id}\n";
		print $hFileSingleFasta "$sequences{$id}\n";
		close $hFileSingleFasta;
		system( $iupredPath . 'iupred seq-' . $id . '.fasta long > iup-' . $id . '-long.txt 2>&1 /dev/null');
		system( $iupredPath . 'iupred seq-' . $id . '.fasta short > iup-' . $id . '-short.txt 2>&1 /dev/null');
		system( $iupredPath . 'iupred seq-' . $id . '.fasta glob > iup-' . $id . '-glob.txt 2>&1 /dev/null');
	}
}
