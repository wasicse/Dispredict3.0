#!/usr/bin/perl 
$LIC = "Licensed to:  jing yan (Academic, university of alberta) located in Edmonton, Canada. \nThis license is for non-commercial use only.
 Please read our LICENSE file for details (http://protein.bio.unipd.it/LICENSE).";
#

# 
#-----------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#   espritz.pl   : Version 1.1
#
#   Description  : v1.0 Disorder based on Sequence and BRNN's
#   		 : v1.1 Disorder based on Sequence and psi-blast multiple sequence alignments with BRNN

#   Author       : Ian Walsh <ian.walsh@bio.unipd.it>
#
#   Date         : v1.0: Jan. 2011
#   		 : v1.1: June 2011
#
#   Arguments	 : In order on command line: 
#   				1. $workdir: the working directory where all there is a list of fasta files
#   				2. $model: 
#					"X" for Xray type predictions.
#					"D" for Disprot type predictions.
#					"N" for NMR type predictions.
#
#					"pX" for Xray type predictions with psi-blast.
#					"pD" for Disprot type predictions with psi-blast.
#					"pN" for NMR type predictions with psi-blast.

#   				3. $sw: 1 maximize Sw measure, 0 use 5% false positive rate
#------------------------------------------------------------------------------------

@timeDataStart = localtime(time);
join(' ', @timeDataStart); 


print "################################################################################################\n";
print "\n\n";
print "$LIC";
print "\nContact silvio.tosatto\@unipd.it for commercial licensing details.";
print "\n\n";
print "################################################################################################\n\n\n";
if ($#ARGV <2 ){
		die "usage: perl esprits.pl \n 
		1. workdir: the working directory where all there is a list of fasta files
		2. model: \n
		\"X\" for Xray type predictions. \n
		\"D\" for Disprot type predictions. \n
		\"N\" for NMR type predictions. \n
		\"pX\" for Xray type predictions with psi-blast. \n
		\"pD\" for Disprot type predictions with psi-blast. \n
		\"pN\" for NMR type predictions with psi-blast. \n
		3. Sw: 1 maximize Sw measure, 0 use 5% false positive rate
		\n\n\n
		If you remain unsure please see EXAMPLES.sh";
}


$workdir = shift @ARGV;
chomp($workdir);
#system("ls $workdir/batch/*.fasta > $workdir/filelist");  depends on kernel
opendir(DIR, "$workdir/");
@files = grep(/\.fasta$/,readdir(DIR));
closedir(DIR);
#$fnamelist = "$workdir/filelist";	#the list of files to predict

$model = shift @ARGV; #disprot or Xray
$sw = shift @ARGV; # use Sw thresholds or 5%FPR
$sw = chop($sw);

if (($sw==1)) {
	if ($model eq "D") {
		$thres = 0.2644;
	}
	elsif ($model eq "X") {
		$thres = 0.0634;
	}
	elsif ($model eq "N") {
		$thres = 0.1803;
	}
	elsif ($model eq "pD") {
		$thres = 0.2644;
	}
	elsif ($model eq "pX") {
		$thres = 0.058;
	}
	elsif ($model eq "pN") {
		$thres = 0.1741;
	}
	else {
		$thres = 0.0634; # default
	}

}
else {
	if ($model eq "D") {
		$thres = 0.5072;
	}
	elsif ($model eq "X") {
		$thres = 0.1434;
	}
	elsif ($model eq "N") {
		$thres = 0.3089;
	}
	elsif ($model eq "pD") {
		$thres = 0.5072;
	}
	elsif ($model eq "pX") {
		$thres = 0.1376;
	}
	elsif ($model eq "pN") {
		$thres = 0.2965;
	}
	else {
		$thres = 0.1434; # default
	}
}


print "working in directory : $workdir \n";
print "model : $model \n";
print "thres : $thres \n";

#open (list, "<$fnamelist");
#@filetext = <list>;
foreach $file (@files) {

	$fname = "$workdir/$file";
	chomp($fname);
	print "$fname\n";

	open (fi, "<$fname");
	@text = <fi>;
	close fi;
	$name = shift @text;
	$name =~s/>//g;
	$name =~s/\s//g;
	$seq = "";
	while(@text){ 
        	$seq  .= shift @text;
	        chomp  $seq;
	}
	$seq =~ y/BOJUZ/CXXXX/;

	if (length($name)>20) {$name = substr($name,0,20);$name.="\n";}

	open(fi2,">$fname.fasta");
	print fi2 "\> $name\n";
	print fi2 "$seq\n";
	close fi2;

	$len =length($seq);
	open (ft, ">$fname.disbin");
	print ft "1\n$len\n$seq\n";
	close ft;

	$JUNK = "$fname.disbin $fname.fasta";

	$tmp = $fname;
	$fname =~ s/\.fasta//g;
	$basename = $fname;
	$fname = $tmp;
	if ($model eq "X") {
		system("./bin/disbinX ./model_definition/ensembleCX $fname.disbin $fname.flatblast $thres\> $basename.espritz");
		print "Finished executing XRAY calpha disorder no psi-blast and threshold=$thres\n";
	}
	elsif ($model eq "D") {
		system("./bin/disbinD ./model_definition/ensembleD $fname.disbin $fname.flatblast $thres\> $basename.espritz");
		print "Finished executing DISPROT disorder no psi-blast and threshold=$thres\n";
	}
	elsif ($model eq "N") {
		system("./bin/disbinN ./model_definition/ensembleN $fname.disbin $fname.flatblast $thres\> $basename.espritz");
		print "Finished executing NMR calpha disorder no psi-blast and threshold=$thres\n";
	}
	elsif ($model eq "pX") {
		system("./align/getAlignments.pl $fname");
		system("./bin/disbin_psi ./model_definition/ensembleCX_psi $fname.disbin $fname.flatblast $thres\> $basename.espritz");
		print "Finished executing XRAY calpha disorder WITH psi-blast and threshold=$thres\n";
	}
	elsif ($model eq "pN") {
		system("./align/getAlignments.pl $fname");
		system("./bin/disbin_psi ./model_definition/ensembleN_psi $fname.disbin $fname.flatblast $thres\> $basename.espritz");
		print "Finished executing NMR disorder WITH psi-blast and threshold=$thres\n";
	}
	elsif ($model eq "pD"){
		system("./align/getAlignments.pl $fname");
		system("./bin/disbin_psi ./model_definition/ensembleD_psi $fname.disbin $fname.flatblast $thres\> $basename.espritz");
		print "Finished executing DISPROT disorder WITH psi-blast and threshold=$thres\n";
	}
	else {
		die "\nUnknown model definition\n";
	}

	system("rm $JUNK");

	# write the prediction in fasta format
 	$PRED = "$basename.espritz";
	open(PRED) or die("Could not open prediction file: $fname.espritz");

	$pred = "";
 	foreach $line (<PRED>) {
		chomp($line);
		@entries = split(/\t/, $line);

		if ($entries[0] eq "D" || $entries[0] eq "O") {
			$pred .= $entries[0];
		}
 	}
   	close (PRED);
	
	open(dis,">$basename.espritz.fasta");
	print dis "\> $name\n";
	print dis "$seq\n";
	print dis "$pred\n";
	close dis;

}
#close list;



@timeDataEnd = localtime(time);
join(' ', @timeDataEnd);

$hrs = $timeDataEnd[2]-$timeDataStart[2];
$min = $timeDataEnd[1]-$timeDataStart[1];
$sec = $timeDataEnd[0]-$timeDataStart[0];
print "\n\nTime for predictions to finish = $hrs hrs : $min mins : $sec secs\n"
