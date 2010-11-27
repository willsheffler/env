#!/usr/bin/perl 
#
# cmd - 6/1/04
# Script to interconvert hydrogen names between three different formats:
# PDB: format in use for pdb files (current cvs rosetta convention)
# IUP: IUPAC format for naming hydrogens 
# ROS: The former rosetta convention for naming hydrogens
#      (This was IUP style for all protons except LEU HB)

# usage: convertHydrogen.pl option
# options are PDB2ROS, PDB2IUP, IUP2ROS, IUP2PDB, ROS2PDB, ROS2IUP

#All pdb files in the current directory are automatically converted
#Newly created files will be placed into a subdirectory "modified_files"

die 
"Converts hydrogen names between different formats.
All pdb files in the current directory are automatically converted.
Newly created files are placed into a subdirectory 'modified_files'

Usage:  convertHydrogen.pl option
   options: PDB2ROS, PDB2IUP, IUP2ROS, IUP2PDB, ROS2PDB or ROS2IUP
       PDB: format in use for PDB files (RCSB standard, current rosetta convention)
       IUP: IUPAC hydrogen names, RCSB field conventions
       ROS: old rosetta convention (IUP style except for LEU HB)
\n" if (@ARGV != 1);
$OPTION = $ARGV[0];

$OPTION_IN = substr($OPTION, 0, 3);
$OPTION_OUT = substr($OPTION, 4, 3);

use vars qw(@AA_NAME,@PDB_PRO,@IUP_PRO,@ROS_PRO,@NUM_PRO,@INPUT_CHECK);

$AA_NAME[1] = "ALA";
$AA_NAME[2] = "CYS";
$AA_NAME[3] = "ASP";
$AA_NAME[4] = "GLU";
$AA_NAME[5] = "PHE";
$AA_NAME[6] = "GLY";
$AA_NAME[7] = "HIS";
$AA_NAME[8] = "ILE";
$AA_NAME[9] = "LYS";
$AA_NAME[10] = "LEU";
$AA_NAME[11] = "MET";
$AA_NAME[12] = "ASN";
$AA_NAME[13] = "PRO";
$AA_NAME[14] = "GLN";
$AA_NAME[15] = "ARG";
$AA_NAME[16] = "SER";
$AA_NAME[17] = "THR";
$AA_NAME[18] = "VAL";
$AA_NAME[19] = "TRP";
$AA_NAME[20] = "TYR";

# side chain proton names indexed by aa#

# PDB names
#CYS
$PDB_PRO[21] = "1HB ";
$PDB_PRO[22] = "2HB ";

#ASP
$PDB_PRO[31] = "1HB ";
$PDB_PRO[32] = "2HB ";

#GLU
$PDB_PRO[41] = "1HB ";
$PDB_PRO[42] = "2HB ";
$PDB_PRO[43] = "1HG ";
$PDB_PRO[44] = "2HG ";

#PHE
$PDB_PRO[51] = "1HB ";
$PDB_PRO[52] = "2HB ";

#GLY
$PDB_PRO[61] = "1HA ";
$PDB_PRO[62] = "2HA ";

#HIS
$PDB_PRO[71] = "1HB ";
$PDB_PRO[72] = "2HB ";

#ILE
$PDB_PRO[81] = "1HG1";
$PDB_PRO[82] = "2HG1";

#LYS
$PDB_PRO[91] = "1HB ";
$PDB_PRO[92] = "2HB ";
$PDB_PRO[93] = "1HG ";
$PDB_PRO[94] = "2HG ";
$PDB_PRO[95] = "1HD ";
$PDB_PRO[96] = "2HD ";
$PDB_PRO[97] = "1HE ";
$PDB_PRO[98] = "2HE ";

#LEU
$PDB_PRO[101] = "1HB ";
$PDB_PRO[102] = "2HB ";

#MET
$PDB_PRO[111] = "1HB ";
$PDB_PRO[112] = "2HB ";
$PDB_PRO[113] = "1HG ";
$PDB_PRO[114] = "2HG ";

#ASN
$PDB_PRO[121] = "1HB ";
$PDB_PRO[122] = "2HB ";

#PRO
$PDB_PRO[131] = "1HD ";
$PDB_PRO[132] = "2HD ";
$PDB_PRO[133] = "1HG ";
$PDB_PRO[134] = "2HG ";
$PDB_PRO[135] = "1HB ";
$PDB_PRO[136] = "2HB ";

#GLN
$PDB_PRO[141] = "1HB ";
$PDB_PRO[142] = "2HB ";

#ARG
$PDB_PRO[151] = "1HB ";
$PDB_PRO[152] = "2HB ";
$PDB_PRO[153] = "1HG ";
$PDB_PRO[154] = "2HG ";
$PDB_PRO[155] = "1HD ";
$PDB_PRO[156] = "2HD ";

#SER
$PDB_PRO[161] = "1HB ";
$PDB_PRO[162] = "2HB ";

#TRP
$PDB_PRO[191] = "1HB ";
$PDB_PRO[192] = "2HB ";

#TYR
$PDB_PRO[201] = "1HB ";
$PDB_PRO[202] = "2HB ";

#
# IUPAC names
#CYS
$IUP_PRO[21] = "2HB ";
$IUP_PRO[22] = "3HB ";

#ASP
$IUP_PRO[31] = "2HB ";
$IUP_PRO[32] = "3HB ";

#GLU
$IUP_PRO[41] = "2HB ";
$IUP_PRO[42] = "3HB ";
$IUP_PRO[43] = "2HG ";
$IUP_PRO[44] = "3HG ";

#PHE
$IUP_PRO[51] = "2HB ";
$IUP_PRO[52] = "3HB ";

#GLY
$IUP_PRO[61] = "2HA ";
$IUP_PRO[62] = "3HA ";

#HIS
$IUP_PRO[71] = "2HB ";
$IUP_PRO[72] = "3HB ";

#ILE
$IUP_PRO[81] = "2HG1";
$IUP_PRO[82] = "3HG1";

#LYS
$IUP_PRO[91] = "2HB ";
$IUP_PRO[92] = "3HB ";
$IUP_PRO[93] = "2HG ";
$IUP_PRO[94] = "3HG ";
$IUP_PRO[95] = "2HD ";
$IUP_PRO[96] = "3HD ";
$IUP_PRO[97] = "2HE ";
$IUP_PRO[98] = "3HE ";

#LEU
$IUP_PRO[101] = "2HB ";
$IUP_PRO[102] = "3HB ";

#MET
$IUP_PRO[111] = "2HB ";
$IUP_PRO[112] = "3HB ";
$IUP_PRO[113] = "2HG ";
$IUP_PRO[114] = "3HG ";

#ASN
$IUP_PRO[121] = "2HB ";
$IUP_PRO[122] = "3HB ";

#PRO
$IUP_PRO[131] = "2HD ";
$IUP_PRO[132] = "3HD ";
$IUP_PRO[133] = "2HG ";
$IUP_PRO[134] = "3HG ";
$IUP_PRO[135] = "2HB ";
$IUP_PRO[136] = "3HB ";

#GLN
$IUP_PRO[141] = "2HB ";
$IUP_PRO[142] = "3HB ";

#ARG
$IUP_PRO[151] = "2HB ";
$IUP_PRO[152] = "3HB ";
$IUP_PRO[153] = "2HG ";
$IUP_PRO[154] = "3HG ";
$IUP_PRO[155] = "2HD ";
$IUP_PRO[156] = "3HD ";

#SER
$IUP_PRO[161] = "2HB ";
$IUP_PRO[162] = "3HB ";

#TRP
$IUP_PRO[191] = "2HB ";
$IUP_PRO[192] = "3HB ";

#TYR
$IUP_PRO[201] = "2HB ";
$IUP_PRO[202] = "3HB ";

#
# Rosetta names (old)
#CYS
$ROS_PRO[21] = "2HB ";
$ROS_PRO[22] = "3HB ";

#ASP
$ROS_PRO[31] = "2HB ";
$ROS_PRO[32] = "3HB ";

#GLU
$ROS_PRO[41] = "2HB ";
$ROS_PRO[42] = "3HB ";
$ROS_PRO[43] = "2HG ";
$ROS_PRO[44] = "3HG ";

#PHE
$ROS_PRO[51] = "2HB ";
$ROS_PRO[52] = "3HB ";

#GLY
$ROS_PRO[61] = "2HA ";
$ROS_PRO[62] = "3HA ";

#HIS
$ROS_PRO[71] = "2HB ";
$ROS_PRO[72] = "3HB ";

#ILE
$ROS_PRO[81] = "2HG1";
$ROS_PRO[82] = "3HG1";

#LYS
$ROS_PRO[91] = "2HB ";
$ROS_PRO[92] = "3HB ";
$ROS_PRO[93] = "2HG ";
$ROS_PRO[94] = "3HG ";
$ROS_PRO[95] = "2HD ";
$ROS_PRO[96] = "3HD ";
$ROS_PRO[97] = "2HE ";
$ROS_PRO[98] = "3HE ";

#LEU
$ROS_PRO[101] = "1HB ";
$ROS_PRO[102] = "2HB ";

#MET
$ROS_PRO[111] = "2HB ";
$ROS_PRO[112] = "3HB ";
$ROS_PRO[113] = "2HG ";
$ROS_PRO[114] = "3HG ";

#ASN
$ROS_PRO[121] = "2HB ";
$ROS_PRO[122] = "3HB ";

#PRO
$ROS_PRO[131] = "2HD ";
$ROS_PRO[132] = "3HD ";
$ROS_PRO[133] = "2HG ";
$ROS_PRO[134] = "3HG ";
$ROS_PRO[135] = "2HB ";
$ROS_PRO[136] = "3HB ";

#GLN
$ROS_PRO[141] = "2HB ";
$ROS_PRO[142] = "3HB ";

#ARG
$ROS_PRO[151] = "2HB ";
$ROS_PRO[152] = "3HB ";
$ROS_PRO[153] = "2HG ";
$ROS_PRO[154] = "3HG ";
$ROS_PRO[155] = "2HD ";
$ROS_PRO[156] = "3HD ";

#SER
$ROS_PRO[161] = "2HB ";
$ROS_PRO[162] = "3HB ";

#TRP
$ROS_PRO[191] = "2HB ";
$ROS_PRO[192] = "3HB ";

#TYR
$ROS_PRO[201] = "2HB ";
$ROS_PRO[202] = "3HB ";

#Number of protons to be modified per amino acid
$NUM_PRO[1] = 0;
$NUM_PRO[2] = 2;
$NUM_PRO[3] = 2;
$NUM_PRO[4] = 4;
$NUM_PRO[5] = 2;
$NUM_PRO[6] = 2;
$NUM_PRO[7] = 2;
$NUM_PRO[8] = 2;
$NUM_PRO[9] = 8;
$NUM_PRO[10] = 2;
$NUM_PRO[11] = 4;
$NUM_PRO[12] = 2;
$NUM_PRO[13] = 6;
$NUM_PRO[14] = 2;
$NUM_PRO[15] = 6;
$NUM_PRO[16] = 2;
$NUM_PRO[17] = 0;
$NUM_PRO[18] = 0;
$NUM_PRO[19] = 2;
$NUM_PRO[20] = 2;


$OUTPUT_DIRECTORY  = "modified_files";
`mkdir $OUTPUT_DIRECTORY` unless (-d $OUTPUT_DIRECTORY);
$INPUT_DIR = ".";
opendir (DIR, $INPUT_DIR) or die "could not open $INPUT_DIR\n";
@FILES = readdir (DIR);
foreach $FILE (@FILES) {
	$_ = $FILE;
	if (/pdb$/) {
		@INPUT_CHECK = 0;
		open (INFILE, $FILE) or die "could not open $FILE\n";
		$OUTFILE = "temp";
		open (OUTFILE, ">$OUTFILE") or die "could not open $OUTFILE\n";
		while (<INFILE>) {
			if (/^ATOM/) {
				$RESTYPE = substr ($_,17,3);
				$ATMNAME = substr ($_,12,4);
				$ATM = substr($ATMNAME,1,1);
				$FRONT = substr($_,0,12);
				$BACK = substr($_,16);
				if ($ATM eq "H") {
					if ($OPTION_IN eq "ROS") {
						$INDEX = ROS_IN($RESTYPE,$ATMNAME);
						}
					elsif ($OPTION_IN eq "PDB") {
						$INDEX = PDB_IN($RESTYPE,$ATMNAME);
						}
					elsif ($OPTION_IN eq "IUP") {
						$INDEX = IUP_IN($RESTYPE,$ATMNAME);
						}
					else {
						die "invalid input option specified\n";
						}
					if ($INDEX == 999) {
						$ATMOUT = $ATMNAME;
						}
					elsif ($OPTION_OUT eq "ROS" ) {
						$ATMOUT = $ROS_PRO[$INDEX];
						}
					elsif ($OPTION_OUT eq "IUP") {
						$ATMOUT = $IUP_PRO[$INDEX];
						}
					elsif ($OPTION_OUT eq "PDB") {
						$ATMOUT = $PDB_PRO[$INDEX];
						}
					else {
						die "invalid output option specified\n";
						}
					printf OUTFILE "%s%s%s", $FRONT, $ATMOUT, $BACK;
					}
				else {
					print OUTFILE "$_";
					}
				}
			else {
				print OUTFILE "$_";
				}
			}
		unless ($INPUT_CHECK[1] == 1 && $INPUT_CHECK[2] == 1) {
			print "WARNING\n";
			print "WARNING you have invoked $OPTION_IN as your input option but the file $FILE does not have the \n";
			print "WARNING characteristic hydrogens for this type of file\n";
			print "WARNING there are several reasons this may be true:\n";
			print "WARNING 1) the file may have no hydrogens to rename\n";
			print "WARNING 2) the protein does not contain at least 1 ASP and 1 LEU\n";
			print "WARNING 3) you have selected the wrong input option for the type of file you are converting\n";
			print "WARNING if #3 is the reason then your modified files are likely to contain errors\n";
			}
		close OUTFILE;
		close INFILE;
 		`mv $OUTFILE $OUTPUT_DIRECTORY/$FILE`;
		}
	}
						
sub PDB_IN {
	my @LIST = @_;
	my $AMINO_ACID = $LIST[0];
	my $HATOM = $LIST[1];
	my $AANUM = 0;
	my $INDEX_L = 999;
	for ($i=1;$i<=20;$i++) {
		if ($AMINO_ACID eq $AA_NAME[$i]) {
			$AANUM = $i;
			last;
			}
		}
	if ($AANUM == 0) {
		print "Abnormal amino acid type $AMINO_ACID found in $FILE\n";
		return $INDEX_L;
		}
	my $AA_INDEX = $AANUM * 10;
	if ($AANUM == 3) {
		if ($HATOM eq $PDB_PRO[31]) {
			$INPUT_CHECK[1] = 1;
			}
		}
	if ($AANUM == 10) {
		if ($HATOM eq $PDB_PRO[101]) {
			$INPUT_CHECK[2] = 1;
			}
		}
	unless ($NUM_PRO[$AANUM] == 0) {
		for ($i=1;$i<=$NUM_PRO[$AANUM];$i++) {
			my $HINDEX = $AA_INDEX + $i;
			if ($HATOM eq $PDB_PRO[$HINDEX]) {
				return $HINDEX;
				}
			}
		}
	return $INDEX_L;
	}

sub IUP_IN {
	my @LIST = @_;
	my $AMINO_ACID = $LIST[0];
	my $HATOM = $LIST[1];
	my $AANUM = 0;
	my $INDEX_L = 999;
	for ($i=1;$i<=20;$i++) {
		if ($AMINO_ACID eq $AA_NAME[$i]) {
			$AANUM = $i;
			last;
			}
		}
	if ($AANUM == 0) {
		print "Abnormal amino acid type $AMINO_ACID found in $FILE\n";
		return $INDEX_L;
		}
	my $AA_INDEX = $AANUM * 10;
	if ($AANUM == 3) {
		if ($HATOM eq $IUP_PRO[32]) {
			$INPUT_CHECK[1] = 1;
			}
		}
	if ($AANUM == 10) {
		if ($HATOM eq $IUP_PRO[102]) {
			$INPUT_CHECK[2] = 1;
			}
		}
	unless ($NUM_PRO[$AANUM] == 0) {
		for ($i=1;$i<=$NUM_PRO[$AANUM];$i++) {
			my $HINDEX = $AA_INDEX + $i;
			if ($HATOM eq $IUP_PRO[$HINDEX]) {
				return $HINDEX;
				}
			}
		}
	return $INDEX_L;
	}

sub ROS_IN {
	my @LIST = @_;
	my $AMINO_ACID = $LIST[0];
	my $HATOM = $LIST[1];
	my $AANUM = 0;
	my $INDEX_L = 999;
	for ($i=1;$i<=20;$i++) {
		if ($AMINO_ACID eq $AA_NAME[$i]) {
			$AANUM = $i;
			last;
			}
		}
	if ($AANUM == 0) {
		print "Abnormal amino acid type $AMINO_ACID found in $FILE\n";
		return $INDEX_L;
		}
	my $AA_INDEX = $AANUM * 10;
	if ($AANUM == 3) {
		if ($HATOM eq $ROS_PRO[32]) {
			$INPUT_CHECK[1] = 1;
			}
		}
	if ($AANUM == 10) {
		if ($HATOM eq $ROS_PRO[101]) {
			$INPUT_CHECK[2] = 1;
			}
		}
	unless ($NUM_PRO[$AANUM] == 0) {
		for ($i=1;$i<=$NUM_PRO[$AANUM];$i++) {
			my $HINDEX = $AA_INDEX + $i;
			if ($HATOM eq $ROS_PRO[$HINDEX]) {
				return $HINDEX;
				}
#			else {
#				print "cant find $HATOM is not $ROS_PRO[$HINDEX] from $HINDEX\n";
#				}
			}
		}
	return $INDEX_L;
	}












