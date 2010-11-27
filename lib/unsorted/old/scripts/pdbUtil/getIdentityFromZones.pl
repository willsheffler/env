#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its
##   disclosure does not constitute publication.  All rights are reserved by
##   University of Washington, the Baker Lab, and Dylan Chivian, except those
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.1.1.1 $
##  $Date: 2003/09/05 01:47:28 $
##  $Author: dylan $
##
###############################################################################
                   

###############################################################################
# conf
###############################################################################

$| = 1;

%res_int_code = ( 'A' => 0,
		  'C' => 1,
		  'D' => 2,
		  'E' => 3,
		  'F' => 4,
		  'G' => 5,
		  'H' => 6,
		  'I' => 7,
		  'K' => 8,
		  'L' => 9,
		  'M' => 10,
		  'N' => 11,
		  'P' => 12,
		  'Q' => 13,
		  'R' => 14,
		  'S' => 15,
		  'T' => 16,
		  'V' => 17,
		  'W' => 18,
		  'Y' => 19
);


$blosum62_table = +[
  [ 4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-2,-1,-1,-1, 1, 0, 0,-3,-2],
  [ 0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2],
  [-2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0,-1,-3,-4,-3],
  [-1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0,-1,-2,-3,-2],
  [-2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3],
  [ 0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3, 0,-2,-2,-2, 0,-2,-3,-2,-3],
  [-2,-3,-1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1,-2,-3,-2, 2],
  [-1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-1, 3,-3,-1],
  [-1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0,-1,-2,-3,-2],
  [-1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-1, 1,-2,-1],
  [-1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1, 1,-1,-1],
  [-2,-3, 1, 0,-3, 0, 1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2],
  [-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2, 7,-1,-2,-1,-1,-2,-4,-3],
  [-1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0,-1,-2,-2,-1],
  [-1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2],
  [ 1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2],
  [ 0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 1, 5, 0,-2,-2],
  [ 0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2, 0, 4,-3,-1],
  [-3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11, 2],
  [-2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7]
];

###############################################################################
# init
###############################################################################

%opts = &getCommandLineOptions ();
$q_fasta_file = $opts{qfastafile};
$p_fasta_file = $opts{pfastafile};
$zones_file   = $opts{zonesfile};

@q_fasta_buf = &fileBufArray ($q_fasta_file);
@p_fasta_buf = &fileBufArray ($p_fasta_file);
@zones_buf   = &fileBufArray ($zones_file);

###############################################################################
# main
###############################################################################

# q_fasta
#
@q_fasta = ();
foreach $line (@q_fasta_buf) {
    next if ($line =~ /^\s*\>/);
    $line =~ s/\s+//g;
    $line = uc $line;
    push (@q_fasta, split (//, $line));
}
$q_len = $#q_fasta + 1;


# p_fasta
#
@p_fasta = ();
foreach $line (@p_fasta_buf) {
    next if ($line =~ /^\s*\>/);
    $line =~ s/\s+//g;
    $line = uc $line;
    push (@p_fasta, split (//, $line));
}
$p_len = $#p_fasta + 1;


# zones
#
$aln_len = 0;
@q2p_map = ();
foreach $line (@zones_buf) {
    next if ($line !~ /^\s*zone\s*(\d+)\s*\-\s*(\d+)\s*\:\s*(\d+)\s*\-\s*(\d+)\s*$/i);
    last if ($line =~ /NO/i);
    $q_start = $1;
    $q_stop  = $2;                                       
    $p_start = $3;
    $p_stop  = $4;
    &abort ("ranges not same length in $zones_file: $line")  if ($q_stop-$q_start != $p_stop-$p_start);

    for ($i=0; $i <= $q_stop-$q_start; ++$i) {
	$q2p_map[$q_start-1+$i] = $p_start-1+$i;
	++$aln_len;
    }
}
print STDERR "no alignment in zones file $zones_file\n"  if ($aln_len == 0);


# get identity and positives
#
$idents = 0;
$positives = 0;
for ($qi=0; $qi <= $#q2p_map; ++$qi) {
    next if (! defined $q2p_map[$qi]);
    $pj = $q2p_map[$qi];
    $q_res = $q_fasta[$qi];
    $p_res = $p_fasta[$pj];
    $q_int_code = $res_int_code{$q_res};
    $p_int_code = $res_int_code{$p_res};
    next if (! defined $q_int_code || ! defined $p_int_code);

    ++$idents if ($q_res eq $p_res);
    ++$positives if ($blosum62_table->[$q_int_code]->[$p_int_code] > 0);
}


# output
#
if ($outfile) {
    open (OUT, '>'.$outfile);
    select (OUT);
}
printf ("IDENT: %-6.4f POS: %-6.4f QUERY_COV: %-6.4f PARENT_COV: %-6.4f\n", 
	$idents/$aln_len,
	$positives/$aln_len,
	$aln_len/$q_len,
	$aln_len/$p_len);
if ($outfile) {
    close (OUT);
    select (STDOUT);
}


# done
exit 0;

###############################################################################
# subs
###############################################################################

# getCommandLineOptions()
#
#  desc: get the command line options
#
#  args: none
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    local $usage = qq{usage: $0 -qfastafile <query_fasta_file> -pfastafile <parent_fasta_file> -zonesfile <zones_file>};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "qfastafile=s", "pfastafile=s", "zonesfile=s");

    # Check for legal invocation
    #
    if (! defined $opts{qfastafile} ||
	! defined $opts{pfastafile} ||
	! defined $opts{zonesfile}
	) {
        print STDERR "$usage\n";
        exit -1;
    }

    # Check for existence
    #
    &checkExistence ('f', $opts{qfastafile});
    &checkExistence ('f', $opts{pfastafile});
    &checkExistence ('f', $opts{zonesfile});

    return %opts;
}

###############################################################################
# util
###############################################################################

sub logMsg {
    local ($msg, $logfile) = @_;

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    else {
	select (STDERR);
    }
    print $msg, "\n";
    if ($logfile) {
        close (LOGFILE);
    }
    select (STDOUT);

    return 'true';
}

sub checkExistence {
    local ($type, $path) = @_;
    if ($type eq 'd') {
	if (! -d $path) { 
            print STDERR "$0: dirnotfound: $path\n";
            exit -3;
	}
    }
    elsif ($type eq 'f') {
	if (! -f $path) {
            print STDERR "$0: filenotfound: $path\n";
            exit -3;
	}
    }
}

sub abort {
    local $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}

sub writeBufToFile {
    ($file, $bufptr) = @_;
    if (! open (FILE, '>'.$file)) {
	&abort ("$0: unable to open file $file for writing");
    }
    print FILE join ("\n", @{$bufptr}), "\n";
    close (FILE);
    return;
}

sub fileBufString {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

sub fileBufArray {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

###############################################################################
# end
###############################################################################
