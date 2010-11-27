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
# init
###############################################################################

$| = 1;

%opts = &getCommandLineOptions ();
$blast_file = $opts{blastfile};
$parentid   = $opts{parentid};

@blast_buf = &fileBufArray ($blast_file);

###############################################################################
# main
###############################################################################

# find place that desired parent came in
#
$place=-1;
for ($i=0, $order_i=1; $i <= $#blast_buf; ++$i) {
    if ($blast_buf[$i] =~ /^Sequences producing significant alignments/) {
	$started = 'true';
	++$i;
	next;
    }
    next if (! $started);
    if ($blast_buf[$i] =~ /^\s*$/) {
	$alignment_start_line = $i+1;
	last;
    }
    ($id, @rest) = split (/\s+/, $blast_buf[$i]);
    $id =~ s/^\>+//;
    if ($id eq $parentid) {
	$place = $order_i;
    }
    ++$order_i;
}
if ($place == -1) {
    print STDERR ("parent '$parentid' not found\n");
    exit 0;
} else {
    print qq{place: $place\n};
}


# parse blast alignment
#
$len_aln = 0;
$equiv_res = 0;
@aln_queryToParent = ();
@aln_parentToQuery = ();
for ($i=$alignment_start_line; $i <= $#blast_buf; ++$i) {
    next if ($blast_buf[$i] !~ /^\>/);
    ($id, @rest) = split (/\s+/, $blast_buf[$i]);
    $id =~ s/^\>+//;
    if ($id eq $parentid) {
	for ($j=$i+1; 
	     $blast_buf[$j] !~ /^\>/ && 
	     $blast_buf[$j] !~ /^\s*Database/; 
	     ++$j) {

	    if ($blast_buf[$j] =~ /^\s*Score/) {
		last if ($score_found);
		$score_found = 'true';
	    }
	    next if ($blast_buf[$j] !~ /^Query/ && $blast_buf[$j] !~ /^Sbjct/);
	    if ($blast_buf[$j] =~ /^Query/) {
		@q_align = ();
		@p_align = ();
		($tag, $q_startres, $align_str, $q_stopres) = split (/\s+/, $blast_buf[$j]);
		@q_align = split (//, $align_str);
	    }
	    elsif ($blast_buf[$j] =~ /^Sbjct/) {
		($tag, $p_startres, $align_str, $p_stopres) = split (/\s+/, $blast_buf[$j]);
		@p_align = split (//, $align_str);

		for ($k=0, $q_i=-1, $p_i=-1; $k <= $#q_align; ++$k) {
		    ++$q_i       if ($q_align[$k] =~ /[A-Z]/);
		    ++$p_i       if ($p_align[$k] =~ /[A-Z]/);
		    if ($q_align[$k] =~ /[A-Z]/ && $p_align[$k] =~ /[A-Z]/) {
			++$len_aln;
			++$equiv_res if ($q_align[$k] eq $p_align[$k]);
			$aln_queryToParent[$q_startres+$q_i-1] = $p_startres+$p_i-1;
			$aln_parentToQuery[$p_startres+$p_i-1] = $q_startres+$q_i-1;
		    }
		}
	    } else {
		&abort ("failure in control flow");
	    }
	}
	last;
    }
}
&abort ("unable to read PSIBLAST alignment for $blast_file\n")  if ($len_aln == 0);


# write identity
#
printf ("identity: %6.4f\n", ($equiv_res/$len_aln));


# find zones
#
for ($qi=0, $start_qi=undef, $last_pj=-1000; 
     $qi <= $#aln_queryToParent+1; 
     ++$qi) {
    $pj = $aln_queryToParent[$qi];
    if (! defined $pj || (defined $pj && $pj != $last_pj+1)) {
	if (defined $start_qi) {
	    ++$start_qi;
	    ++$start_pj;
	    ++$last_qi;
	    ++$last_pj;
	    if ($reverse) {
		#print qq{zone $start_pj-$last_pj:$start_qi-$last_qi\n};
		printf (qq{zone %4d-%-4d:%4d-%-4d\n}, $start_pj, $last_pj, $start_qi, $last_qi);
	    } else {
		#print qq{zone $start_qi-$last_qi:$start_pj-$last_pj\n};
		printf (qq{zone %4d-%-4d:%4d-%-4d\n}, $start_qi, $last_qi, $start_pj, $last_pj);
	    }
	}
	if (defined $pj) {
	    $start_qi = $qi;
	    $start_pj = $pj;
	} else {
	    $start_qi = undef;
	    $start_pj = undef;
	}
    }
    $last_qi = $qi;
    $last_pj = $pj  if (defined $pj);
}

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
    local $usage = qq{usage: $0
\t-blastfile   <blastfile>
\t-parentid    <parentid>
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "blastfile=s", "parentid=s");

    # Check for legal invocation
    #
    if (! defined $opts{blastfile} ||
	! defined $opts{parentid}
	) {
        print STDERR "$usage\n";
        exit -1;
    }

    # Check for existence
    #
    &checkExistence ('f', $opts{blastfile});

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
