#!/usr/bin/perl

###############################################################################
# init
###############################################################################
$| = 1;

%opts = &getCommandLineOptions ();
$a_zones_file = $opts{azones};
$b_zones_file = $opts{bzones};
$length_from  = $opts{lengthfrom};

@a_zones = &fileBufArray ($a_zones_file);
@b_zones = &fileBufArray ($b_zones_file);

###############################################################################
# main
###############################################################################

# a_zones
#
$len_a = 0;
$max_a = -1;
@a_q2p_map = ();
foreach $line (@a_zones) {
    next if ($line !~ /^\s*zone\s*(\d+)\s*\-\s*(\d+)\s*\:\s*(\d+)\s*\-\s*(\d+)\s*$/i);
    last if ($line =~ /NO/i);
    $q_start = $1;
    $q_stop  = $2;                                       
    $p_start = $3;
    $p_stop  = $4;
    &abort ("ranges not same length in $a_zones_file: $line")  if ($q_stop-$q_start != $p_stop-$p_start);

    for ($i=0; $i <= $q_stop-$q_start; ++$i) {
	$a_q2p_map[$q_start-1+$i] = $p_start-1+$i;
	++$len_a;
	$max_a = $q_stop-1;
    }
}
print STDERR "no alignment in a_zones \n"  if ($max_a < 0);


# b_zones
#
$len_b = 0;
$max_b = -1;
@b_queryToParent = ();
foreach $line (@b_zones) {
    next if ($line !~ /^\s*zone\s*(\d+)\s*\-\s*(\d+)\s*\:\s*(\d+)\s*\-\s*(\d+)\s*$/i);
    last if ($line =~ /NO/i);
    $q_start = $1;
    $q_stop  = $2;                                       
    $p_start = $3;
    $p_stop  = $4;
    &abort ("ranges not same length in $b_zones_file: $line")  if ($q_stop-$q_start != $p_stop-$p_start);

    for ($i=0; $i <= $q_stop-$q_start; ++$i) {
	$b_q2p_map[$q_start-1+$i] = $p_start-1+$i;
	++$len_b;
	$max_b = $q_stop-1;
    }
}
print STDERR "no alignment in b_zones \n"  if ($max_b < 0);


# get values
#
$shift0 = 0;
$shift1 = 0;
$shift2 = 0;
$shift3 = 0;
$shift4 = 0;
$length = 0;
for ($i=0; $i <= (($max_a > $max_b)?$max_a:$max_b); ++$i) {
    next if (! defined $a_q2p_map[$i] && ! defined $b_q2p_map[$i]);
    ++$length  if (! defined $length_from || $length_from eq 'a|b');
    ++$length  if (defined $a_q2p_map[$i] && $length_from eq 'a');
    ++$length  if (defined $b_q2p_map[$i] && $length_from eq 'b');
    next if (! defined $a_q2p_map[$i] || ! defined $b_q2p_map[$i]);
    ++$shift0  if (abs($a_q2p_map[$i] - $b_q2p_map[$i]) == 0);
    ++$shift1  if (abs($a_q2p_map[$i] - $b_q2p_map[$i]) <= 1);
    ++$shift2  if (abs($a_q2p_map[$i] - $b_q2p_map[$i]) <= 2);
    ++$shift3  if (abs($a_q2p_map[$i] - $b_q2p_map[$i]) <= 3);
    ++$shift4  if (abs($a_q2p_map[$i] - $b_q2p_map[$i]) <= 4);
}
$shift0 /= $length;
$shift1 /= $length;
$shift2 /= $length;
$shift3 /= $length;
$shift4 /= $length;


# output
#
printf ("length: %5d+%-5d sh0: %6.4f sh1: %6.4f sh2: %6.4f sh3: %6.4f sh4: %6.4f\n",
	$len_a, $len_b,
	$shift0, $shift1, $shift2, $shift3, $shift4);


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
    local $usage = qq{usage: $0
\t -azones     <a_zones_file>
\t -bzones     <b_zones_file>
\t[-lengthfrom <a/b/a|b>]          (def: a|b (i.e. the union))
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "azones=s", "bzones=s", "lengthfrom=s");

    # Check for legal invocation
    #
    if (! defined $opts{azones} ||
	! defined $opts{bzones}
	) {
        print STDERR "$usage\n";
        exit -1;
    }

    # Defaults
    $opts{lengthfrom} = 'a|b'  if (! defined $opts{lengthfrom});

    # Check for existence
    #
    &checkExistence ('f', $opts{azones});
    &checkExistence ('f', $opts{bzones});

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
