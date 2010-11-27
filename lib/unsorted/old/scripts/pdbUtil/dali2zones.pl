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
$query_id     = $opts{queryid};
$parent_id    = $opts{parentid};
$query_fasta  = $opts{queryfasta};
$parent_fasta = $opts{parentfasta};
$DALI_file    = $opts{dalifile};
$shift0       = $opts{shift0};
$shift1       = $opts{shift1};
$reverse      = $opts{reverse};

($query_id, $query_base, $query_chain, $query_dom)     = &parseID ($query_id);
($parent_id, $parent_base, $parent_chain, $parent_dom) = &parseID ($parent_id);

@DALI             = &fileBufArray ($DALI_file);

###############################################################################
# main
###############################################################################

# query fasta
if ($query_fasta) {
    @query_fasta_buf  = &fileBufArray ($query_fasta);
    foreach $line (@query_fasta_buf) {
	if ($line !~ /^\>/ && $line !~ /^\s*$/) {
	    @chars = split (/\s*/, $line);
	    shift (@chars)  if ($chars[0] eq '');
	    pop (@chars)  if ($chars[$#chars] eq '');
	    push (@query_fasta, @chars);
	}
    }
}

# parent fasta
if ($parent_fasta) {
    @parent_fasta_buf = &fileBufArray ($parent_fasta);
    foreach $line (@parent_fasta_buf) {
	if ($line !~ /^\>/ && $line !~ /^\s*$/) {
	    @chars = split (/\s*/, $line);
	    shift (@chars)  if ($chars[0] eq '');
	    pop (@chars)  if ($chars[$#chars] eq '');
	    push (@parent_fasta, @chars);
	}
    }
}


# parse DALI alignment
#
$len_aln = 0;
$equiv_res = 0;
@aln_queryToParent = ();
@aln_parentToQuery = ();

$started = undef;
foreach $line (@DALI) {
    if ($line =~ /^  NR. STRID1 STRID2      STRID1 \<=\> STRID2/) {
	$started = 'true';
	next;
    }
    next if (! $started);
    last if ($line =~ /^\s*$/);
	
# no longer including chain because dali is such a pain about it!
#	$query_match = ($query_chain eq '_') ? $query_id 
#	                                     : "$query_id-$query_chain";
#	$parent_match = ($parent_chain eq '_') ? $parent_id 
#	                                     : "$parent_id-$parent_chain";
    $line =~ s/^\s+|\s+$//g;
    @dali_dat = split (/\s+/, $line);
#   next if ($dali_dat[1] ne $query_match || $dali_dat[2] ne $parent_match)
    next if ((lc $dali_dat[1]) ne $query_base || (lc $dali_dat[2]) ne $parent_base)
	;
    $found = 'true';

    $query_start = substr ($line, 52, 4) - 1;
    $query_stop  = substr ($line, 64, 4) - 1;
    $parent_start = substr ($line, 78, 4) - 1;
    $parent_stop  = substr ($line, 90, 4) - 1;

    $query_start  += $shift0;
    $query_stop   += $shift0;
    $parent_start += $shift1;
    $parent_stop  += $shift1;

    for ($i=0; $i <= $query_stop - $query_start; ++$i) {
	$aln_queryToParent[$query_start+$i] = $parent_start+$i;
	++$len_aln;
	++$equiv_res  if ($query_fasta && $parent_fasta &&
			  $query_fasta[$query_start+$i] eq $parent_fasta[$parent_start+$i]);
    }
    for ($i=0; $i <= $parent_stop - $parent_start; ++$i) {
	$aln_parentToQuery[$parent_start+$i] = $query_start+$i;
    }
}

if (! $found) {
    print "NO DALI ALIGNMENT\n";
    exit 0;
}


# write identity
#
if ($query_fasta && $parent_fasta) {
    printf ("identity: %6.4f\n", ($equiv_res/$len_aln));
}


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


# parseID()
#
sub parseID {
    my $id = shift;
    my ($base, $chain, $domain) = (undef, undef, undef);

    $id =~ s/^(\w\w\w\w)$/$1.'_'/e;
    $id =~ s/^(\w\w\w\w\w)$/$1.'_'/e;
    $id =~ s/^(\w\w\w\w)(\w)(\w)$/(lc $1).(uc $2).(uc $3)/e;
    $base   = lc $1;
    $chain  = uc $2;
    $domain = uc $3;

    return ($id, $base, $chain, $domain);
}


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
\t -queryid      <queryid>
\t -parentid     <parentid>
\t -dalifile     <dalifile>
\t[-queryfasta   <queryfasta>]
\t[-parentfasta  <parentfasta>]
\t[-shift0       <shift0>]
\t[-shift1       <shift1>]
\t[-reverse      <T/F>]
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "queryid=s", "parentid=s", "queryfasta=s", "parentfasta=s", "dalifile=s", "shift0=s", "shift1=s", "reverse=s");

    # Check for legal invocation
    #
    if (! defined $opts{queryid} ||
	! defined $opts{parentid} ||
	! defined $opts{dalifile}
	) {
        print STDERR "$usage\n";
        exit -1;
    }

    # Check for existence
    #
    &checkExist ('f', $opts{dalifile});
    &checkExist ('f', $opts{queryfasta})  if (defined $opts{queryfasta});
    &checkExist ('f', $opts{parentfasta}) if (defined $opts{parentfasta});

    # defaults
    #
    $opts{shift0} = 0   if (! defined $opts{shift0});
    $opts{shift1} = 0   if (! defined $opts{shift1});


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

sub checkExist {
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
