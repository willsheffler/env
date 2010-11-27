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
$CE_file     = $opts{cefile};
$shift0      = $opts{shift0};
$shift1      = $opts{shift1};
$reverse     = $opts{reverse};

@CE = &fileBufArray ($CE_file);

###############################################################################
# main
###############################################################################

# parse CE alignment
#
$len_aln = 0;
$equiv_res = 0;
@aln_queryToParent = ();
@aln_parentToQuery = ();
foreach $line (@CE) {
    next if ($line !~ /^([A-Z])\s+(\d+)\s+\w+\s+\w+\s+([A-Z])\s+(\d+)\s+\w+\s+\w+$/);
    $query_res      = $1;
    $query_res_num  = $2;
    $parent_res     = $3;                                       
    $parent_res_num = $4;                                       

    $query_res_num  += $shift0;
    $parent_res_num += $shift1;

    &abort ("unreal query_res_num: $query_res_num")   if ($query_res_num  < 1);
    &abort ("unreal parent_res_num: $parent_res_num") if ($parent_res_num < 1);

    $aln_queryToParent[$query_res_num-1] = $parent_res_num-1;
    $aln_parentToQuery[$parent_res_num-1] = $query_res_num-1;

    ++$len_aln;
    ++$equiv_res  if ($query_res eq $parent_res);
}
&abort ("unable to read CE alignment for $CE_file\n")  if ($len_aln == 0);


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
\t -cefile   <cefile>
\t[-shift0   <shift0>]
\t[-shift1   <shift1>]
\t[-reverse  <T/F>]
};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "cefile=s", "shift0=s", "shift1=s", "reverse=s");

    # Check for legal invocation
    #
    if (! defined $opts{cefile}) {
        print STDERR "$usage\n";
        exit -1;
    }

    # Check for existence
    #
    &checkExistence ('f', $opts{cefile});

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
