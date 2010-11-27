#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.3 $
##  $Date: 2004/05/24 19:45:13 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering

$man_method_str  = 
qq{
ROSETTA provides both ab initio and
comparative models of protein domains. It
uses the ROSETTA fragment insertion method
[Simons et al. J Mol Biol 1997;268:209-225].
Comparative models are built from structures
detected by PSI-BLAST, FFAS03, or 3DJury-A1
and aligned by the K*SYNC alignment method.
Loop regions are assembled from fragments and 
optimized to fit the aligned template structure.
SS predictions are accomplished by the JUFO
ss prediction method utilizing ROSETTA decoys.
JUFO_3D:  artificial neural net
  input:  amino acid property profile, psipred
          scoring matrix and low resolution 3D
          structural model(s)
 output:  three-state probabilities
};

$auto_method_str =
qq{
ROBETTA provides both ab initio and
comparative models of protein domains. It
uses the ROSETTA fragment insertion method
[Simons et al. J Mol Biol 1997;268:209-225].
Comparative models are built from structures
detected by PSI-BLAST, FFAS03, or 3DJury-A1
and aligned by the K*SYNC alignment method.
Loop regions are assembled from fragments and 
optimized to fit the aligned template structure.
SS predictions are accomplished by the JUFO
ss prediction method utilizing ROSETTA decoys.
JUFO_3D:  artificial neural net
  input:  amino acid property profile, psipred
          scoring matrix and low resolution 3D
          structural model(s)
 output:  three-state probabilities
The procedure is fully automated.
};

###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
$conf             = $opts{conf};
$target_name      = $opts{targetname};
$model_num        = $opts{modelnum};
$cutdefs_file     = $opts{cutdefs};
$domainfiles_list = $opts{domainfileslist};
$parents_list     = $opts{parentslist};
$out_file         = $opts{outfile};

$outbuf = '';
@domainfiles = split (/,/, $domainfiles_list);

###############################################################################
# main
###############################################################################

# get group id
$group_id = ($conf =~ /^M/i) 
              ? '4793-7654-3148'   # MANUAL: BAKER (CASP-6)
              #: '6750-1618-6576';  # AUTO:   BAKER-ROBETTA (CASP-6)
              : 'BAKER-ROBETTA';  # AUTO:   BAKER-ROBETTA (CASP-6)


# get method desc
$method_str = ($conf =~ /^M/i) 
              ? $man_method_str    # MANUAL: BAKER
              : $auto_method_str;  # AUTO:   BAKER-ROBETTA
$method_str =~ s/^\s+|\s+$//g;


# read cutdefs
@cutdefs = &fileBufArray ($cutdefs_file);
for ($line_i=1; $line_i <= $#cutdefs; ++$line_i) {
    $domain_i = $line_i - 1;
    $cutdefs[$line_i] =~ s!^\s+|\s+$!!g;
    ($q_beg[$domain_i], 
     $q_end[$domain_i], 
     $q_len[$domain_i], 
     $m_beg[$domain_i], 
     $m_end[$domain_i], 
     $m_len[$domain_i], 
     $p_beg[$domain_i], 
     $p_end[$domain_i],
     $p_id[$domain_i],
     $confidence[$domain_i],
     $seq[$domain_i]) = split (/\s+/, $cutdefs[$line_i]);
     
     $parents[$domain_i] = (! $parents_list && $p_id[$domain_i] !~ /^n\/?a$/i)
                             ? $p_id[$domain_i]
                             : 'N/A';
}


# get parents
if ($parents_list) {
    @cmdline_parents = split (/,/, $parents_list);
    for ($domain_i=0; $domain_i <= $#cmdline_parents; ++$domain_i) {
	if ($cmdline_parents[$domain_i]) {
	    @these_parents = split (/:/, $cmdline_parents[$domain_i]);
	    for ($parent_i=0; $parent_i <= $#these_parents; ++$parent_i) {
		$these_parents[$parent_i] =~ s/^(\w\w\w\w)(\w).*$/$1.'_'.$2/e;
		$these_parents[$parent_i] =~ s/__$//;
	    }
	    $parents[$domain_i] = join (' ', @these_parents);
	}
    }
}


# start outbuf
$outbuf = qq{PFRMAT SS\nTARGET $target_name\nAUTHOR $group_id\n};
foreach $line (split (/\n/, $method_str)) {
    $outbuf .= qq{METHOD $line\n};
}            
$outbuf .= qq{MODEL  $model_num\n};


# assemble ss preds
for ($domain_i=0; $domain_i <= $#domainfiles; ++$domain_i) {
#    $outbuf .= qq{PARENT $parents[$domain_i]\n};  

    $global_shift = $q_beg[$domain_i];
    $local_shift  = $q_beg[$domain_i] - $m_beg[$domain_i] + 1; 
    $start_res    = $q_beg[$domain_i] - $m_beg[$domain_i] + 1;
    $stop_res     = $q_end[$domain_i] - $m_beg[$domain_i] + 1; 

    # body
    $ss_buf = &bigFileBufArray ($domainfiles[$domain_i]);
    $started = undef;
    $model_resnum = 0;
    for ($i=0; $i <= $#{$ss_buf}; ++$i) {
	last if ($ss_buf->[$i] =~ /^END/);
	if ($ss_buf->[$i] =~ /^MODEL/) {
	    $started = 'true';
	    next;
	}
	if ($ss_buf->[$i] =~ /^[a-zA-Z] [HEC] [\d\.]+$/) {
	    $started = 'true';
	}
	next if (! $started);
	
	++$model_resnum;
	($res, $ss, $confidence) = split (/\s+/, $ss_buf->[$i]);
	    
	# adjust residue numbering if necessary, and skip extra coordinates
	next if ($model_resnum < $start_res || $model_resnum > $stop_res);

	# append
	$outbuf .= sprintf ("%1s %1s %4.2f\n", $res, $ss, $confidence);
    }
    #$outbuf .= qq{TER\n};
}
# add HETATMS?
$outbuf .= qq{END\n};


# output
if ($out_file) {
    open (OUT, '>'.$out_file);
    select (OUT);
}
print $outbuf;
if ($out_file) {
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
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{usage: $0
\t -targetname      <target_name>
\t -modelnum        <model_num>
\t -cutdefs         <cutdefs_file>
\t -domainfileslist <domainfile1,domainfile2,...>
\t[-conf            <MANUAL/AUTOMATIC>]               (def: MANUAL)
\t[-parentslist     <parent1,parent2a:parent2b,...>]  (def: reads cutdefs)
\t[-outfile         <out_file>]                       (def: STDOUT)
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "targetname=s",
		 "modelnum=i",
		 "cutdefs=s",
		 "domainfileslist=s",
		 "conf=s",
		 "parentslist=s",
		 "outfile=s"
		 );

    # Check for legal invocation
    if (! defined $opts{targetname} ||
	! defined $opts{modelnum} ||
	! defined $opts{cutdefs} ||
	! defined $opts{domainfileslist}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }

    # Defaults
    $opts{conf} = 'MANUAL'  if (! defined $opts{conf});

    # Existence checks
    &checkExist ('f', $opts{cutdefs});

    return %opts;
}

###############################################################################
# util
###############################################################################

# insertSortIndexList ()
#
sub insertSortIndexList {
    my ($val_list, $direction) = @_;
    my $index_list = +[];
    my ($index, $val, $i, $i2, $assigned);

    $index_list->[0] = 0;
    for ($index=1; $index <= $#{$val_list}; ++$index) {
        $assigned = undef;
        $val = $val_list->[$index];
        for ($i=0; $i <= $#{$index_list}; ++$i) {
            if ($direction eq 'decreasing') {
                if ($val > $val_list->[$index_list->[$i]]) {
                    for ($i2=$#{$index_list}; $i2 >= $i; --$i2) {
                        $index_list->[$i2+1] = $index_list->[$i2];
                    }
                    $index_list->[$i] = $index;
                    $assigned = 'true';
                    last;
                }
            }
            else {
                if ($val < $val_list->[$index_list->[$i]]) {
                    for ($i2=$#{$index_list}; $i2 >= $i; --$i2) {
                        $index_list->[$i2+1] = $index_list->[$i2];
                    }
                    $index_list->[$i] = $index;
                    $assigned = 'true';
                    last;
                }
            }
        }
        $index_list->[$#{$index_list}+1] = $index  if (! $assigned);
    }
    return $index_list;
}

# readFiles
#
sub readFiles {
    my ($dir, $fullpath_flag) = shift;
    my $inode;
    my @inodes = ();
    my @files = ();
    
    opendir (DIR, $dir);
    @inodes = readdir (DIR);
    closedir (DIR);
    foreach $inode (@inodes) {
	next if (! -f "$dir/$inode");
	next if ($inode =~ /^\./);
	push (@files, ($fullpath_flag) ? "$dir/$inode" : "$inode");
    }
    return @files;
}

# createDir
#
sub createDir {
    my $dir = shift;
    if (! -d $dir && (system (qq{mkdir -p $dir}) != 0)) {
	print STDERR "$0: unable to mkdir -p $dir\n";
	exit -2;
    }
    return $dir;
}

# copyFile
#
sub copyFile {
    my ($src, $dst) = @_;
    if (system (qq{cp $src $dst}) != 0) {
	print STDERR "$0: unable to cp $src $dst\n";
	exit -2;
    }
    return $dst;
}

# zip
#
sub zip {
    my $file = shift;
    if ($file =~ /^\.Z/ || $file =~ /\.gz/) {
	print STDERR "$0: ABORT: already a zipped file $file\n";
	exit -2;
    }
    if (system (qq{gzip -9 $file}) != 0) {
	print STDERR "$0: unable to gzip -9 $file\n";
	exit -2;
    }
    $file .= ".gz";
    return $file;
}

# unzip
#
sub unzip {
    my $file = shift;
    if ($file !~ /^\.Z/ && $file !~ /\.gz/) {
	print STDERR "$0: ABORT: not a zipped file $file\n";
	exit -2;
    }
    if (system (qq{gzip -d $file}) != 0) {
	print STDERR "$0: unable to gzip -d $file\n";
	exit -2;
    }
    $file =~ s/\.Z$|\.gz$//;
    return $file;
}

# remove
#
sub remove {
    my $inode = shift;
    if (system (qq{rm -rf $inode}) != 0) {
	print STDERR "$0: unable to rm -rf $inode\n";
	exit -2;
    }
    return $inode;
}
     
# runCmd
#
sub runCmd {
    my $cmd = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print "[$date]:$0: $cmd\n" if ($debug);
    &abort ("failure running cmd: $cmd\n") if (system ($cmd) != 0);
    #if (system ($cmd) != 0) {
    #    print STDERR ("failure running cmd: $cmd\n");
    #    return -2;
    #}
    return 0;
}

# logMsg()
#
sub logMsg {
    my ($msg, $logfile) = @_;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    else {
	select (STDERR);
    }
    print "[$date]:$0: $msg\n";
    if ($logfile) {
        close (LOGFILE);
    }
    select (STDOUT);

    return 'true';
}

# checkExist()
#
sub checkExist {
    my ($type, $path) = @_;
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
	elsif (! -s $path) {
            print STDERR "$0: emptyfile: $path\n";
            exit -3;
	}
    }
}

# abort()
#
sub abort {
    my $msg = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  chomp $date;
    print STDERR "[$date]:$0: $msg\n";
    exit -2;
}

# writeBufToFile()
#
sub writeBufToFile {
    my ($file, $bufptr) = @_;
    if (! open (FILE, '>'.$file)) {
	&abort ("$0: unable to open file $file for writing");
    }
    print FILE join ("\n", @{$bufptr}), "\n";
    close (FILE);
    return;
}

# fileBufString()
#
sub fileBufString {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

# fileBufArray()
#
sub fileBufArray {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

# bigFileBufArray()
#
sub bigFileBufArray {
    my $file = shift;
    my $buf = +[];
    if ($file =~ /\.gz|\.Z/) {
        if (! open (FILE, "gzip -dc $file |")) {
            &abort ("$0: unable to open file $file for gzip -dc");
        }
    }
    elsif (! open (FILE, $file)) {
        &abort ("$0: unable to open file $file for reading");
    }
    while (<FILE>) {
        chomp;
        push (@$buf, $_);
    }
    close (FILE);
    return $buf;
}     

###############################################################################
# end
1;                                                     # in case it's a package
###############################################################################
