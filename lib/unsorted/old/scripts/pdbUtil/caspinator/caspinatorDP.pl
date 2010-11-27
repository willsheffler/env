#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its 
##   disclosure does not constitute publication.  All rights are reserved by 
##   University of Washington, the Baker Lab, and Dylan Chivian, except those 
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.11 $
##  $Date: 2004/06/08 22:31:57 $
##  $Author: dylan $
##
###############################################################################


###############################################################################
# conf
###############################################################################

$| = 1;                                              # disable stdout buffering

# paths
if (! defined $ENV{'BAKER_HOME'}) {
    &abort ("must define \$BAKER_HOME in environment");
}
if (! defined $ENV{'SCRATCH_DIR'}) {
    &abort ("must define \$SCRATCH_DIR in environment");
}
$src_dir = $ENV{'BAKER_HOME'}."/src";
$dat_dir = $ENV{'BAKER_HOME'}."/dat";
$tmp_dir = (-d $ENV{'SCRATCH_DIR'}."/tmp") ? $ENV{'SCRATCH_DIR'}."/tmp" : "/tmp";

$slicePdb                = "$src_dir/pdbUtil/slicePdb.pl";
$getTaylorParseConsensus = "$src_dir/mstraUtil/getTaylorParseConsensus.pl";

# vars
$min_ai_domain_to_parse = 160;


# CASP group config
#
$man_group_id  = '4793-7654-3148';  # MANUAL: BAKER (CASP-6)
#$auto_group_id = '2692-3823-6511';  # AUTO:   BAKER-ROBETTA-GINZU (CASP-6)
$auto_group_id = 'BAKER-ROBETTA-GINZU';  # AUTO:   BAKER-ROBETTA-GINZU (CASP-6)


# Methods description config
#
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

Domain Predictions by GINZU method, based on
scanning for homologs with structures, followed
by HMM-PFAM and then by PSI-BLAST multiple
sequence clusters.  ROBETTA model subsequently
parsed with consensus variant of Taylor's
method.
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
The procedure is fully automated.

Domain Predictions by GINZU method, based on
scanning for homologs with structures, followed
by HMM-PFAM and then by PSI-BLAST multiple
sequence clusters.  ROBETTA model subsequently
parsed with consensus variant of Taylor's
method.
};

###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
$fasta_file        = $opts{fastafile};
$conf              = $opts{conf};
$target_name       = $opts{targetname};
$cutdefs_file      = $opts{cutdefs};
$domainfiles_list  = $opts{domainfileslist};
$taylor_cons_parse = $opts{taylorconsparse};
$parse_ai          = $opts{parseai};
$out_file          = $opts{outfile};
$inet_host         = $opts{inethost};
#$model_num         = $opts{modelnum};
$model_num         = 1;

@domainfiles = split (/,/, $domainfiles_list);
foreach $domainfile (@domainfiles) {
     &checkExist ('f', $domainfile);
}


# out path
$out_dir = $fasta_file;
$out_dir =~ s!^(.*)/[^/]+$!$1!;
$out_dir = '.'  if ($out_dir eq $fasta_file);
if ($out_dir !~ /^\//) {
     $cwd = `pwd`;  chomp $cwd;
     $out_dir = $cwd.'/'.$out_dir;
}

# zero out structs
$outbuf = '';
$outrows = ();

###############################################################################
# main
###############################################################################

# get group id
$group_id = ($conf =~ /^M/i) 
              ? $man_group_id
              : $auto_group_id;


# get method desc
$method_str = ($conf =~ /^M/i) 
              ? $man_method_str
              : $auto_method_str;
$method_str =~ s/^\s+|\s+$//g;


# clean target name
$target_name =~ s/^t(\d\d\d\d)$/T$1/;
$target_name =~ s/^t(\d\d\d)_?$/T0$1/i;


# read fasta
@fasta = ();
foreach $line (&fileBufArray ($fasta_file)) {
    next if ($line =~ /^\s*>/);
    $line =~ s/\s+//g;
    push (@fasta, split (//, $line));
}
$outrows[0] = join ('', @fasta);


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
     $p_src[$domain_i],
     $seq[$domain_i]) = split (/\s+/, $cutdefs[$line_i]);

     $parents[$domain_i] = (! $parents_list && 
			    $p_id[$domain_i] !~ /^n\/?a$/i &&
			    $p_src[$domain_i] !~ /^pfam$/i)
                            ? $p_id[$domain_i]
                             : 'N/A';

     # until rich's confidence function is here
     $confs[$domain_i] = (! $parents_list && 
			  $p_id[$domain_i] !~ /^n\/?a$/i &&
			  $p_src[$domain_i] !~ /^pfam$/i)
                           ? sprintf ("%-4.2f", 1-(1/(1+$confidence[$domain_i])))
                           : '0.50';

    # retain pfam family in p_src 
    if ($p_src[$domain_i] =~ /^pfam$/) {
	$p_src[$domain_i] .= ":$p_id[$domain_i]";
    }
}
# define linkers
for ($domain_i=1; $domain_i <= $#m_beg; ++$domain_i) {
    for ($res_n=$m_beg[$domain_i]; $res_n <= $m_end[$domain_i-1]; ++$res_n) {
	$linker[$res_n] = 'true';
    }
}


# start outbuf
#
$outbuf = qq{PFRMAT DP\nTARGET $target_name\nAUTHOR $group_id\n};
foreach $line (split (/\n/, $method_str)) {
    $outbuf .= qq{METHOD    $line\n};
}            
$outbuf .= qq{MODEL     $model_num\n};


# add domain predictions to output
#
$domain_num = 0;
for ($domain_i=0; $domain_i <= $#q_beg; ++$domain_i) {
    
    $outbuf .= qq{REMARK PARENT $parents[$domain_i]\n};
    $outbuf .= qq{REMARK PARSRC $p_src[$domain_i]\n};  
    $outbuf .= qq{REMARK SCORE  $confs[$domain_i]\n};  


    # either straight from cuts file
    #    (ai domains that too short)
    #
    if (! $taylor_cons_parse ||
	($parents[$domain_i] =~ /^n\/?a$/i && 
	 (! $parse_ai ||
	 $q_end[$domain_i] - $q_beg[$domain_i] - 1 < $min_ai_domain_to_parse
	  )
	 )
	) {

	++$domain_num;
	for ($res_i=$q_beg[$domain_i]-1; $res_i <= $q_end[$domain_i]-1; ++$res_i) {
	    $res_n = $res_i + 1;	
	    $confidence = $confs[$domain_i];
	    $confidence -= 0.25  if ($linker[$res_n]);
	    $outbuf .= join ("\t", $res_n, $fasta[$res_i], $domain_num, $confidence) . "\n";
	    $outrows[1] .= &base36 ($domain_num);
	}
    }

    # or cut body to consensus taylor parse each piece
    #
    else {
    
	@domain_num_mask = ();
        #$pdb_buf = &bigFileBufArray ($domainfiles[0]);
    
	$q_range = $q_beg[$domain_i].'-'.$q_end[$domain_i];
	$domain_dir = "$out_dir/model_dom_parse/t000_.$q_range";
	system (qq{mkdir -p $domain_dir})  if (! -f $domain_dir);

	$domain_slice_file      = "$domain_dir/t000_.pdb";
	&runCmd (qq{$slicePdb $domainfiles[0] _ $q_beg[$domain_i] $q_end[$domain_i] $domain_slice_file FALSE RENUM SEQUENTIAL});

	chdir $domain_dir;

	$domain_parse_info_file = "$domain_dir/t000_.taylor_cons.doms";
	if (! -s $domain_parse_info_file) {
	    $inet_host_opt = (defined $inet_host) ? "-inethost $inet_host" : "";
	    &runCmd (qq{$getTaylorParseConsensus -id t000_ -pdbfile $domain_slice_file -outdir $domain_dir $inet_host_opt});
	}
	
	foreach $domain_parse_info_line (&fileBufArray ($domain_parse_info_file)) {
	    next if ($domain_parse_info_line =~ /^DOMAIN/i);
	    ++$domain_num;
	    ($domain_id, $coord_range, $seq_range) = split (/\s+/, $domain_parse_info_line);
	    $seq_range =~ s/[a-zA-Z]//g;
	    @ranges = split (/,/, $seq_range);
	    foreach $range (@ranges) {
		($start_n, $stop_n) = split (/:/, $range);
		for ($dom_res_i=0; $dom_res_i <= $stop_n-$start_n; ++$dom_res_i) {
		    $domain_num_mask[$start_n-1+$dom_res_i] = $domain_num;
		}
	    }
	}

	# add to outbuf
	#
	for ($dom_res_i=0; $dom_res_i <= $q_end[$domain_i]-$q_beg[$domain_i]; ++$dom_res_i) {
	    $res_i = $q_beg[$domain_i]-1+$dom_res_i;
	    $res_n = $res_i + 1;
	    $confidence = $confs[$domain_i];
	    $confidence -= 0.25  if ($linker[$res_n]);
	    $outbuf .= join ("\t", $res_n, $fasta[$res_i], $domain_num_mask[$dom_res_i], $confidence) . "\n";
	    $outrows[1] .= &base36 ($domain_num_mask[$dom_res_i]);
	}

	chdir $out_dir;
	#&runCmd (qq{rm -rf $domain_dir});
    }
}
$outbuf .= qq{END\n};


# output
#
if ($out_file) {
    open (OUT, '>'.$out_file);
    select (OUT);
}
print $outbuf;
if ($out_file) {
    close (OUT);
    select (STDOUT);
}
print STDERR join ("\n", 
    "================================================",
    "= CaspinatorDP Ginzu+Robetta domain assignment =", 
    "================================================",
    @outrows)
    . "\n";


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
\t -fastafile        <fasta_file>
\t -targetname       <target_name>
\t -cutdefs          <cutdefs_file>
\t[-domainfileslist  <domainfile1,domainfile2,...>]  (def: undef)
\t[-taylorconsparse  <TRUE/FALSE>]                   (def: FALSE)
\t[-parseai          <TRUE/FALSE>]                   (def: FALSE)
\t[-conf             <MANUAL/AUTOMATIC>]             (def: MANUAL)
\t[-inethost         <inet_host>]                    (def: undef)
\t[-outfile          <out_file>]                     (def: STDOUT)
};
#\t -modelnum        <model_num>

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "fastafile=s",
		 "targetname=s",
		 "cutdefs=s",
		 "domainfileslist=s",
		 "taylorconsparse=s",
		 "parseai=s",
		 "conf=s",
		 "inethost=s",
		 "outfile=s"
		 );
#		 "modelnum=i",

    # Check for legal invocation
    if (! defined $opts{fastafile} ||
	! defined $opts{targetname} ||
#	! defined $opts{modelnum} ||
	! defined $opts{cutdefs} 
        ) {
        print STDERR "$usage\n";
        exit -1;
    }

    # Defaults
    $opts{conf}            = 'MANUAL'  if (! defined $opts{conf});
    $opts{taylorconsparse} = undef     if ($opts{taylorconsparse} =~ /^F/i);
    $opts{parseai}         = undef     if ($opts{parseai} =~ /^F/i);

    # complex checks
    if (defined $opts{taylorconsparse} && ! defined $opts{domainfileslist}) {
	&abort ("must provide model with -domainfileslist if you wish to -taylorconsparse TRUE");
    }

    # Existence checks
    &checkExist ('f', $opts{fastafile});
    &checkExist ('f', $opts{cutdefs});

    return %opts;
}

###############################################################################
# util
###############################################################################

# base36 ()
#
sub base36 {
    my $d = shift;
    return $d  if ($d =~ /^\d$/);
    return chr ($d - 10 + ord ('A')); 
}


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
