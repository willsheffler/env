#!/usr/bin/perl
##
## Copyright 2002, University of Washington, the Baker Lab, and Dylan Chivian.
##   This document contains private and confidential information and its
##   disclosure does not constitute publication.  All rights are reserved by
##   University of Washington, the Baker Lab, and Dylan Chivian, except those
##   specifically granted by license.
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 1.2 $
##  $Date: 2004/07/15 23:04:57 $
##  $Author: dylan $
##
###############################################################################

###############################################################################
# base conf
###############################################################################
                   
# general
$| = 1;                                              # disable stdout buffering
$debug = 1;                                             # chatter while running

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

###############################################################################
# init
###############################################################################

if ($#ARGV != 1) {
    print STDERR "usage: $0 <query.pdb> <parent.pdb>\n";
    exit -1;
}
$query_pdb  = $ARGV[0];
$parent_pdb = $ARGV[1];

$query_id      = $query_pdb;
$query_id      =~ s!^(.*?)([^/]{4})([^/])[^/]*$!$2!;
$q_dir         = $1;
$q_ch          = $3;
$q_dir         = './'  if ($q_dir eq '');
$q_dir         =~ s!/$!!;
$query_fullid  = "$query_id$q_ch";

$parent_id     = $parent_pdb;
$parent_id     =~ s!^(.*?)([^/]{4})([^/])[^/]*$!$2!;
$p_dir         = $1;
$p_ch          = $3;
$p_dir         = './'  if ($p_dir eq '');
$p_dir         =~ s!/$!!;
$parent_fullid = "$parent_id$p_ch";

###############################################################################
# conf
###############################################################################

$DaliLite              = "$src_dir/shareware/Dali/DaliLite.dylan";
$ENV{DALI_SERVER_HOME} = $ENV{'SCRATCH_DIR'}."/DaliLite";
$ENV{DALIDATDIR}       = "$p_dir/DAT/";

$removeMissingDensity = "$src_dir/pdbUtil/removeMissingDensity.pl";
$removeChain          = "$src_dir/pdbUtil/removeChain.pl";

###############################################################################
# main
###############################################################################

# prepare dali format pdbs 
# (no -1 occupancy atoms, no chains, what the hell... no header!)
#
$query_pdb_dali = "$q_dir/$query_id.pdb";
$parent_pdb_dali = "$p_dir/$parent_id.pdb";
system ("$removeMissingDensity -p $query_pdb | grep '^ATOM\\|^HETATM' | sed 's/^\\(.\\{21\\}\\)./\\1 /' > $query_pdb_dali");
system ("$removeMissingDensity -p $parent_pdb | grep '^ATOM\\|^HETATM' | sed 's/^\\(.\\{21\\}\\)./\\1 /' > $parent_pdb_dali");


# run Dali
#
system ("$DaliLite -readbrk $query_pdb_dali $query_id > /dev/null 2>&1");
system ("$DaliLite -readbrk $parent_pdb_dali $parent_id > /dev/null 2>&1");
system ("$DaliLite -align $query_id $parent_id > /dev/null 2>&1");
system ("$DaliLite -fssp $query_id $query_id.dccp $query_fullid-$parent_fullid.dali > /dev/null 2>&1");

exit 0;

###############################################################################
# end
###############################################################################
