#!/usr/bin/env perl 

#==================================================================================#
#                                   crud.pl                                        #
#==================================================================================#
#                                                                                  #
# This file is part of the AIRSS structure prediction package.                     #
#                                                                                  #
# AIRSS is free software; you can redistribute it and/or modify it under the terms #
# of the GNU General Public License version 2 as published by the Free Software    #
# Foundation.                                                                      #
#                                                                                  #
# This program is distributed in the hope that it will be useful, but WITHOUT ANY  #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  #
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.        #           
#                                                                                  #
# You should have received a copy of the GNU General Public License along with this#
# program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street,#                   
# Fifth Floor, Boston, MA  02110-1301, USA.                                        #
#                                                                                  #
#----------------------------------------------------------------------------------#
# The CASTEP Run Daemon (CRUD)                                                     #
#----------------------------------------------------------------------------------#
# Written by Chris Pickard, Copyright (c) 2005-2020                                #
#----------------------------------------------------------------------------------#
#                                                                                  #
#==================================================================================#


use strict;
use Getopt::Long;
use File::Copy;
use List::Util qw ( shuffle );
use Time::HiRes qw ( sleep );

sub usage {
  printf STDERR "Usage: crud.pl [-launch] [-exec] [-mpinp xx] [-repose] [-ramble] [-keep] [-nostop] [-cycle] Run many Castep computations\n";
  printf STDERR "       -exec           Use this executable\n";
  printf STDERR "       -launch         Use this parallel launcher\n";
  printf STDERR "       -mpinp          Number of cores per mpi Castep (0)\n";
  printf STDERR "       -repose         Use repose (0)\n";
  printf STDERR "       -ramble         Use ramble (0)\n";
  printf STDERR "       -keep           Keep all output files (0)\n";
  printf STDERR "       -nostop         Keep the script running (0)\n";
  printf STDERR "       -cycle          Retry failed runs (0)\n";
  exit();
}

my ($opt_launch,$opt_exec,$opt_mpinp,$opt_repose,$opt_ramble,$opt_keep,$opt_nostop,$opt_cycle,$opt_help) = ("","",0,0,0,0,0,0);

my $commandline = (join " ", @ARGV);

&GetOptions("launch=s"     => \$opt_launch,
	    "exec=s"       => \$opt_exec,
	    "mpinp=n"      => \$opt_mpinp,
	    "repose"       => \$opt_repose,
	    "ramble"       => \$opt_ramble,
	    "keep"         => \$opt_keep,
	    "nostop"       => \$opt_nostop,
	    "cycle"        => \$opt_cycle,
       "h|help"       => \$opt_help)|| &usage();

if ($opt_help) {
  &usage();
}
;

do {
  
  my $gotfile=1;
  my $file="";
  my @files = <hopper/*-*.res>;
  if (scalar @files == 0) { sleep 1 ; if ( $opt_nostop == 0) { exit } };
  foreach (shuffle(@files)) {
    $file = $_;
    my $stat=system("mv $file .");
    if ( $stat == 0) {
      print $file."\n";
      $gotfile = 0;
      last;
    }
  }  

  if ($gotfile == 0) {

    my @tmp = split('hopper/',$file);@tmp=split('\.',$tmp[1]);my $seed=$tmp[0];@tmp = split('-',$seed);my $root=$tmp[0];
    
    system("(cabal res cell < $seed.res ; sed -e '/^%BLOCK [Ll][Aa][Tt]*/, /^%ENDBLOCK [Ll][Aa][Tt]*/d' $root.cell | sed -e '/^%BLOCK [Pp][Oo][Ss]*/, /^%ENDBLOCK [Pp][Oo][Ss]*/d') > $seed.cell");
     
    my $executable= 'castep'; 

    if ( $opt_exec ne '' ) {$executable=$opt_exec};

    my $launcher='mpirun -np ';

    if ( $opt_launch ne '' ) {$launcher=$opt_launch};
    
    if ( $opt_mpinp > 0 ) { $executable= '"'.$launcher.$opt_mpinp.' '.$executable.'"' }
    
    #system("ln -sf $root.param $seed.param ; eval $executable $seed");
    system("cp $root.param $seed.param");
    system("cp $root.par $seed.par");
    
    if( ($opt_repose > 0) or ($opt_ramble > 0) ) {system("cp $root.ddp $seed.ddp 2>/dev/null ; cp $root.eddp $seed.eddp 2>/dev/null")} 
    
    # Ensure that the spin in the cell and param are consistent

    my $spintot=`grep -v "#SPIN=" $seed.cell | grep SPIN= | awk 'BEGIN {FS="SPIN="};{ sum += \$2 } END {printf "%10.3f",sum}'`;

    open  PARAMFILE, "$seed.param" or die $!;
    my @paramdata = <PARAMFILE>;
    close PARAMFILE;
      
    open  PARAMFILE, ">$seed.param" or die $!;
    foreach (@paramdata) {
      my @vec = split(' ',$_);
      if ( lc $vec[0] ne "spin" ) {
	print PARAMFILE $_;
      }
      ;
    }
    print PARAMFILE "spin : ".$spintot;
    close PARAMFILE;
      
    if ($opt_mpinp > 0 ) {
      if($opt_repose > 0) {
	      system("(repose_relax repose $opt_mpinp $seed)")
      } elsif($opt_ramble > 0){
      	system("(repose_relax ramble $opt_mpinp $seed)")
      } else {
      	system("eval $executable $seed");
      }
    } else {
      if($opt_repose > 0) {
      	system("repose_relax repose 1 $seed")
      }elsif($opt_ramble > 0){
         system("repose_relax ramble 1 $seed")  
      } else {
      	system("$executable $seed");
      }
    }
      
    if (-e $seed."-out.cell") {
	
      if (-e $seed.".dome_bin") {
	system("(echo 'task : dos';echo 'broadening : linear'; echo 'compute_band_gap : true')>$seed.odi;optados $seed;rm $seed.linear.dat;sed 's/legend 0.85, 0.8/legend 0.2, 0.8/g' $seed.linear.agr | sed 's/Electronic Density of States/$seed/g' > $seed.dos.agr ; (echo '\@target G0.s1'; echo '\@type xy' ; echo 'EfD 0'; echo 'EfD 5' ; echo '&' ;echo ' @    xaxis ticklabel font 4';echo ' @    yaxis ticklabel font 4') >> $seed.dos.agr");
      }
	
      # Construct the res file

      system("castep2res $seed > $seed.res");

      # Record the command line used
    
      system("sed -i 's/REM COMMAND_LINE/REM cmdline: crud.pl $commandline /g' $seed.res");
	
      # Add data to DOS plot
	
      if (-e $seed.".dos.agr") {
	my $efd='';
	if (-e $seed.".odo") {
	  chomp($efd = `grep EfD $seed.odo | tail -1 | awk '{print \$7}'`);
	}
	system("mv $seed.dos.agr $seed.dos.tmp.agr ; sed 's/EfD/$efd/g' $seed.dos.tmp.agr | sed 's/Generated by OptaDOS//g' > $seed.dos.agr ; rm $seed.dos.tmp.agr ")
      }
	  
    } else {
      chomp(my $err_count=`grep -c -m 1 electronic_minimisation $seed.*.err`);
      if ( (($err_count)&&($opt_cycle)) > 0 ) {
	system("mv $seed.* hopper")
      } else {
	system("mkdir -p bad_castep ; mv $seed.* bad_castep")
      }
    }
	
    if (! $opt_keep) {
      my @delete_files = <$seed.*>;
      foreach (@delete_files) {
	my @tmp = split('\.',$_); my $type=$tmp[1];
	if (($type ne "res")&&($type ne "cif")&&($type ne "magres")&&($type ne "castep")&&($type ne "odo")&&($type ne "dos")&&($type ne "den_fmt")) {
	  unlink($_);
	}	
      }
      unlink($seed."-out.cell");
    }
      

    system("mkdir -p good_castep ; mv -f $seed.* good_castep ;  mv -f $seed-out.cell good_castep 2> /dev/null")
	
    
  }
    
} until( -e "FINISH");
  
  
