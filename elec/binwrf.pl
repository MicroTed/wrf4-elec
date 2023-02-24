#!/usr/bin/perl
# a perl program to read in a flash summary file and report the
#   ICs and CGs in time bins of whatever width. (hardwired value of intervalm, below,
#   normally set to 1 or 60 minute bins, but can be any integer value)
#
# Usage:
#  binwrf.pl rsl.out.0000
#  binwrf.pl rsl.out.0000 > output.txt
#
# Looks for the lines:
# $infile = 
#open(INPUT, "$run.hlgttyphsll") or die "file not found";
#open(INPUT, "linegrab") or die "file not found";
use Time::Local;

#open(INPUT, "$run.log") or die "file not found";

#open(OUT, ">$run.flashbin") or die "cannot create file";
 open(OUT, '>-');
$cr = chr(13);  # newline character for Mac
$lf = chr(10);  # newline character for Unix
#$/ = $lf;  # set record separator to Unix newline (default is chr(13) )
# foreach $file (@ARGV) {
#   open INPUT, "<$file";
#   $next = <INPUT>;
#   chomp($next);
#    @parts = split( / +/, $next);
#    print OUT ("$next \n");
#    }
# foreach ...
    $intervalminv = 1;
#    $intervalm = 60 ; # 1/$intervalminv;   # time discretization (usually 1 minute)
    $intervalm = 1 ; # 1/$intervalminv;   # time discretization (usually 1 minute)
    $intervals = $intervalm*60;  
#    print ("intervalminv,intervalm,intervals = $intervalminv, $intervalm, $intervals\n");
    $dxdy = 1; #1000*1000;
    $istop = 0;
    $binshift = 0;
    $tries = 0;
    $bintry = 0;
    $newbin = 1;
    $curbin = 0;
    $ic = 0;
    $cgn = 0;
    $cgp = 0;
    $icx = 0;
    $cgnx = 0;
    $cgpx = 0;
    $ictot = 0;
    $cgntot = 0;
    $cgptot = 0;
    $icdistot = 0;
    $icchargetot = 0;
    $cgchargetot = 0;
    $cgpchargetot = 0;
    $cgpchargetot = 0;
    $cgpcharge = 0;
    $cgncharge = 0;
    $triestot = 0;
    $addtime = 0;
    $addstep = 0;
    $nonindtotp = 0;
    $nonindtotn = 0;
    $indtot = 0;
    $nonindp = 0;
    $nonindn = 0;
    $ind = 0;
    $ghsn = 0;
    $ghsp = 0;
    $ghin = 0;
    $ghip = 0;
    $netchg = 0;
    $netchgp = 0;
    $netchgn = 0;
    $netchgcld = 0;
    $netchgpcld = 0;
    $netchgncld = 0;
    $irestart = 0;
    $dv = 1.0; 
    $irst2 = 0;
    $readfirst = 0;
    $dt = 4;
    $com2trmm = 0;
    $startbin = -1;
    $times = 0;
    $firsttime = 0;
    $yearadd = 0;
    $cwmass = 0;
    $rwmass = 0;
    $icemass = 0;
    $snowmass = 0;
    $grmass = 0;
    $hlmass = 0;
    $wvol5 = 0;
    $wvol10 = 0;
    $filecount = 0;

foreach $file (@ARGV) {

    open INPUT, "<$file";
    $start = 0;
    while ( (<INPUT>)  ) {
    if ( /^STOP HERE/ ) {
      goto "foo";
      }
      
    if ( $firsttime == 0 ) {
      $n = 0;
      if ( /^Timing for Writing/ ) {
        $firsttime = 1;
       chomp();
       @parts = split(/ +/);
       $n = 3;
#        print( "$parts[3] \n");
#        $wrfout = $parts[3];
#        @parts2 = split(/_/,$wrfout);
#        print( "parts2 = @parts2\n");
#        @yymmdd = split(/-/, $parts2[2]);
#        print( "yymmdd = @yymmdd\n");
#        @hhmmss = split(/:/, $parts2[3]);
#        print( "hhmmss = @hhmmss\n");
#        $time = timelocal($hhmmss[2],$hhmmss[1],$hhmmss[0],$yymmdd[2],$yymmdd[1],$yymmdd[0]);
#        print( "time = $time\n");

       
       } elsif ( /^ RESTART run: / ) {
        $firsttime = 1;
       chomp();
       @parts = split(/ +/);
       $n = 4;
# 	 print( "$parts[3] \n");
# 	 $wrfout = $parts[4];
# 	 @parts2 = split(/_/,$wrfout);
# 	 print( "parts2 = @parts2\n");
# 	 @yymmdd = split(/-/, $parts2[2]);
# 	 print( "yymmdd = @yymmdd\n");
# 	 @hhmmss = split(/:/, $parts2[3]);
# 	 print( "hhmmss = @hhmmss\n");
# 	 $timestart = timelocal($hhmmss[2],$hhmmss[1],$hhmmss[0],$yymmdd[2],$yymmdd[1],$yymmdd[0]);
# 	 print( "timestart = $timestart\n");
        
       }
       if ( $n > 0 ) {
     #  print( "parts = @parts\n");
       $wrfout = $parts[$n];
       @parts2 = split(/_/,$wrfout);
     #  print( "parts2 = @parts2\n");
       @yymmdd = split(/-/, $parts2[2]);
     #  print( "yymmdd = @yymmdd\n");
       @hhmmss = split(/:/, $parts2[3]);
     #  print( "hhmmss = @hhmmss\n");
       if ($yymmdd[0] > 1900) { 
         $yearadd = -1900; # - $yymmdd[0]; 
      #   print("yearadd = $yearadd ",$yymmdd[0]+$yearadd,"\n");
         }
       $timestart = timelocal($hhmmss[2],$hhmmss[1],$hhmmss[0],$yymmdd[2],$yymmdd[1]-1,$yymmdd[0]+$yearadd);
      # print( "timestart = $timestart\n");
       $time1 = $timestart;
       }
    
     }
    if ( /^ TIME STEP  / ) {
      chomp();
      @parts = split(/ +/);
      $dt = $parts[4];
#      print("dt = $dt\n");
#      $next = <INPUT>;
#      chomp($next);
#      @parts = split(/ +/, $next);
       }
    if ( /^ NX / ) {
      chomp();
      @parts = split(/ +/);
      $dx = $parts[8];
      $next = <INPUT>;
      $next = <INPUT>;
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $dy = $parts[8];
      $next = <INPUT>;
      $next = <INPUT>;
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $dz = $parts[8];
 #     print("dx = $dx\n");
 #     print("dy = $dy\n");
 #     print("dz = $dz\n");
      $dv = 1.0; #$dx*$dy*$dz;
      }
    if ( / SIMULATION STOP  TIME/ && $start == 0 ) {
      chomp();
      @parts = split( / +/ );
      $nstop = $parts[$#parts];
#      print("tstop = $nstop\n");
#      print("nstop=$nstop\n");
      }
    if ( /^  COMMAS:   RESTART/ ) {
      $irestart = 1;
      $irst2 = 1;
     }


    if ( /^Timing for main:/ ) {
    @parts = split(/ +/);
#    print( "$parts[3] \n");
    $wrfout = $parts[4];
    @parts2 = split(/_/,$wrfout);
#    print( "parts2 = @parts2\n");
    @yymmdd = split(/-/, $parts2[0]);
#    print( "yymmdd = @yymmdd\n");
    @hhmmss = split(/:/, $parts2[1]);
#    print( "hhmmss = @hhmmss\n");
     # range for month is 0-11!!!!
    $time2 = timelocal($hhmmss[2],$hhmmss[1],$hhmmss[0],$yymmdd[2],$yymmdd[1]-1,$yymmdd[0]+$yearadd);
#    print( "time2 = $time2\n");
     $dt = $time2 - $time1;
     $time1 = $time2;
     $times = $time2 - $timestart; # $times + $dt;
#     print( "time = $times\n" );
     $min = int(($times-$dt)/$intervals);
     $newbin = $min + 1;
     #print("min,curbin,newbin = $min,$curbin,$newbin\n");
    
     }
     
    if ( /^END OF STEP: nstep, time/ ) {
      @parts = split(/ +/);
      $nsteptst = $parts[6];
#      print ("end of step: $nsteptst, $nstop, $min, $curbin, $binshift\n");
      if ( $nsteptst == $nstop && $curbin == 1 ) {
       $min = $min+1;
       $newbin = $min - $binshift;
       }
#       if ( $irst2 == 1 ) {  print ("HERE I AM: $newbin, $curbin\n"); }
     }
    if ( /^NSTEPXX,/ ) {
     $next = <INPUT>;
     chomp($next);
     @parts = split(/ +/, $next);
     $nstep = $parts[1];
     $nstop = $parts[3];
     $times = $parts[4];
#     $min = int(($times-$dt)/60.0);
     $min = int(($times-$dt)/$intervals);
#      if ( $irst2 == 1 ) {
#       print ("binshift,curbin,min = $binshift, $curbin, $min\n");
#       }
     if ($binshift == 0) {
      $binshift = $min - 1;
      $curbin = $min - $binshift;
#       if ( $irst2 >= 0 ) {  print ("HERE I AM 2: $newbin, $curbin\n"); }
#      @min = ($min);
     }
     $newbin = $min - $binshift;
      if ( $irst2 == 1 ) {
 #      print ("newbin,curbin,min,binshift = $newbin, $curbin,$min,$binshift\n");
       }
 #     if ( $times >= 10790 ) {
 #     print OUT ("$times $min $curbin $newbin $cgp $ic \n");
 #     }
     }

#    if (  /^ Integration Done/ ) {
    if (  /SUCCESS COMPLETE WRF/ ) {
     #    print("end of file\n" );
         $min = int(($times)/$intervals); 
         $tmp = (($times)/$intervals);
  #       print ("check Integration Done: times,min = $times,$min,$tmp,$newbin,$curbin \n");
  # second part of 'if' test added to account for "off" restart times that
  # can result from EnKF timing (oldobsmethod)
      if ( $newbin == $curbin &&  int(($times)/$intervals) == (($times)/$intervals)) { 
         $newbin = $curbin + 1; 
         $min = int(($times)/$intervals); 
         $tmp = (($times)/$intervals);
  #       print ("Integration Done: times,min = $times,$min,$tmp \n");
         }
     }
    if (  $newbin > $curbin) {

     $bintime[$curbin] = $min*$intervalm;
#     print OUT ("curbin, bintime, min, newbin = $curbin, $bintime[$curbin], $min, $times, $newbin\n");

     $ictot = $ictot + $ic;
     $ic[$curbin] = $ic;
     $icpermin[$curbin] = $ic/$intervalm;
     $ic = 0;
     
     $cgptot = $cgptot + $cgp;
     $cgp[$curbin] = $cgp;
     $cgp = 0;
     
     $cgntot = $cgntot + $cgn;
     $cgn[$curbin] = $cgn;
     $cgn = 0;
     
     $triestot = $triestot + $tries;
     $tries[$curbin] = $tries;
     $tries = 0;
     
     $icdistot = $icdistot + $icdis; # *$dv;
     $icdis[$curbin] = $icdis; # *$dv;
     $icdis = 0;

     $icchargetot = $icchargetot + $iccharge; #*$dv;
     $iccharge[$curbin] = $iccharge; #*$dv;
     $iccharge = 0;
     
     $cgpchargetot = $cgpchargetot + $cgpcharge; # *$dv;
     $cgpcharge[$curbin] = $cgpcharge; # *$dv;
     $cgpcharge = 0;

     $cgnchargetot = $cgnchargetot + $cgncharge; # *$dv;
     $cgncharge[$curbin] = $cgncharge; # *$dv;
     $cgncharge = 0;
     
     $nonindtotn = $nonindtotn + $nonindn*$dt; # *$dv
     $nonindn[$curbin] = $nonindn*$dt; # *$dv
#      print ("nonindn: $nonindn[$curbin], $nonindn ($curbin)\n");
     $nonindn = 0;

     $nonindtotp = $nonindtotp + $nonindp*$dt; # *$dv
     $nonindp[$curbin] = $nonindp*$dt;  # *$dv
     $nonindp = 0;
     
     $ghsn[$curbin] = $ghsn*$dt; # *$dv
     $ghsn = 0;

     $ghsp[$curbin] = $ghsp*$dt;  # *$dv
     $ghsp = 0;

     $ghin[$curbin] = $ghin*$dt;  # *$dv
     $ghin = 0;

     $ghip[$curbin] = $ghip*$dt; # *$dv
     $ghip = 0;

     $indtotn = $indtotn + $indn*$dv*$dt;
     $indn[$curbin] = $indn*$dv*$dt;
     $indn = 0;

     $indtotp = $indtotp + $indp*$dv*$dt;
     $indp[$curbin] = $indp*$dv*$dt;
     $indp = 0;
     
     $wmax[$curbin] = $wmax;
     
     $cwmass[$curbin] = $cwmass/$intervals;
     $cwmass = 0;

     $rwmass[$curbin] = $rwmass/$intervals;
     $rwmass = 0;

     $icemass[$curbin] = $icemass/$intervals;
     $icemass = 0;

     $snowmass[$curbin] = $snowmass/$intervals;
     $snowmass = 0;

     $grmass[$curbin] = $grmass/$intervals;
     $grmass = 0;

     $hlmass[$curbin] = $hlmass/$intervals;
     $hlmass = 0;

     $wvol5[$curbin] = $wvol5/$intervals;
     $wvol5 = 0;

     $wvol10[$curbin] = $wvol10/$intervals;
     $wvol10 = 0;
     
#     if ( $intervalm == 1 ) {
#     $wmax[$curbin] = $wmax;
#     } else {
#     $wmax[$curbin] = $wmaxtot/$intervals;
#     $wmaxtot = 0;
#     }

     $efield[$curbin] = $efield;

     $efield2[$curbin] = $efield2;
     $efield2 = 0;

     $netchg[$curbin] = $netchg;
     $netchgp[$curbin] = $netchgp;
     $netchgn[$curbin] = $netchgn;
#     print ("$min, $curbin, $times \n");

     $netchgcld[$curbin] = $netchgcld;
     $netchgpcld[$curbin] = $netchgpcld;
     $netchgncld[$curbin] = $netchgncld;
#       if ( $irst2 == 1 ) {  print ("HERE I AM 3: $newbin, $curbin\n"); }
     $curbin = $newbin;
      }
     
#    if ( /^wmax/ ) {
# microphysics_driver: GLOBAL w max/min =    41.61107      -15.55726
    if ( /^ microphysics_driver/ ) {
      chomp();
      @parts = split( / +/ );
      $wmax = $parts[6];
      $wmaxtot = $wmaxtot + $wmax*$dt;
      }
# postlight pos/neg/net charge (C):   1.29770E+03  -1.36479E+03  -6.70915E+01
    if ( /^postlight/ ) {
      chomp();
      @parts = split( / +/ );
      $netchgp = $parts[4];
      $netchgn = $parts[5];
      $netchg = $parts[6];
      }

# storm postlight pos/neg/net charge (C):   1.29770E+03  -1.36479E+03  -6.70915E+01
#    if ( /^storm postlight/ ) {
#pos/neg0:    3.34135E+03  -2.27763E+03   1.06372E+03
    if ( /^pos\/neg0/ ) {
      chomp();
      @parts = split( / +/ );
      $netchgp = $parts[1];
      $netchgn = $parts[2];
      $netchg = $parts[3];
      }
#premic pos/neg/net charge (C):   1.44624008E+02  -1.43335462E+02   1.28854578E+00
    if ( /^premic/ && $com2trmm == 1 ) {
      chomp();
      @parts = split( / +/ );
      $netchgp = $parts[4];
      $netchgn = $parts[5];
      $netchg = $parts[6];
      
#      $ionchg2 = $ionchg2 + $netchg - $netchg1;
#      $ionchgcum2 = $ionchgcum2 + $netchg - $netchg1;
#      $ionerrtot = $ionerrtot + $netchg - $netchg1 - $scrchg0;
      }
    if ( /^IEMAXX, IEMAXY, IEMAXZ,/ ) {
      chomp();
      @parts = split( / +/ );
      $efield = $parts[$#parts];
        if ( $readfirst == 1 ) {
          $readfirst = 0;
          $efield2 = Max($efield2, $efield);
        }
      }

#cwmass1,2,tot =     1.33540E+08    5.23609E+07    1.85901E+08
    if ( /^cwmass1/ ) {
      @parts = split( / +/ );
      $cwmass = $cwmass + $parts[4]*$dt;
    }
    if ( /^rwmass1/ ) {
      @parts = split( / +/ );
      $rwmass = $rwmass + $parts[4]*$dt;
    }
    if ( /^icemass1/ ) {
      @parts = split( / +/ );
      $icemass = $icemass + $parts[4]*$dt;
    }
    if ( /^swmass1/ ) {
      @parts = split( / +/ );
      $snowmass = $snowmass + $parts[4]*$dt;
    }
    if ( /^grmass/ ) {
      @parts = split( / +/ );
      $grmass = $grmass + $parts[4]*$dt;
    }
    if ( /^hlmass/ ) {
      @parts = split( / +/ );
      $hlmass = $hlmass + $parts[4]*$dt;
    }
    if ( /^wvol5/ ) {
      @parts = split( / +/ );
      $wvol5 = $wvol5 + $parts[2]*$dt;
      $wvol10 = $wvol10 + $parts[3]*$dt;
    }
#    $rwmass = 0;
#    $icemass = 0;
#    $snowmass = 0;
#    $grmass = 0;
#    $hlmass = 0;
#    $wvol5 = 0;
#    $wvol10 = 0;

#    if ( /^Charge fallout/ ) {
    if ( /^chgiona1/ ) {
        $readfirst = 1;
      }
    if ( /^ EFIELD2--E-MAX/ ) {
      chomp();
      @parts = split( / +/ );
#      $efield2 = $parts[$#parts];
      }
   if ( $nstep <= $nstop  )  {
    if ( /ctswin,ctswip/ ) {
      chomp();
      @parts = split(/ +|,/);
      $nonindn = $nonindn + $parts[3]*$dxdy;
      $nonindp = $nonindp + $parts[5]*$dxdy;

      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +|,/,$next);
      $nonindn = $nonindn + $parts[3]*$dxdy;
      $ghsn = $ghsn + $parts[3]*$dxdy;
      $nonindp = $nonindp + $parts[5]*$dxdy;
      $ghsp = $ghsp + $parts[5]*$dxdy;

      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +|,/,$next);
      $nonindn = $nonindn + $parts[3]*$dxdy;
      $ghin = $ghin + $parts[3]*$dxdy;
      $nonindp = $nonindp + $parts[5]*$dxdy;
      $ghip = $ghip + $parts[5]*$dxdy;
#      print ("ghin/ghip: $ghin, $ghip\n");
#      print ("nonindn/nonindp: $nonindn, $nonindp\n");

      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +|,/,$next);
      $indn = $indn + $parts[3]*$dxdy;
      $indp = $indp + $parts[5]*$dxdy;
      }
    if ( /^ COM2TRMM: IC flashes,/ ) {
     $com2trmm = 1;
     @parts = split(/ +/ );
     $ic  = $parts[7];
     $cgn = $parts[8];
     $cgp = $parts[9];
     }
    #if ( /^Try number/ ) {
    if ( /^ RETRY/ ) {
     $tries = $tries + 1;
     $bintry = $bintry + 1;
        if ( $type == 1 ) {
           $iccharge = $iccharge - $pos - $neg;
           $ic = $ic - 1;
           $icx = $icx + 1;
           if ( abs($pos) < abs($neg) ) {
             $icdis = $icdis - abs($pos);
            } else {
             $icdis = $icdis - abs($neg); 
            }
         } elsif ( $type == 2 ) {
           $cgpcharge = $cgpcharge - $pos - $neg;
           $cgp = $cgp - 1;
           $cgpx = $cgpx + 1;
         } elsif ( $type == 3 ) {
           $cgncharge = $cgncharge - $pos - $neg;
           $cgn = $cgn - 1;
           $cgnx = $cgnx + 1;
            }
      }
    if ( /^ WARNING: Lightning increased the total energy/ ) {
     $tries = $tries + 1;
     $bintry = $bintry + 1;
    }
    if ( /IC DISCHARGE/ ) {
      $ic = $ic + 1;
      $type = 1;
   #   print OUT ("ic = $ic time = $times \n");
     }
    if ( /DISCHARGE IS POSITIVE/ ) {
      $cgp = $cgp + 1;
   #   print OUT ("cgp = $cgp time = $times \n");
      $type = 2;
     }
    if ( /DISCHARGE IS NEGATIVE/ ) {
      $type = 3;
      $cgn = $cgn + 1;
     }
#    if (/^check totals/ ) {
#      @parts = split(/ +/ );
#      $pos = $parts[2];
#      $neg = $parts[3];
#      if ( ($pos < 0.0001) && (($neg) > -0.0001) ) {
#         $pos = $pos*$dv;
#         $neg = $neg*$dv;
# #       print OUT ( "This is odd: pos = $pos, neg = $neg\n");
#         }
#         if ( $type == 1 ) {
#            $iccharge = $iccharge + $pos + $neg;
#            if ( abs($pos) < abs($neg) ) {
#              $icdis = $icdis + abs($pos);
#             } else {
#              $icdis = $icdis + abs($neg); 
#             }
#          } elsif ( $type == 2 ) {
#            $cgpcharge = $cgpcharge + $pos + $neg;
#   #         print OUT ("cgp : $pos $neg\n");
#          } elsif ( $type == 3 ) {
#            $cgncharge = $cgncharge + $pos + $neg;
#             }
#      }
    if (/START LIGHT/) {
      $pos = 0;
      $neg = 0;
     }
    if (/^ADJ COULOMBS SCNETP/ ) {
     @parts = split(/ +/ );
     $pos = $parts[3];
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +|,/,$next);
     $neg = $parts[3];
        if ( $type == 1 ) {
           $iccharge = $iccharge + $pos + $neg;
           if ( abs($pos) < abs($neg) ) {
             $icdis = $icdis + abs($pos);
            } else {
             $icdis = $icdis + abs($neg); 
            }
         } elsif ( $type == 2 ) {
           $cgpcharge = $cgpcharge + $pos + $neg;
  #         print OUT ("cgp : $pos $neg\n");
         } elsif ( $type == 3 ) {
           $cgncharge = $cgncharge + $pos + $neg;
            }
     }
    if (/^ WARNING: Lightning increased the total energy/) {
        if ( $type == 1 ) {
           $iccharge = $iccharge - $pos - $neg;
           $ic = $ic - 1;
           $icx = $icx + 1;
           if ( abs($pos) < abs($neg) ) {
             $icdis = $icdis - abs($pos);
            } else {
             $icdis = $icdis - abs($neg); 
            }
         } elsif ( $type == 2 ) {
           $cgpcharge = $cgpcharge - $pos - $neg;
           $cgp = $cgp - 1;
           $cgpx = $cgpx + 1;
         } elsif ( $type == 3 ) {
           $cgncharge = $cgncharge - $pos - $neg;
           $cgn = $cgn - 1;
           $cgnx = $cgnx + 1;
            }
     }    
     } # if ( $nstep < $nstop - 1)  
    } # while
   } # foreach

  foo: 
   
    if ( $icx > 0 || $cgpx > 0 || $cgnx > 0 ) {}
    print OUT ( "Removed $icx ICs, $cgpx +CGs, $cgnx -CGs\n" );
    
    $ictot = 0; 
    $cgptot = 0;  
    $cgntot = 0;
    $triestot = 0;
    $icdistot = 0;
    $cgpchargetot = 0;
    $cgnchargetot = 0;
    
#    print ("curbin = $curbin; newbin = $newbin\n");
    
   for ($i = 1; $i < $curbin; $i++ ) {
     $ictot = $ictot + $ic[$i];
     $cgptot = $cgptot + $cgp[$i];
     $cgntot = $cgntot + $cgn[$i];
     $triestot = $triestot + $tries[$i];
     $icdistot = $icdistot + $icdis[$i];
     $cgpchargetot = $cgpchargetot + $cgpcharge[$i];
     $cgnchargetot = $cgnchargetot + $cgncharge[$i];
    }
    
    print OUT ( "Totals: $ictot, $cgptot, $cgntot, $triestot, ",
           "$icdistot ", $icdistot/Max($ictot,1), " $cgpchargetot $cgnchargetot \n" ); 
    print OUT ( "bin ICs CGPs CGNs Tries ICDIS  chg/ic   CGPchg   CGNchg   Netchg ",
                "   NONIn    NONIp    INDn     INDp   wmax  Emax Emax2 NONItot cwmass rainmass icemass snowmass grmass hailmass wvol5 wvol10\n" );
   for ($i = 1; $i < $curbin; $i++ ) {
    if ( $intervalminv == 1 ) {
#    printf OUT ( "%8.3f %3d %3d %3d %3d %8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f %  8.3f %5.2f %e %e\n",
    printf OUT ( "%03d %3d %3d %3d %3d %8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f %  8.3f %5.2f %e %e %e %e %e %e %e %e %e %e %e\n",
      $bintime[$i],  $ic[$i],  $cgp[$i],  $cgn[$i],  
      $tries[$i], $icdis[$i], $icdis[$i]/Max(1,$ic[$i]), $cgpcharge[$i],
      $cgncharge[$i], $netchg[$i],
      $nonindn[$i], $nonindp[$i], $indn[$i], 
      $indp[$i], $wmax[$i], $efield[$i] ,$efield2[$i], $nonindp[$i]-$nonindn[$i], $cwmass[$i], $rwmass[$i], $icemass[$i], 
      $snowmass[$i], $grmass[$i], $hlmass[$i], $wvol5[$i], $wvol10[$i]  );
      } else {

     if ( $bintime[$i] - int($bintime[$i]) == 0  ) {
     printf OUT ( "%8.1f %3d %3d %3d %3d %8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f %  8.3f %5.2f %e %e\n",
      $bintime[$i],  $ic[$i]+$ic[$i-1],  $cgp[$i]+$cgp[$i-1],  $cgn[$i]+ $cgn[$i-1],  
      $tries[$i]+$tries[$i-1], $icdis[$i]+$icdis[$i-1], ($icdis[$i]+$icdis[$i-1])/Max(1,$ic[$i]+$ic[$i-1]), $cgpcharge[$i]+$cgpcharge[$i-1],
      $cgncharge[$i]+$cgncharge[$i-1], $netchg[$i],
      $nonindn[$i]+$nonindn[$i-1], $nonindp[$i]+$nonindp[$i-1], $indn[$i]+$indn[$i-1], 
      $indp[$i]+$indp[$i-1], $wmax[$i], $efield[$i] ,$efield2[$i], $grmass[$i] );
      }
      }

#    print OUT ( $i+$binshift+1,"  $ic[$i]  $cgp[$i]  $cgn[$i] ", 
#      "$tries[$i] $icdis[$i] $iccharge[$i] $cgcharge[$i]\n" );
    }

    sub Max {
        my ($max,$tem);
        $max = shift(@_);
        foreach $tem (@_) {
            $max = $tem if $max < $tem;
        }
        return $max;
    }
