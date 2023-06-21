#!/usr/bin/perl -w

#
# this test program goes through $sourcefile and figures
# out what all of the included header files are.
#

$sourcefile = "Grid_FinishFFT.C";

open(SOURCE,"<$sourcefile") || die "can't open $sourcefile\n";

while($linein = <SOURCE>){
  chomp($linein);

  # does the line start with #include?
  if($linein =~ /^\#include/) {
    print "$linein\n";
    # if so, we need to strip off the quotation marks
    # or < > stuff.
    if($linein =~ /\"$/) {
      @bob = split(/\"/,$linein);
      print "\t$bob[1]\n";
    } elsif($linein =~ /\>/){
      @bob = split(/\<|\>/,$linein);
      print "\tdoh @bob[1]\n";
    }
    
  }

}

close(SOURCE);
