#!/usr/bin/perl -w

#
#  This test code takes a given function
#  and checks all of the other files in the
#  list of source files to see if it is referred
#  to by it, making sure to exclude files that have
#  been commented out by // (in C++) or 'c' or '!' (in fortran)
#  This function also makes sure to NOT check in the
#  function's own source file!
#  It also makes sure to only count files once!
#

$sourcefile = "star_maker2.src";
$functionname = "star_maker2";
$filenames = "sourcefiles.txt";
$lastfile = NULL;

open(FILES,"<$filenames") || die "can't open $filenames\n";

while($linein = <FILES>) {
  chomp($linein);
  
  open(SOURCEFILE,"$linein") || die "can't open $linein\n";

  while($sourcelinein = <SOURCEFILE>) {
    if($sourcelinein =~ /$functionname/ && $linein ne $sourcefile) {

      if($sourcelinein =~ /$\// && !($sourcelinein =~ /^c/) 
				     && !($sourcelinein =~ /^!/)  ) {
	if( $lastfile ne $linein){
	  print "file $linein contains reference to $functionname!!!\n";
	}
	$lastfile = $linein;
      }
    }
  }

  close(SOURCEFILE);


}

close(FILES);
