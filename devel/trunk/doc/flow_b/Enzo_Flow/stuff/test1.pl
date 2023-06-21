#!/usr/bin/perl -w

#
#  this test function parses a text file with a bunch of
#  file names in it, figures out what kind of file it is
#  (ie, C++, fortran source files or header files of some
#  sort) and then, if it's a source file, gets the name
#  of the function.  This means stripping off .C or .src,
#  and, if it's a Grid class object, the Grid_ part.
#

open(SOURCEFILE,"<sourcefiles.txt") || die "can't open sourcefiles.txt";

while($linein = <SOURCEFILE>){
chomp($linein);
#print "$linein";


if($linein =~ /def$/) {
  print "$linein \t\tis a fortran .def file!\n"; 
}

if($linein =~ /h$/) {
  print "$linein \t\tis a C++ header file!\n"; 
}

if($linein =~ /C$/) {
  print "$linein \t\tis a C++ source file!\n"; 
  @chunks = split(/\./,$linein);
  print "\t function name is:  $chunks[0]\n";
  if($linein =~ /^Grid/) {
    @gridchunks = split(/Grid_/,$chunks[0]);
    print "\t\tGrid split off gives us: $gridchunks[1]\n";

  }

}

if($linein =~ /src$/) {
  print "$linein \t\tis a fortran source file!\n"; 
  @chunks = split(/\./,$linein);
  print "\t function name is:  $chunks[0]\n";
}


}



close(SOURCEFILE);

