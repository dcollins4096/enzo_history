#!/usr/bin/perl -w
#
#  enzocodeflow.pl
#
#  This script reads in the user-defined file below, which is
#  supposed to contain the names of all of the Enzo source files
#  and gets us:
#    1)  what type of file it is (C++, f77, or header file)
#    2)  the file name with all extensions and Grid_ stripped off
#        (ie, the function name)
#    3)  the file name with no extensions
#    4)  the names of all of the functions which CALL that function
#    5)  the names of all of the functions which are called by that
#        function.
#    6)  the names of all header files that are included by that function
#
#  Assumptions:
#    I)  The file $fileoffilenames contains a list of files, not
#        necessarily in alphabetical order.  If any of these are NOT
#        source files, Bad Things could happen.
#    II) we ignore all HDF 4 files - hdf 4 is included in enzo
#        only for backwards compatibility and is not supported.
#    III) All source files (or at least a copy of them) are in
#        the directory $htmldir/src .  This is for linking
#        purposes (the #include statements in the web pages).
#        Apache 1.3 seems to not function so well with 
#        a #include "../../foo/blah.C" type instruction.
#        An alternate solution to this is to use a softlink to
#        wherever your source directory is.
#    IV) Assumes that all of the code is enzo code, ie, in the 
#        Enzo code format.
#
#  Notes:
#    The code takes about 10 minutes to run on a 1 ghz pentium III.
#    Most of the slowness is due to disk access times.
#    This isn't perfect - if particular pieces of the code don't 
#    go along with the standard enzo naming convention of the file
#    having the same name as the function (discounting Grid_, HDF5,
#    etc.) then it's not going to be included.  I'll work on that
#    when I get the time.  But hey, if you're advanced enough to
#    notice that, what the heck are you doing reading this stuff?
#
#  Written by bwo (bwoshea@cosmos.ucsd.edu), October 2003
#  This script is made available using the GNP General Public
#  License (http://www.gnu.org/licenses/gpl.html)
#
# TO-DO:
#
# 1. WriteDataHierarchy and EvolveLevel are called recursively - maybe
#    make sure that all of the stuff that's called recursively will be
#    recognized.
#
# 2. The current version doesn't take into account the HDF4/HDF5 compile-
#    time crap.  Also, the 'called from' stuff doesn't ignore comments, so
#    if a routine is mentioned in the comments the browser will include it in
#    the 'calls this function' list.  Perhaps search for function foo as "foo("
#    to eliminate that.
#
# 3. if #2 is done, we'll have to fix the CheckForOverlap files by hand, as it
#    accepts pointers to some of the functions it calls.  Those files are:
#    AddOverlappingParticleMassField, CopyPotentialField, CopyZonesFromGrid, and
#    CopyOverlappingMassField.
#
#
#------------------------- USER DEFINES ----------------------------
$fileoffilenames = "sourcefiles_enzocode.txt";  # file name with all of the
# source file names in it - include the path if necessary
# These files should be all of the .C, .h, .src and .def files!

#
$sourcedir = "src";  # relative path to current dir - no slashes (/)
# if it's the current dir put '.'

$htmldir = "html";  # relative path to current dir - no slashes (/)
# if it's the current dir put '.'

$date = "16 October 2003"; # date for html files

$debug = 0;  # 1 = on, 0 = off.  This controls the verbosity
# of output to the screen.  Even with it off, the script is
# relatively verbose.  With it ON, however, you might want to
# redirect the output into a log file.
#
#-------------------------------------------------------------------

# this forces stuff to be written out immediately, instead of
# buffered.  This is in lieu of writing things to stderr, which
# would make capturing log info a pain.
use IO::Handle;
STDOUT ->autoflush(1);


# declare arrays!
my @filenames; # file names
my @filetype; #file types
my @filename_stripped;  #file name with all Grid_ and exts.
# stripped off of it

my @filename_noextension;  # file name with only extensions stripped
# this is different from @filename_stripped only for C++ Grid class
# functions

my @calledbyfunction;  # 2D array with file names with extensions
# stripped of the functions which CALL this function

my @calledbyfunction_files; # associated file name with no
# extensions

my @functioncalls;  # 2D array with file names with exts. stripped
# of the functions which are called by this function
my @functioncalls_files;  # associated file names (no ext.)
my @headerfile_calls; # 2D arrays with file names - .h and .def kept
my @headerfile_files_stripped; # file names, stripped

# read in list of source files and add to our array
open(FILENAMES,"<$fileoffilenames") || die "can't open $fileoffilenames";

while( $linein = <FILENAMES> ){
  chomp($linein);
  @filenames = (@filenames,$linein);
}

close(FILENAMES);

$numsrcfiles = @filenames;

print "\nThere are $numsrcfiles source files to read through.\n";
print "Now would be a good time to go get some coffee.\n";
if($debug == 1){ print "last line is $filenames[$numsrcfiles-1]\n"; }


# figure out what kind of source file, function name,
# stripped file name, file name w/o extension
print "Getting information from file names...\n";
get_file_info();

# for each file, figure out what functions it calls and store them
print "Getting call tree for each function (this may take a while)...\n";
get_functions_it_calls();

# for each function, figure out what functions call it and store them.
print "Inverting call tree...\n";
get_functions_that_call();


# for each function, print out a web page!

# if it's a Grid:: file, make sure to append
# Grid_ to the beginning of the html file name.
# make sure that if it's not the main function, it at least
# is called by a single file - if it isn't, complain!

print "Creating web pages (this may take a while too)...\n";
make_web_page();

print "Doing some error checking and reporting...\n";
check_calls_and_called();

print "\n\nDone!  Exiting...\n\n";


#
# -------------------- SUBROUTINES -----------------------
#

#--------------------------------------------------------------
# loops over all source files and figures out what kind of
# file it is - C++, f77, C++ header or f77 .def file.  Then
# get the type, the stripped filename, and the filename with
# no extension (the latter two are only different in the case of
# Grid and ExteriorBoundary class files, and also those that
# have HDF5 in the name) and put them into arrays.
#
sub get_file_info {

  # loop over all source files 
  for ($i = 0; $i < $numsrcfiles; $i++){
    # fortran .def file?
    if($filenames[$i] =~ /def$/) {
      @filetype = (@filetype,"f77 .def file");
      @chunks = split(/\./,$filenames[$i]);
      @filename_stripped = (@filename_stripped,$chunks[0]);
      @filename_noextension = (@filename_noextension,$chunks[0]);
    }

    # C++ header file?
    if($filenames[$i] =~ /h$/) {
      @filetype = (@filetype,"C++ .h file");
      @chunks = split(/\./,$filenames[$i]);
      @filename_stripped = (@filename_stripped,$chunks[0]);
      @filename_noextension = (@filename_noextension,$chunks[0]);
    }

    # C++ source (non-header) file?
    # if so, figure out what kind of file it is
    if($filenames[$i] =~ /C$/) {
      @filetype = (@filetype,"C++");
      @chunks = split(/\./,$filenames[$i]);

      # is this a member of the Grid class?
      if($filenames[$i] =~ /^Grid/) {
	@gridchunks = split(/Grid_/,$chunks[0]);
	# is this file an HDF5-based file?
	# if so, strip off the HDF5.
	if($gridchunks[1] =~ /HDF5$/){
	  @hdfchunks = split(/HDF5/,$gridchunks[1]);

	  @filename_stripped = (@filename_stripped,$hdfchunks[0]);
	  @filename_noextension = (@filename_noextension,$chunks[0]);
	} else {
	  @filename_stripped = (@filename_stripped,$gridchunks[1]);
	  @filename_noextension = (@filename_noextension,$chunks[0]);
	}

	# is this a member of the ExteriorBoundary class?
      } elsif($filenames[$i] =~ /^ExternalBoundary_/){
	@ebchunks = split(/ExternalBoundary_/,$chunks[0]);

	# hdf 5-based file?  If so, strip HDF5 off the file name too
	if($ebchunks[1] =~ /HDF5$/){
	  @hdfchunks = split(/HDF5/,$ebchunks[1]);

	  @filename_stripped = (@filename_stripped,$hdfchunks[0]);
	  @filename_noextension = (@filename_noextension,$chunks[0]);

	} else {
	  @filename_stripped = (@filename_stripped,$ebchunks[1]);
	  @filename_noextension = (@filename_noextension,$chunks[0]);
	}

	# if it's not a Grid:: or ExteriorBoundary:: file, is
	# it an HDF 5 file?
      } elsif($chunks[0] =~ /HDF5$/){

	@hdfchunks = split(/HDF5/,$gridchunks[1]);	
	@filename_stripped = (@filename_stripped,$hdfchunks[0]);
	@filename_noextension = (@filename_noextension,$chunks[0]);

	# if it isn't any of those, do nothing exciting.
      } else {
	@filename_stripped = (@filename_stripped,$chunks[0]);
	@filename_noextension = (@filename_noextension,$chunks[0]);
      }

    }

    # fortran source (non-header file)?
    if($filenames[$i] =~ /src$/) {
      @filetype = (@filetype,"f77");
      @chunks = split(/\./,$filenames[$i]);
      @filename_stripped = (@filename_stripped,$chunks[0]);
      @filename_noextension = (@filename_noextension,$chunks[0]);
    }

    if($debug==1) {
      print "$filenames[$i] $filetype[$i] $filename_stripped[$i] $filename_noextension[$i]\n";
    }

  }
}


#--------------------------------------------------------------------
#  Go through the list of files.  Parse each line in each C++ and f77
#  source file (skipping header files) and see if that line has a call
#  to any of the functions or header files in the list.
sub get_functions_it_calls{

  # loop over all source files 
  for ($i = 0; $i < $numsrcfiles; $i++){

    if($debug == 1){ 
      print "\t Reading file $filenames[$i]\n"; 
    } else {
      print ".";
    }

    # skip over .h and .def files since they don't call anything (or
    # at least not anything we really care about)
    if($filetype[$i] eq "f77 .def file" || $filetype[$i] eq "C++ .h file"){
      if($debug ==1){ print "$filenames[$i] is a .h or .def file, skipping.\n"; }
      next;
    }

    # open .src and .C files.
    open(SOURCEFILE,"<$sourcedir/$filenames[$i]") || die "can't open $filenames[$i]\n";

    # go through entire file, line-by-line
    while($sourcelinein = <SOURCEFILE>){
      chomp($sourcelinein);

      # Is this line a comment?  If so, just skip it.
      if($sourcelinein =~ /^\// || $sourcelinein =~ /^c/ ||
	 $sourcelinein =~ /^!/ ) {
	next;
      }

      # actually go through all of the stripped file names
      # (ie, function names) and figure out what we've got.
      for($j = 0; $j < $numsrcfiles; $j++){

	# we don't want to count ourselves!
	if($i == $j || $filenames[$j] eq "main.C"){
	  next;
	}


	# does the stripped file name match this line?  If so, 
	# we need to figure out what exactly it IS.  ie, is it
	# header or a function?
	if($sourcelinein =~ /$filename_stripped[$j]/  ) {

	  # is it a header file?
	  if($sourcelinein =~ /^\#include/){
	    if($debug ==1){ print "$filenames[$i] includes the header file $filenames[$j]\n"; }

	    # write to list of header files included,also file name w/no extensions

	    push @{ $headerfile_calls[$i] }, $filenames[$j];
	    push @{ $headerfile_files_stripped[$i] }, $filenames[$j];

	    if($debug ==1){ 
	      print "$filenames[$i]:  added  $filenames[$j] and $filenames[$j] to arrays.\n"; 
	    }

	  } else { # if not, it's a function of some sort! (subject to caveat...)

	    # takes care of all references that are actually just part of another
	    # file name
	    if($filetype[$j] eq "f77 .def file" || $filetype[$j] eq "C++ .h file"){
	      if($debug == 1){
		print "$filenames[$i]: we don't want to include $filename_stripped[$j]";
		print "in $filenames[$j] - header or slip-through!\n";
	      }
	      next;
	    }

	    if($debug == 1){
	      print "$filenames[$i] includes $filename_stripped[$j] in $filenames[$j]\n";
	    }

	    # write to list of functions called
	    push @{ $functioncalls[$i] }, $filename_stripped[$j];
	    push @{ $functioncalls_files[$i] }, $filenames[$j];

	    if($debug == 1){
	      print "$filenames[$i]:  added  $filename_stripped[$j] and $filenames[$j] to arrays.\n";
	    }


	  } #  } else { # if not, it's a function of som

	} # if($sourcelinein =~ /$filename_stripped[$j]/  ) {
	
      } # for($j ...
    }  # while($sourcelinein..


    close(SOURCEFILE);

    # print out all of the function and header files which are called by this function,
    # making sure to only print out things ONCE!
    if($debug == 1){

      print "Summary:\n";

      if($headerfile_calls[$i]){
	$loopsize = @{ $headerfile_calls[$i] };
	if($loopsize > 0){
	  print "header files included by $filenames[$i]: \n";
	  for($k = 0; $k < $loopsize; $k++){
	    if($k > 0) {
	      $printthis = 1;
	
	      for($l=0; $l < $k; $l++){
		if($headerfile_calls[$i][$k] eq $headerfile_calls[$i][$l]){
		  $printthis = 0;
		}
	      }
	      if($printthis == 1){
		print "\t$headerfile_calls[$i][$k] \n";
	      }
	
	    } else {
	      print "\t$headerfile_calls[$i][$k] \n";
	    }
	  }
       
	  print "\n";
	}
      }
      
      
      if($functioncalls[$i]){
	$loopsize = @{ $functioncalls[$i] };
	if($loopsize > 0){
	  print "functions called  by $filenames[$i]: \n";
	  for($k = 0; $k < $loopsize; $k++){
	    if($k > 0) {
	      $printthis = 1;
	      
	      for($l=0; $l < $k; $l++){
		if($functioncalls[$i][$k] eq $functioncalls[$i][$l]){
		  $printthis = 0;
		}
	      }
	      if($printthis == 1){
		print "\t$functioncalls[$i][$k] \n";
	      }
	      
	    } else {
	      print "\t$functioncalls[$i][$k] \n";
	    }
	  }
	
	  print "\n";
	}
      }
    }

  } # for($i = 0....

}


#--------------------------------------------------------------
#  This routine basically inverts get_functions_it_calls.  For
#  each function we're interested in, go through the list of
#  @functioncalls and figure out which routines call that
#  particular function or header file.
#
sub get_functions_that_call{

  # loop over all source files.
  for($i = 0; $i <  $numsrcfiles; $i++){
    # is this a header file?  If so, do one thing...  If not, it's a 
    # source file so we treat it differently.
    if($filetype[$i] eq "f77 .def file" || $filetype[$i] eq "C++ .h file"){

      # for file $i, go through all other arrays and figure out if it's
      # called from those files.
      for($j = 0; $j < $numsrcfiles; $j++){
	if($i == $j) { next; }

	# are there any headerfile calls?
	if($headerfile_calls[$j]){
	  $loopsize = @{ $headerfile_calls[$j] };
	  for($k = 0; $k < $loopsize; $k++){
	    if($headerfile_calls[$j][$k] eq $filenames[$i]){
	      if($debug == 1){ print "$filenames[$i] called by $filename_stripped[$j]\n"; }

	      push @{ $calledbyfunction[$i] }, $filename_stripped[$j];
	      push @{ $calledbyfunction_files[$i] }, $filenames[$j];

	    }
	  }
	}
      }

    } else {  # must be a source file

      # for file $i, go through all other arrays and figure out if it's
      # called from those files.
      for($j = 0; $j < $numsrcfiles; $j++){
	if($i == $j) { next; }

	if($functioncalls[$j]){
	  $loopsize = @{ $functioncalls[$j] };
	  for($k = 0; $k < $loopsize; $k++){
	    if($functioncalls[$j][$k] eq $filename_stripped[$i]){
	      if($debug == 1){ print "$filenames[$i] called by $filename_stripped[$j]\n"; }

	      push @{ $calledbyfunction[$i] }, $filename_stripped[$j];
	      push @{ $calledbyfunction_files[$i] }, $filenames[$j];

	    }
	  }
	}
      }
    } # else { # must be a source file

  }  # for($i = 0; $i <  $numsrcfiles; $i++){

}


#--------------------------------------------------------------
#  Put together the set of web pages, in whatever directory
#  $htmldir is in.  This makes some attempt to make
#  maintainable and extensible pages, which is really fairly
#  simple to do since we're not asking for a lot out of our
#  web pages here.
sub make_web_page{

  # loop over source code files
  for($i = 0; $i < $numsrcfiles; $i++){

    # get html file name
    $htmlfile = $htmldir . "/" .$filenames[$i] . ".html";
    if($debug == 1){ 
      print "\t writing file $htmlfile\n"; 
    } else {
      print ".";
    }

    open(HTML,">$htmlfile") || die "can't open $htmlfile\n";

    print HTML "<html>\n<head>\n<title>$filename_stripped[$i] ($filetype[$i])</title>\n</head>\n<body>\n\n\n";
    print HTML "<center><h1>$filename_stripped[$i] ($filetype[$i])</h1></center>\n\n";
    #    print HTML "<h2>Inputs/Outputs</h2>
    #<p>This routine takes as inputs</p>
    #<p>This routine modifies </p>\n\n";
    print HTML "<h2>Included files</h2>\n
<p>This function includes the following files:<p>
\n";


    # add header file includes here, making sure to only print each one once
    # loop over all headerfile calls, checking to make sure that each one is
    # only printed out once.  No attempt is made to alphabetize.
    if($headerfile_calls[$i]){
      print HTML "<p>\n";
      $loopsize = @{ $headerfile_calls[$i] };
      if($loopsize > 0){
	#print "header files included by $filenames[$i]: \n";
	for($k = 0; $k < $loopsize; $k++){
	  if($k > 0) {
	    $printthis = 1;

	    for($l=0; $l < $k; $l++){
	      if($headerfile_calls[$i][$k] eq $headerfile_calls[$i][$l]){
		$printthis = 0;
	      }
	    }
	    if($printthis == 1){
	      $headerhtml = $headerfile_files_stripped[$i][$k] . ".html";
	      print HTML "<a href=\"$headerhtml\">$headerfile_calls[$i][$k]</a><br>\n";
	    }

	  } else {
	    $headerhtml = $headerfile_files_stripped[$i][$k] . ".html";
	    print HTML "<a href=\"$headerhtml\">$headerfile_calls[$i][$k]</a><br>\n";
	    #print "$headerfile_calls[$i][$k] \n";
	  }
	}

	print HTML "</p>\n";
      }
    } else {
      print HTML "<p>No header files were included in this file.</p>\n";
    }

    #  Add function calls here.  Loop over list of functions called
    #  from this sub function, making sure to only print each one
    #  once (since a lot of them call things multiple times).  Also,
    #  make sure to print out which class the function belongs to,
    #  if it does belong to a class.
    print HTML "<h2>Function Calls</h2>\n";
    if($functioncalls[$i]){
      print HTML "<p>This function calls the following functions:</p>\n\n";
      $loopsize = @{ $functioncalls[$i] };
      if($loopsize > 0){
	#print "functions called  by $filenames[$i]: \n";
	for($k = 0; $k < $loopsize; $k++){
	  if($k > 0) {
	    $printthis = 1;

	    for($l=0; $l < $k; $l++){
	      if($functioncalls[$i][$k] eq $functioncalls[$i][$l]){
		$printthis = 0;
	      }
	    }
	    if($printthis == 1){
	      $funcname_html = $functioncalls_files[$i][$k] . ".html";
	      #print HTML "<a href=\"$funcname_html\">$functioncalls[$i][$k]</a><br>\n";
	      if($functioncalls_files[$i][$k] =~ /^Grid_/) {
		print HTML "<a href=\"$funcname_html\">Grid::$functioncalls[$i][$k]</a><br>\n";
	      } elsif ( $functioncalls_files[$i][$k] =~ /^ExternalBoundary_/){
	        print HTML "<a href=\"$funcname_html\">ExternalBoundary::$functioncalls[$i][$k]</a><br>\n";
	      } else {
		print HTML "<a href=\"$funcname_html\">$functioncalls[$i][$k]</a><br>\n";
	      }
	    }

	  } else {
	      $funcname_html = $functioncalls_files[$i][$k] . ".html";
	      #print HTML "<a href=\"$funcname_html\">$functioncalls[$i][$k]</a><br>\n";
	      if($functioncalls_files[$i][$k] =~ /^Grid_/) {
		print HTML "<a href=\"$funcname_html\">Grid::$functioncalls[$i][$k]</a><br>\n";
	      } elsif ( $functioncalls_files[$i][$k] =~ /^ExternalBoundary_/){
	        print HTML "<a href=\"$funcname_html\">ExternalBoundary::$functioncalls[$i][$k]</a><br>\n";
	      } else {
		print HTML "<a href=\"$funcname_html\">$functioncalls[$i][$k]</a><br>\n";
	      }
	  }
	}

	print HTML "\n";
      }
    } else {
      print HTML "<p>This function called no subfunctions.</p>\n";
    }


    #  Put list of functions that this function is called by,
    #  only printing out each one once and making sure to get
    #  class calls right (etc etc, see above)
    print HTML "<h2>Called by</h2>\n";
    if($calledbyfunction[$i]){
      print HTML "<p>This function is called by the following functions:</p>\n\n";
      $loopsize = @{ $calledbyfunction[$i] };
      if($loopsize > 0){
	#print "functions called  by $filenames[$i]: \n";
	for($k = 0; $k < $loopsize; $k++){
	  if($k > 0) {
	    $printthis = 1;

	    for($l=0; $l < $k; $l++){
	      if($calledbyfunction[$i][$k] eq $calledbyfunction[$i][$l]){
		$printthis = 0;
	      }
	    }
	    if($printthis == 1){
	      $funcname_html = $calledbyfunction_files[$i][$k] . ".html";
	      if($calledbyfunction_files[$i][$k] =~ /^Grid_/) {
		print HTML "<a href=\"$funcname_html\">Grid::$calledbyfunction[$i][$k]</a><br>\n";
	      } elsif ( $calledbyfunction_files[$i][$k] =~ /^ExternalBoundary_/){
	        print HTML "<a href=\"$funcname_html\">ExternalBoundary::$calledbyfunction[$i][$k]</a><br>\n";
	      } else {
		print HTML "<a href=\"$funcname_html\">$calledbyfunction[$i][$k]</a><br>\n";
	      }

	    }

	  } else {
	      $funcname_html = $calledbyfunction_files[$i][$k] . ".html";
	      
	      if($calledbyfunction_files[$i][$k] =~ /^Grid_/) {
		print HTML "<a href=\"$funcname_html\">Grid::$calledbyfunction[$i][$k]</a><br>\n";
	      } elsif ( $calledbyfunction_files[$i][$k] =~ /^ExternalBoundary_/){
	        print HTML "<a href=\"$funcname_html\">ExternalBoundary::$calledbyfunction[$i][$k]</a><br>\n";
	      } else {
		print HTML "<a href=\"$funcname_html\">$calledbyfunction[$i][$k]</a><br>\n";
	      }
	    }
	  
	  
	  print HTML "\n";
	}
      } else {
	print HTML "<p>This function is not called by any other function.</p>\n";
      }
    }

#  print the #include for the source file, and also some misc stuff.
#  make sure to set the date and who wrote it so they can hunt me
#  down and complain when something's not 100% accurate.
$sourcefilename = "src/" . $filenames[$i];

    print HTML "
<h2>Notes</h2>

<p>&nbsp;</p>
<p>&nbsp;</p>

<center><h2>Source Code</h2></center>

<center><hr size=4 width =\"90\%\"></center>
<pre><!--#include file=\"$sourcefilename\"--></pre>
<center><hr size=4 width =\"90\%\"></center>

<center><p><b>This web page was automagically created by <br>
<i>enzocodeflow</i> on $date.<br>
Email <a href=\"mailto:bwoshea\@cosmos.ucsd.edu\">bwo</a>
with questions or comments.</b></p></center>
</body></html>\n";

    close(HTML);

  } # for($i = 0; $i < $numsrcfiles; $i++){
}


#-----------------------------------------------------------
# This function basically just prints out which files aren't
# called by anything (except for main.C, which shouldn't ever
# be called).  It's for my own error-checking purposes, really.
# Put it all into a html file called "not_called.html" in the
# html directory.
sub check_calls_and_called {

  $htmlfile = $htmldir . "/" . "not_called.html";
  open(NOTCALLED,">$htmlfile") || die "can't open $htmlfile";
  print NOTCALLED "
<html>
<head>
<title>uncalled routines/files</title>
</head>
<body>
<p><center><h1>uncalled routines/files</h1></center></p>
<p>&nbsp;</p>
<p>";

  for($i = 0; $i < $numsrcfiles; $i++){

    if( @{ $calledbyfunction[$i] } > 0){
      $calledby = @{ $calledbyfunction[$i] };
    } else {
      $calledby = 0;
      print "******* $filenames[$i] wasn't called by anybody!  $calledby *********\n";
      $uncalledhtml = $filenames[$i] . ".html";
      print NOTCALLED "<a href=\"$uncalledhtml\">$filenames[$i]</a><br>\n";
    }

  }

print NOTCALLED "
</p>
</body>
</html>";

close(NOTCALLED);

}
