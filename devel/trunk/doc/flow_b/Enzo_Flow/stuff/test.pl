#!/usr/bin/perl -w

#
#  this test function writes a web page, with a whole bunch of fairly
#  neat stuff in it.
#
#

$htmlfile = "test.html";
$sourcefile = "bob.src";
$sourcefile_stripped = "bob";
$c_or_f = "f77";
$date = "13 September 2003";

open(HTML,">$htmlfile") || die "can't open $htmlfile\n";
print HTML "<html>\n<head>\n<title>$sourcefile_stripped ($c_or_f)</title>\n</head>\n<body>\n\n\n";

print HTML "<h1>$sourcefile_stripped ($c_or_f)</h1>\n\n";

print HTML "<h2>Inputs/Outputs</h2>
<p>This routine takes as inputs</p>
<p>This routine modifies </p>

<h2>Function Calls</h2>

";

#
# add function calls here!
#

print HTML "<h2>Header Files</h2>

";

#
# add header files here!
#

print HTML "<h2>Notes</h2>
<p> </p>

<h2>Source Code</h2>

<center><hr size=4 width =\"90\%\"></center>
<pre><!--#include file=\"$sourcefile\"--></pre>
<center><hr size=4 width =\"90\%\"></center>

<center><p><b>Last modified by 
<a href=\"mailto:bwoshea\@cosmos.ucsd.edu\">bwo</a><br>
on $date</b></p></center>";




print HTML "</body></html>\n";
close(HTML);
