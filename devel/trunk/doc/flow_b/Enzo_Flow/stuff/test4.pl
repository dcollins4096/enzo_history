#!/usr/bin/perl -w

#
# load all of the source file names into memory!
#

$fileoffilenames = "sourcefiles.txt";
$bobfile = "bob";

@filenames = qw();
@calledby = qw();
@itcalls = qw();

open(FILENAMES,"<$fileoffilenames") || die "can't open $fileoffilenames";

while($linein = <FILENAMES>){
  chomp($linein);
  @filenames = (@filenames,$linein);
  @calledby = (@calledby,$linein);
  @itcalls = (@itcalls,$linein);
}

close(FILENAMES);

$sizeof = @filenames;
$sizeofcalledby = @calledby;
$sizeofitcalls = @itcalls;

print "this thing has a size of $sizeof\n";
print "$sizeofcalledby $sizeofitcalls\n";

push @{ $calledby[3] }, $bobfile;


$sizeof = @filenames;
$sizeofcalledby = @calledby;
$sizeofitcalls = @itcalls;

print "this thing has a size of $sizeof\n";
print "$sizeofcalledby $sizeofitcalls\n";
print "$calledby[3][0]\n";

$blah = $calledby[3][0];
$blah1 = @{$calledby[3]};

print "$blah $blah1 blah1\n";
