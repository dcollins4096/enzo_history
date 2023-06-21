#! /usr/bin/perl

$target = shift @ARGV;
$destination = shift @ARGV;
@files = @ARGV;

foreach $file (@files) {
    foreach $thisFile (glob($file)) {
	my $changes = 0;
	open (IN,"<$thisFile");
	my @lines = <IN>;
	close (IN);
	foreach $line (@lines) {
	    if ($line =~ /$target/) {
		$line =~ s/$target/$destination/g;
		$changes++;
	    }
	}
	open (OUT,">$thisFile");
	print OUT @lines;
	close (OUT);
	print "Made $changes replacements in $thisFile.\n" if ($changes);
    }
}
