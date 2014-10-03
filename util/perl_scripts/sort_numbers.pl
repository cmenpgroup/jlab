#!/usr/bin/perl -w

$file = $ARGV[0];

my @numbers;

$i=0;


open(LIST, "$file") || die "\n$0: Could not open $file.\n";
while (defined($line=<LIST>)) {
    chomp($line);

    $numbers[$i]=$line;
    $i++;
}
close(LIST);

my @sorted_numbers = sort { $a <=> $b } @numbers;

open(OUTFILE,"> $file-sorted");

for (my $j=0; $j<$i; $j++) {
    print OUTFILE "$sorted_numbers[$j]\n";
}
close(OUTFILE);
