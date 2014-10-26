use Test::More tests => 3;

BEGIN {
    use_ok( 'Statistics::Sequences::Runs' ) || print "Bail out!\n";
}

# usually need to know the frequency of each element from loaded/given data:

my @ari = (1, 1, 0, 1);
my @frq = Statistics::Sequences::Runs::_calc_run_Ns(\@ari);
@frq = sort {$b <=> $a } @frq; # ensure have the max and then min frequencies
ok(equal($frq[0], 3), "element 1 frequency: observed = $frq[0]\texpected = 3");
ok(equal($frq[1], 1), "element 2 frequency: observed = $frq[1]\texpected = 1");

sub equal {
    return 1 if $_[0] == $_[1];
    return 0;
}