#!perl -T

use Test::More tests => 6;
use constant EPS => 1e-3;
use Array::Compare;

BEGIN {
    use_ok( 'Statistics::Sequences::Runs' ) || print "Bail out!\n";
}

my $runs = Statistics::Sequences::Runs->new();
my $val;

my %ref = (
    swed_1943_1 => {
        observed => 5,
        expected => 9,
        z_value => -2.29,
        p_value => .010973, # based on deviation ratio with ccorr & tails = 1
        p_exact => .0183512,  # SW 1943 p. 71 Table 1
        trials => [5, 20],
        data => [qw/H H H H H H H H H H D H D D D D H H H H H H H H H/],
    },
    swed_1943_2 => {
        observed => 15,
        trials => [20, 20],
        p_exact => .0379982, # SW 1943 p. 82 Table 1
        p_value => .039036, # based on deviation ratio with ccorr & tails = 1
    },
    swed_1943_3 => {
        data => [qw/E O E E O E E E O E E E O E O E/],
        observed => 11,
        expected => 7.88,
        trials => [5, 11],
        p_exact => .0576923, # SW 1943 p. 82 Table 1
    },
);

# test that observed value is less than expectation:
$runs->load(swed => $ref{'swed_1943_1'}->{'data'});
$val = $runs->p_value(exact => -1, tails => 1, ccorr => 1); # this should be below .05, according to Siegal, 1956, Table F;
ok(equal($val, $ref{'swed_1943_1'}->{'p_exact'}) , "exact p_value observed $val != $ref{'swed_1943_1'}->{'p_exact'}");

$val = $runs->p_value(trials => $ref{'swed_1943_2'}->{'trials'}, observed => $ref{'swed_1943_2'}->{'observed'}, exact => -1, tails => 1, ccorr => 1); # this should be below .05, according to Siegal, 1956, Table F;
ok( equal($val, $ref{'swed_1943_2'}->{'p_exact'}) , "exact p_value observed $val != $ref{'swed_1943_2'}->{'p_exact'}");

# test that observed value is greater than expectation:
$val = $runs->p_value(trials => [5, 11], observed => 11, exact => 1, precision_p => 7);
ok( equal($val, $ref{'swed_1943_3'}->{'p_exact'}) , "exact p_value observed $val != $ref{'swed_1943_3'}->{'p_exact'}");

# set hypothesis by the observed deviation:
$runs->load($ref{'swed_1943_1'}->{'data'});
$val = $runs->test(exact => 2); # should be same effect as exact => -1 as rco < rce
ok(equal($val, $ref{'swed_1943_1'}->{'p_exact'}) , "exact p_value observed $val != $ref{'swed_1943_1'}->{'p_exact'}");

$runs->load($ref{'swed_1943_3'}->{'data'});
$val = $runs->test(exact => 2); # should be same effect as exact => 1 as rco > rce
ok(equal($val, $ref{'swed_1943_3'}->{'p_exact'}) , "exact p_value observed $val != $ref{'swed_1943_3'}->{'p_exact'}");

sub equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}
