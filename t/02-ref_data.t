#!perl -T

use Test::More tests => 28;
use constant EPS => 1e-2;
use Array::Compare;

BEGIN {
    use_ok( 'Statistics::Sequences::Runs' ) || print "Bail out!\n";
}

my $runs = Statistics::Sequences::Runs->new();

my %refdat = (
    sw2 => {
        observed  => 4, p_value => 0.0099, data => [qw//],
    },
    siegal_p142 => {
        observed => 6,
        expected => 12.586206896551724137931034482759,
        variance => 103152/23548, # 4.38
        z_value => -2.9,
        data => [qw/c c c c c c c c c c c c c c e c c c c c c e c e e e e e e/],
        ztest_ok => 1, # Siegal's rule
    },
    swed_1943 => {
        observed => 5,
        expected => 9,
        z_value => -2.29,
        p_value => .010973,
        p_exact => .0183512,
        data => [qw/H H H H H H H H H H D H D D D D H H H H H H H H H/],
    },
    reiter => { # sample from reiter1.com
        data => [qw/m m w w m w m m w w m w m w m w m m w m m w m w m/],
        observed => 19,
        expected => 13.32,
        stdev => 2.41,
        variance => 2.41**2,
        z_value => 2.15,
        p_value => 1 - .984,
    },
);
my @raw_data = ();
my $val;

eval { $runs->load($refdat{'siegal_p142'}->{'data'});};
ok(!$@, do {chomp $@; "Data load failed: $@";});
$val = $runs->observed();
ok(equal($val, $refdat{'siegal_p142'}->{'observed'}), "runcount_observed  $val != $refdat{'siegal_p142'}->{'observed'}");
$val = $runs->expected();
ok(equal($val, $refdat{'siegal_p142'}->{'expected'}), "runcount_expected  $val != $refdat{'siegal_p142'}->{'expected'}");
$val = $runs->variance();
ok(equal($val, $refdat{'siegal_p142'}->{'variance'}), "runcount_variance  $val != $refdat{'siegal_p142'}->{'variance'}");
my $stdev = sqrt($val);
$val = $runs->stdev();
ok(equal($val, $stdev), "runcount stdev observed  $val != $stdev");
my $obsdev = $refdat{'siegal_p142'}->{'observed'} - $refdat{'siegal_p142'}->{'expected'}; 
$val = $runs->obsdev();
ok(equal($val, $obsdev), "runcount obsdev observed  $val != $obsdev");
$val = $runs->z_value(ccorr => 1);
ok(equal($val, $refdat{'siegal_p142'}->{'z_value'}), "runcount_z_value  $val != $refdat{'siegal_p142'}->{'z_value'}");
$val = $runs->ztest_ok();
ok(equal($val, $refdat{'siegal_p142'}->{'ztest_ok'}), "run ztest_ok  $val != $refdat{'siegal_p142'}->{'ztest_ok'}");

eval { $runs->load($refdat{'swed_1943'}->{'data'});};
ok(!$@, do {chomp $@; "Data load failed: $@";});
$val = $runs->observed();
ok(equal($val, $refdat{'swed_1943'}->{'observed'}), "runcount_observed  $val != $refdat{'swed_1943'}->{'observed'}");
$val = $runs->expected();
ok(equal($val, $refdat{'swed_1943'}->{'expected'}), "runcount_expected  $val != $refdat{'swed_1943'}->{'expected'}");
$val = $runs->z_value();
ok(equal($val, $refdat{'swed_1943'}->{'z_value'}), "runcount_z_value  $val != $refdat{'swed_1943'}->{'z_value'}");
$val = $runs->test('stat' => 'runs', tails => 1);
ok(equal($val, $refdat{'swed_1943'}->{'p_value'}), "$val = $refdat{'swed_1943'}->{'p_value'}");
$runs->unload();

# Data from Swed & Eisenhart after pooling:
@raw_data = (1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0);
eval {$runs->load(@raw_data);};
ok(!$@);

$val = $runs->observed();
ok(equal($val, $refdat{'sw2'}->{'observed'}), "runcount_observed  $val != $refdat{'sw2'}->{'observed'}");

$val = $runs->test(tails => 1, ccorr => 1);
ok(equal($val, $refdat{'sw2'}->{'p_value'}), "runcount_p_value  $val != $refdat{'sw2'}->{'p_value'}");

# reiter1.com example:
eval { $runs->load($refdat{'reiter'}->{'data'});};
ok(!$@, do {chomp $@; "Data load failed: $@";});
$val = $runs->observed();
ok(equal($val, $refdat{'reiter'}->{'observed'}), "runcount_observed  $val != $refdat{'reiter'}->{'observed'}");
$val = $runs->expected();
ok(equal($val, $refdat{'reiter'}->{'expected'}), "runcount_expected  $val != $refdat{'reiter'}->{'expected'}");
$val = $runs->z_value(ccorr => 1);
ok(equal($val, $refdat{'reiter'}->{'z_value'}), "runcount_z_value  $val != $refdat{'reiter'}->{'z_value'}");
$val = $runs->test('stat' => 'runs', tails => 1);
ok(equal($val, $refdat{'reiter'}->{'p_value'}), "$val = $refdat{'reiter'}->{'p_value'}");

#$val = $runs->test(exact => 1, tails => 1);
#ok(equal($val, $refdat{'reiter'}->{'p_value'}), "runcount_p_value  $val != $refdat{'reiter'}->{'p_value'}");

# using labelled data:
eval { $runs->load(reiterdata => $refdat{'reiter'}->{'data'});};
ok(!$@, do {chomp $@; "Data load failed: $@";});
$val = $runs->observed(label => 'reiterdata');
ok(equal($val, $refdat{'reiter'}->{'observed'}), "runcount_observed  $val != $refdat{'reiter'}->{'observed'}");

# trying to access by index:
eval { $runs->load(['somedata'], $refdat{'reiter'}->{'data'});};
ok(!$@);
$val = $runs->observed(index => 1);
ok(equal($val, $refdat{'reiter'}->{'observed'}), "runcount_observed  $val != $refdat{'reiter'}->{'observed'}");

# a direct "load"
$val = $runs->observed(data => $refdat{'reiter'}->{'data'});
ok(equal($val, $refdat{'reiter'}->{'observed'}), "runcount_observed  $val != $refdat{'reiter'}->{'observed'}");
# -- even works all the way to deviation p_value?
$val = $runs->p_value(data => $refdat{'reiter'}->{'data'}, exact => 0, ccorr => 1, tails => 1);
ok(equal($val, $refdat{'reiter'}->{'p_value'}), "runcount_p_value  $val != $refdat{'reiter'}->{'p_value'}");

sub equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}
