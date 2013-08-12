#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Statistics::Sequences::Runs', 0.12 ) || print "Bail out!\n";
}

diag( "Testing Statistics::Sequences::Runs $Statistics::Sequences::Runs::VERSION, Perl $], $^X" );
