package Statistics::Sequences::Runs;

use 5.008008;
use strict;
use warnings;
use Carp qw(carp croak);
use vars qw($VERSION @ISA);
use Statistics::Sequences 0.10;
#@ISA = qw(Statistics::Sequences);
use Moose;
extends 'Statistics::Sequences';
use List::AllUtils qw(sum uniq);
use Number::Misc 'is_even';
use Statistics::Zed 0.072;
our $zed = Statistics::Zed->new();

$VERSION = '0.10';

=pod

=head1 NAME

Statistics::Sequences::Runs - observed, expected and variance counts, deviation and combinatorial tests, of Wald-type runs

=head1 SYNOPSIS

 use strict;
 use Statistics::Sequences::Runs 0.10; # methods/args here are not compatible with earlier versions
 my $runs = Statistics::Sequences::Runs->new();
 $runs->load(qw/1 0 0 0 1 1 0 1 1 0 0 1 0 0 1 1 1 1 0 1/); # dichotomous sequence (any values); or send as "data => $aref" with each stat call
 my $val = $runs->observed(); # other methods include: expected(), variance(), obsdev() and stdev()
 $val = $runs->expected(trials => [11, 9]); # by-passing need for data; also works with other methods except observed()
 $val = $runs->z_value(tails => 1, ccorr => 1); # or want an array & get back both z- and p-value
 $val = $runs->z_value(trials => [11, 9], observed => 11, tails => 1, ccorr => 1); # by-pass need for data; also works with p_value()
 $val = $runs->p_value(exact => 1, tails => 1); # assuming data are loaded; alias: test()
 my $href = $runs->stats_hash(values => {observed => 1, p_value => 1}, exact => 1); # include any other stat-method as needed
 $runs->dump(values => {observed => 1, expected => 1, p_value => 1}, exact => 1, flag => 1, precision_s => 3, precision_p => 7);
 # prints: observed = 11.000, expected = 10.900, p_value = 0.5700167

=head1 DESCRIPTION

A Wald-type run is a sequence of identical events on 1 or more consecutive trials. For example, in a signal-detection test with sequence over time of match (H) and miss (M) events like H-H-M-H-M-M-M-M-H, there are 5 runs: 3 Hs, 2 Ms. This number of runs over time between two events can be compared with the number expected to occur by chance over the number of trials, and variance has also been worked out (see REFERENCES).

More runs than expected ("negative serial dependence") can denote irregularity, instability, mixing up of alternatives - lack of prejudice, even. Fewer runs than expected ("positive serial dependence") can denote cohesion, insulation, isolation of alternatives - loneliness, even. Both can indicate sequential dependency: either negative (a bias to produce too many alternations), or positive (a bias to produce too many repetitions).

The distribution of runs is asymptotically normal, and a deviation-based test of extra-chance occurrence when at least one alternative has more than 20 occurrences (Siegal rule), or both event occurrences exceed 10 (Kelly, 1982), is conventionally considered reliable; otherwise, consider exact test option.

Have non-dichotomous, continuous or multinomial data? See L<Statistics::Data::Dichotomize> for how to prepare non-dichotomous data, whether numerical or made up of categorical events, for test of runs.

=head1 METHODS

=head2 new

 $runs = Statistics::Sequences::Runs->new();

Returns a new Runs object. Expects/accepts no arguments but the classname.

=head2 load

 $runs->load(@data); # anonymously
 $runs->load(\@data);
 $runs->load('sample1' => \@data); # labelled whatever

Loads data anonymously or by name - see L<load|Statistics::Data/load, load_data> in the Statistics::Data manpage for details on the various ways data can be loaded and then retrieved (more than shown here).

After the load, the data are L<read|Statistics::Data/read, read_data, get_data> to ensure that they contain only two unique elements - if not, carp occurs and 0 rather than 1 is returned. 

Alternatively, skip this action; data don't always have to be loaded to use the stats methods here. To get the observed number of runs, data of course have to be loaded, but other stats can be got if given the observed count - otherwise, they too depend on data having been loaded.

Every load unloads all previous loads and any additions to them.

sub load {
    my $self = shift;
    $self->SUPER::load(@_);
    my $data = $self->read(@_);
    my $nuniq = scalar(uniq(@{$data}));
    if ($nuniq > 2) {
        carp __PACKAGE__, ' More than two elements were found in the data: ' . join(' ', uniq(@$data));
        return 0;
    }
    else {
        return 1;
    }
}

=head2 add, read, unload

See L<Statistics::Data> for these additional operations on data that have been loaded.

=head2 observed, runcount_observed, rco

 $v = $runs->observed(); # use the first data loaded anonymously
 $v = $runs->observed(index => 1); # ... or give the required "index" for the loaded data
 $v = $runs->observed(label => 'mysequence'); # ... or its "label" value
 $v = $runs->observed(data => [1, 0, 1, 1]); # ... or just give the data now

Returns the observed number of runs in the loaded or given data.

=cut

sub observed {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my $data = ref $args->{'data'} ? $args->{'data'} : $self->read($args);
    my $num = scalar(@$data);
    my ($rco, $i) = (0, 0);
    foreach (; $i < $num; $i++) {
        $rco++ if !$i || ( $data->[$i] ne $data->[$i - 1] );
    }
    return $rco;
}
*runcount_observed = \&observed;
*rco = \&observed;

=head2 expected, runcount_expected, rce

 $v = $runs->expected(); # or specify loaded data by "index" or "label", or give it as "data" - see observed()
 $v = $runs->expected(data => [qw/blah bing blah blah blah/]); # use these data
 $v = $runs->expected(trials => [12, 7]); # don't use actual data; calculate from these two Ns

Returns the expected number of runs across the loaded data. Expectation is given as follows: 

=for html <p>&nbsp;&nbsp;<i>E[R]</i> = ( (2<i>n</i><sub>1</sub><i>n</i><sub>2</sub>) / (<i>n</i><sub>1</sub> + <i>n</i><sub>2</sub>) ) + 1</p>

where I<n>(I<i)> is the number of observations of each element in the data.

=cut 

sub expected {
   my $self = shift;
   my $args = ref $_[0] ? $_[0] : {@_};
   my ($n1, $n2) = _get_run_Ns($self, $args);
   my $sum = $n1 + $n2;
   return $sum ? ( ( 2 * $n1 * $n2 ) /  $sum ) + 1 : undef;
}
*rce = \&expected;
*runcount_expected = \&expected;

=head2 variance, runcount_variance, rcv

 $v = $runs->variance(); # use data already loaded - anonymously; or specify its "label" or "index" - see observed()
 $v = $runs->variance(data => [qw/blah bing blah blah blah/]); # use these data
 $v = $runs->variance(trials => [5, 12]); # use these trial numbers - not any particular sequence of data

Returns the variance in the number of runs for the given data.

=for html <p>&nbsp;&nbsp;<i>V[R]</i> = ( (2<i>n</i><sub>1</sub><i>n</i><sub>2</sub>)([2<i>n</i><sub>1</sub><i>n</i><sub>2</sub>] &ndash; [<i>n</i><sub>1</sub> + <i>n</i><sub>2</sub>]) ) / ( ((<i>n</i><sub>1</sub> + <i>n</i><sub>2</sub>)<sup>2</sup>)((<i>n</i><sub>1</sub> + <i>n</i><sub>2</sub>) &ndash; 1) ) </p>

defined as above for L<runcount_expected|expected, runcount_expected, rce>.

The data to test can already have been L<load|load>ed, or you send it directly as a flat referenced array keyed as B<data>.

=cut

sub variance {
   my $self = shift;
   my $args = ref $_[0] ? $_[0] : {@_};
   my ($n1, $n2) = _get_run_Ns($self, $args);
   my $sum = $n1 + $n2;
   return $sum < 2 ? 1 : (( 2 * $n1 * $n2 * ( ( 2 * $n1 * $n2 ) - $sum) ) 
            /
           ( ( $sum**2 ) * ( $sum - 1 ) ));
}
*rcv = \&variance;
*runcount_variance = \&variance;

=head2 obsdev, observed_deviation

 $v = $runs->obsdev(); # use data already loaded - anonymously; or specify its "label" or "index" - see observed()
 $v = $runs->obsdev(data => [qw/blah bing blah blah blah/]); # use these data

Returns the deviation of (difference between) observed and expected runs for the loaded/given sequence (I<O> - I<E>). 

=cut

sub obsdev {
    return observed(@_) - expected(@_);
}
*observed_deviation = \&obsdev;

=head2 stdev, standard_deviation

 $v = $runs->stdev(); # use data already loaded - anonymously; or specify its "label" or "index" - see observed()
 $v = $runs->stdev(data => [qw/blah bing blah blah blah/]);

Returns square-root of the variance.

=cut

sub stdev {
    return sqrt(variance(@_));
}
*standard_deviation = \&stdev;

=head2 z_value, runcount_zscore, rzs, zscore

 $v = $runs->z_value(ccorr => 1); # use data already loaded - anonymously; or specify its "label" or "index" - see observed()
 $v = $runs->z_value(data => $aref, ccorr => 1);
 ($zvalue, $pvalue) = $runs->z_value(data => $aref, ccorr => 1, tails => 2); # same but wanting an array, get the p-value too

Returns the zscore from a test of runcount deviation, taking the runcount expected away from that observed and dividing by the root expected runcount variance, by default with a continuity correction to expectation. Called wanting an array, returns the z-value with its I<p>-value for the tails (1 or 2) given.

The data to test can already have been L<load|load>ed, or sent directly as an aref keyed as B<data>.

Other options are B<precision_s> (for the z_value) and B<precision_p> (for the p_value).

=cut

sub z_value {
   my $self = shift;
   my $args = ref $_[0] ? shift : {@_};
   my $ccorr = defined $args->{'ccorr'} ? $args->{'ccorr'} : 1;
   my $tails = $args->{'tails'} || 2;
   my $precision_s = $args->{'precision_s'};
   my $precision_p = $args->{'precision_p'};
   my %stat_args = ();
   if (ref $args->{'trials'}) { # allow stats to be calculated from Ns
       $stat_args{'trials'} = $args->{'trials'};
   }
   else {
       $stat_args{'data'} = ref $args->{'data'} ? $args->{'data'} : $self->read($args);
   }
   my $rco = defined $args->{'observed'} ? $args->{'observed'} : $self->rco(\%stat_args); # giving rco the data
   my ($zval, $pval) = $zed->zscore(
        observed => $rco,
        expected => $self->rce(\%stat_args),
        variance => $self->rcv(\%stat_args),
        ccorr => $ccorr,
        tails => $tails,
        precision_s => $precision_s, 
        precision_p => $precision_p,
     );
    return wantarray ? ($zval, $pval) : $zval;
}
*rzs = \&z_value;
*runcount_zscore = \&z_value;
*zscore = \&z_value;

=head2 p_value, test, runs_test, rct

 $p = $runs->p_value(); # using loaded data and default args
 $p = $runs->p_value(ccorr => 0|1, tails => 1|2); # normal-approximation based on loaded data
 $p = $runs->p_value(exact => 0|1|-1); # if not zero: calc combinatorially for either upper- or lower-tail test
 $p = $runs->p_value(data => [1, 0, 1, 1, 0], exact => 1); #  using given data (by-passing load and read)
 $p = $runs->p_value(trials => [12, 12], observed => 8); # without using data, specifying N-per-event and run-count

Returns the probability of getting the observed number of runs or a smaller number given the number of each of the two events. By default, a large sample is assumed, and the probability is obtained from the normalized deviation, as given by the L<zscore|Statistics::Sequences:Runs/zscore> method.

If the option B<exact> is defined and does not equal zero, then the probability is worked out combinatorially, as per Swed & Eisenhart (1943), Eq. 1, p. 66 (see also Siegal, 1956, Eqs. 6.12a and 6.12b, p. 138). By default, this is a one-tailed test, testing the hypotheses that there are either too many or too few runs relative to chance expectation. Which hypothesis can be specified by giving B<exact> => 1 (upper-tail test that observed runs >= expected runs), or to -1 (lower-tail test that observed runs < expected runs). Alternatively, if B<exact> is set to any other value (e.g., 'auto'), the "correct" hypothesis is tested based on the expected value returned by the L<expected|Statistics::Sequences::Runs/expected, runcount_expected, rce> method. Setting B<tails> => 2 simply doubles the one-tailed I<p>-value from any of these tests. Output from these tests has been checked against the tables and examples in Swed & Eisenhart (given to 7 decimal places), and found to agree.

The option B<precision_p> gives the returned I<p>-value to so many decimal places.

=cut

sub p_value {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my $pval;
    if (defined $args->{'exact'} and $args->{'exact'} ne '0') {
        $pval = _test_exact($self, $args); # from whatever tail
        $pval = sprintf('%.' . $args->{'precision_p'} . 'f', $pval) if $args->{'precision_p'};
    }
    else {
        my @vals = $self->zscore($args);
        $pval = $vals[1];
    }
    return $pval;
}
*test = \&p_value;
*runs_test = \*p_value;
*rct = \*p_value;

=head2 ztest_ok

Returns true for the loaded sequence if its constituent sample numbers are sufficient for their expected runs to be normally approximated - using Siegal's (1956, p. 140) rule - ok if I<either> of the two I<N>s are greater than 20.

=cut

sub ztest_ok {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my ($n1, $n2) = _get_run_Ns($self, $args);
    #my $sum = $n1 + $n2;
    my $retval = $n1 > 20 || $n2 > 20 ? 1 : 0; # Siegal's rule (p. 140) - ok if either of the two Ns are greater than 20
    return $retval;
}

=head2 stats_hash

 $href = $runs->stats_hash(values => {observed => 1, expected => 1, variance => 1, z_value => 1, p_value => 1}, exact => 0, ccorr => 1);

Returns a hashref for the counts and stats as specified in its "values" argument, and with any options for calculating them (e.g., exact for p_value). See L<Statistics::Sequences/stats_hash> for details. If calling via a "runs" object, the option "stat => 'runs'" is not needed (unlike when using the parent "sequences" object).

=head2 dump

 $runs->dump(values => { observed => 1, variance => 1, p_value => 1}, exact => 1, flag => 1,  precision_s => 3); # among other options

Print Runs-test results to STDOUT. See L<Statistics::Sequences/dump> for details of what stats to dump (default is observed() and p_value()). Optionally also give the data directly.

=cut

sub dump {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    $args->{'stat'} = 'runs';
    $self->SUPER::dump($args);
}

=head2 dump_data

 $runs->dump_data(delim => "\n"); # print whatevers loaded (or specify by label, index, or as "data") 

See L<Statistics::Sequences/dump_data> for details.

=cut

sub _test_exact {
    my ($self, $args) = @_;
    
    # Number of runs per event to use in calculation:
    my ($n1, $n2) = _get_run_Ns($self, $args);
    my $sum = $n1 + $n2;
    my $n1m = $n1 - 1;
    my $n2m = $n2 - 1;
    
    # Calc appropriate rco value based on hypothesis being tested:
    my $rco = defined $args->{'observed'} ? $args->{'observed'} : $self->rco($args);
    my $correct_h = ( $rco - $self->rce($args) >= 0 ) ? 1 : -1; # lower or upper tail test?
    my $wanted_h = $args->{'exact'} eq '1' ? 1 : $args->{'exact'} eq '-1' ? -1 : $correct_h; 
    $rco-- if $wanted_h == 1; # test Hypothesis that rco is greater than expectation
    
    # Calc the probability: from Swed & Eisenhart 1943, p. 66, Eq 1:
    my ($psum, $pval, $i, $k) = (0);
    for ($i = 2; $i <= $rco; $i++) {
        if (is_even($i)) {
            $k = $i / 2;
            $psum += 2 * _choose($n1m, $k -1 ) * _choose($n2m, $k - 1 );
        }
        else {
            $k = ( $i + 1 ) / 2;
            $psum += (
                _choose($n1m, $k - 1) * _choose($n2m, $k - 2)
                +
                _choose($n1m, $k - 2) * _choose($n2m, $k - 1)
            );
        }
    }
    $pval = ( 1 / _choose($sum, $n1) ) * $psum;

    # Form probability according to hypothesis tested:
    $pval = 1 - $pval if $wanted_h == 1; # i.e., if upper-tail test for rco >= $rce
    $pval *= 2 if $args->{'tails'} and $args->{'tails'} == 2;

    return $pval;

    sub _choose { # from Orwant et al., p. 573
        my ($n, $k) = @_;
        my ($res, $j) = (1, 1);
        return 0 if $k > $n || $k < 0;
        $k = ($n - $k) if ($n - $k) < $k;
        while ($j <= $k) {
            $res *= $n--;
            $res /= $j++;
        }
        return $res;
    }
}

sub _get_run_Ns {
   my ($self, $args) = @_;
   return @{$args->{'trials'}} if ref $args->{'trials'};
   return _calc_run_Ns( ref $args->{'data'} ? $args->{'data'} : $self->read($args));
}

sub _calc_run_Ns {
    # find  frequency of all unique elements in array without knowing what elements are - expecting only two but 1 is ok
    my ($data, @vals, %states) = (shift);
    $states{$_}++ for @$data; # init a hash keying each element and its frequency:
    @vals = values %states;# get the two (?) values
    push(@vals, 0) if scalar @vals < 2;
    return @vals;
}

__END__

=head1 EXAMPLE

=head2 Seating at the diner

Swed and Eisenhart (1943) list the occupied (O) and empty (E) seats in a row at a lunch counter. Have people taken up their seats on a random basis?

 use Statistics::Sequences::Runs;
 my $runs = Statistics::Sequences::Runs->new();
 my @seating = (qw/E O E E O E E E O E E E O E O E/); # data already form a single sequence with dichotomous observations
 $runs->dump(data => \@seating, exact => 0, tails => 1);

Suggesting some non-random basis for people taking their seats, this prints:

 observed = 11, p_value = 0.054834

But these data would fail Siegal's rule (L<ztest_ok|Statistics::Sequences::Runs/ztest_ok> = 0) (neither state has 20 observations). So just check exact probability of the hypothesis that the observed deviation is greater than zero (1-tailed):

 $runs->dump(data => \@seating, values => {'p_value'}, exact => 1, tails => 1);

This prints a I<p>-value of .0576923 (so the normal approximation seems good in any case).

These data are also used in an example of testing for L<Vnomes|Statistics::Sequences::Vnomes/EXAMPLE>.

=head2 Runs in multinomial matching

In a single run of a classic ESP test, there are 25 trials, each composed of a randomly generated event (typically, one of 5 possible geometric figures), and a human-generated event arbitrarily drawn from the same pool of alternatives. Tests of the match between the random and human data are typically for number of matches observed versus expected. The I<runs> of matches and misses can be tested by dichotomizing the data on the basis of the L<match|Statistics::Data::Dichotomize/match> of the random "targets" with the human "responses", as described by Kelly (1982):

 use Statistics::Sequences::Runs;
 use Statistics::Data::Dichotomize;
 my @targets = (qw/p c p w s p r w p c r c r s s s s r w p r w c w c/);
 my @responses = (qw/p c s c s s p r w r w c c s s r w s w p c r w p r/);

 # Test for runs of matches between targets and responses:
 my $runs = Statistics::Sequences::Runs->new();
 my $ddat = Statistics::Data::Dichotomize->new();
 $runs->load($ddat->match(data => [\@targets, \@responses]));
 $runs->dump_data(delim => ' '); # have a look at the match sequence; prints "1 1 0 0 1 0 0 0 0 0 0 1 0 1 1 0 0 0 1 1 0 0 0 0 0\n"
 print "Probability of these many runs vs expectation: ", $runs->test(), "\n"; # 0.51436
 # or test for runs in matching when responses are matched to targets one trial behind:
 print $runs->test(data => $ddat->match(data => [\@targets, \@responses], lag => -1)), "\n"; # 0.73766
 
=head1 REFERENCES

These papers provide the implemented algorithms and/or the sample data used in tests. 

Kelly, E. F. (1982). On grouping of hits in some exceptional psi performers. I<Journal of the American Society for Psychical Research>, I<76>, 101-142.

Siegal, S. (1956). I<Nonparametric statistics for the behavioral sciences>. New York, NY, US: McGraw-Hill.

Swed, F., & Eisenhart, C. (1943). Tables for testing randomness of grouping in a sequence of alternatives. I<Annals of Mathematical Statistics>, I<14>, 66-87.

Wald, A., & Wolfowitz, J. (1940). On a test whether two samples are from the same population. I<Annals of Mathematical Statistics>, I<11>, 147-162.

Wolfowitz, J. (1943). On the theory of runs with some applications to quality control. I<Annals of Mathematical Statistics>, I<14>, 280-288. 

=head1 SEE ALSO

L<Statistics::Sequences|Statistics::Sequences> for other tests of sequences, and for sharing data between these tests.

=head1 REVISION HISTORY

See CHANGES in installation dist for revisions.

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2013 Roderick Garton

rgarton AT cpan DOT org

This program is free software. It may be used, redistributed and/or modified under the same terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=back

=head1 DISCLAIMER

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=head1 END

This ends documentation of a Perl implementation of the Wald-Walfowitz Runs test for randomness and group differences within a sequence.

=cut
