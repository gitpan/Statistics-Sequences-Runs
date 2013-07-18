use 5.006;
use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::Sequences::Runs',
    license             => 'perl',
    dist_author         => q{Roderick Garton <rgarton@cpan.org>},
    dist_version_from   => 'lib/Statistics/Sequences/Runs.pm',
    build_requires => {
        'Test::More' => 0,
        'Statistics::Zed' => '0.072',
        'Statistics::Sequences' => '0.11',
        #'Moose' => '2.0604',
        'Number::Misc' => '1.0',
        'List::AllUtils' => '0.03',
    },
    configure_requires => {
        'Module::Build' => '0.36',
    },
    requires => {
        'perl' => 5.006,
    },
    add_to_cleanup      => [ 'Statistics-Sequences-Runs-0.11-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
