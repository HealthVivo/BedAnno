# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl BedAnno.t'

our $data;
our $extradb;
our $config;

BEGIN {
    unless ( grep /blib/, @INC ) {
        chdir 't' if -d 't';
        unshift @INC, '../plugins' if -d '../plugins';
        unshift @INC, '../lib'     if -d '../lib';
        $data    = '../data';
        $extradb = '../db';
	$config  = '../config';
    }
}
$data    ||= "data";
$extradb ||= "db";
$config  ||= "config";
my %opts = (
    db    => "$data/test_db.bed.gz",
    tr    => "$data/test.fas.gz",
    trans => "$data/trans.list",
    genome => "$data/hg19_chM.fa.rz",
    quiet => 1,
);
my %bare_opts = %opts;

use Test::Most;
BEGIN { use_ok('BedAnno') }

my $bare_beda = BedAnno->new(%bare_opts);

my $genomicWalker_on_ins_rep = bless(
    {
        '+' => {
            'ba'  => 'AGTAGTAGT',
            'bal' => 9,
            'bp'  => 14415,
            'br'  => '',
            'brl' => 0
        },
        '-' => {
            'ba'  => 'AGTAGTAGT',
            'bal' => 9,
            'bp'  => 14412,
            'br'  => '',
            'brl' => 0
        },
        'a'      => 'AGTAGTAGTAGT',
        'al'     => 12,
        'alt'    => 'AGTAGTAGTAGT',
        'alt_cn' => 4,
        'altlen' => 12,
        'chr'    => '1',
        'end'    => 14415,
        'guess'  => 'delins',
        'imp'    => 'rep',
        'p'      => 14412,
        'pos'    => 14412,
        'r'      => 'AGT',
        'ref'    => 'AGT',
        'ref_cn' => 1,
        'reflen' => 3,
        'rep'    => 'AGT',
        'replen' => 3,
        'rl'     => 3,
        'sm'     => 2
    },
    'BedAnno::Var'
);


my $genomicWalker_on_span_exon_intron = bless( {
   '+' => {
     'ba' => 'TCCAGGAAGCCT',
     'bal' => 12,
     'bp' => 55248992,
     'br' => '',
     'brl' => 0
   },
   '-' => {
     'ba' => 'TCCAGGAAGCCT',
     'bal' => 12,
     'bp' => 55248980,
     'br' => '',
     'brl' => 0
   },
   'a' => 'TCCAGGAAGCCTTCCAGGAAGCCT',
   'al' => 24,
   'alt' => 'TCCAGGAAGCCTTCCAGGAAGCCT',
   'alt_cn' => 2,
   'altlen' => 24,
   'chr' => '7',
   'end' => 55248992,
   'guess' => 'delins',
   'imp' => 'rep',
   'p' => 55248980,
   'pos' => 55248980,
   'r' => 'TCCAGGAAGCCT',
   'ref' => 'TCCAGGAAGCCT',
   'ref_cn' => 1,
   'reflen' => 12,
   'rep' => 'TCCAGGAAGCCT',
   'replen' => 12,
   'rl' => 12,
   'sm' => 2
 }, 'BedAnno::Var' );

test_genomic_walker(
    "genomic insertion walker span exon intron region",
    $genomicWalker_on_span_exon_intron,
    "7", 55248980, 55248980, "", "TCCAGGAAGCCT"
);

test_genomic_walker(
    "genomic insertion walker for new rep element insertion",
    $genomicWalker_on_ins_rep,
    "1", 14412, 14412, "", "AGTAGTAGT"
);

my $genomicWalker_on_span_exon_intron_varname = 'NM_005228.3: c.2284-5_2290dupTCCAGGAAGCCT (p.A763_Y764insFQEA)';
my ( $t_anno, $noneed ) = $bare_beda->varanno($genomicWalker_on_span_exon_intron);
ok ( $t_anno->{var}->{varName} eq $genomicWalker_on_span_exon_intron_varname.
     "for [genomicWalker_on_span_exon_intron_varname]"
) or explain "The anno info: ", $t_anno;

$bare_beda->DESTROY();
undef $bare_beda;

done_testing();

exit 0;

sub test_genomic_walker {
    my ( $tag, $expect, @args ) = @_;
    my $rvar = BedAnno::Var->new(@args);
    my @walkingRst = $bare_beda->genomicWalker($rvar);
    my $rwalkingRst = $walkingRst[0];
    if ( !is_deeply( $rwalkingRst, $expect, "for [ $tag ]" ) ) {
        explain "The genomic walking result are: ", $rwalkingRst;
    }
}

