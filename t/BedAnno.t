# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl BedAnno.t'

our $data;
our $extradb;

BEGIN {
    unless ( grep /blib/, @INC ) {
        chdir 't' if -d 't';
        unshift @INC, '../plugins' if -d '../plugins';
        unshift @INC, '../lib'     if -d '../lib';
        $data    = '../data';
        $extradb = '../db';
    }
}
$data    ||= "data";
$extradb ||= "db";
my %opts = (
    db    => "$data/test_db.bed.gz",
    tr    => "$data/test.fas.gz",
    batch => 1,
);

if ( -e $extradb and -r $extradb ) {
    if ( -e "$extradb/cytoBand/cytoBand_hg19_grch37.txt.gz" ) {
        $opts{cytoBand} = "$extradb/cytoBand/cytoBand_hg19_grch37.txt.gz";
    }
    if ( -e "$extradb/dbsnp/snp137.bed.gz" ) {
        $opts{dbSNP} = "$extradb/dbsnp/snp137.bed.gz";
    }
    if ( -e "$extradb/tgp/tgp_phaseI_version3_hg19.dbdump.bed.gz" ) {
        $opts{tgp} = "$extradb/tgp/tgp_phaseI_version3_hg19.dbdump.bed.gz";
    }
    if ( -e "$extradb/CG/CG54_20130709_stats_filtered_ct0.tsv.gz" ) {
        $opts{cg54} = "$extradb/CG/CG54_20130709_stats_filtered_ct0.tsv.gz";
    }
    if ( -e "$extradb/NHLBI/ESP6500SI-V2-SSA137.NHLBI.bed.rmanchor.uniq.gz" ) {
        $opts{esp6500} =
          "$extradb/NHLBI/ESP6500SI-V2-SSA137.NHLBI.bed.rmanchor.uniq.gz";
    }
    if ( -e "$extradb/pfam/Pfam-A-ncbi_2012-12-21.bed.gz" ) {
        $opts{pfam} = "$extradb/pfam/Pfam-A-ncbi_2012-12-21.bed.gz";
    }
    if ( -e "$extradb/predictions/predictDB_for_anno104.tab.gz" ) {
        $opts{prediction} = "$extradb/predictions/predictDB_for_anno104.tab.gz";
    }
    if ( -e "$extradb/phyloP/phyloP_all3class_combin_2013-09-25.bed.gz" ) {
        $opts{phyloP} =
          "$extradb/phyloP/phyloP_all3class_combin_2013-09-25.bed.gz";
    }
}

use Test::Most;
BEGIN { use_ok('BedAnno') }

my $beda = BedAnno->new(%opts);

#explain "The database are:", $beda;
# ncRNA part
my $crawler_input = {
    chr   => 1,
    begin => 14409,
    end   => 14410,
    referenceSequence   => 'c',
    variantSequence   => 'a',
};

my $crawler_input2 = {
    chr   => 1,
    begin => 14410,
    referenceSequence   => "C",
    variantSequence   => "CGAATAGCTA",
};
my $crawler_input3 = {
    chr	  => "1",
    begin => 14410,
    end   => 14417,
    referenceSequence => "TAGATCG",
    variantSequence => undef,
};

my $crawler_input4 = {
    chr	  => "1",
    begin => "14410",
    end	  => "14410",
    variantSequence => "GAATAGCTA",
};

my $snv_parse = bless(
    {
        'alt'    => 'A',
        'altlen' => 1,
        'chr'    => '1',
        'end'    => 14410,
        'guess'  => 'snv',
        'imp'    => 'snv',
        'pos'    => 14409,
        'ref'    => 'C',
        'reflen' => 1,
        'sm'     => 1
    },
    'BedAnno::Var'
);

my $insert_parse = bless(
    {
        'alt'    => 'GAATAGCTA',
        'altlen' => 9,
        'chr'    => '1',
        'end'    => 14410,
        'guess'  => 'ins',
        'imp'    => 'ins',
        'pos'    => 14410,
        'ref'    => '',
        'reflen' => 0,
        'sm'     => 0
    },
    'BedAnno::Var'
);

my $del_parse = bless(
    {
        'alt'    => '',
        'altlen' => 0,
        'chr'    => '1',
        'end'    => 14417,
        'guess'  => 'del',
        'imp'    => 'del',
        'pos'    => 14410,
        'ref'    => 'TAGATCG',
        'reflen' => 7,
        'sm'     => 2
    },
    'BedAnno::Var'
);

my $delins_parse = bless(
    {
        'alt'    => 'GT',
        'altlen' => 2,
        'chr'    => '1',
        'end'    => 14414,
        'guess'  => 'delins',
        'imp'    => 'delins',
        'pos'    => 14409,
        'ref'    => 'CTAGA',
        'reflen' => 5,
        'sm'     => 2
    },
    'BedAnno::Var'
);

my $rep_parse = bless(
    {
        '+' => {
            'ba'  => 'TAGTAGTAG',
            'bal' => 9,
            'bp'  => 14413,
            'br'  => '',
            'brl' => 0
        },
        '-' => {
            'ba'  => 'TAGTAGTAG',
            'bal' => 9,
            'bp'  => 14410,
            'br'  => '',
            'brl' => 0
        },
        'a'      => 'TAGTAGTAGTAG',
        'al'     => 12,
        'alt'    => 'TAGTAGTAGTAGA',
        'alt_cn' => 4,
        'altlen' => 13,
        'chr'    => '1',
        'end'    => 14414,
        'guess'  => 'delins',
        'imp'    => 'rep',
        'p'      => 14410,
        'pos'    => 14410,
        'r'      => 'TAG',
        'ref'    => 'TAGA',
        'ref_cn' => 1,
        'reflen' => 4,
        'rep'    => 'TAG',
        'replen' => 3,
        'rl'     => 3,
        'sm'     => 2
    },
    'BedAnno::Var'
);

my $subs_parse = bless(
    {
        'alt'      => 'TGTGT',
        'altlen'   => 5,
        'chr'      => '1',
        'end'      => 14414,
        'guess'    => 'delins',
        'imp'      => 'delins',
        'pos'      => 14409,
        'ref'      => 'CATGA',
        'reflen'   => 5,
        'sep_snvs' => [ 14410, 14411, 14414 ],
        'sm'       => 3
    },
    'BedAnno::Var'
);

my $no_call_parse = bless(
    {
        'alt'    => '?',
        'chr'    => '1',
        'end'    => 14410,
        'guess'  => 'no-call',
        'imp'    => 'ins',
        'pos'    => 14410,
        'ref'    => '=',
        'reflen' => 0,
        'sm'     => 0
    },
    'BedAnno::Var'
);

my $no_call_edge_parse = bless(
    {
        'reflen'    => 10000,
        'imp'       => 'delins',
        'pos'       => 0,
        'guess' => 'no-call',
        'alt'   => '?',
        'ref'   => '=',
        'chr'   => '1',
        'sm'       => 3,
        'end'      => 10000,
    },
    'BedAnno::Var'
);

# For coding RNA test
my $cds_rna_delins = bless(
    {
        'alt'    => 'A',
        'altlen' => 1,
        'chr'    => '8',
        'end'    => 24811066,
        'guess'  => 'snv',
        'imp'    => 'snv',
        'pos'    => 24811065,
        'ref'    => 'G',
        'reflen' => 1,
        'sm'     => 1
    },
    'BedAnno::Var'
);
my $cds_no_change = bless(
    {
        'alt'    => '',
        'altlen' => 0,
        'chr'    => '8',
        'end'    => 24811065,
        'guess'  => 'del',
        'imp'    => 'del',
        'pos'    => 24811064,
        'ref'    => 'G',
        'reflen' => 1,
        'sm'     => 1
    },
    'BedAnno::Var'
);
my $cds_snv = bless(
    {
        'alt'    => 'T',
        'altlen' => 1,
        'chr'    => '8',
        'end'    => 24811067,
        'guess'  => 'snv',
        'imp'    => 'snv',
        'pos'    => 24811066,
        'ref'    => 'G',
        'reflen' => 1,
        'sm'     => 1
    },
    'BedAnno::Var'
);

my $cds_del = bless(
    {
        'alt'    => '',
        'altlen' => 0,
        'chr'    => '8',
        'end'    => 24811067,
        'guess'  => 'del',
        'imp'    => 'del',
        'pos'    => 24811066,
        'ref'    => 'G',
        'reflen' => 1,
        'sm'     => 1
    },
    'BedAnno::Var'
);

my $large_del = bless(
    {
        'alt'    => '',
        'altlen' => 0,
        'chr'    => '1',
        'end'    => 100327884,
        'guess'  => 'del',
        'imp'    => 'del',
        'pos'    => 100327863,
        'ref'    => 'TGGTGCTGATAATCATGTGCT',
        'reflen' => 21,
        'sm'     => 2
    },
    'BedAnno::Var'
);

my $cds_ins = bless(
    {
        'alt'    => 'GGG',
        'altlen' => 3,
        'chr'    => '8',
        'end'    => 24811067,
        'guess'  => 'ins',
        'imp'    => 'ins',
        'pos'    => 24811067,
        'ref'    => '',
        'reflen' => 0,
        'sm'     => 0
    },
    'BedAnno::Var'
);
my $cds_rep = bless(
    {
        '+' => {
            'ba'  => '',
            'bal' => 0,
            'bp'  => 49395685,
            'br'  => 'GCCGCC',
            'brl' => 6
        },
        '-' => {
            'ba'  => '',
            'bal' => 0,
            'bp'  => 49395673,
            'br'  => 'GCCGCC',
            'brl' => 6
        },
        'a'      => 'GCCGCCGCCGCC',
        'al'     => 12,
        'alt'    => 'GCCGCCGCCGCC',
        'alt_cn' => 4,
        'altlen' => 12,
        'chr'    => '3',
        'end'    => 49395691,
        'guess'  => 'delins',
        'imp'    => 'rep',
        'p'      => 49395673,
        'pos'    => 49395673,
        'r'      => 'GCCGCCGCCGCCGCCGCC',
        'ref'    => 'GCCGCCGCCGCCGCCGCC',
        'ref_cn' => 6,
        'reflen' => 18,
        'rep'    => 'GCC',
        'replen' => 3,
        'rl'     => 18,
        'sm'     => 2
    },
    'BedAnno::Var'
);

my $cds_delins = bless(
    {
        'alt'      => 'TC',
        'altlen'   => 2,
        'chr'      => '1',
        'end'      => 865534,
        'guess'    => 'delins',
        'imp'      => 'delins',
        'pos'      => 865532,
        'ref'      => 'AG',
        'reflen'   => 2,
        'sep_snvs' => [ 865533, 865534 ],
        'sm'       => 3
    },
    'BedAnno::Var'
);
my $cds_no_call = bless(
    {
        'alt'    => 'N',
        'altlen' => 1,
        'chr'    => '1',
        'end'    => 861320,
        'guess'  => 'no-call',
        'imp'    => 'snv',
        'pos'    => 861319,
        'ref'    => 'T',
        'reflen' => 1,
        'sm'     => 1
    },
    'BedAnno::Var'
);

my $cds_rna_snv1 = bless(
    {
        'alt'    => 'T',
        'altlen' => 1,
        'chr'    => '3',
        'end'    => 49395709,
        'guess'  => 'snv',
        'imp'    => 'snv',
        'pos'    => 49395708,
        'ref'    => 'C',
        'reflen' => 1,
        'sm'     => 1
    },
    'BedAnno::Var'
);
my $cds_rna_snv2 = bless(
    {
        'alt'    => 'A',
        'altlen' => 1,
        'chr'    => '3',
        'end'    => 49395705,
        'guess'  => 'snv',
        'imp'    => 'snv',
        'pos'    => 49395704,
        'ref'    => 'C',
        'reflen' => 1,
        'sm'     => 1
    },
    'BedAnno::Var'
);

my $snv_varanno = bless(
    {
        'trInfo' => {
            'NR_024540.1' => {
                'c'             => 'n.1721G>T',
                'cdsBegin'      => '',
                'cdsEnd'        => '',
                'ei_Begin'      => 'EX11E',
                'ei_End'        => 'EX11E',
                'exin'          => 'EX11E',
                'exonIndex'     => '11',
                'func'          => 'ncRNA',
                'funcSO'        => '',
                'funcSOname'    => 'unknown',
                'geneId'        => '653635',
                'geneSym'       => 'WASH7P',
                'genepart'      => 'ncRNA',
                'genepartIndex' => '11',
                'genepartSO'    => 'SO:0000655',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1722,
                    'r'    => 'R11E'
                },
                'preStart' => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1720,
                    'r'    => 'R11E'
                },
		'primaryTag'=> 'Y',
                'prot'      => '',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'R11E',
                'r_Begin'   => 'R11E',
                'r_End'     => 'R11E',
                'rnaBegin'  => 1721,
                'rnaEnd'    => 1721,
                'strd'      => '-',
                'trAlt'     => 'T',
                'trRef'     => 'G',
                'trRefComp' => {
                    'EX11E' => 1
                }
            }
        },
        'var' => bless(
            {
                'alt'       => 'A',
                'altlen'    => 1,
                'chr'       => '1',
                'end'       => 14410,
                'gHGVS'     => 'g.14409C>A',
                'guess'     => 'snv',
                'imp'       => 'snv',
                'pos'       => 14409,
                'ref'       => 'C',
                'refbuild'  => 'GRCh37',
                'reflen'    => 1,
                'sm'        => 1,
                'varTypeSO' => 'SO:0001483'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $ins_varanno = bless(
    {
        'trInfo' => {
            'NR_024540.1' => {
                'c'             => 'n.1720_1721insTAGCTATTC',
                'cdsBegin'      => '',
                'cdsEnd'        => '',
                'ei_Begin'      => 'EX11E',
                'ei_End'        => 'EX11E',
                'exin'          => 'EX11E',
                'exonIndex'     => '11',
                'func'          => 'ncRNA',
                'funcSO'        => '',
                'funcSOname'    => 'unknown',
                'geneId'        => '653635',
                'geneSym'       => 'WASH7P',
                'genepart'      => 'ncRNA',
                'genepartIndex' => '11',
                'genepartSO'    => 'SO:0000655',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1721,
                    'r'    => 'R11E'
                },
                'preStart' => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1720,
                    'r'    => 'R11E'
                },
		'primaryTag'=> 'Y',
                'prot'      => '',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'R11E',
                'r_Begin'   => 'R11E',
                'r_End'     => 'R11E',
                'rnaBegin'  => 1721,
                'rnaEnd'    => 1720,
                'strd'      => '-',
                'trAlt'     => 'TAGCTATTC',
                'trRef'     => '',
                'trRefComp' => {
                    'EX11E' => 0
                }
            }
        },
        'var' => bless(
            {
                'alt'       => 'GAATAGCTA',
                'altlen'    => 9,
                'chr'       => '1',
                'end'       => 14410,
                'gHGVS'     => 'g.14410_14411insGAATAGCTA',
                'guess'     => 'ins',
                'imp'       => 'ins',
                'pos'       => 14410,
                'ref'       => '',
                'refbuild'  => 'GRCh37',
                'reflen'    => 0,
                'sm'        => 0,
                'varTypeSO' => 'SO:0000667'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $del_varanno = bless(
    {
        'trInfo' => {
            'NR_024540.1' => {
                'c'             => 'n.1714_1720del',
                'cdsBegin'      => '',
                'cdsEnd'        => '',
                'ei_Begin'      => 'EX11E',
                'ei_End'        => 'EX11E',
                'exin'          => 'EX11E',
                'exonIndex'     => '11',
                'func'          => 'ncRNA',
                'funcSO'        => '',
                'funcSOname'    => 'unknown',
                'geneId'        => '653635',
                'geneSym'       => 'WASH7P',
                'genepart'      => 'ncRNA',
                'genepartIndex' => '11',
                'genepartSO'    => 'SO:0000655',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1721,
                    'r'    => 'R11E'
                },
                'preStart' => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1713,
                    'r'    => 'R11E'
                },
		'primaryTag'=> 'Y',
                'prot'      => '',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'R11E',
                'r_Begin'   => 'R11E',
                'r_End'     => 'R11E',
                'rnaBegin'  => 1714,
                'rnaEnd'    => 1720,
                'strd'      => '-',
                'trAlt'     => '',
                'trRef'     => 'GAACTGA',
                'trRefComp' => {
                    'EX11E' => 7
                }
            }
        },
        'var' => bless(
            {
                'alt'       => '',
                'altlen'    => 0,
                'chr'       => '1',
                'end'       => 14417,
                'gHGVS'     => 'g.14411_14417delTAGATCG',
                'guess'     => 'del',
                'imp'       => 'del',
                'pos'       => 14410,
                'ref'       => 'TAGATCG',
                'refbuild'  => 'GRCh37',
                'reflen'    => 7,
                'sm'        => 2,
                'varTypeSO' => 'SO:0000159'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $rep_varanno = bless(
    {
        'trInfo' => {
            'NR_024540.1' => {
                'c'             => 'n.1718_1719delinsCTACTACTACT',
                'cdsBegin'      => '',
                'cdsEnd'        => '',
                'ei_Begin'      => 'EX11E',
                'ei_End'        => 'EX11E',
                'exin'          => 'EX11E',
                'exonIndex'     => '11',
                'func'          => 'ncRNA',
                'funcSO'        => '',
                'funcSOname'    => 'unknown',
                'geneId'        => '653635',
                'geneSym'       => 'WASH7P',
                'genepart'      => 'ncRNA',
                'genepartIndex' => '11',
                'genepartSO'    => 'SO:0000655',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1721,
                    'r'    => 'R11E'
                },
                'preStart' => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1717,
                    'r'    => 'R11E'
                },
		'primaryTag'=> 'Y',
                'prot'      => '',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'R11E',
                'r_Begin'   => 'R11E',
                'r_End'     => 'R11E',
                'rnaBegin'  => 1718,
                'rnaEnd'    => 1720,
                'strd'      => '-',
                'trAlt'     => 'CTACTACTACTA',
                'trRef'     => 'TGA',
                'trRefComp' => {
                    'EX11E' => 3
                }
            }
        },
        'var' => bless(
            {
                '+' => {
                    'ba'  => 'TAGTAGTAG',
                    'bal' => 9,
                    'bp'  => 14413,
                    'br'  => '',
                    'brl' => 0
                },
                '-' => {
                    'ba'  => 'TAGTAGTAG',
                    'bal' => 9,
                    'bp'  => 14410,
                    'br'  => '',
                    'brl' => 0
                },
                'a'         => 'TAGTAGTAGTAG',
                'al'        => 12,
                'alt'       => 'TAGTAGTAGTAGA',
                'alt_cn'    => 4,
                'altlen'    => 13,
                'chr'       => '1',
                'end'       => 14414,
                'gHGVS'     => 'g.14411TAG[1>4]',
                'guess'     => 'delins',
                'imp'       => 'rep',
                'p'         => 14410,
                'pos'       => 14410,
                'r'         => 'TAG',
                'ref'       => 'TAGA',
                'ref_cn'    => 1,
                'refbuild'  => 'GRCh37',
                'reflen'    => 4,
                'rep'       => 'TAG',
                'replen'    => 3,
                'rl'        => 3,
                'sm'        => 2,
                'varTypeSO' => 'SO:1000032'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $delins_varanno = bless(
    {
        'trInfo' => {
            'NR_024540.1' => {
                'c'             => 'n.1717_1721delinsAC',
                'cdsBegin'      => '',
                'cdsEnd'        => '',
                'ei_Begin'      => 'EX11E',
                'ei_End'        => 'EX11E',
                'exin'          => 'EX11E',
                'exonIndex'     => '11',
                'func'          => 'ncRNA',
                'funcSO'        => '',
                'funcSOname'    => 'unknown',
                'geneId'        => '653635',
                'geneSym'       => 'WASH7P',
                'genepart'      => 'ncRNA',
                'genepartIndex' => '11',
                'genepartSO'    => 'SO:0000655',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1722,
                    'r'    => 'R11E'
                },
                'preStart' => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1716,
                    'r'    => 'R11E'
                },
		'primaryTag'=> 'Y',
                'prot'      => '',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'R11E',
                'r_Begin'   => 'R11E',
                'r_End'     => 'R11E',
                'rnaBegin'  => 1717,
                'rnaEnd'    => 1721,
                'strd'      => '-',
                'trAlt'     => 'AC',
                'trRef'     => 'CTGAG',
                'trRefComp' => {
                    'EX11E' => 5
                }
            }
        },
        'var' => bless(
            {
                'alt'       => 'GT',
                'altlen'    => 2,
                'chr'       => '1',
                'end'       => 14414,
                'gHGVS'     => 'g.14410_14414delCTAGAinsGT',
                'guess'     => 'delins',
                'imp'       => 'delins',
                'pos'       => 14409,
                'ref'       => 'CTAGA',
                'refbuild'  => 'GRCh37',
                'reflen'    => 5,
                'sm'        => 2,
                'varTypeSO' => 'SO:1000032'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $subs_varanno = bless(
    {
        'trInfo' => {
            'NR_024540.1' => {
                'c'             => 'n.1717_1721delinsACACA',
                'cdsBegin'      => '',
                'cdsEnd'        => '',
                'ei_Begin'      => 'EX11E',
                'ei_End'        => 'EX11E',
                'exin'          => 'EX11E',
                'exonIndex'     => '11',
                'func'          => 'ncRNA',
                'funcSO'        => '',
                'funcSOname'    => 'unknown',
                'geneId'        => '653635',
                'geneSym'       => 'WASH7P',
                'genepart'      => 'ncRNA',
                'genepartIndex' => '11',
                'genepartSO'    => 'SO:0000655',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1722,
                    'r'    => 'R11E'
                },
                'preStart' => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1716,
                    'r'    => 'R11E'
                },
		'primaryTag'=> 'Y',
                'prot'      => '',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'R11E',
                'r_Begin'   => 'R11E',
                'r_End'     => 'R11E',
                'rnaBegin'  => 1717,
                'rnaEnd'    => 1721,
                'strd'      => '-',
                'trAlt'     => 'ACACA',
                'trRef'     => 'CTGAG',
                'trRefComp' => {
                    'EX11E' => 5
                }
            }
        },
        'var' => bless(
            {
                'alt'       => 'TGTGT',
                'altlen'    => 5,
                'chr'       => '1',
                'end'       => 14414,
                'gHGVS'     => 'g.14410_14414delCATGAinsTGTGT',
                'guess'     => 'delins',
                'imp'       => 'delins',
                'pos'       => 14409,
                'ref'       => 'CATGA',
                'refbuild'  => 'GRCh37',
                'reflen'    => 5,
                'sep_snvs'  => [ 14410, 14411, 14414 ],
                'sm'        => 3,
                'varTypeSO' => 'SO:1000032'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $no_call_varanno = bless(
    {
        'trInfo' => {
            'NR_024540.1' => {
                'c'             => 'n.1720_1721ins?',
                'cdsBegin'      => '',
                'cdsEnd'        => '',
                'ei_Begin'      => 'EX11E',
                'ei_End'        => 'EX11E',
                'exin'          => 'EX11E',
                'exonIndex'     => '11',
                'func'          => 'unknown-no-call',
                'funcSO'        => '',
                'funcSOname'    => 'unknown-no-call',
                'geneId'        => '653635',
                'geneSym'       => 'WASH7P',
                'genepart'      => 'ncRNA',
                'genepartIndex' => '11',
                'genepartSO'    => 'SO:0000655',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1721,
                    'r'    => 'R11E'
                },
                'preStart' => {
                    'cDot' => '',
                    'exin' => 'EX11E',
                    'nDot' => 1720,
                    'r'    => 'R11E'
                },
		'primaryTag'=> 'Y',
                'prot'      => '',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'R11E',
                'r_Begin'   => 'R11E',
                'r_End'     => 'R11E',
                'rnaBegin'  => 1721,
                'rnaEnd'    => 1720,
                'strd'      => '-',
                'trAlt'     => '?',
                'trRef'     => '=',
                'trRefComp' => {
                    'EX11E' => 0
                }
            }
        },
        'var' => bless(
            {
                'alt'       => '?',
                'chr'       => '1',
                'end'       => 14410,
                'gHGVS'     => 'g.14410_14411ins?',
                'guess'     => 'no-call',
                'imp'       => 'ins',
                'pos'       => 14410,
                'ref'       => '=',
                'refbuild'  => 'GRCh37',
                'reflen'    => 0,
                'sm'        => 0,
                'varTypeSO' => 'no-call'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $no_call_edge_varanno = bless(
    {
        'trInfo' => {
            'NR_046018.2' => {
                'r'             => 'PROM',
                'protBegin'     => '',
                'cdsBegin'      => '',
                'genepartIndex' => 0,
                'ei_End'        => '.',
                'trRef'         => '=',
                'protEnd'       => '',
                'cdsEnd'        => '',
                'exonIndex'     => '.',
                'trAlt'         => '?',
                'r_Begin'       => 'PROM',
                'intronIndex'   => '.',
                'exin'          => '.',
                'prot'          => '',
                'strd'          => '+',
                'trRefComp'     => {
                    'P0' => [ 0, 10000 ]
                },
                'geneSym' => 'DDX11L1',
                'r_End'   => 'PROM',
                'postEnd' => {
                    'r'    => 'PROM',
                    'cDot' => '',
                    'nDot' => -1873,
                    'exin' => '.'
                },
                'c'          => 'n.-11873_-1874delins?',
                'geneId'     => '100287102',
                'rnaEnd'     => -1874,
                'genepart'   => 'promoter',
                'func'       => 'unknown-no-call',
                'rnaBegin'   => -11873,
                'funcSO'     => '',
                'ei_Begin'   => '.',
                'genepartSO' => 'SO:0000167',
                'funcSOname' => 'unknown-no-call',
                'primaryTag' => 'Y'
            }
        },
        'var' => bless(
            {
                'refbuild'  => 'GRCh37',
                'varTypeSO' => 'no-call',
                'reflen'    => 10000,
                'imp'       => 'delins',
                'pos'       => 0,
                'guess' => 'no-call',
                'gHGVS' => 'g.1_10000delins?',
                'alt'   => '?',
                'ref'   => '=',
                'chr'   => '1',
                'sm'       => 3,
                'end'      => 10000,
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);


my $cds_rna_snv1_anno = {
    'c'             => 'c.3G>A',
    'cdsBegin'      => '3',
    'cdsEnd'        => '3',
    'ei_Begin'      => 'EX1',
    'ei_End'        => 'EX1',
    'exin'          => 'EX1',
    'exonIndex'     => '1',
    'func'          => 'init-loss',
    'funcSO'        => 'SO:0001582',
    'funcSOname'    => 'initiator_codon_variant',
    'geneId'        => '2876',
    'geneSym'       => 'GPX1',
    'genepart'      => 'CDS',
    'genepartIndex' => '1',
    'genepartSO'    => 'SO:0000316',
    'intronIndex'   => '.',
    'p'             => 'p.0?',
    'postEnd'       => {
        'cDot' => '4',
        'exin' => 'EX1',
        'nDot' => 84,
        'r'    => 'C1'
    },
    'prAlt'    => 'I',
    'prRef'    => 'M',
    'preStart' => {
        'cDot' => '2',
        'exin' => 'EX1',
        'nDot' => 82,
        'r'    => 'C1'
    },
    'primaryTag' => 'Y',
    'prot'       => 'NP_000572.2',
    'protBegin'  => 1,
    'protEnd'    => 1,
    'r'          => 'C1',
    'r_Begin'    => 'C1',
    'r_End'      => 'C1',
    'rnaBegin'   => 83,
    'rnaEnd'     => 83,
    'strd'       => '-',
    'trAlt'      => 'A',
    'trRef'      => 'G',
    'trRefComp'  => {
        'EX1' => 1
    }
};

my $cds_rna_snv2_anno = {
    'c'             => 'c.7G>T',
    'cc'            => 'GCT=>TCT',
    'cdsBegin'      => '7',
    'cdsEnd'        => '7',
    'ei_Begin'      => 'EX1',
    'ei_End'        => 'EX1',
    'exin'          => 'EX1',
    'exonIndex'     => '1',
    'func'          => 'missense',
    'funcSO'        => 'SO:0001583',
    'funcSOname'    => 'missense_variant',
    'geneId'        => '2876',
    'geneSym'       => 'GPX1',
    'genepart'      => 'CDS',
    'genepartIndex' => '1',
    'genepartSO'    => 'SO:0000316',
    'intronIndex'   => '.',
    'p'             => 'p.A3S',
    'polar'         => 'NP=>P0',
    'postEnd'       => {
        'cDot' => '8',
        'exin' => 'EX1',
        'nDot' => 88,
        'r'    => 'C1'
    },
    'prAlt'    => 'S',
    'prRef'    => 'A',
    'preStart' => {
        'cDot' => '6',
        'exin' => 'EX1',
        'nDot' => 86,
        'r'    => 'C1'
    },
    'primaryTag' => 'Y',
    'prot'       => 'NP_000572.2',
    'protBegin'  => 3,
    'protEnd'    => 3,
    'r'          => 'C1',
    'r_Begin'    => 'C1',
    'r_End'      => 'C1',
    'rnaBegin'   => 87,
    'rnaEnd'     => 87,
    'strd'       => '-',
    'trAlt'      => 'T',
    'trRef'      => 'G',
    'trRefComp'  => {
        'EX1' => 1
    }
};

my $cds_rna_delins_anno = bless(
    {
        'trInfo' => {
            'NM_006158.3' => {
                'c'             => 'c.1412_1413insT',
                'cdsBegin'      => '1413',
                'cdsEnd'        => '1413',
                'ei_Begin'      => 'EX3',
                'ei_End'        => 'EX3',
                'exin'          => 'EX3',
                'exonIndex'     => '3',
                'func'          => 'frameshift',
                'funcSO'        => 'SO:0001589',
                'funcSOname'    => 'frameshift_variant',
                'geneId'        => '4747',
                'geneSym'       => 'NEFL',
                'genepart'      => 'CDS',
                'genepartIndex' => '3',
                'genepartSO'    => 'SO:0000316',
                'intronIndex'   => '.',
                'p'             => 'p.S472fs*2',
		'prAlt'		=> 'PL*',
		'prRef'		=> 'PSEGEAEEEEKDKEEAEEEEAAEEEEAAKEESEEAKEEEEGGEGEEGEETKEAEEEEKKVEGAGEEQAAKKKD*',
                'postEnd'       => {
                    'cDot' => '1414',
                    'exin' => 'EX3',
                    'nDot' => '1516',
                    'r'    => 'C3'
                },
                'preStart' => {
                    'cDot' => '1412',
                    'exin' => 'EX3',
                    'nDot' => 1514,
                    'r'    => 'C3'
                },
		'primaryTag'=> 'N',
                'prot'      => 'NP_006149.2',
                'protBegin' => 471,
                'protEnd'   => 544,
                'r'         => 'C3',
                'r_Begin'   => 'C3',
                'r_End'     => 'C3',
                'rnaBegin'  => 1515,
                'rnaEnd'    => '1515',
                'strd'      => '-',
                'trAlt'     => 'TC',
                'trRef'     => 'C',
                'trRefComp' => {
                    'EX3' => 1
                }
            }
        },
        'var' => bless(
            {
                'alt'       => 'A',
                'altlen'    => 1,
                'chr'       => '8',
                'end'       => 24811066,
                'gHGVS'     => 'g.24811065G>A',
                'guess'     => 'snv',
                'imp'       => 'snv',
                'pos'       => 24811065,
                'ref'       => 'G',
                'refbuild'  => 'GRCh37',
                'reflen'    => 1,
                'sm'        => 1,
                'varTypeSO' => 'SO:0001483'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);
my $cds_no_change_anno = bless(
    {
        'trInfo' => {
            'NM_006158.3' => {
                'c'             => 'c.=',
                'cdsBegin'      => '1414',
                'cdsEnd'        => '1413',
                'ei_Begin'      => 'EX3',
                'ei_End'        => 'EX3',
                'exin'          => 'EX3',
                'exonIndex'     => '3',
                'func'          => 'no-change',
                'funcSO'        => '',
                'funcSOname'    => 'no-change',
                'geneId'        => '4747',
                'geneSym'       => 'NEFL',
                'genepart'      => 'CDS',
                'genepartIndex' => '3',
                'genepartSO'    => 'SO:0000316',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '1414',
                    'exin' => 'EX3',
                    'nDot' => 1516,
                    'r'    => 'C3'
                },
                'preStart' => {
                    'cDot' => '1413',
                    'exin' => 'EX3',
                    'nDot' => 1515,
                    'r'    => 'C3'
                },
		'primaryTag'=> 'N',
                'prot'      => 'NP_006149.2',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'C3',
                'r_Begin'   => 'C3',
                'r_End'     => 'C3',
                'rnaBegin'  => '1516',
                'rnaEnd'    => '1515',
                'strd'      => '-',
                'trAlt'     => '',
                'trRef'     => '',
                'trRefComp' => {
                    'EX3' => 0
                }
            }
        },
        'var' => bless(
            {
                'alt'       => '',
                'altlen'    => 0,
                'chr'       => '8',
                'end'       => 24811065,
                'gHGVS'     => 'g.24811065delG',
                'guess'     => 'del',
                'imp'       => 'del',
                'pos'       => 24811064,
                'ref'       => 'G',
                'refbuild'  => 'GRCh37',
                'reflen'    => 1,
                'sm'        => 1,
                'varTypeSO' => 'SO:0000159'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $cds_snv_anno = bless(
    {
        'trInfo' => {
            'NM_006158.3' => {
                'c'             => 'c.1412C>A',
                'cc'            => 'CCC=>CAC',
                'cdsBegin'      => '1412',
                'cdsEnd'        => '1412',
                'ei_Begin'      => 'EX3',
                'ei_End'        => 'EX3',
                'exin'          => 'EX3',
                'exonIndex'     => '3',
                'func'          => 'missense',
                'funcSO'        => 'SO:0001583',
                'funcSOname'    => 'missense_variant',
                'geneId'        => '4747',
                'geneSym'       => 'NEFL',
                'genepart'      => 'CDS',
                'genepartIndex' => '3',
                'genepartSO'    => 'SO:0000316',
                'intronIndex'   => '.',
                'p'             => 'p.P471H',
                'polar'         => 'NP=>P+',
                'postEnd'       => {
                    'cDot' => '1413',
                    'exin' => 'EX3',
                    'nDot' => 1515,
                    'r'    => 'C3'
                },
                'preStart' => {
                    'cDot' => '1411',
                    'exin' => 'EX3',
                    'nDot' => 1513,
                    'r'    => 'C3'
                },
		'primaryTag'=> 'N',
		'prAlt'	    => 'H',
		'prRef'	    => 'P',
                'prot'      => 'NP_006149.2',
                'protBegin' => 471,
                'protEnd'   => 471,
                'r'         => 'C3',
                'r_Begin'   => 'C3',
                'r_End'     => 'C3',
                'rnaBegin'  => 1514,
                'rnaEnd'    => 1514,
                'strd'      => '-',
                'trAlt'     => 'A',
                'trRef'     => 'C',
                'trRefComp' => {
                    'EX3' => 1
                }
            }
        },
        'var' => bless(
            {
                'alt'       => 'T',
                'altlen'    => 1,
                'chr'       => '8',
                'end'       => 24811067,
                'gHGVS'     => 'g.24811066G>T',
                'guess'     => 'snv',
                'imp'       => 'snv',
                'pos'       => 24811066,
                'ref'       => 'G',
                'refbuild'  => 'GRCh37',
                'reflen'    => 1,
                'sm'        => 1,
                'varTypeSO' => 'SO:0001483'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);
my $cds_del_anno = bless(
    {
        'trInfo' => {
            'NM_006158.3' => {
                'c'             => 'c.1412delC',
                'cdsBegin'      => '1412',
                'cdsEnd'        => '1412',
                'ei_Begin'      => 'EX3',
                'ei_End'        => 'EX3',
                'exin'          => 'EX3',
                'exonIndex'     => '3',
                'func'          => 'frameshift',
                'funcSO'        => 'SO:0001589',
                'funcSOname'    => 'frameshift_variant',
                'geneId'        => '4747',
                'geneSym'       => 'NEFL',
                'genepart'      => 'CDS',
                'genepartIndex' => '3',
                'genepartSO'    => 'SO:0000316',
                'intronIndex'   => '.',
                'p'             => 'p.S472fs*78',
                'postEnd'       => {
                    'cDot' => '1413',
                    'exin' => 'EX3',
                    'nDot' => 1515,
                    'r'    => 'C3'
                },
                'preStart' => {
                    'cDot' => '1411',
                    'exin' => 'EX3',
                    'nDot' => 1513,
                    'r'    => 'C3'
                },
		'prAlt' => 'PLKEKPRRRRRTRKRPRKRRQLKRKKLPRKSLKKQKKKKKEVKVKKERKPKKLKRRRRKLKVLGRNKQLRRKIEPPFP*',
		'prRef' => 'PSEGEAEEEEKDKEEAEEEEAAEEEEAAKEESEEAKEEEEGGEGEEGEETKEAEEEEKKVEGAGEEQAAKKKD*',
		'primaryTag'=> 'N',
                'prot'      => 'NP_006149.2',
                'protBegin' => 471,
                'protEnd'   => 544,
                'r'         => 'C3',
                'r_Begin'   => 'C3',
                'r_End'     => 'C3',
                'rnaBegin'  => 1514,
                'rnaEnd'    => 1514,
                'strd'      => '-',
                'trAlt'     => '',
                'trRef'     => 'C',
                'trRefComp' => {
                    'EX3' => 1
                }
            }
        },
        'var' => bless(
            {
                'alt'       => '',
                'altlen'    => 0,
                'chr'       => '8',
                'end'       => 24811067,
                'gHGVS'     => 'g.24811067delG',
                'guess'     => 'del',
                'imp'       => 'del',
                'pos'       => 24811066,
                'ref'       => 'G',
                'refbuild'  => 'GRCh37',
                'reflen'    => 1,
                'sm'        => 1,
                'varTypeSO' => 'SO:0000159'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $large_del_anno = {
    'c'             => 'c.345_365del',
    'cdsBegin'      => '345',
    'cdsEnd'        => '365',
    'ei_Begin'      => 'EX4',
    'ei_End'        => 'EX4',
    'exin'          => 'EX4',
    'exonIndex'     => '4',
    'func'          => 'cds-del',
    'funcSO'        => 'SO:0001822',
    'funcSOname'    => 'inframe_deletion',
    'geneId'        => '178',
    'geneSym'       => 'AGL',
    'genepart'      => 'CDS',
    'genepartIndex' => '4',
    'genepartSO'    => 'SO:0000316',
    'intronIndex'   => '.',
    'p'             => 'p.G116_L122del',
    'postEnd'       => {
        'cDot' => '366',
        'exin' => 'EX4',
        'nDot' => 766,
        'r'    => 'C3'
    },
    'prAlt'    => 'V',
    'prRef'    => 'VGADNHVL',
    'preStart' => {
        'cDot' => '344',
        'exin' => 'EX4',
        'nDot' => 744,
        'r'    => 'C3'
    },
    'primaryTag' => 'Y',
    'prot'       => 'NP_000633.2',
    'protBegin'  => 115,
    'protEnd'    => 122,
    'r'          => 'C3',
    'r_Begin'    => 'C3',
    'r_End'      => 'C3',
    'rnaBegin'   => 745,
    'rnaEnd'     => 765,
    'strd'       => '+',
    'trAlt'      => '',
    'trRef'      => 'TGGTGCTGATAATCATGTGCT',
    'trRefComp'  => {
        'EX4' => 21
    }
};

my $cds_ins_anno = bless(
    {
        'trInfo' => {
            'NM_006158.3' => {
                'c'             => 'c.1411_1412insCCC',
                'cdsBegin'      => '1412',
                'cdsEnd'        => '1411',
                'ei_Begin'      => 'EX3',
                'ei_End'        => 'EX3',
                'exin'          => 'EX3',
                'exonIndex'     => '3',
                'func'          => 'cds-ins',
                'funcSO'        => 'SO:0001821',
                'funcSOname'    => 'inframe_insertion',
                'geneId'        => '4747',
                'geneSym'       => 'NEFL',
                'genepart'      => 'CDS',
                'genepartIndex' => '3',
                'genepartSO'    => 'SO:0000316',
                'intronIndex'   => '.',
                'p'             => 'p.P471dup',
                'postEnd'       => {
                    'cDot' => '1412',
                    'exin' => 'EX3',
                    'nDot' => 1514,
                    'r'    => 'C3'
                },
                'preStart' => {
                    'cDot' => '1411',
                    'exin' => 'EX3',
                    'nDot' => 1513,
                    'r'    => 'C3'
                },
		'primaryTag'=> 'N',
		'prAlt'	    => 'PP',
		'prRef'	    => 'P',
                'prot'      => 'NP_006149.2',
                'protBegin' => 471,
                'protEnd'   => 471,
                'r'         => 'C3',
                'r_Begin'   => 'C3',
                'r_End'     => 'C3',
                'rnaBegin'  => 1514,
                'rnaEnd'    => 1513,
                'strd'      => '-',
                'trAlt'     => 'CCC',
                'trRef'     => '',
                'trRefComp' => {
                    'EX3' => 0
                }
            }
        },
        'var' => bless(
            {
                'alt'       => 'GGG',
                'altlen'    => 3,
                'chr'       => '8',
                'end'       => 24811067,
                'gHGVS'     => 'g.24811067_24811068insGGG',
                'guess'     => 'ins',
                'imp'       => 'ins',
                'pos'       => 24811067,
                'ref'       => '',
                'refbuild'  => 'GRCh37',
                'reflen'    => 0,
                'sm'        => 0,
                'varTypeSO' => 'SO:0000667'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);
my $cds_rep_anno = bless(
    {
        'trInfo' => {
            'NM_000581.2' => {
                'c'             => 'c.21GGC[6>4]',
                'cdsBegin'      => '21',
                'cdsEnd'        => '38',
                'ei_Begin'      => 'EX1',
                'ei_End'        => 'EX1',
                'exin'          => 'EX1',
                'exonIndex'     => '1',
                'func'          => 'cds-del',
                'funcSO'        => 'SO:0001822',
                'funcSOname'    => 'inframe_deletion',
                'geneId'        => '2876',
                'geneSym'       => 'GPX1',
                'genepart'      => 'CDS',
                'genepartIndex' => '1',
                'genepartSO'    => 'SO:0000316',
                'intronIndex'   => '.',
                'p'             => 'p.A7[7>5]',
                'postEnd'       => {
                    'cDot' => '39',
                    'exin' => 'EX1',
                    'nDot' => 119,
                    'r'    => 'C1'
                },
                'preStart' => {
                    'cDot' => '20',
                    'exin' => 'EX1',
                    'nDot' => 100,
                    'r'    => 'C1'
                },
		'primaryTag'=> 'Y',
		'prAlt'	    => 'AAAAA',
		'prRef'	    => 'AAAAAAA',
                'prot'      => 'NP_000572.2',
                'protBegin' => 7,
                'protEnd'   => 13,
                'r'         => 'C1',
                'r_Begin'   => 'C1',
                'r_End'     => 'C1',
                'rnaBegin'  => 101,
                'rnaEnd'    => 118,
                'strd'      => '-',
                'trAlt'     => 'GGCGGCGGCGGC',
                'trRef'     => 'GGCGGCGGCGGCGGCGGC',
                'trRefComp' => {
                    'EX1' => 18
                }
            },
            'NM_201397.1' => {
                'c'             => 'c.21GGC[6>4]',
                'cdsBegin'      => '21',
                'cdsEnd'        => '38',
                'ei_Begin'      => 'EX1E',
                'ei_End'        => 'EX1E',
                'exin'          => 'EX1E',
                'exonIndex'     => '1',
                'func'          => 'cds-del',
                'funcSO'        => 'SO:0001822',
                'funcSOname'    => 'inframe_deletion',
                'geneId'        => '2876',
                'geneSym'       => 'GPX1',
                'genepart'      => 'CDS',
                'genepartIndex' => '1',
                'genepartSO'    => 'SO:0000316',
                'intronIndex'   => '.',
                'p'             => 'p.A7[7>5]',
                'postEnd'       => {
                    'cDot' => '39',
                    'exin' => 'EX1E',
                    'nDot' => 119,
                    'r'    => 'C1E'
                },
                'preStart' => {
                    'cDot' => '20',
                    'exin' => 'EX1E',
                    'nDot' => 100,
                    'r'    => 'C1E'
                },
		'primaryTag'=> 'N',
		'prAlt'	    => 'AAAAA',
		'prRef'	    => 'AAAAAAA',
                'prot'      => 'NP_958799.1',
                'protBegin' => 7,
                'protEnd'   => 13,
                'r'         => 'C1E',
                'r_Begin'   => 'C1E',
                'r_End'     => 'C1E',
                'rnaBegin'  => 101,
                'rnaEnd'    => 118,
                'strd'      => '-',
                'trAlt'     => 'GGCGGCGGCGGC',
                'trRef'     => 'GGCGGCGGCGGCGGCGGC',
                'trRefComp' => {
                    'EX1E' => 18
                }
            }
        },
        'var' => bless(
            {
                '+' => {
                    'ba'  => '',
                    'bal' => 0,
                    'bp'  => 49395685,
                    'br'  => 'GCCGCC',
                    'brl' => 6
                },
                '-' => {
                    'ba'  => '',
                    'bal' => 0,
                    'bp'  => 49395673,
                    'br'  => 'GCCGCC',
                    'brl' => 6
                },
                'a'         => 'GCCGCCGCCGCC',
                'al'        => 12,
                'alt'       => 'GCCGCCGCCGCC',
                'alt_cn'    => 4,
                'altlen'    => 12,
                'chr'       => '3',
                'end'       => 49395691,
                'gHGVS'     => 'g.49395674GCC[6>4]',
                'guess'     => 'delins',
                'imp'       => 'rep',
                'p'         => 49395673,
                'pos'       => 49395673,
                'r'         => 'GCCGCCGCCGCCGCCGCC',
                'ref'       => 'GCCGCCGCCGCCGCCGCC',
                'ref_cn'    => 6,
                'refbuild'  => 'GRCh37',
                'reflen'    => 18,
                'rep'       => 'GCC',
                'replen'    => 3,
                'rl'        => 18,
                'sm'        => 2,
                'varTypeSO' => 'SO:1000032'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);
my $cds_delins_anno = bless(
    {
        'trInfo' => {
            'NM_152486.2' => {
                'c'             => 'c.73-2_73-1delinsTC',
                'cdsBegin'      => '73-2',
                'cdsEnd'        => '73-1',
                'ei_Begin'      => 'IVS2',
                'ei_End'        => 'IVS2',
                'exin'          => 'IVS2',
                'exonIndex'     => '.',
                'func'          => 'splice-3',
                'funcSO'        => '',
                'funcSOname'    => 'unknown-likely-deleterious',
                'geneId'        => '148398',
                'geneSym'       => 'SAMD11',
                'genepart'      => 'three_prime_cis_splice_site',
                'genepartIndex' => '2',
                'genepartSO'    => 'SO:0000164',
                'intronIndex'   => '2',
                'postEnd'       => {
                    'cDot' => '73',
                    'exin' => 'EX3',
                    'nDot' => 153,
                    'r'    => 'C2'
                },
                'preStart' => {
                    'cDot' => '73-3',
                    'exin' => 'IVS2',
                    'nDot' => '153-3',
                    'r'    => 'IC1'
                },
		'primaryTag'=> 'Y',
                'prot'      => 'NP_689699.2',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'AC1',
                'r_Begin'   => 'AC1',
                'r_End'     => 'AC1',
                'rnaBegin'  => '153-2',
                'rnaEnd'    => '153-1',
                'strd'      => '+',
                'trAlt'     => 'TC',
                'trRef'     => 'AG',
                'trRefComp' => {
                    'IVS2' => [ 0, 2 ]
                }
            }
        },
        'var' => bless(
            {
                'alt'       => 'TC',
                'altlen'    => 2,
                'chr'       => '1',
                'end'       => 865534,
                'gHGVS'     => 'g.865533_865534delAGinsTC',
                'guess'     => 'delins',
                'imp'       => 'delins',
                'pos'       => 865532,
                'ref'       => 'AG',
                'refbuild'  => 'GRCh37',
                'reflen'    => 2,
                'sep_snvs'  => [ 865533, 865534 ],
                'sm'        => 3,
                'varTypeSO' => 'SO:1000032'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $cds_no_call_anno = bless(
    {
        'trInfo' => {
            'NM_152486.2' => {
                'c'             => 'c.-2T>N',
                'cdsBegin'      => '-2',
                'cdsEnd'        => '-2',
                'ei_Begin'      => 'EX2',
                'ei_End'        => 'EX2',
                'exin'          => 'EX2',
                'exonIndex'     => '2',
                'func'          => 'utr-5',
                'funcSO'        => '',
                'funcSOname'    => 'unknown',
                'geneId'        => '148398',
                'geneSym'       => 'SAMD11',
                'genepart'      => 'five_prime_UTR',
                'genepartIndex' => '2',
                'genepartSO'    => 'SO:0000204',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '-1',
                    'exin' => 'EX2',
                    'nDot' => 80,
                    'r'    => '5U1'
                },
                'preStart' => {
                    'cDot' => '-3',
                    'exin' => 'EX2',
                    'nDot' => 78,
                    'r'    => '5U1'
                },
		'primaryTag'=> 'Y',
                'prot'      => 'NP_689699.2',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => '5U1',
                'r_Begin'   => '5U1',
                'r_End'     => '5U1',
                'rnaBegin'  => 79,
                'rnaEnd'    => 79,
                'strd'      => '+',
                'trAlt'     => 'N',
                'trRef'     => 'T',
                'trRefComp' => {
                    'EX2' => 1
                }
            },
            'NR_026874.1' => {
                'c'             => 'n.-6503A>N',
                'cdsBegin'      => '',
                'cdsEnd'        => '',
                'ei_Begin'      => '.',
                'ei_End'        => '.',
                'exin'          => '.',
                'exonIndex'     => '.',
                'func'          => 'promoter',
                'funcSO'        => '',
                'funcSOname'    => 'unknown',
                'geneId'        => '100130417',
                'geneSym'       => 'LOC100130417',
                'genepart'      => 'promoter',
                'genepartIndex' => 0,
                'genepartSO'    => 'SO:0000167',
                'intronIndex'   => '.',
                'postEnd'       => {
                    'cDot' => '',
                    'exin' => '.',
                    'nDot' => -6502,
                    'r'    => 'PROM'
                },
                'preStart' => {
                    'cDot' => '',
                    'exin' => '.',
                    'nDot' => -6504,
                    'r'    => 'PROM'
                },
		'primaryTag'=> 'Y',
                'prot'      => '',
                'protBegin' => '',
                'protEnd'   => '',
                'r'         => 'PROM',
                'r_Begin'   => 'PROM',
                'r_End'     => 'PROM',
                'rnaBegin'  => -6503,
                'rnaEnd'    => -6503,
                'strd'      => '-',
                'trAlt'     => 'N',
                'trRef'     => 'A',
                'trRefComp' => {
                    'P0' => [ 0, 1 ]
                }
            }
        },
        'var' => bless(
            {
                'alt'       => 'N',
                'altlen'    => 1,
                'chr'       => '1',
                'end'       => 861320,
                'gHGVS'     => 'g.861319T>N',
                'guess'     => 'no-call',
                'imp'       => 'snv',
                'pos'       => 861319,
                'ref'       => 'T',
                'refbuild'  => 'GRCh37',
                'reflen'    => 1,
                'sm'        => 1,
                'varTypeSO' => 'no-call'
            },
            'BedAnno::Var'
        )
    },
    'BedAnno::Anno'
);

my $no_call_edge_ins_parse = bless(
    {
        'reflen'    => 0,
        'imp'       => 'ins',
        'pos'       => 6526167,
        'guess'     => 'no-call',
        'alt'       => '?',
        'ref'       => '=',
        'chr'       => '1',
        'sm'        => 0,
        'end'       => 6526167,
    },
    'BedAnno::Var'
);

my $no_call_edge_ins_trinfo = {
    'r'         => '5U1E',
    'protBegin' => '',
    'cdsBegin'  => '1',
    'preStart'  => {
        'r'    => '5U1E',
        'cDot' => '-1',
        'nDot' => 88,
        'exin' => 'EX1'
    },
    'genepartIndex' => '1',
    'ei_End'        => 'EX1',
    'trRef'         => '=',
    'protEnd'       => '',
    'cdsEnd'        => '-1',
    'exonIndex'     => '1',
    'trAlt'         => '?',
    'r_Begin'       => 'C1',
    'intronIndex'   => '.',
    'exin'          => 'EX1',
    'prot'          => 'NP_683866.1',
    'strd'          => '-',
    'trRefComp'     => {
        'EX1' => 0
    },
    'geneSym' => 'TNFRSF25',
    'r_End'   => '5U1E',
    'postEnd' => {
        'r'    => 'C1',
        'cDot' => '1',
        'nDot' => 89,
        'exin' => 'EX1'
    },
    'c'          => 'c.-1_1ins?',
    'geneId'     => '8718',
    'rnaEnd'     => 88,
    'genepart'   => 'five_prime_UTR',
    'func'       => 'unknown-no-call',
    'rnaBegin'   => 89,
    'funcSO'     => '',
    'genepartSO' => 'SO:0000204',
    'ei_Begin'   => 'EX1',
    'funcSOname' => 'unknown-no-call',
    'primaryTag' => 'Y'
};

test_parse_var( "crawler_snv_parse",           $snv_parse,    $crawler_input );
test_parse_var( "crawler_vcf_insert_parse",    $insert_parse, $crawler_input2 );
test_parse_var( "crawler_var_undef_del_parse", $del_parse,    $crawler_input3 );
test_parse_var( "crawler_no_ref_ins_parse",    $insert_parse, $crawler_input4 );
test_parse_var( "snv_parse", $snv_parse, "chr1", 14410, "C", "A" );
test_parse_var( "insert_parse", $insert_parse, "chr1", 14410, "C",
    "CGAATAGCTA" );
test_parse_var( "del_parse", $del_parse, "chr1", 14410, "CTAGATCG", "C" );
test_parse_var( "rep_parse", $rep_parse, "chr1", 14410, "CTAGA",
    "CTAGTAGTAGTAGA" );
test_parse_var( "delins_parse", $delins_parse, "chr1", 14410, "CTAGA", "GT" );
test_parse_var( "subs_parse", $subs_parse, "chr1", 14410, "CATGA", "TGTGT" );
test_parse_var( "no_call_parse", $no_call_parse, "chr1", 14410, 14410, "=",
    '?' );
test_parse_var( "no_call_edge_parse", $no_call_edge_parse, "chr1", 0, 10000,
    "=", "?" );
test_parse_var( "no_call_edge_ins_parse", $no_call_edge_ins_parse, "chr1",
    6526167, 6526167, "=", "?" );


test_varanno( "snv_varanno",     $snv_varanno,     $snv_parse );
test_varanno( "ins_varanno",     $ins_varanno,     $insert_parse );
test_varanno( "del_varanno",     $del_varanno,     $del_parse );
test_varanno( "rep_varanno",     $rep_varanno,     $rep_parse );
test_varanno( "delins_varanno",  $delins_varanno,  $delins_parse );
test_varanno( "subs_varanno",    $subs_varanno,    $subs_parse );
test_varanno( "no_call_varanno", $no_call_varanno, $no_call_parse );
test_varanno( "cds_rna_delins_anno", $cds_rna_delins_anno, $cds_rna_delins );
test_varanno( "cds_no_change_anno",  $cds_no_change_anno,  $cds_no_change );
test_varanno( "cds_snv_anno",        $cds_snv_anno,        $cds_snv );
test_varanno( "cds_del_anno",        $cds_del_anno,        $cds_del );
test_varanno( "cds_ins_anno",        $cds_ins_anno,        $cds_ins );
test_varanno( "cds_rep_anno",        $cds_rep_anno,        $cds_rep );
test_varanno( "cds_delins_anno",     $cds_delins_anno,     $cds_delins );
test_varanno( "cds_no_call_anno",    $cds_no_call_anno,    $cds_no_call );

test_varanno( "no_call_edge_varanno", $no_call_edge_varanno,
    $no_call_edge_parse );
test_ok ( "no_call_edge_ins_trinfo", $no_call_edge_ins_trinfo, 'NM_148965.1',
   "chr1", 6526167, 6526167, "=", "?" ); 
test_ok ( "cds_rna_snv1_anno", $cds_rna_snv1_anno, 'NM_000581.2',
   "chr3", 49395708, 49395709, "C", "T" );
test_ok ( "cds_rna_snv2_anno", $cds_rna_snv2_anno, 'NM_000581.2',
   "chr3", 49395704, 49395705, "C", "A" );
test_ok ( "large_del_anno", $large_del_anno, 'NM_000642.2',
   "chr1", 100327863, 100327884, "TGGTGCTGATAATCATGTGCT", "" );

done_testing();
exit 0;

sub test_ok {
    my ( $tag, $expect, $tid, @args ) = @_;
    my $ranno = $beda->anno(@args);
    if (
        !is_deeply( $ranno->{trInfo}->{$tid}, $expect, "for [ $tag ]")
      )
    {
        explain "The anno infor are: ", $ranno;
    }
}

sub test_parse_var {
    my ( $tag, $expect, @args ) = @_;
    my $rvar = BedAnno::Var->new(@args);
    if ( !is_deeply( $rvar, $expect, "for [ $tag ]" ) ) {
        explain "The var infomations are: ", $rvar;
    }
}

sub test_varanno {
    my ( $tag, $expect, $vara )   = @_;
    my ( $vAnno,  $noneed ) = $beda->varanno($vara);
    if ( !is_deeply( $vAnno, $expect, "for [ $tag ]" ) ) {
        explain "The anno infomations are: ", $vAnno;
    }
}

# test if the tabix is available
sub test_tabix {
    my ($file) = @_;
    my $cmd;
    $cmd = "tabix -f -p bed $file";
    system($cmd);
    if ( !is( $?, 0, "tabix .. $cmd" ) ) {
        exit 1;
    }
}
