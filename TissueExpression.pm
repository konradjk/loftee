=head1 CONTACT                                                                                                       

 Konrad Karczewski <konradjkarczewski@gmail.com>
 
=cut

=head1 NAME

 TissueExpression

=head1 SYNOPSIS

 mv TissueExpression.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin TissueExpression

=head1 DESCRIPTION

 A VEP plugin that overlays GTEx data on transcripts.
 Requires DBD::SQLite (>=1.4.2).

=cut

package TissueExpression;

use strict;
use warnings;

our $debug = 0;
our $ddebug = 0;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use DBI;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub get_header_info {
    return {
        TissueExpression => "GTEx data"
    };
}

sub feature_types {
    return ['Transcript'];
}

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);
    
    #$self->{has_cache} = 1;
    $self->{expressed_cutoff} = $self->params->[0] || 1;
    $self->{database} = DBI->connect("dbi:SQLite:dbname=gtex.db", "", "") or die "Cannot find gtex.db\n";
    
    return $self;
}

sub run {
    my ($self, $transcript_variation_allele) = @_;
    
    my $transcript = $transcript_variation_allele->transcript_variation->transcript;
    
    my $transcript_tissue;
    if (exists($transcript->{expression_cache})) {
        #print "Got expression cache!\n" if $ddebug;
        $transcript_tissue = $transcript->{expression_cache};
    } else {
        my $sth = $self->{database}->prepare("SELECT * FROM tissues WHERE transcript = ?");
        $sth->execute($transcript->stable_id());
        
        my @tissue_entries = ();
        while (my $entry = $sth->fetchrow_hashref) {
            $entry->{tissue} =~ s/ /_/g;
            push(@tissue_entries, $entry->{tissue} . ":" . $entry->{expression});
        }
        $transcript_tissue = join(",", @tissue_entries);
        
        print "Tissues: " . $transcript_tissue . "\n" if $ddebug;
        $transcript->{expression_cache} = $transcript_tissue;
    }

    return { TissueExpression => $transcript_tissue };
}

1;