=head1 CONTACT                                                                                                       

 Konrad Karczewski <konradjkarczewski@gmail.com>
 
=cut

=head1 NAME

 LoF

=head1 SYNOPSIS

 mv LoF.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin LoF

=head1 DESCRIPTION

 A VEP plugin that determines loss-of-function variation.

=cut

package LoF;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub get_header_info {
    return {
        LoF => "Loss-of-function annotation (HC = High Confidence)"
    };
}

sub feature_types {
    return ['Transcript'];
}

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);
    
    $self->{filter_position} = $self->params->[0] || 0.05;
    $self->{min_intron_size} = $self->params->[1] || 15;
    
    return $self;
}

our $debug = 0;
our $ddebug = 0;

sub run {
    my ($self, $tva) = @_;
    
    my $transcript = $tva->transcript;
    my $transcript_variation = $tva->transcript_variation;
    my $variation_feature = $tva->variation_feature;
    my $base_variation_feature = $tva->base_variation_feature;
    
    #unless ($variation_feature->variation_name eq 'rs139717535') { return {} }
    
    my @consequences = map { $_->SO_term } @{ $tva->get_all_OverlapConsequences };
    
    my $output;
    
    # Filter in
    unless ($transcript->biotype eq "protein_coding") {
        return {}
    }
    
    unless ("stop_gained" ~~ @consequences || "splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences || "frameshift_variant" ~~ @consequences) {
        return {}
    }
    
    # Filter out
    if ("stop_gained" ~~ @consequences || "frameshift_variant" ~~ @consequences){
        if (check_for_exon_annotation_errors($transcript_variation)) {
            print "Exon annotation Error\n" if $debug;
            return { LoF => "LOF=LC", }
        }
        if (check_for_single_exon($transcript_variation)) {
            print "Single Exon\n" if $debug;
            return { LoF => "LOF=LC", }
        }
        if (check_surrounding_introns($transcript_variation, $self->{min_intron_size})) {
            print "Small surrounding introns\n" if $debug;
            return { LoF => "LOF=LC", }
        }
        if (check_position($transcript_variation, $self->{filter_position})) {
            print "Last X\n" if $debug;
            return { LoF => "LOF=LC", }
        }
    }
    
    if ("splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences) {
        # Currently in this order intentionally (each one before is checking undefined)
        if (check_for_intron_annotation_errors($transcript_variation)) {
            print "Intron annotation error\n" if $debug;
            return { LoF => "LOF=LC", }
        }
        if (check_intron_size($transcript_variation, $self->{min_intron_size})) {
            print "Small intron\n" if $debug;
            return { LoF => "LOF=LC", }
        }
        if (check_for_non_canonical_intron_motif($transcript_variation)) {
            print "Not a canonical motif\n" if $debug;
            return { LoF => "LOF=LC", }
        }
    }
    
    #if (check_for_ancestral_allele($variation_feature)) {
    #    print "Ancestral allele\n" if $debug;
    #    return { LoF => "LOF=LC", }
    #}
    
    return { LoF => "LOF=HC", }
}

# Global functions
sub small_intron {
    my $transcript_variation = shift;
    my $intron_number = shift;
    my $min_intron_size = shift;
    my @gene_introns = @{$transcript_variation->transcript->get_all_Introns()};
    return (length($gene_introns[$intron_number]->seq) < $min_intron_size)
}

sub intron_motif_start {
    my $transcript_variation = shift;
    my $intron_number = shift;
    my @gene_introns = @{$transcript_variation->transcript->get_all_Introns()};
    return (substr($gene_introns[$intron_number]->seq, 0, 2) ne 'GT')
}

sub intron_motif_end {
    my $transcript_variation = shift;
    my $intron_number = shift;
    my @gene_introns = @{$transcript_variation->transcript->get_all_Introns()};
    return (substr($gene_introns[$intron_number]->seq, length($gene_introns[$intron_number]->seq) - 2, 2) ne 'AG')
}

# Stop-gain and frameshift annotations
sub check_for_exon_annotation_errors {
    my $transcript_variation = shift;
    return (!defined($transcript_variation->exon_number))
}

sub check_position {
    my ($transcript_variation, $cutoff) = @_;
    # Need to resolve this: CDS should be multiple of 3 and isn't...
    my $transcript_cds_length = $transcript_variation->transcript->cdna_coding_end - $transcript_variation->transcript->cdna_coding_start + 2;
    my $variant_cds_position = $transcript_variation->cdna_start;
    print $transcript_variation->transcript->stable_id . ": " . $variant_cds_position . " out of " . $transcript_cds_length . "\n" if $ddebug;
    print "CDS not multiple of 3\n" if ($debug && $transcript_cds_length % 3);
    return ($variant_cds_position/$transcript_cds_length >= 1-$cutoff);
}

sub check_for_single_exon {
    my $transcript_variation = shift;
    my @exons = split /\//, ($transcript_variation->exon_number);
    return ($exons[1] == 1)
}

sub check_surrounding_introns {
    my $transcript_variation = shift;
    my $min_intron_size = shift;
    my ($exon_number, $total_exons) = split /\//, ($transcript_variation->exon_number);
    $exon_number--;
    
    if ($exon_number == 0) {
        return (small_intron($transcript_variation, $exon_number, $min_intron_size) ||
               intron_motif_start($transcript_variation, $exon_number))
    } elsif ($exon_number == $total_exons - 1) {
        return (small_intron($transcript_variation, $exon_number - 1, $min_intron_size) ||
                intron_motif_end($transcript_variation, $exon_number - 1))
    } else {
        return (small_intron($transcript_variation, $exon_number, $min_intron_size) ||
                small_intron($transcript_variation, $exon_number - 1, $min_intron_size) ||
                intron_motif_start($transcript_variation, $exon_number) ||
                intron_motif_end($transcript_variation, $exon_number - 1))
    }
}

# Splicing annotations
sub check_for_intron_annotation_errors {
    # Currently only checking if intron number is defined (is this necessary?)
    my $transcript_variation = shift;
    return (!defined($transcript_variation->intron_number))
}

sub check_intron_size {
    my $transcript_variation = shift;
    my $min_intron_size = shift;
    
    my ($intron_number, $total_introns) = split /\//, ($transcript_variation->intron_number);
    $intron_number--;
    
    return small_intron($transcript_variation, $intron_number, $min_intron_size)
}

sub check_for_non_canonical_intron_motif {
    my $transcript_variation = shift;
    my ($intron_number, $total_introns) = split /\//, ($transcript_variation->intron_number);
    $intron_number--;
    return (intron_motif_start($transcript_variation, $intron_number) || intron_motif_end($transcript_variation, $intron_number))
}

#sub check_for_ancestral_allele {
#    my $variation_feature = shift;
#
#    my $reg = "Bio::EnsEMBL::Registry";
#    
#    my $mlss_adaptor = $reg->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet") or die "Failed to connect to compara database\n";
#    my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("EPO", "primates");
#    my $aln_adaptor = $reg->get_adaptor('Multi', 'compara', 'GenomicAlignBlock') or die "Failed to fetch conservation adaptor\n";
#    
#    my $species = "Homo sapiens";
#    my $chr = '1';
#    my $start = 879431;
#    my $end = 879431;
#    #my $chr = $variation_feature->seq_region_name();
#    #my $start = $variation_feature->seq_region_start();
#    #my $end = $variation_feature->seq_region_end();
#    
#    my $query_slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, "core", "Slice");
#    my $query_slice = $query_slice_adaptor->fetch_by_region("chromosome", $chr, $start, $end);
#    
#    my $genomic_align_blocks = $aln_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $query_slice);
#    
#    my $all_aligns;
#    
#    foreach my $this_genomic_align_block (@$genomic_align_blocks) {
#        my $simple_align = $this_genomic_align_block->get_SimpleAlign;
#        print $simple_align;
#    }
#    
#    return 0;
#}

1;