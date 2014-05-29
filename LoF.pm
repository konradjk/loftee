=head1 CONTACT                                                                                                       

 Konrad Karczewski <konradjkarczewski@gmail.com>
 
=cut

=head1 NAME

 LoF

=head1 SYNOPSIS

 mv LoF.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin LoF
 perl variant_effect_predictor.pl -i variations.vcf --plugin LoF,filter_position:0.05,...

=head1 DESCRIPTION

 A VEP plugin that filters loss-of-function variation.

=cut

package LoF;

use strict;
use warnings;

our $debug = 0;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub get_header_info {
    return {
        LoF => "Loss-of-function annotation (HC = High Confidence; LC = Low Confidence)",
        LoF_filter => "Reason for LoF not being HC",
        LoF_flags => "Possible warning flags for LoF"
    };
}

sub feature_types {
    return ['Transcript'];
}

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);
    
    foreach my $parameter (@{$self->params}) {
        my @param = split /:/, $parameter;
        if (scalar @param == 2) {
            $self->{$param[0]} = $param[1];
        }
    }
    
    $self->{filter_position} = $self->{filter_position} || 0.05;
    $self->{min_intron_size} = $self->{min_intron_size} || 15;
    $self->{fast_length_calculation} = $self->{fast_length_calculation} || 'fast';
    $self->{human_ancestor_fa} = $self->{human_ancestor_fa} || 'human_ancestor.fa.rz';
    $self->{check_complete_cds} = $self->{check_complete_cds} || 'false';
    
    if ($debug) {
        print "Read LOFTEE parameters\n";
        while (my ($key, $value) = each(%$self)) {
            print $key . " : " . $value . "\n";
        }
    }
    
    return $self;
}

sub run {
    my ($self, $transcript_variation_allele) = @_;
    
    my $transcript_variation = $transcript_variation_allele->transcript_variation;
    my $variation_feature = $transcript_variation_allele->variation_feature;
    
    my @consequences = map { $_->SO_term } @{ $transcript_variation_allele->get_all_OverlapConsequences };
    
    my @filters = ();
    my @flags = ();
    
    # Filter in
    unless ($transcript_variation->transcript->biotype eq "protein_coding") {
        return {};
    }
    unless ("stop_gained" ~~ @consequences || "splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences || "frameshift_variant" ~~ @consequences) {
        return {};
    }
    
    my $confidence = 'HC';
    
    # Filter out
    if ("stop_gained" ~~ @consequences || "frameshift_variant" ~~ @consequences){
        push(@filters, 'END_TRUNC') if (check_position($transcript_variation, $self->{filter_position}, $self->{fast_length_calculation}));
        if (check_for_exon_annotation_errors($transcript_variation)) {
            push(@filters, 'EXON_INTRON_UNDEF');
        } elsif (check_for_single_exon($transcript_variation)) {
            push(@flags, 'SINGLE_EXON');
        } else {
            if (lc($self->{check_complete_cds}) eq 'true') {
                push(@filters, 'INCOMPLETE_CDS') if (check_incomplete_cds($transcript_variation));
            }
            push(@filters, 'NON_CAN_SPLICE_SURR') if (check_surrounding_introns($transcript_variation, $self->{min_intron_size}));
        }
    }

    if ("splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences) {
        if (check_for_intron_annotation_errors($transcript_variation)) {
            push(@filters, 'EXON_INTRON_UNDEF');
        } else {
            push(@filters, 'SMALL_INTRON') if (check_intron_size($transcript_variation, $self->{min_intron_size}));
            push(@filters, 'NON_CAN_SPLICE') if (check_for_non_canonical_intron_motif($transcript_variation));
            if ("splice_acceptor_variant" ~~ @consequences) {
                push(@flags, 'NAGNAG_SITE') if (check_nagnag_variant($variation_feature));
            }
        }
    }
    
    if (lc($self->{human_ancestor_fa}) ne 'false') {
        push(@filters, 'ANC_ALLELE') if (check_for_ancestral_allele($transcript_variation_allele, $self->{human_ancestor_fa}));
    }
    
    if ($confidence eq 'HC' && scalar @filters > 0) {
        $confidence = 'LC';
    }
    
    return { LoF => $confidence, LoF_filter => join(',', @filters), LoF_flags => join(',', @flags) };
}

# Global functions
sub small_intron {
    my $transcript_variation = shift;
    my $intron_number = shift;
    my $min_intron_size = shift;
    my @gene_introns = @{$transcript_variation->transcript->get_all_Introns()};
    
    return ($gene_introns[$intron_number]->length < $min_intron_size);
}

sub intron_motif_start {
    my ($transcript_variation, $intron_number) = @_;
    
    my $transcript = $transcript_variation->transcript;
    my @gene_introns = @{$transcript->get_all_Introns()};
    
    # Cache intron sequence
    unless (exists($transcript->{intron_cache}->{$intron_number})) {
        $transcript->{intron_cache}->{$intron_number} = $gene_introns[$intron_number]->seq;
    }
    my $sequence = $transcript->{intron_cache}->{$intron_number};
    
    return (substr($sequence, 0, 2) ne 'GT');
}

sub intron_motif_end {
    my ($transcript_variation, $intron_number) = @_;
    
    my $transcript = $transcript_variation->transcript;
    my @gene_introns = @{$transcript->get_all_Introns()};
    
    # Cache intron sequence
    unless (exists($transcript->{intron_cache}->{$intron_number})) {
        $transcript->{intron_cache}->{$intron_number} = $gene_introns[$intron_number]->seq;
    }
    my $sequence = $transcript->{intron_cache}->{$intron_number};
    
    return (substr($sequence, length($sequence) - 2, 2) ne 'AG')
}

sub get_cds_length_fast {
    my $transcript = shift;
    
    my $transcript_cds_length = $transcript->cdna_coding_end - $transcript->cdna_coding_start + 1;
    return $transcript_cds_length;
}

sub get_cds_length {
    my $transcript = shift;
    
    # Cache CDS sequence
    unless (exists($transcript->{cds_seq_cache})) {
        $transcript->{cds_seq_cache} = $transcript->translateable_seq;
    }
    my $transcript_cds_length = length($transcript->{cds_seq_cache});
    return $transcript_cds_length;
}

# Stop-gain and frameshift annotations
sub check_incomplete_cds {
    my $transcript_variation = shift;
    
    my $transcript = $transcript_variation->transcript;
    my $start_annotation = $transcript->get_all_Attributes('cds_start_NF');
    my $end_annotation = $transcript->get_all_Attributes('cds_end_NF');
    
    return (defined($start_annotation) || defined($end_annotation));
}

sub check_for_exon_annotation_errors {
    my $transcript_variation = shift;
    return (!defined($transcript_variation->exon_number))
}

sub check_position {
    my ($transcript_variation, $cutoff, $speed) = @_;
    
    # 2 ways to get length: fast and approximate, or slow and accurate
    my $transcript_cds_length;
    if ($speed eq 'fast') {
        $transcript_cds_length = get_cds_length_fast($transcript_variation->transcript);
    } else {
        $transcript_cds_length = get_cds_length($transcript_variation->transcript);
    }
    my $variant_cds_position = $transcript_variation->cdna_end;
    
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
    
    # Check for small introns and GT..AG motif
    # Only next intron if in first exon, only previous intron if in last exon, otherwise both previous and next
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
sub check_nagnag_variant {
    my $variation_feature = shift;
    
    # Cache splice sites
    unless (exists($variation_feature->{splice_context_seq_cache})) {
        $variation_feature->{splice_context_seq_cache} = uc($variation_feature->feature_Slice->expand(4, 4)->seq);
    }
    my $sequence = $variation_feature->{splice_context_seq_cache};
    
    # Only consider NAGNAG sites for SNPs for now
    if (length($sequence) == 9) {
        return ($sequence =~ m/AG.AG/);
    } else {
        return 0;
    }
}
sub check_for_intron_annotation_errors {
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

sub check_for_ancestral_allele {
    my $transcript_variation_allele = shift;
    my $human_ancestor_location = shift;
    my $variation_feature = $transcript_variation_allele->variation_feature;
    my $aff_allele = $transcript_variation_allele->variation_feature_seq;
    
    # Get ancestral allele from human_ancestor.fa.rz
    my $region = $variation_feature->seq_region_name() . ":" . $variation_feature->seq_region_start() . '-' . $variation_feature->seq_region_end();
    my $faidx = `samtools faidx $human_ancestor_location $region`;
    my @lines = split(/\n/, $faidx);
    shift @lines;
    my $ancestral_allele = uc(join('', @lines));
    
    return ($ancestral_allele eq $aff_allele)
}

1;