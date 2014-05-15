=head1 CONTACT                                                                                                       

 Konrad Karczewski <konradjkarczewski@gmail.com>
 
=cut

=head1 NAME

 LoF

=head1 SYNOPSIS

 mv LoF.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin LoF
 perl variant_effect_predictor.pl -i variations.vcf --plugin LoF,filter_position:0.05,min_intron_size:15

=head1 DESCRIPTION

 A VEP plugin that filters loss-of-function variation.

=cut

package LoF;

use strict;
use warnings;

our $debug = 1;
our $ddebug = 0;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub get_header_info {
    return {
        LoF => "Loss-of-function annotation (HC = High Confidence)",
        LoF_filter => "Reason for LoF not being HC"
    };
}

sub feature_types {
    return ['Transcript'];
}

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);
    #$self->{has_cache} = 1;
    
    foreach my $parameter (@{$self->params}) {
        my @param = split /:/, $parameter;
        if (scalar @param == 2) {
            $self->{$param[0]} = $param[1];
        }
    }
    
    $self->{filter_position} = $self->{filter_position} || 0.05;
    $self->{min_intron_size} = $self->{min_intron_size} || 15;
    $self->{fast_length_calculation} = $self->{fast_length_calculation} || 'fast';
    
    if ($debug) {
        print "Read parameters\n";
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
    
    # Filter in
    unless ($transcript_variation->transcript->biotype eq "protein_coding") {
        return {};
    }
    unless ("stop_gained" ~~ @consequences || "splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences || "frameshift_variant" ~~ @consequences) {
        return {};
    }
    
    my $confidence = 'HC';
    my @filters = ();
    
    # Filter out
    if ("stop_gained" ~~ @consequences || "frameshift_variant" ~~ @consequences){
        push(@filters, 'END_TRUNC') if (check_position($transcript_variation, $self->{filter_position}, $self->{fast_length_calculation}));
        if (check_for_exon_annotation_errors($transcript_variation)) {
            push(@filters, 'EXON_INTRON_UNDEF');
        } elsif (check_for_single_exon($transcript_variation)) {
            push(@filters, 'SINGLE_EXON');
        } else {
            #push(@filters, 'INCOMPLETE_CDS') if (check_incomplete_cds($transcript_variation));
            push(@filters, 'NON_CAN_SPLICE_SURR') if (check_surrounding_introns($transcript_variation, $self->{min_intron_size}));
        }
    }

    if ("splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences) {
        if (check_for_intron_annotation_errors($transcript_variation)) {
            push(@filters, 'EXON_INTRON_UNDEF');
        } else {
            push(@filters, 'SMALL_INTRON') if (check_intron_size($transcript_variation, $self->{min_intron_size}));
            push(@filters, 'NON_CAN_SPLICE') if (check_for_non_canonical_intron_motif($transcript_variation));
            push(@filters, 'NAGNAG_SITE') if (check_nagnag_variant($variation_feature));
        }
    }
    
    push(@filters, 'ANC_ALLELE') if (check_for_ancestral_allele($transcript_variation_allele));
    
    if ($confidence eq 'HC' && scalar @filters > 0) {
        $confidence = 'LC';
    }
    
    return { LoF => $confidence, LoF_filter => join(',', @filters) };
}

# Global functions
sub small_intron {
    my $transcript_variation = shift;
    my $intron_number = shift;
    my $min_intron_size = shift;
    my @gene_introns = @{$transcript_variation->transcript->get_all_Introns()};
    unless(defined($gene_introns[$intron_number])) {
        print "Issue with " . $transcript_variation->transcript->display_id . ": intron " . $intron_number . "\n";
        print $transcript_variation->intron_number . "\n";
    }
    
    if ($gene_introns[$intron_number]->length < $min_intron_size) {
        print "Small intron (" . $intron_number . ")\n" if $debug;
    }
    #print "Got length: " . $gene_introns[$intron_number]->length . "\n" if $ddebug;
    return ($gene_introns[$intron_number]->length < $min_intron_size);
}

sub intron_motif_start {
    my ($transcript_variation, $intron_number) = @_;
    print "Checking start motif\n" if $ddebug;
    
    my $transcript = $transcript_variation->transcript;
    my @gene_introns = @{$transcript->get_all_Introns()};
    my $sequence;
    if (exists($transcript->{intron_cache}->{$intron_number})) {
        print "Got intron cache!\n" if $ddebug;
        $sequence = $transcript->{intron_cache}->{$intron_number};
    } else {
        $sequence = $gene_introns[$intron_number]->seq;
        $transcript->{intron_cache}->{$intron_number} = $sequence;
    }
    if (substr($sequence, 0, 2) ne 'GT') {
        print "\nIssue with " . $transcript_variation->variation_feature->variation_name . " in " . $transcript_variation->transcript->display_id . "\n";
        print "Intron " . $intron_number . " (length: " . length($sequence) . "), seq is: " . $sequence . "\n";
    }
    print "Got start motif: Intron " . $intron_number . " (length: " . length($sequence) . ")\n" if $ddebug;
    return (substr($sequence, 0, 2) ne 'GT');
}

sub intron_motif_end {
    my ($transcript_variation, $intron_number) = @_;
    print "Checking end motif\n" if $ddebug;
    
    my $transcript = $transcript_variation->transcript;
    my @gene_introns = @{$transcript->get_all_Introns()};
    my $sequence;
    if (exists($transcript->{intron_cache}->{$intron_number})) {
        print "Got intron cache!\n" if $ddebug;
        $sequence = $transcript->{intron_cache}->{$intron_number};
    } else {
        $sequence = $gene_introns[$intron_number]->seq;
        $transcript->{intron_cache}->{$intron_number} = $sequence;
    }
    print "Got end motif: Intron " . $intron_number . " (length: " . length($sequence) . ")\n" if $ddebug;
    return (substr($sequence, length($sequence) - 2, 2) ne 'AG')
}

sub get_cds_length_fast {
    my $transcript = shift;
    
    my $transcript_cds_length = $transcript->cdna_coding_end - $transcript->cdna_coding_start + 1;
    
    print ('Tx (' . $transcript->display_id . ') length is: ' . ($transcript_cds_length) . "\n") if $ddebug;
    print "CDS length (" . $transcript_cds_length . ") not multiple of 3 (" . $transcript->display_id . ")\n" if ($debug && $transcript_cds_length % 3);
    
    return $transcript_cds_length;
}

sub get_cds_length {
    my $transcript = shift;
    
    my $transcript_cds_length;
    if (exists($transcript->{seq_cache})) {
        print "Got seq cache\n" if $ddebug;
        $transcript_cds_length = length($transcript->{seq_cache});
    } else {
        my $seq = $transcript->translateable_seq;
        $transcript_cds_length = length($seq);
        $transcript->{seq_cache} = $seq;
    }
    return $transcript_cds_length;
}

# Stop-gain and frameshift annotations
sub check_incomplete_cds {
    my $transcript_variation = shift;
    
    my $transcript = $transcript_variation->transcript;
    my $start_annotation = $transcript->get_all_Attributes('cds_start_NF');
    my $end_annotation = $transcript->get_all_Attributes('cds_end_NF');
    if (defined($start_annotation)) {
        foreach my $annot (@{$start_annotation}) {
            print "Start annotation: " . $annot->code . "\n";
        }
    }
    if (defined($end_annotation)) {
        foreach my $annot (@{$end_annotation}) {
            print "End annotation: " . $annot->code .  "\n";
        }
    }
    return (defined($start_annotation) || defined($end_annotation));
}

sub check_for_exon_annotation_errors {
    my $transcript_variation = shift;
    return (!defined($transcript_variation->exon_number))
}

sub check_position {
    my ($transcript_variation, $cutoff, $speed) = @_;
    
    my $transcript_cds_length;
    if ($speed eq 'fast') {
        $transcript_cds_length = get_cds_length_fast($transcript_variation->transcript);
    } else {
        $transcript_cds_length = get_cds_length($transcript_variation->transcript);
    }
    my $variant_cds_position = $transcript_variation->cdna_end;
    
    print $transcript_variation->transcript->stable_id . ": " . $variant_cds_position . " out of " . $transcript_cds_length . "\n" if $ddebug;
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
    
    print "Checking exon " . $exon_number . " (out of " . ($total_exons - 1) . " exons) \n" if $ddebug;
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
    
    my $seq;
    if (exists($variation_feature->{splice_context_seq_cache})) {
        print "Got splice context seq cache\n" if $ddebug;
        $seq = $variation_feature->{splice_context_seq_cache};
    } else {
        $seq = uc($variation_feature->feature_Slice->expand(4, 4)->seq);
        print "Got splice context seq: " . $seq . "\n" if $ddebug;
        $variation_feature->{splice_context_seq_cache} = $seq;
    }
    return (length($seq) == 9 && $seq =~ m/AG.AG/);
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

sub check_for_ancestral_allele_faidx {
    my $transcript_variation_allele = shift;
    my $variation_feature = $transcript_variation_allele->variation_feature;
    
    my $aff_allele = $transcript_variation_allele->variation_feature_seq;
    print "Allele is " . $aff_allele . "\n" if $ddebug;
    
    my $region = $variation_feature->seq_region_name() . ":" . $variation_feature->seq_region_start() . '-' . $variation_feature->seq_region_end();
    my $faidx = `samtools faidx human_ancestor.fa.rz $region`;
    my @lines = split(/\n/, $faidx);
    shift @lines;
    my $ancestral_allele = uc(join('', @lines));
    print $variation_feature->variation_name . "; Allele is: " . $aff_allele . "; Ancestral allele is: " . $ancestral_allele . "\n" if $ddebug;
    
    return ($ancestral_allele eq $aff_allele)
}

sub check_for_ancestral_allele_var_api {
    my $transcript_variation_allele = shift;
    my $variation_feature = $transcript_variation_allele->variation_feature;
    
    my $aff_allele = $transcript_variation_allele->variation_feature_seq;
    print "Allele is " . $aff_allele . "\n" if $ddebug;
    
    # This doesn't seem to work
    my $ancestral_allele = $variation_feature->variation->ancestral_allele;
    
    print $variation_feature->variation_name . "; Allele is: " . $aff_allele . "; Ancestral allele is: " . $ancestral_allele . "\n" if $ddebug;
    
    return ($ancestral_allele eq $aff_allele)
}

sub check_for_ancestral_allele {
    my $transcript_variation_allele = shift;
    return check_for_ancestral_allele_faidx($transcript_variation_allele)
}


1;