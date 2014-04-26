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

 A VEP plugin that determines loss-of-function variation.

=cut

package LoF;

use strict;
use warnings;

our $debug = 0;
our $ddebug = 0;
our $slow = 1;

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
    $self->{has_cache} = 1;
    
    foreach my $parameter (@{$self->params}) {
        my @param = split /:/, $parameter;
        if (scalar @param == 2) {
            $self->{$param[0]} = $param[1];
        }
    }
    
    $self->{filter_position} = $self->{filter_position} || 0.05;
    $self->{min_intron_size} = $self->{min_intron_size} || 15;
    
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
    
    print $transcript_variation->transcript->stable_id() . ": " . $transcript_variation->transcript->description() . "\n";
    my @consequences = map { $_->SO_term } @{ $transcript_variation_allele->get_all_OverlapConsequences };
    
    my $confidence = 'HC';
    my @filters = ();
    
    # Filter in
    unless ($transcript_variation->transcript->biotype eq "protein_coding") {
        return {}
    }
    unless ("stop_gained" ~~ @consequences || "splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences || "frameshift_variant" ~~ @consequences) {
        return {}
    }
    
    # Filter out
    if ("stop_gained" ~~ @consequences || "frameshift_variant" ~~ @consequences){
        push(@filters, 'EXON_INTRON_UNDEF') if (check_for_exon_annotation_errors($transcript_variation));
        push(@filters, 'SINGLE_EXON') if (check_for_single_exon($transcript_variation));
        push(@filters, 'END_TRUNC') if (check_position($transcript_variation, $self->{filter_position}));
        push(@filters, 'NON_CAN_SPLICE_SURR') if (check_surrounding_introns($transcript_variation, $self->{min_intron_size}));
    }

    if ("splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences) {
        push(@filters, 'EXON_INTRON_UNDEF') if (check_for_intron_annotation_errors($transcript_variation));
        push(@filters, 'SMALL_INTRON') if (check_intron_size($transcript_variation, $self->{min_intron_size}));
        push(@filters, 'NON_CAN_SPLICE') if (check_for_non_canonical_intron_motif($transcript_variation));
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
    if ($gene_introns[$intron_number]->length < $min_intron_size) {
        print "Small intron (" . $intron_number . ")\n" if $debug;
    }
    #print "Got length: " . $gene_introns[$intron_number]->length . "\n" if $ddebug;
    return ($gene_introns[$intron_number]->length < $min_intron_size);
}

sub intron_motif_start {
    my ($transcript_variation, $intron_number) = @_;
    print "Checking start motif\n" if $ddebug;
    my @gene_introns = @{$transcript_variation->transcript->get_all_Introns()};
    my $sequence = $gene_introns[$intron_number]->seq;
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
    my @gene_introns = @{$transcript_variation->transcript->get_all_Introns()};
    my $sequence = $gene_introns[$intron_number]->seq;
    print "Got end motif: Intron " . $intron_number . " (length: " . length($sequence) . ")\n" if $ddebug;
    return (substr($sequence, length($sequence) - 2, 2) ne 'AG')
}

sub get_cds_length {
    my $transcript_variation = shift;
    
    my $transcript_cds_length = $transcript_variation->transcript->cdna_coding_end - $transcript_variation->transcript->cdna_coding_start + 1;
    
    print ('Tx (' . $transcript_variation->transcript->display_id . ') length is: ' . ($transcript_variation->transcript->cdna_coding_end - $transcript_variation->transcript->cdna_coding_start + 1) . "\n") if $ddebug;
    print "CDS length (" . $transcript_cds_length . ") not multiple of 3 (" . $transcript_variation->transcript->display_id . ")\n" if ($debug && $transcript_cds_length % 3);
    
    return $transcript_cds_length;
}

# Stop-gain and frameshift annotations
sub check_for_exon_annotation_errors {
    my $transcript_variation = shift;
    return (!defined($transcript_variation->exon_number))
}

sub check_position {
    my ($transcript_variation, $cutoff) = @_;
    
    my $transcript_cds_length = get_cds_length($transcript_variation);
    #my $variant_cds_position = $transcript_variation->cdna_start;
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
    
    # This doesn't seem to work for some reason
    my $ancestral_allele = $variation_feature->variation->ancestral_allele;
    
    print $variation_feature->variation_name . "; Allele is: " . $aff_allele . "; Ancestral allele is: " . $ancestral_allele . "\n" if $ddebug;
    
    return ($ancestral_allele eq $aff_allele)
}

sub check_for_ancestral_allele {
    my $transcript_variation_allele = shift;
    return check_for_ancestral_allele_faidx($transcript_variation_allele)
}

#sub check_for_ancestral_allele_api {
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