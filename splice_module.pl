use strict;
require "score5.pl";
require "score3.pl";

sub check_5UTR_splice {
    my $transcript_variation = shift;
    my $variation_feature = shift;
    my $strand = $variation_feature->{strand};

    # Start is always less than or equal to end regardless of the orientation of the slice (true for introns as well).
    my $slice = $variation_feature->feature_Slice();

    # filter if variant occurs upstream of first intron
    if ($strand == 1) {
        my $cds_start = $transcript_variation->transcript->{coding_region_start};
        return ($slice->end() < $cds_start);
    } elsif ($strand == -1) {
        my $cds_start = $transcript_variation->transcript->{coding_region_end};
        return ($slice->start() > $cds_start);
    }
}

sub check_3UTR_splice {
    my $transcript_variation = shift;
    my $variation_feature = shift;
    my $strand = $variation_feature->{strand};

    # Start is always less than or equal to end regardless of the orientation of the slice (true for introns as well).
    my $slice = $variation_feature->feature_Slice();
    
    # filter if variant occurs downstream of last intron (i.e. not actually a splice donor site)
    if ($strand == 1) {
        my $cds_end = $transcript_variation->transcript->{coding_region_end};
        return ($slice->start() > $cds_end);
    } elsif ($strand == -1) {
        my $cds_end = $transcript_variation->transcript->{coding_region_start};
        return ($slice->end() < $cds_end);
    }
}

# Determine whether two given intervals overlap eachother. Returns boolean.
sub overlap {
    my ($a, $b, $c, $d) = @_[0..3];
    if ($a < $c) {
        return ($b >= $c);
    } elsif ($a > $c) {
        return ($a <= $d);
    } else {
        return 1;
    }
}

sub check_if_in_splice_site {
    my $transcript_variation = shift;
    my $variation_feature = shift;
    my $allele_string = shift;

    my $slice = $variation_feature->feature_Slice();
    my $strand = $transcript_variation->transcript->strand();
    my ($ref, $alt) = split /\//, ($allele_string);
    my $insertion = ($ref eq "-"); # slice start > slice end by 1

    my $donor = 0;
    my $acceptor = 0;
    if ($transcript_variation->exon_number) {
        my @exons = @{ $transcript_variation->transcript->get_all_Exons };
        my ($exon_num, $number_of_exons) = split /\//, ($transcript_variation->exon_number);
        my $exon = $exons[$exon_num - 1];
        if ($strand == 1) {
            if ($insertion) {
                $donor = ($exon->{end} - 2 <= $slice->{end}) && ($exon->{end} + 6 >= $slice->{start});
                $acceptor = ($exon->{start} - 20 <= $slice->{end}) && ($exon->{start} + 2 >= $slice->{start});
            } else {
                $donor = overlap($slice->{start}, $slice->{end}, $exon->{end} - 2, $exon->{end} + 6);
                $acceptor = overlap($slice->{start}, $slice->{end}, $exon->{start} - 20, $exon->{start} + 2);
            }
        } elsif ($strand == -1) {
            if ($insertion) {
                $donor = ($exon->{start} - 6 <= $slice->{end}) && ($exon->{start} + 2 >= $slice->{start});
                $acceptor = ($exon->{end} - 2 <= $slice->{end}) && ($exon->{end} + 20 >= $slice->{start});
            } else {
                $donor = overlap($slice->{start}, $slice->{end}, $exon->{start} - 6, $exon->{start} + 2); 
                $acceptor = overlap($slice->{start}, $slice->{end}, $exon->{end} - 2, $exon->{end} + 20); 
            }
        }
    } 

    elsif ($transcript_variation->intron_number) {
        my @introns = @{ $transcript_variation->transcript->get_all_Introns };
        my ($intron_num, $number_of_introns) = split /\//, ($transcript_variation->intron_number);
        my $intron = $introns[$intron_num - 1];
        if ($strand == 1) {
            if ($insertion) {
                $donor = ($intron->{start} - 3 <= $slice->{end}) && ($intron->{start} + 5 >= $slice->{start});
                $acceptor = ($intron->{end} - 19 <= $slice->{end}) && ($intron->{end} + 3 >= $slice->{start});
            } else {
                $donor = overlap($slice->{start}, $slice->{end}, $intron->{start} - 3, $intron->{start} + 5);
                $acceptor = overlap($slice->{start}, $slice->{end}, $intron->{end} - 19, $intron->{end} + 3);
            }
        } elsif ($strand == -1) {
            if ($insertion) {
                $donor = ($intron->{end} - 5 <= $slice->{end}) && ($intron->{end} + 3 >= $slice->{start});
                $acceptor = ($intron->{start} - 3 <= $slice->{end}) && ($intron->{start} + 19 >= $slice->{start});
            } else {
                $donor = overlap($slice->{start}, $slice->{end}, $intron->{end} - 5, $intron->{end} + 3);
                $acceptor = overlap($slice->{start}, $slice->{end}, $intron->{start} - 3, $intron->{start} + 19);
            }
        }
    }
    else {
        # INSERTION OCCURRING RIGHT A INTRON-EXON JUNCTION
        return (2, 2);
    }
    return ($donor, $acceptor);
}

sub get_exon {
    my $transcript_variation = $_[0];
    my $delta = defined $_[1] ? $_[1] : 0;
    my ($exon_num, $number_of_exons) = split /\//, ($transcript_variation->exon_number);
    my $exon_idx = $exon_num - 1 + $delta;
    my @exons = @{ $transcript_variation->transcript->get_all_Exons };
    my $exon = $exons[$exon_idx];
    return ($exon, $exon_idx, $number_of_exons);
}

sub get_intron {
    my $transcript_variation = shift;
    my $delta = defined $_[0] ? $_[0] : 'default';
    my ($intron_num, $number_of_introns) = split /\//, ($transcript_variation->intron_number);
    my $intron_idx = $intron_num - 1 + $delta;
    my @introns = @{ $transcript_variation->transcript->get_all_Introns };
    my $intron = $introns[$intron_idx];
    return ($intron, $intron_idx, $number_of_introns);
}

sub truncate_utr_from_exon {
    my $transcript_variation = shift;
    my $exon = shift;

    my $strand = $transcript_variation->transcript->strand();
    # if exon is first/last exon, adjust it's start/end position to be position of start/stop codon
    my $start_codon_pos = 0;
    my $stop_codon_pos = 0;
    my $exon_start = $exon->{start};
    my $exon_end = $exon->{end};
    # adjust exon start and end to agree with location of start/stop codons
    if ($strand == 1) {
        $start_codon_pos = $transcript_variation->transcript->{coding_region_start};
        $stop_codon_pos = $transcript_variation->transcript->{coding_region_end};
        # exon spans start codon
        if ($exon_start < $start_codon_pos && $exon_end > $start_codon_pos) {
            $exon_start = $start_codon_pos;
        } 
        # exon spans stop codon
        if ($exon_end > $stop_codon_pos && $exon_start < $stop_codon_pos) {
            $exon_end = $stop_codon_pos;
        }
    } elsif ($strand == -1) {
        $start_codon_pos = $transcript_variation->transcript->{coding_region_end};
        $stop_codon_pos = $transcript_variation->transcript->{coding_region_start};
        # exon spans start codon
        if ($exon_end > $start_codon_pos && $exon_start < $start_codon_pos) {
            $exon_end = $start_codon_pos;
        } 
        # exon spans stop codon
        if ($exon_start < $stop_codon_pos && $exon_end > $stop_codon_pos) {
            $exon_start = $stop_codon_pos;
        }
    }
    return ($exon_start, $exon_end);
}

# for fully exonic and fully intronic variants:
    # mut_junc_idx is defined as wherever the exon ends
# for junction-spanning variants:
    # mut_junc_idx is defined as the nearest location which preserves both reading frame AND length of exon
sub get_ruler_cis {
    my $transcript_variation = shift;
    my $variation_feature = shift;
    my $donor = shift;
    my $allele_string = shift;
    my %ruler = ();
    
    $ruler{donor} = $donor;
    my $strand = $transcript_variation->transcript->strand();
    # exonic
    my $exon, my $intron;
    if ($transcript_variation->exon_number) {
        ($exon, my $exon_idx) = get_exon($transcript_variation);
        my @introns = @{ $transcript_variation->transcript->get_all_Introns };
        $intron = ($donor) ? $introns[$exon_idx] : $introns[$exon_idx - 1];
        $ruler{exon_idx} = $exon_idx; # need for determining reading frame later
    # intronic
    } elsif ($transcript_variation->intron_number) {
        ($intron, my $intron_idx) = get_intron($transcript_variation);
        my @exons = @{ $transcript_variation->transcript->get_all_Exons };
        $exon = ($donor) ? $exons[$intron_idx] : $exons[$intron_idx + 1];
        $ruler{exon_idx} = ($donor) ? $intron_idx : $intron_idx + 1; # need for determining reading frame later
    }

    $ruler{original_intron_length} = $intron->{end} - $intron->{start} + 1;
    $ruler{intronic} = ($transcript_variation->intron_number) ? 1 : 0;
    $ruler{exonic} = ($transcript_variation->exon_number) ? 1: 0;

    my $slice = $variation_feature->feature_Slice();
    my ($ref, $alt) = split /\//, ($allele_string);
    # change slice bounds for insertions
    if ($ref eq "-") {
        if ($strand == 1) {
            $slice->{end} = $slice->{start};
        } else {
            $slice->{start} = $slice->{end};
        }
    }

    my $left_bound, my $right_bound, my $dist_to_mut_junc, my $dist_to_ref_junc;
    my $insertion = ($ref eq "-");
    my $deletion = ($alt eq "-");
    if ($donor) {
        if ($strand == 1) {
            $left_bound = $exon->{start};
            $right_bound = $intron->{end} - 10;
            if ($ruler{intronic}) {
                $dist_to_mut_junc = $intron->{start} - $slice->{start};
            } else {
                $dist_to_mut_junc = $intron->{start} - $slice->{end};
                $dist_to_mut_junc = $dist_to_mut_junc + length $alt if ($insertion);
                $dist_to_mut_junc = $dist_to_mut_junc - 1 if ($deletion);
            }
            $dist_to_ref_junc = $intron->{start} - $slice->{start};
        } elsif ($strand == -1) {
            $left_bound = $intron->{start} + 10;
            $right_bound = $exon->{end};
            if ($ruler{intronic}) {
                $dist_to_mut_junc = $slice->{end} - $intron->{end};
            } else {
                $dist_to_mut_junc = $slice->{start} - $intron->{end};
                $dist_to_mut_junc = $dist_to_mut_junc + length $alt if ($insertion);
                $dist_to_mut_junc = $dist_to_mut_junc - 1 if ($deletion);
            }
            $dist_to_ref_junc = $slice->{end} - $intron->{end};
        }
    } else {
        if ($strand == 1) {
            $left_bound = $intron->{start} + 10;
            $right_bound = $exon->{end};
            if ($ruler{intronic}) {
                $dist_to_mut_junc = $exon->{start} - $slice->{end};
                $dist_to_mut_junc = $dist_to_mut_junc + length $alt if ($insertion);
                $dist_to_mut_junc = $dist_to_mut_junc - 1 if ($deletion);
            } else {
                $dist_to_mut_junc = $exon->{start} - $slice->{start};
            }
            $dist_to_ref_junc = $exon->{start} - $slice->{start};
        } elsif ($strand == -1) {            
            $left_bound = $exon->{start};
            $right_bound = $intron->{end} - 10;
            if ($ruler{intronic}) {
                $dist_to_mut_junc = $slice->{start} - $exon->{end};
                $dist_to_mut_junc = $dist_to_mut_junc + length $alt if ($insertion);
                $dist_to_mut_junc = $dist_to_mut_junc - 1 if ($deletion);
            } else {
                $dist_to_mut_junc = $slice->{end} - $exon->{end};
            }
            $dist_to_ref_junc = $slice->{end} - $exon->{end};
        }
    }
    my $left_flank = $slice->{start} - $left_bound;
    my $right_flank = $right_bound - $slice->{end};  
    $slice->{strand} = 1; # weird thing here, but necessary because function assumes slice expansion is positive direction    
    my $seq = uc($slice->expand($left_flank, $right_flank)->seq());
    $seq = ($strand == -1) ? reverse_complement($seq)->seq() : $seq;
    $ruler{reference_seq} = $seq;
    if ($strand == -1) {
        my $old = $right_flank;
        $right_flank = $left_flank;
        $left_flank = $old;
    }
    $ruler{mut_junc_idx} = $left_flank + $dist_to_mut_junc;
    $ruler{ref_junc_idx} = $left_flank + $dist_to_ref_junc;

    # print ">$slice->{start}-$slice->{end} : $ref -> $alt\n";
    # insertion
    if ($ref eq "-") {
        $alt = ($strand == -1) ? reverse_complement($alt)->seq() : $alt;
        $ruler{mutated_seq} = (substr $seq, 0, $left_flank) . $alt . (substr $seq, $left_flank);
        $ruler{indel_size} = length $alt;
    # deletion
    } elsif ($alt eq "-") {
        $ruler{mutated_seq} = (substr $seq, 0, $left_flank) . (substr $seq, ($left_flank + length $ref));
        $ruler{indel_size} = (length $ref) * -1;
    # snp
    } else {
        $alt = ($strand == -1) ? reverse_complement($alt)->seq() : $alt;
        $ruler{mutated_seq} = (substr $seq, 0, $left_flank) . $alt . (substr $seq, $left_flank + 1);
        $ruler{indel_size} = 0;
    }
    return \%ruler;
}

sub get_consensus_ref {
    my $ruler_ref = shift;
    my $donor = shift;
    my $dist = shift;
    my %ruler = %$ruler_ref;
    my $consensus;
    if ($donor) {
        $consensus = substr $ruler{reference_seq}, $ruler{ref_junc_idx} - 3 + $dist, 9;
    } else {
        $consensus = substr $ruler{reference_seq}, $ruler{ref_junc_idx} - 20 + $dist, 23;
    } 
}

sub scan_for_donor_rescue {
    my $ruler_ref = shift;
    my %ruler = %$ruler_ref;
    my $seq = $ruler{mutated_seq};
    my $l = length $seq;
    my $a = $ruler{mut_junc_idx} - 3;
    my $b = $ruler{mut_junc_idx} + 3;
    my $authentic = get_consensus_ref($ruler_ref, 1, 0);
    my $authentic_mes = mes_donor($authentic);

    for (my $i=3; $i < $l - 5; $i++) {
        my $dist = $i - $ruler{mut_junc_idx};
        # my $donor = get_consensus_alt($ruler_ref, 1, $dist);
        my $donor = substr $ruler{mutated_seq}, $i - 3, 9;
        my $e_delta = ($ruler{intronic} == 1) ? $dist : $dist + $ruler{indel_size}; # net change in exonic bases
        my $inframe = (($e_delta % 3) == 0) ? 1 : 0;
        my $i_delta = ($l - $i + 10) - $ruler{original_intron_length};
        my $mes = mes_donor($donor);
        # check if rescue
        if (($mes - $authentic_mes) > -7.2 && $e_delta > -16 && $e_delta < 7 && $inframe) {
            return 1;
        } 
    }
    return 0;
}

sub scan_for_acceptor_rescue {
    my $ruler_ref = shift;
    my %ruler = %$ruler_ref;
    my $seq = $ruler{mutated_seq};
    my $l = length $seq;
    my $a = $ruler{mut_junc_idx} - 3;
    my $b = $ruler{mut_junc_idx} + 3;
    my $authentic = get_consensus_ref($ruler_ref, 0, 0);
    my $authentic_mes = mes_acceptor($authentic);
    # print "authentic=$authentic, authentic_mes=$authentic_mes\n";

    for (my $i=20; $i < $l - 2; $i++) {
        my $dist = $i - $ruler{mut_junc_idx};
        # my $acceptor = get_consensus_alt($ruler_ref, 0, $dist);
        my $acceptor = substr $ruler{mutated_seq}, $i - 20, 23;
        my $e_delta = ($ruler{intronic} == 1) ? $dist : $dist - $ruler{indel_size}; # net change in exonic bases
        my $inframe = (($e_delta % 3) == 0) ? 1 : 0;
        my $i_delta = $i + 10 - $ruler{original_intron_length};
        my $mes = mes_acceptor($acceptor);
        # check if rescue
        if (($mes - $authentic_mes) > -5 && $e_delta > -13 && $e_delta < 7 && $inframe) {
            return 1;
        } 
    }
    return 0;
}

# parameters based on maximizing sensitivity while ensuring a fpr < 0.03
sub check_for_splice_disruption {
    my $ruler_ref = shift;
    my $donor = shift;
    my %ruler = %$ruler_ref;

    my $flag = 'fine';
    my $i = $ruler{mut_junc_idx};
    if ($donor) {
        my $disrupted = ($donor) ? substr $ruler{mutated_seq}, $i - 3, 9 : substr $ruler{mutated_seq}, $i - 20, 23;
        my $authentic = get_consensus_ref($ruler_ref, 1, 0);
        my $mes_old = mes_donor($authentic);
        my $mes_new = mes_donor($disrupted);
        $flag = 'DISRUPTED_DONOR' if ($mes_new - $mes_old) < -4.52;
    } else {
        my $disrupted = substr $ruler{mutated_seq}, $i - 20, 23;
        my $authentic = get_consensus_ref($ruler_ref, 0, 0);
        my $mes_old = mes_acceptor($authentic);
        my $mes_new = mes_acceptor($disrupted);
        $flag = 'DISRUPTED_ACCEPTOR' if ($mes_new - $mes_old) < -6.04;
    }
    return $flag; 
}


sub scan_for_donor {
    my $ruler_ref = shift;
    my %ruler = %$ruler_ref;
    my $seq = $ruler{mutated_seq};
    my $l = length $seq;
    my $a = $ruler{mut_junc_idx} - 3;
    my $b = $ruler{mut_junc_idx} + 3;
    my $authentic = get_consensus_ref($ruler_ref, 1, 0);
    my $authentic_mes = mes_donor($authentic);

    my @candidate_mes = ();
    my @candidate_dist = ();
    for (my $i=3; $i < $l - 5; $i++) {
        my $dist = $i - $ruler{mut_junc_idx};
        # my $donor = get_consensus_alt($ruler_ref, 1, $dist);
        my $donor = substr $ruler{mutated_seq}, $i - 3, 9;
        my $e_delta = ($ruler{intronic} == 1) ? $dist : $dist + $ruler{indel_size}; # net change in exonic bases
        my $inframe = (($e_delta % 3) == 0) ? 1 : 0;
        my $i_delta = ($l - $i + 10) - $ruler{original_intron_length};
        my $mes = mes_donor($donor);
        push(@candidate_mes, $mes);
        push(@candidate_dist, $dist);
    }
    return (join('_', @candidate_mes), join('_', @candidate_dist), $authentic_mes);
}

# return scalar original MES, array of new MES, array of corresponding distances
sub scan_for_acceptor {
    my $ruler_ref = shift;
    my %ruler = %$ruler_ref;
    my $seq = $ruler{mutated_seq};
    my $l = length $seq;
    my $a = $ruler{mut_junc_idx} - 3;
    my $b = $ruler{mut_junc_idx} + 3;
    my $authentic = get_consensus_ref($ruler_ref, 0, 0);
    my $authentic_mes = mes_acceptor($authentic);
    
    my @candidate_mes = ();
    my @candidate_dist = ();
    for (my $i=20; $i < $l - 2; $i++) {
        my $dist = $i - $ruler{mut_junc_idx};
        my $acceptor = substr $ruler{mutated_seq}, $i - 20, 23;
        my $e_delta = ($ruler{intronic} == 1) ? $dist : $dist - $ruler{indel_size}; # net change in exonic bases
        my $inframe = (($e_delta % 3) == 0) ? 1 : 0;
        my $i_delta = $i + 10 - $ruler{original_intron_length};
        my $mes = mes_acceptor($acceptor);
        push(@candidate_mes, $mes);
        push(@candidate_dist, $dist);
    }
    return (join('_', @candidate_mes), join('_', @candidate_dist), $authentic_mes);
}


1;