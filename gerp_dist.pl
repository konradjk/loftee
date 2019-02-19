use strict;
use Bio::DB::BigWig 'binMean';

sub get_gerp_weighted_dist {
    my ($tr, $pos, $bigwig_file) = @_[0..2];

    # collect some variables
    my $chr = $tr->seq_region_name();
    my @exons = @{ $tr->get_all_Exons };
    my $strand = $tr->strand();
    my $transcript_id = $tr->{stable_id};
    my $number_of_exons = scalar @exons;

    # determine boundaries of CDS sequence
    my ($stop_codon_pos, $start_codon_pos);
    if ($strand == 1) {
        $stop_codon_pos = $tr->{coding_region_end};
        $start_codon_pos = $tr->{coding_region_start};
    } elsif ($strand == -1) {
        $stop_codon_pos = $tr->{coding_region_start};
        $start_codon_pos = $tr->{coding_region_end};
    }

    # get distance to from variant to stop codon, weighted by GERP
    my $weighted_dist = 0;
    my $dist = 0;
    for (my $i=0; $i <= $number_of_exons - 1; $i++) {
        my $current_exon = $exons[$i];
        # skip exons upstream of variant
        if ($strand == -1) {
            next if $pos < $current_exon->start;
        } else {
            next if $pos > $current_exon->end;
        }

        # determine if last exon by checking if exon spans stop codon position
        my $last_exon = 0;
        if ($strand == 1) {
            $last_exon = ($current_exon->start < $stop_codon_pos) && ($current_exon->end >= $stop_codon_pos);
        } else {
            $last_exon = ($current_exon->end > $stop_codon_pos) && ($current_exon->start <= $stop_codon_pos);
        }

        # get contribution of current exon to total weighted distance
        my ($start, $end, $wd);
        my $in_affected_exon = ($pos >= $current_exon->start) && ($pos <= $current_exon->end);
        if ($last_exon) {
            if ($in_affected_exon) {
                $start = $pos;
            } else {
                $start = ($strand == 1) ? $current_exon->start : $current_exon->end;
            }
            $end = $stop_codon_pos;
        } elsif ($in_affected_exon) {
            $start = $pos;
            $end = ($strand == 1) ? $current_exon->{end} : $current_exon->{start};
        } else {
            # my $exon_num = $i + 1;
            # $wd = get_exon_gerp($transcript_id, $exon_num, $cons_db);
            $start = $current_exon->start;
            $end = $current_exon->end;
        }
        $wd = get_interval_gerp($chr, $start, $end, $bigwig_file);
        $weighted_dist = $weighted_dist + $wd;
        $dist = $dist + (abs ($end - $start));
    }
    return ($weighted_dist, $dist);
}

sub get_interval_gerp {
    my ($chrom, $a, $b, $bigwig_file) = @_[0..3];
    my $wig = Bio::DB::BigWig->new(-bigwig=>$bigwig_file);
    my @feats = $wig->features(-type=>'summary', -seq_id=>$chrom, -start=>$a, -end=>$b);
    my $total = 0;
    for my $c (@feats) {
        $total += $c->length * binMean($c->score());
    }
    return ($total);
}

1;
