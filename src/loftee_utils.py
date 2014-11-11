__author__ = 'konradjk'

# Note that this is the current as of v77 with 2 included for backwards compatibility (VEP <= 75)
csq_order = ["transcript_ablation",
"splice_donor_variant",
"splice_acceptor_variant",
"stop_gained",
"frameshift_variant",
"stop_lost",
"initiator_codon_variant",
"transcript_amplification",
"inframe_insertion",
"inframe_deletion",
"missense_variant",
"splice_region_variant",
"incomplete_terminal_codon_variant",
"stop_retained_variant",
"synonymous_variant",
"coding_sequence_variant",
"mature_miRNA_variant",
"5_prime_UTR_variant",
"3_prime_UTR_variant",
"non_coding_transcript_exon_variant",
"non_coding_exon_variant",  # deprecated
"intron_variant",
"NMD_transcript_variant",
"non_coding_transcript_variant",
"nc_transcript_variant",  # deprecated
"upstream_gene_variant",
"downstream_gene_variant",
"TFBS_ablation",
"TFBS_amplification",
"TF_binding_site_variant",
"regulatory_region_ablation",
"regulatory_region_amplification",
"regulatory_region_variant",
"feature_elongation",
"feature_truncation",
"intergenic_variant",
""]
csq_order_dict = dict(zip(csq_order, range(len(csq_order))))
rev_csq_order_dict = dict(zip(range(len(csq_order)), csq_order))


def csq_max_compare_lists(ann1, ann2):
    return rev_csq_order_dict[min(csq_max_list(ann1), csq_max_list(ann2))]


def csq_max_vep(ann_list):
    return rev_csq_order_dict[csq_max_list(ann_list.split('&'))]


def csq_max(ann_list):
    if len(ann_list) == 0: return ''
    return rev_csq_order_dict[csq_max_list(ann_list)]


def csq_max_list(ann_list):
    return min([csq_order_dict[ann] for ann in ann_list])