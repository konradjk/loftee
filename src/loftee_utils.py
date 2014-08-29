__author__ = 'konradjk'

csq_order = ["transcript_ablation",
"splice_donor_variant",
"splice_acceptor_variant",
"stop_gained",
"frameshift_variant",
"stop_lost",
"initiator_codon_variant",
"inframe_insertion",
"inframe_deletion",
"missense_variant",
"transcript_amplification",
"splice_region_variant",
"incomplete_terminal_codon_variant",
"synonymous_variant",
"stop_retained_variant",
"coding_sequence_variant",
"mature_miRNA_variant",
"5_prime_UTR_variant",
"3_prime_UTR_variant",
"non_coding_exon_variant",
"nc_transcript_variant",
"intron_variant",
"NMD_transcript_variant",
"upstream_gene_variant",
"downstream_gene_variant",
"TFBS_ablation",
"TFBS_amplification",
"TF_binding_site_variant",
"regulatory_region_variant",
"regulatory_region_ablation",
"regulatory_region_amplification",
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