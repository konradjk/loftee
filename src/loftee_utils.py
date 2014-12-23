__author__ = 'konradjk'
from operator import itemgetter
import re

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


def worst_csq_index(csq_list):
    """
    Input list of consequences (e.g. ['frameshift_variant', 'missense_variant'])
    Return index of the worst annotation (In this case, index of 'frameshift_variant', so 4)
    Works well with csqs = 'non_coding_exon_variant&nc_transcript_variant' by worst_csq_index(csqs.split('&'))

    :param annnotation:
    :return most_severe_consequence_index:
    """
    return min([csq_order_dict[ann] for ann in csq_list])


def worst_csq_from_list(csq_list):
    """
    Input list of consequences (e.g. ['frameshift_variant', 'missense_variant'])
    Return the worst annotation (In this case, 'frameshift_variant')
    Works well with csqs = 'non_coding_exon_variant&nc_transcript_variant' by worst_csq_from_list(csqs.split('&'))

    :param annnotation:
    :return most_severe_consequence:
    """
    return rev_csq_order_dict[worst_csq_index(csq_list)]


def worst_csq_from_csq(csq):
    """
    Input possibly &-filled csq string (e.g. 'non_coding_exon_variant&nc_transcript_variant')
    Return the worst annotation (In this case, 'non_coding_exon_variant')

    :param consequence:
    :return most_severe_consequence:
    """
    return rev_csq_order_dict[worst_csq_index(csq.split('&'))]


def order_vep_by_csq(annotation_list):
    output = sorted(annotation_list, cmp=lambda x, y: compare_two_consequences(x, y), key=itemgetter('Consequence'))
    for ann in output:
        ann['major_consequence'] = worst_csq_from_csq(ann['Consequence'])
    return output


def worst_csq_with_vep(annotation_list):
    """
    Takes list of VEP annotations [{'Consequence': 'frameshift', Feature: 'ENST'}, ...]
    Returns most severe annotation (as full VEP annotation [{'Consequence': 'frameshift', Feature: 'ENST'}])
    Also tacks on worst consequence for that annotation (i.e. worst_csq_from_csq)
    :param annotation_list:
    :return worst_annotation:
    """
    if len(annotation_list) == 0: return None
    worst = annotation_list[0]
    for annotation in annotation_list:
        if compare_two_consequences(annotation['Consequence'], worst['Consequence']) < 0:
            worst = annotation
        elif compare_two_consequences(annotation['Consequence'], worst['Consequence']) == 0 and annotation['CANONICAL'] == 'YES':
            worst = annotation
    worst['major_consequence'] = worst_csq_from_csq(worst['Consequence'])
    return worst


def compare_two_consequences(csq1, csq2):
    if csq_order_dict[worst_csq_from_csq(csq1)] < csq_order_dict[worst_csq_from_csq(csq2)]:
        return -1
    elif csq_order_dict[worst_csq_from_csq(csq1)] == csq_order_dict[worst_csq_from_csq(csq2)]:
        return 0
    return 1


def simplify_polyphen(polyphen_list):
    """
    Takes list of polyphen score/label pairs (e.g. ['probably_damaging(0.968)', 'benign(0.402)'])
    Returns worst (worst label and highest score) - in this case, 'probably_damaging(0.968)'
    """
    max_score = 0.0
    max_label = 'benign'
    for polyphen in polyphen_list:
        label, score = polyphen.rstrip(')').split('(')
        if float(score) > max_score:
            max_score = float(score)
            max_label = label
    return max_score, max_label


def simplify_sift(sift_list):
    """
    Takes list of SIFT score/label pairs (e.g. ['tolerated(0.26)', 'deleterious(0)'])
    Returns worst (worst label and highest score) - in this case, 'deleterious(0)'
    """
    max_score = 1.0
    max_label = 'tolerated'
    for sift in sift_list:
        label, score = sift.rstrip(')').split('(')
        if float(score) < max_score:
            max_score = float(score)
            max_label = label
    return max_score, max_label

POLYPHEN_SIFT_REGEX = re.compile('^[a-z]+\([0-9\.]+\)$', re.IGNORECASE)


def simplify_polyphen_sift(input_list, type):
    if len(input_list) == 0 or not all([POLYPHEN_SIFT_REGEX.match(x) for x in input_list]):
        return None
    if type.lower() == 'polyphen':
        return simplify_polyphen(input_list)
    elif type.lower() == 'sift':
        return simplify_sift(input_list)
    raise Exception('Type is not polyphen or sift')