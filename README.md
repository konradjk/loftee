# LOFTEE (Loss-Of-Function Transcript Effect Estimator)

## Loss-of-function pipeline inspired by MacArthur et al., 2012.

A VEP plugin to identify LoF (loss-of-function) variation.

Currently assesses variants that are:

-   Stop-gained
-   Splice site disrupting
-   Frameshift variants

For stop-gained and frameshift variants, LOFTEE removes:

-   Variants that are in the last X% of the transcript (default = 5%)
-   Variants in genes with only a single exon
-   Variants that land in an exon with non-canonical splice sites around it (i.e. intron does not start with GT and end with AG)

For splice-site variants, LOFTEE removes:

-   Variants in small introns (default = 15 bp
-   Variants that fall in an intron with a non-canonical splice site (i.e. intron does not start with GT and end with AG).
-   Variants in NAGNAG sites (acceptor sites rescued by in-frame acceptor site)

For all variants, LOFTEE removes:

-   Variants where the variant is the ancestral state (shared with chimps)

## Requirements

-   VEP
-   Ancestral sequence (human_ancestor.fa.rz)

## Usage

Basic usage:

    perl variant_effect_predictor.pl [--other options to VEP] --plugin LoF

Advanced usage:

    perl variant_effect_predictor.pl [--other options to VEP] --plugin LoF,filter_position:0.05

Options:

-   `filter_position`

Position in transcript where a variant should be filtered. Default is 0.05, corresponding to last 5% of transcript.

-   `min_intron_size`

Minimum intron size, below which a variant should be filtered.

-   `fast_length_calculation`

The Ensembl API can be used to calculate transcript length in two different methods: one approximate (fast; usually within 3 bp of correct length) and one perfect (slow). Default: fast.