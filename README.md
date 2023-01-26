# LOFTEE (Loss-Of-Function Transcript Effect Estimator)

## Loss-of-function pipeline (inspired by MacArthur et al., 2012, published in Karczewski et al., 2020).

A VEP plugin to identify LoF (loss-of-function) variation.

Currently assesses variants that are:

-   Stop-gained
-   Splice site disrupting
-   Frameshift variants

Note: the master branch does not work with GRCh38. Please use the grch38 branch.

### Filters

LOFTEE implements a set of filters to deem a LoF as "low-confidence" (LC). Variants that pass these filters are labeled as "high-confidence" (HC).

For stop-gained and frameshift variants, LOFTEE removes:

-   Variants that are near the end of the transcript (based on the 50 bp rule, with modifications as described in [Karczewski et al., 2019 supplement](https://www.biorxiv.org/content/biorxiv/early/2019/01/30/531210/DC1/embed/media-1.pdf?download=true))
-   Variants that land in an exon with non-canonical splice sites around it (i.e. intron does not start with GT and end with AG)

For splice-site variants, LOFTEE removes:

-   Variants that only affect splicing of UTRs
-   Variants that are not predicted to affect a donor site (GC -> GT)
-   Variants where MaxEntScan does not predict an effect on splicing
-   Variants that are "rescued" by nearby, in-frame splice sites (max_scan_distance determines distance from original splice site where rescue splice sites can occur; default = 15 bp)
-   Variants in small introns (min_intron_size; default = 15 bp; only relevant to older versions of Gencode)

For all variants, LOFTEE removes:

-   Variants where the purported LoF allele is the ancestral state (across primates)
-   Variants in incomplete transcripts (only relevant to older versions of Gencode)

### Flags

LOFTEE implements a series of flags in addition to the above filters. Flagged variants should be treated with caution, particularly when doing genome-wide scans of LoF variation. However, they largely relate to the properties of individual transcripts or exons, so domain knowledge of a given gene will typically outperform these flags.

For stop-gained and frameshift variants, LOFTEE flags:

-   Variants in genes with only a single exon
-   Variants in exons that do not have the evolutionary signature of a protein-coding gene based on PhyloCSF
-   Variants where no exon number is indicated (apparently because the variant overlaps an intron)

For splice-site variants, LOFTEE flags:

-   Variants in NAGNAG sites (acceptor sites rescued by in-frame acceptor site)
-   Variants that fall in an intron with a non-canonical splice site (i.e. intron does not start with GT and end with AG).

### Predictions of splice-altering variants

LOFTEE also makes predictions of other splice (OS) variants that may cause LoF by disrupting normal splicing patterns.

For variants that occur in the extended (but not essential) splice sites, LOFTEE uses logistic regression models to predict whether the splice site is significantly disrupted. 

LOFTEE also uses an SVM model to predict variants that cause LoF by creating de novo donor splice sites leading to a frameshift.

## Requirements

-   VEP
-   Perl >= 5.10.1
-   Ancestral sequence (human\_ancestor.fa[.gz|.rz])
-   Samtools (must be on path)
-   PhyloCSF database (phylocsf.sql) for conservation filters

## Usage

LOFTEE is easiest run when cloned from Github and passed to VEP using `--dir_plugins` (or move all files in the directory into `~/.vep/Plugins/`).

Basic usage:

    perl variant_effect_predictor.pl [--other options to VEP] --plugin LoF,loftee_path:/path/to/loftee --dir_plugins /path/to/loftee

Pass additional options to LOFTEE by:

    perl variant_effect_predictor.pl [--other options to VEP] --plugin LoF,loftee_path:/path/to/loftee,human_ancestor_fa:/path/to/human_ancestor.fa.gz

Options:

-   `loftee_path`

Path to loftee directory. Default is the current working directory. **Note: Your PERL5LIB should also contain this path.**

-   `min_intron_size`

Minimum intron size, below which a variant should be filtered.

-   `fast_length_calculation`

The Ensembl API can be used to calculate transcript length in two different methods: one approximate (fast; usually within 3 bp of correct length) and one perfect (slow). Default: fast.

-   `human_ancestor_fa`

Location of human\_ancestor.fa file (need associated tabix index file),
available for download here (for samtools 0.1.19 and older):
[http://www.broadinstitute.org/~konradk/loftee/human\_ancestor.fa.rz](http://www.broadinstitute.org/~konradk/loftee/human_ancestor.fa.rz)
and [http://www.broadinstitute.org/~konradk/loftee/human\_ancestor.fa.rz.fai](http://www.broadinstitute.org/~konradk/loftee/human_ancestor.fa.rz.fai).
Courtesy of Javier Herrero and the 1000 Genomes Project
(source:
[ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/)). samtools 1.x
uses bgzipped inputs for samtools faidx and downloads are available here:
[https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz](https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz),
[https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai](https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai),
[https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.gzi](https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.gzi).
If this flag is set to 'false', the ancestral allele will not be checked and filtered.

-   `conservation_file`

The required SQL database (gzip) can be downloaded [here](https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz).
Alternatively, this can be loaded into MySQL by downloading the source file [here](https://www.broadinstitute.org/~konradk/loftee/phylocsf_data.tsv.gz)
and loaded into MySQL with the schema available [here](https://www.broadinstitute.org/~konradk/loftee/phylocsf_data_schema.sql). This route requires an additional load of the GERP [base](https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/GERP_scores.final.sorted.txt.gz) and [exon](https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/GERP_scores.exons.txt.gz) files into the same database of `gerp_bases` and `gerp_exons` repsectively.
You will then need to create a \[loftee\] entry in your `~/.my.cnf` (creating one if it does not exist) that looks like:

<pre>
[loftee]
host=your_mysql_host
user=your_mysql_user
password=your_mysql_pass
database=your_mysql_db
</pre>

-   `check_complete_cds`

The Ensembl API contains a "Complete CDS" annotation that indicates that a start and stop codon has been identified for this transcript.
This flag unfortunately requires Ensembl database access, and thus, severely decreases performance and is disabled by default.

-   `get_splice_features` 

Flag indicating whether or not to write splice prediction features to LoF_info field. Default: 1.

-   `donor_disruption_cutoff` 

The minimum cutoff on DONOR_DISRUPTION_PROB (computed from logistic regression model) used to predict a DONOR_DISRUPTION LoF. Default: 0.98.

-   `acceptor_disruption_cutoff` 

The minimum cutoff on ACCEPTOR_DISRUPTION_PROB (computed from logistic regression model) used to predict a ACCEPTOR_DISRUPTION LoF. Default: 0.99.

-   `donor_disruption_mes_cutoff` 

If no conservation_file is specified, then LOFTEE cannot use the logistic regression model to compute DONOR_DISRUPTION_PROB. Instead, it will predict donor disruption using only the impact of the variant on the splice site’s MES score. In this case, donor_disruption_mes_cutoff is the minimum cutoff used to predict DONOR_DISRUPTION. Default: 6 (i.e. the variant must lower the MES score of the splice site by at least 6 to activate DONOR_DISRUPTION).

-   `acceptor_disruption_mes_cutoff` 

Ditto for variants affecting the acceptor site. Default: 7.

-   `max_scan_distance` 

The maximum distance (in bp) from the disrupted donor or acceptor splice site where LOFTEE will look for "rescue" splice sites. Default: 15.

-   `donor_rescue_cutoff` 

The minimum cutoff on RESCUE_DONOR_MES (i.e. the highest MES score out of all in-frame donor splice sites within max_scan_distance bp of the original splice site) used to activate the RESUCE_DONOR filter. Default: 8.5.

-   `acceptor_rescue_cutoff` 

The minimum cutoff on RESCUE_ACCEPTOR_MES used to activate the RESCUE_ACCEPTOR filter. Default: 8.5.

-   `exonic_denovo_only` 

If this flag is set to true, LOFTEE will only look for de novo donor splice sites occuring in the exon. Default: 1.

-   `weak_donor_cutoff` 

Minimum MES of the annotated donor site for LOFTEE to consider any potential de novo donor alternatives. This is necessary because instances of annotated sites with very low MES scores lead to the false prediction of many de novo donor-creating variants. Default: -4.

-   `max_denovo_donor_distance`

The maximum distance from the original donor splice site where LOFTEE will look for de novo donor splice sites. Default: 200.

-   `denovo_donor_cutoff`

The minimum cutoff on DE_NOVO_DONOR_PROB (computed from SVM model) used to predict a DE_NOVO_DONOR LoF. Default: 0.995.

## Output

The output is the standard VEP output, or standard VEP VCF if `--vcf` is passed to VEP.
For those unfamiliar with VEP's VCF output, the annotations are written to the CSQ attribute of the INFO field.
Here, a comma-separated list of consequences, corresponding to each transcript-(alternate)allele pair, is written with each entry as a pipe-delimited set of annotations.
With more alleles and transcripts (and especially with the `--everything` flag), this will inevitably make for some very long INFO fields that are difficult to parse by eye.

See `src/read_vep_vcf.py` for a barebones example of a parsing script, or the section below on Parsing the VEP/LoF VCF for some tips and tricks.

From VEP, a VCF line may look like:

<pre>
1       1178848 rs115005664     G       A       1000.0   PASS   AC=1;AF=0.5;AN=2;CSQ=A|ENSG00000184163|ENST00000468365|Transcript|non_coding_exon_variant&nc_transcript_variant|445|||||||-1|||,A|ENSG00000184163|ENST00000462849|Transcript|upstream_gene_variant|||||||5|-1|||,A|ENSG00000184163|ENST00000486627|Transcript|downstream_gene_variant|||||||513|-1|||,A|ENSG00000184163|ENST00000330388|Transcript|stop_gained|648|616|206|Q/*|Cag/Tag|||-1|||HC,A|ENSG00000184163|ENST00000478606|Transcript|upstream_gene_variant|||||||310|-1|||`
</pre>

This is comprehensive, but the crucial information is in the `CSQ=` part, so here we have the line split up by allele-transcript pair:

<pre>
A|ENSG00000184163|ENST00000468365|Transcript|non_coding_exon_variant&nc_transcript_variant|445|||||||-1|||,
A|ENSG00000184163|ENST00000462849|Transcript|upstream_gene_variant|||||||5|-1|||,
A|ENSG00000184163|ENST00000486627|Transcript|downstream_gene_variant|||||||513|-1|||,
A|ENSG00000184163|ENST00000330388|Transcript|stop_gained|648|616|206|Q/*|Cag/Tag|||-1|||HC,
A|ENSG00000184163|ENST00000478606|Transcript|upstream_gene_variant|||||||310|-1|||
</pre>

The overall format of a VCF is described on the [VCF Specification Page](http://samtools.github.io/hts-specs/VCFv4.1.pdf).
The key to parsing this section is in the header line added by VEP.

<pre>
##INFO=&lt;ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|LoF_flags|LoF_filter|LoF"&gt;
</pre>

This line contains the corresponding mappings to these fields after `Format:`:

<pre>
Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|LoF_flags|LoF_filter|LoF
</pre>

### Parsing the VEP/LoF VCF

One approach that simplifies extracting data from the INFO column is to convert the text representation to a dictionary, where each annotation entry is a key-value pair.
Since CSQ entries in the INFO field can contain multiple allele-transcript pairs, we will need to make a list of these dictionaries -- one dictionary per pair.

In Python, the header line line can be read with:

<pre>
vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
</pre>

To access the annotations, the VCF record can then be read with:
<pre># Split VCF line
fields = vcf_line.split('\t')

# Split INFO field by semicolons (using lookahead regular expressions due to VEP introducing semi-colons into the INFO field in some scenarios)
# This creates a dictionary with key-value pairs of info field.
info_field = dict([(x.split('=', 1)) for x in re.split(';(?=\w)', fields[7]) if x.find('=') > -1])

# For instance, info_field['AF'] would return the allele frequency of that variant.

# Pull together the VEP field names from before with CSQ attribute from INFO field (which are pipe-delimited).
annotations = [dict(zip(vep_field_names, x.split('|'))) for x in info_field['CSQ'].split(',')]
</pre>

`annotations` is now a list of dictionaries like so:

<pre>
[{'Allele': 'A',
  'Amino_acids': '',
  'CDS_position': '',
  'Codons': '',
  'Consequence': 'non_coding_exon_variant&nc_transcript_variant',
  'DISTANCE': '',
  'Existing_variation': '',
  'Feature': 'ENST00000468365',
  'Feature_type': 'Transcript',
  'Gene': 'ENSG00000184163',
  'LoF': '',
  'LoF_filter': '',
  'LoF_flags': '',
  'Protein_position': '',
  'STRAND': '-1',
  'cDNA_position': '445'},
  ...
  {'Allele': 'A',
  'Amino_acids': 'Q/*',
  'CDS_position': '616',
  'Codons': 'Cag/Tag',
  'Consequence': 'stop_gained',
  'DISTANCE': '',
  'Existing_variation': '',
  'Feature': 'ENST00000330388',
  'Feature_type': 'Transcript',
  'Gene': 'ENSG00000184163',
  'LoF': 'HC',
  'LoF_filter': '',
  'LoF_flags': '',
  'Protein_position': '206',
  'STRAND': '-1',
  'cDNA_position': '648'}, ...]
</pre>

CSQ entries in the INFO field for a given variant can now be accessed easily.
For example, to get the LoF_filter field for the first allele transcript pair use: `annotations[0]['LoF_filter']`.

LoF annotations can be extracted as:

<pre>
lof_annotations = [x for x in annotations if x['LoF'] == 'HC']
</pre>

`HC` refers to high-confidence LoF variants (i.e. does not fail any filters). `LC` denotes low-confidence, failing at least one filter, which are written to the `LoF_filter` field.

Possible values for the `LoF_filter` field are:

- END_TRUNC
The LoF variant falls in the last `filter_position` of the transcript (default = 0.05).

- INCOMPLETE_CDS
The LoF falls in a transcript whose start or stop codons are not known.

- EXON\_INTRON_UNDEF
The LoF falls in a transcript whose exon/intron boundaries are undefined in the EnsEMBL API.

- SMALL_INTRON
The LoF falls in a splice site of a small (biologically unlikely; default < 15 bp) intron.

- ANC_ALLELE
The alternate allele of the LoF reverts the sequence back to the ancestral state.

- NON_DONOR_DISRUPTING
An essential splice donor variant’s DISRUPTION_PROB fails to exceed the donor_disruption_cutoff.

- NON_ACCEPTOR_DISRUPTING
An essential splice acceptor variant’s DISRUPTION_PROB fails to exceed the acceptor_disruption_cutoff.

- RESCUE_DONOR
A splice donor-disrupting variant (either essential or extended with sufficient DONOR_DISRUPTION_PROB) is rescued by an alternative splice site (less than max_scan_distance bp away) with an MES score above donor_rescue_cutoff. The variant in question, which was formerly determined to disrupt an existing splice site, gets downgraded to an LC LoF.

- RESCUE_ACCEPTOR
Ditto for splice acceptor-disrupting variants.

- GC_TO_GT_DONOR
Essential donor splice variant creates a more canonical splice site (strengthening the site if anything, thus unlikely to disrupt splicing).

- 5UTR_SPLICE and 3UTR_SPLICE
Essential splice variant LoF occurs in the UTR of the transcript.

Possible values for the `Lof_flags` field are:

- SINGLE_EXON
The LoF falls in a single exon transcript.

- NAGNAG_SITE
The LoF is a splice variant that falls into a NAGNAG sequence, which may indicate a frame-restoring splice site.

- PHYLOCSF_WEAK
The LoF falls in an exon that does not exhibit a pattern of conservation typical of a protein-coding exon.

- PHYLOCSF_UNLIKELY_ORF
The LoF falls in an exon that exhibits a pattern of conservation typical of a protein-coding exon, but the reading frame is likely offset.

- NON\_CAN_SPLICE
The LoF is a splice variant that falls in a non-canonical splice site (not GT..AG).

Special thanks to Monkol Lek for the initial implementation of the software and developing many of these filters.

## Citing LOFTEE

To cite loftee, please use the citation: [Karczewski et al., 2020](https://www.nature.com/articles/s41586-020-2308-7).
