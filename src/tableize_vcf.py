__author__ = 'konradjk'

import argparse
import gzip
import re
import sys
from loftee_utils import *
import copy
try:
    from minimal_representation import get_minimal_representation
except ImportError, e:
    get_minimal_representation = None
    print >> sys.stderr, "WARNING: Did not find minimal_representation. Outputting raw positions."

def main(args):
    # Read parameters
    f = gzip.open(args.vcf) if args.vcf.endswith('.gz') else open(args.vcf)
    if args.output is None:
        args.output = '.table'.join(args.output.rsplit('.vcf', 1))
    if args.output == args.vcf:
        print >> sys.stderr, "VCF filename has no '.vcf' and no output file name was provided. Exiting."
        sys.exit(1)
    g = gzip.open(args.output, 'w') if args.output.endswith('.gz') else open(args.output, 'w')

    desired_info = [] if args.info is None else args.info.split(',')
    desired_vep_info = [] if args.vep_info is None else args.vep_info.split(',')
    if args.simplify:
        args.max_csq = True
        args.lof_only = True
        args.collapse_annotations = True

    missing_string = '\N' if args.mysql else 'NA'

    header = None
    vep_field_names = None
    info_from_header = {}
    started = False

    output_header = 'CHROM\tPOS\tREF\tALT\t'
    if args.include_id: output_header += 'ID\t'
    if not args.omit_filter: output_header += 'FILTER\t'

    for line in f:
        line = line.strip()

        # Reading header lines to get VEP and individual arrays
        if line.startswith('#'):
            line = line.lstrip('#')
            if line.startswith('INFO=<ID='):
                try:
                    header_metadata = dict([x.split('=', 1) for x in line.split('<')[1].split('>')[0].split(',', 3) if '=' in x])
                except Exception, e:
                    print >> sys.stderr, "Malformed header line: %s" % line
                    sys.exit(1)
                info_from_header[header_metadata['ID']] = header_metadata
            if 'ID=CSQ' in line:
                vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
                vep_info_from_header = dict(zip(vep_field_names, range(len(vep_field_names))))
            if line.startswith('CHROM'):
                header = line.split()
                header = dict(zip(header, range(len(header))))
                if args.options:
                    print >> sys.stderr, "######### OPTIONS FOR INFO #########"
                    for info in info_from_header:
                        print >> sys.stderr, '%s\t%s' % (info, info_from_header[info]['Description'])
                    if vep_field_names is not None:
                        print >> sys.stderr, "######### OPTIONS FOR VEP_INFO #########"
                        print >> sys.stderr, '\n'.join(vep_field_names)
                    sys.exit(0)
            continue

        if len(desired_vep_info) > 0:
            if vep_field_names is None:
                print >> sys.stderr, "VEP info requested, but VCF file does not have a VEP header line. Exiting."
                sys.exit(1)
            if 'ALLELE_NUM' not in vep_info_from_header:
                print >> sys.stderr, "VEP output does not have ALLELE_NUM which is required for extraction. Please re-run VEP with --allele_number. Exiting."
                sys.exit(1)
        if header is None:
            print >> sys.stderr, "VCF file does not have a header line (CHROM POS etc.). Exiting."
            sys.exit(1)

        if not started:
            # Allowing entries even if not found in the header line, with some caveats
            original_desired_info = copy.deepcopy(desired_info)
            desired_info = []
            any_missing = 0
            for info in original_desired_info:
                if info in info_from_header:
                    print >> sys.stderr, 'SUCCESS: Found %s: %s' % (info, info_from_header[info]['Description'])
                    desired_info.append(info)
                else:
                    matches = 0
                    for header_record in info_from_header:
                        if re.search('^%s$' % info, header_record):
                            matches += 1
                            desired_info.append(header_record)
                            print >> sys.stderr, 'SUCCESS: Found %s (matching %s): %s' % (header_record, info, info_from_header[header_record]['Description'])
                    if not matches:
                        print >> sys.stderr, 'WARNING: Did not find %s in header.' % info
                        any_missing += 1
                        desired_info.append(info)

            # Only allowing entries in VEP header.
            original_desired_vep_info = copy.deepcopy(desired_vep_info)
            desired_vep_info = []
            for info in original_desired_vep_info:
                if info in vep_info_from_header:
                    desired_vep_info.append(info)
                    print >> sys.stderr, 'SUCCESS: Found %s' % info
                else:
                    print >> sys.stderr, 'WARNING: Did not find %s in VEP header. Not including from here on out.' % info

            # Warnings/errors for missing data
            if any_missing: print >> sys.stderr, 'WARNING: At least one INFO line requested was not found. Continuing, but results may be off.'
            if len(desired_info) + len(desired_vep_info) == 0:
                print >> sys.stderr, 'No fields left in requested info/VEP info. Exiting.'
                sys.exit(1)
            if args.lof_only and 'LoF' not in desired_vep_info:
                print >> sys.stderr, '--lof_only was used, but no LoF tag found in VEP field. Exiting.'
                sys.exit(1)

            # Ready to go.
            output_header += '\t'.join(desired_info)
            if len(desired_vep_info) > 0: output_header += '\t' + '\t'.join(desired_vep_info)
            print >> g, output_header
            started = True

        # Pull out annotation info from INFO and ALT fields
        fields = line.split('\t')
        info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[header['INFO']])])

        # Only get VEP info if requested
        if len(desired_vep_info) > 0:
            if 'CSQ' in info_field:
                # if statement here is a fix for VEP's occasional introduction of a semi-colon into the CSQ.
                # Can be removed once that is completely fixed.
                annotations = [dict(zip(vep_field_names, x.split('|'))) for x in info_field['CSQ'].split(',') if len(vep_field_names) == len(x.split('|'))]
                if args.lof_only: annotations = [x for x in annotations if x['LoF'] == 'HC']
            else:
                annotations = []
            if args.lof_only and len(annotations) == 0: continue
        alts = fields[header['ALT']].split(',')

        # Default is split line into all alternate alleles
        if not args.preserve_multiallelic:
            for index, alt in enumerate(alts):
                if get_minimal_representation is not None:
                    new_pos, new_ref, new_alt = get_minimal_representation(fields[header['POS']], fields[header['REF']], alt)
                    output = [fields[header['CHROM']], str(new_pos), new_ref, new_alt]
                else:
                    output = [fields[header['CHROM']], fields[header['POS']], fields[header['REF']], alt]
                if args.include_id: output.append(fields[header['ID']])
                if not args.omit_filter: output.append(fields[header['FILTER']])

                # Get info and VEP info
                for info in desired_info:
                    if info in info_field:
                        if info in info_from_header and 'Number' in info_from_header[info] and info_from_header[info]['Number'] == 'A':
                            output.append(info_field[info].split(',')[index])
                        else:
                            output.append(info_field[info])
                    else:
                        output.append(missing_string)
                if len(desired_vep_info) > 0:
                    # Filter to this allele
                    this_alt_annotations = [x for x in annotations if int(x['ALLELE_NUM']) - 1 == index]
                    if args.lof_only and len(this_alt_annotations) == 0: continue
                    for info in desired_vep_info:
                        this_alt_vep_info = [x[info] for x in this_alt_annotations if x[info] != '']

                        # Process options
                        if args.max_csq and info == 'Consequence': this_alt_vep_info = [csq_max_vep(x) for x in this_alt_vep_info]
                        if args.simplify_gtex and info == 'TissueExpression':
                            # Converting from tissue1:value1&tissue2:value2 to [tissue1, tissue2]
                            this_alt_vep_info = set([y.split(':')[0] for x in this_alt_vep_info for y in x.split('&')])
                        if args.collapse_annotations:
                            this_alt_vep_info = set(this_alt_vep_info)
                            # Collapse consequence further
                            if args.max_csq and info == 'Consequence': this_alt_vep_info = [csq_max(this_alt_vep_info)]

                        annotation_output = ','.join(this_alt_vep_info)
                        if annotation_output == '': annotation_output = missing_string
                        output.append(annotation_output)
                print >> g, '\t'.join(output)
        else:
            output = [fields[header['CHROM']], fields[header['POS']], fields[header['REF']], fields[header['ALT']]]
            if args.include_id: output.append(fields[header['ID']])
            if not args.omit_filter: output.append(fields[header['FILTER']])
            for info in desired_info:
                if info in info_field:
                    output.append(info_field[info])
                else:
                    output.append(missing_string)
            if len(desired_vep_info) > 0:
                for info in desired_vep_info:
                    this_vep_info = [x[info] for x in annotations]

                    # Process options
                    if args.max_csq and info == 'Consequence': this_vep_info = [csq_max(x) for x in this_vep_info]
                    if args.simplify_gtex and info == 'TissueExpression':
                        this_vep_info = set([y.split(':')[0] for x in this_vep_info for y in x.split('&')])
                    if args.collapse_annotations:
                        this_vep_info = set(this_vep_info)
                        if args.max_csq and info == 'Consequence': this_vep_info = [csq_max(this_vep_info)]

                    annotation_output = ','.join(this_vep_info)
                    if annotation_output == '': annotation_output = missing_string
                    output.append(annotation_output)
            print >> g, '\t'.join(output)

    f.close()

if __name__ == '__main__':
    INFO = '''Parses VCF to extract data from INFO field and CSQ (from VEP) inside INFO field.
By default, splits VCF record into one allele per line and creates R/MySQL/etc readable table.'''
    parser = argparse.ArgumentParser(description=INFO)

    parser.add_argument('--vcf', '--input', '-i', help='Input VCF file; may be gzipped', required=True)
    parser.add_argument('--output', '-o', help='Output table file (default=input{-.vcf}.table[.gz]); may be gzipped')
    parser.add_argument('--info', help='Comma separated list of INFO fields to extract')
    parser.add_argument('--omit_filter', action='store_true', help='Omit FILTER field from output')
    parser.add_argument('--include_id', action='store_true', help='Include ID field in output')
    parser.add_argument('--vep_info', help='Comma separated list of CSQ sub-fields to extract')
    parser.add_argument('--lof_only', action='store_true', help='Limit output to HC LoF')
    parser.add_argument('--max_csq', action='store_true', help='Max Consequence for each annotation')
    parser.add_argument('--collapse_annotations', action='store_true', help='Collapse identical annotations')
    parser.add_argument('--simplify', action='store_true', help='Alias for --lof_only --max_csq --collapse_annotations')
    parser.add_argument('--simplify_gtex', action='store_true', help='Simplify GTEx info (only print expressed tissues, not expression values)')
    parser.add_argument('--preserve_multiallelic', '-p', action='store_true', help='Preserve multi-allelic records as one line')
    parser.add_argument('--options', action='store_true', help='Print possible info and vep_info options (from header) and exit')
    parser.add_argument('--mysql', action='store_true', help='Uses \N for missing data for easy reading into MySQL (default = NA)')
    args = parser.parse_args()
    main(args)