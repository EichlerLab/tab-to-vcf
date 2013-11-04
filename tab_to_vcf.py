#!/bin/env python
"""
Convert the given tab-delimited document into a VCF 4.1 document for annotation with Seattle Seq.
"""
import argparse
import csv
import vcf
from vcf.model import _Record


TEMPLATE_VCF_FILE = "template-4.1.vcf"
VCF_TO_FIELDS = (
    ("#CHROM", "Chrom"),
    ("POS", "Pos(hg19)"),
    ("ID", "Unique id"),
    ("REF", "Ref"),
    ("ALT", "Allele"),
    ("QUAL", "QUAL"),
    ("FILTER", "FILTER")
)
POSITION_INDEX=1
ALT_INDEX=4


def tab_to_vcf(input_file, output_file):
    """
    Convert tab-delimited file to VCF.

    Support for the fixed VCF fields: #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO

    PyVCF's _Record class requires the following arguments:

    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample_indexes
    """
    with open(input_file, "r") as input_fh:
        reader = csv.DictReader(input_fh, delimiter="\t")

        with open(TEMPLATE_VCF_FILE, "r") as template_fh:
            vcf_reader = vcf.Reader(template_fh)

            with open(output_file, "w") as output_fh:
                vcf_writer = vcf.Writer(output_fh, vcf_reader, lineterminator='\n')

                for row in reader:
                    args = [row.get(tab_field, ".")
                            for vcf_field, tab_field in VCF_TO_FIELDS]

                    # Convert position to an integer.
                    args[POSITION_INDEX] = int(args[POSITION_INDEX])

                    # Convert alternate allele scalar to a list.
                    args[ALT_INDEX] = [args[ALT_INDEX]]

                    # Add empty entries for INFO, FORMAT, and sample_indexes.
                    args.extend([{}, ".", []])

                    record = _Record(*args)
                    vcf_writer.write_record(record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="tab-delimited input")
    parser.add_argument("output_file", help="VCF 4.1 output")
    parser.add_argument("reference_file", help="reference assembly for variants in a single FASTA file")
    args = parser.parse_args()

    tab_to_vcf(args.input_file, args.output_file, args.reference_file)
