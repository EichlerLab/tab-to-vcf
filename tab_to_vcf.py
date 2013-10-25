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
    ("FILTER", "FILTER"),
    ("INFO", "INFO")
)


def tab_to_vcf(input_file, output_file):
    """
    Convert tab-delimited file to VCF.

    Support for the fixed VCF fields: #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
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
                    args[1] = int(args[1])
                    args.append([])
                    record = _Record(*args)
                    vcf_writer.write_record(record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="tab-delimited input")
    parser.add_argument("output_file", help="VCF 4.1 output")
    args = parser.parse_args()

    tab_to_vcf(args.input_file, args.output_file)
