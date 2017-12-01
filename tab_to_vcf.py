#!/bin/env python
"""
Convert the given tab-delimited document into a VCF 4.1 document for annotation with Seattle Seq.
"""
import argparse
import csv
import os
import yaml
from fastahack import FastaHack
import vcf
from vcf.model import _Record


TEMPLATE_VCF_FILE = os.path.join(os.path.dirname(__file__), "template-4.1.vcf")

# each pair here represents a REF==>ALT mapping
# For example, IUPAC "R":
#   REF      ALT
#   A   -->  G
#   G   -->  A

IUPAC_CODES = {
    "R": {"A":"G", "G":"A"},
    "Y": {"C":"T", "T":"C"},
    "S": {"G":"C", "C":"G"},
    "W": {"A":"T", "T":"A"},
    "K": {"G":"T", "T":"G"},
    "M": {"A":"C", "C":"A"},
}


VCF_COLUMN_ORDER = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]

CHROMOSOME_INDEX = 0
POSITION_INDEX = 1
REF_INDEX = 3
ALT_INDEX = 4

def get_sequence(reference_dict, chrom, position):
    position = int(position)
    return reference_dict["%s:%s-%s" % (chrom, position, position)].upper()


def gatk_indel_to_vcf(vcf_row, reference_dict):
    """
    Convert the given indel from GATK format to VCF 4.1 standard. For example,
    the following lines should be converted from:

    2       60689253        1720    *       +G      .       .       .       .
    21      38877833        1721    *       -C      .       .       .       .
    21      47958429        1722    *       +CTGGTCT        .       .       .       .

    to:

    2       60689253        1720    A       AG      .       .       .       .
    21      38877833        1721    GC      G       .       .       .       .
    21      47958429        1722    A       ACTGGTCT        .       .       .       .

    >>> reference_dict = FastaHack("human_g1k_v37.fasta")
    >>> gatk_indel_to_vcf(['2', 60689253, '1720', '*', '+G', '.', '.', '.', '.'], reference_dict)
    ['2', 60689253, '1720', 'A', 'AG', '.', '.', '.', '.']
    >>> gatk_indel_to_vcf(['21', 38877833, '1721', '*', '-C', '.', '.', '.', '.'], reference_dict)
    ['21', 38877833, '1721', 'GC', 'G', '.', '.', '.', '.']
    >>> gatk_indel_to_vcf(['21', 47958429, '1722', '*', '+CTGGTCT', '.', '.', '.', '.'], reference_dict)
    ['21', 47958429, '1722', 'A', 'ACTGGTCT', '.', '.', '.', '.']
    """
    # Load the base at the given position.
    reference_base = get_sequence(reference_dict, vcf_row[CHROMOSOME_INDEX], vcf_row[POSITION_INDEX])

    # Create a new reference allele based on the event type (the position's base
    # for insertions, the position base plus the deleted base(s) for deletions).

    # Create a new alternate allele based on the event type (the position's base
    # plus the inserted base(s) for insertions, the position's base for
    # deletions).
    if vcf_row[ALT_INDEX].startswith("+"):
        vcf_row[REF_INDEX] = reference_base
        vcf_row[ALT_INDEX] = vcf_row[ALT_INDEX].replace("+", reference_base)
    elif vcf_row[ALT_INDEX].startswith("-"):
        vcf_row[REF_INDEX] = "%s%s" % (reference_base, vcf_row[ALT_INDEX].lstrip("-"))
        vcf_row[ALT_INDEX] = reference_base

    return vcf_row


def _convert_iupac(vcf_row):
    """
    Convert a REF/ALT genotype, where the ALT is an IUPAC code
    Does not support triallelic iupac codes. Errors print warning and return input
    >>> convert_iupac(['21', 47958429, '1722', 'T', 'W', '.', '.', '.', '.'])
    ['21', 47958429, '1722', 'T', 'A', '.', '.', '.', '.']
    """
    if vcf_row[ALT_INDEX] in IUPAC_CODES.keys():
        try:
            vcf_row[ALT_INDEX] = IUPAC_CODES[vcf_row[ALT_INDEX]][vcf_row[REF_INDEX]]
        except KeyError:
            print "WARNING: could not convert IUPAC code (triallelic, malformed or other?) for row:"
            print "\t".join(vcf_row)
        finally:
            return vcf_row
    else:
        return vcf_row

def tab_to_vcf(input_file, output_file, reference_file, columns, info_fields, convert_iupac=False):
    """
    Convert tab-delimited file to VCF.

    Support for the fixed VCF fields: #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO

    PyVCF's _Record class requires the following arguments:

    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample_indexes

    convert_iupac (bool) : When present, convert IUPAC codes to the non-reference allele.
        This is only possible for when the reference and IUPAC-determined alternates share 
        at least one allele. Tri-allelic conversion is not supported and will emit a warning.
        IUPAC codes: http://www.bioinformatics.org/sms/iupac.html
    """
    reference_dict = FastaHack(reference_file)

    with open(input_file, "r") as input_fh:
        reader = csv.DictReader(input_fh, delimiter="\t")

        with open(TEMPLATE_VCF_FILE, "r") as template_fh:
            vcf_reader = vcf.Reader(template_fh)

            with open(output_file, "w") as output_fh:
                vcf_writer = vcf.Writer(output_fh, vcf_reader, lineterminator='\n')

                for row in reader:
                    
                    args = [row.get(columns.get(f,None), ".") for f in VCF_COLUMN_ORDER]
                    # Convert position to an integer.
                    args[POSITION_INDEX] = int(args[POSITION_INDEX])

                    # Convert indels from GATK to VCF format.
                    if args[ALT_INDEX].startswith(("+", "-")) and not "/" in args[ALT_INDEX]:
                        args = gatk_indel_to_vcf(args, reference_dict)

                    # Optionally convert IUPAC code
                    if convert_iupac:
                        args = _convert_iupac(args)

                    # Convert alternate allele scalar to a list.
                    args[ALT_INDEX] = [args[ALT_INDEX]]

                    # Convert info fields
                    if info_fields:
                        INFO = {}
                        for vcf_field,tab_field in info_fields.items():
                            if tab_field in row:
                                INFO[vcf_field] = row[tab_field]
                    else:
                        INFO = {}
                    # Add empty entries for INFO, FORMAT, and sample_indexes.
                    args.extend([INFO, ".", []])

                    record = _Record(*args)
                    vcf_writer.write_record(record)

def load_config(config_file):
    y = yaml.load(open(config_file))
    return y

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="tab-delimited input")
    parser.add_argument("output_file", help="VCF 4.1 output")
    parser.add_argument("--config", help="config.yaml configuration file", required=False)
    parser.add_argument("--reference_file", help="reference assembly for variants in a single FASTA file")
    parser.add_argument("--convert-iupac", help="Convert IUPAC codes to alternate allele only", 
        required=False, default=False, action="store_true")
    parser.add_argument("--info-fields", help="input:ouput (comma separated) mapping for INFO field", 
        required=False, default=None)

    args = parser.parse_args()
    if args.config:
        config = load_config(args.config)
        args.convert_iupac = config["format"]["convert_iupac"]
        args.reference_file = config["format"]["reference"]
        COLUMNS = config["columns"]
        INFO = config["info"]
    else:
        COLUMNS = {
            "CHROM": "Chrom",
            "POS": "Pos(hg19)",
            "ID": "Unique id",
            "REF": "Ref",
            "ALT": "Allele",
            "QUAL": "QUAL",
            "FILTER": "FILTER",
        }
        if args.info_fields:
            INFO_FIELDS = {}
            mapping = args.info_fields.split(",")
            for m in mapping:
                k,v = m.strip(" ").split(":")
                INFO_FIELDS[v] = k

    tab_to_vcf(args.input_file, args.output_file, args.reference_file, columns=COLUMNS, info_fields=INFO,
        convert_iupac=args.convert_iupac)
