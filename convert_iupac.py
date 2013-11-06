import pandas as pd
from cStringIO import StringIO


IUPAC_CODES = {
    "R": {"A":"G", "G":"A"},
    "Y": {"C":"T", "T":"C"},
    "S": {"G":"C", "C":"G"},
    "W": {"A":"T", "T":"A"},
    "K": {"G":"T", "T":"G"},
    "M": {"A":"C", "C":"A"},
}

def read_vcf(vcf_filename, columns=None):
    columns = columns or ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","DATA"]
    s = StringIO()
    vcf_header_lines = ""
    with open(vcf_filename) as f:
        for line in f:
            if line.startswith('#'):
                vcf_header_lines += line
            else:
                s.write(line)
    s.seek(0)
    df = pd.read_csv(s, sep="\t",names=columns)
    return df, vcf_header_lines

def convert_iupac(row):
    if row["ALT"] not in IUPAC_CODES.keys():
        return row["ALT"]
    else:
        return IUPAC_CODES[row["ALT"]][row["REF"]]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="tab-delimited input")
    parser.add_argument("output_file", help="VCF 4.1 output")
    
    vcf, header = read_vcf(args.input_file,
        columns=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"])
    vcf["ALT"] = map(lambda x: convert_iupac(x[1]), vcf.iterrows())

    vcf.to_csv(args.output_file, sep="\t", index=False)
