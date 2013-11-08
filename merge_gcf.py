import pandas as pd
import argparse
from cStringIO import StringIO
import difflib

def read_vcf(vcf_filename, columns=None):
    columns = None
    s = StringIO()
    vcf_header_lines = ""
    with open(vcf_filename) as f:
        for line in f:
            if line.startswith('#'):
        	if line.startswith('#'):
     	           if line.startswith('#CHROM'):
                       columns = line.lstrip("#").split()
	        vcf_header_lines += line
            else:
                s.write(line)
    s.seek(0)
    df = pd.read_csv(s, sep="\t",names=columns)
    return df, vcf_header_lines, columns

def extract_info_field(info_field):
    info_field = info_field.rstrip(";").split(";") or []
    out = {}
    for i in info_field:
        k,v = i.split("=")
	out[k] = v
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_a")
    parser.add_argument("vcf_b")
    parser.add_argument("outfile")
    parser.add_argument("--merge-keys", required=False, nargs="?", 
        default=["CHROM","POS","SAMPLES"],
        help="List of VCF column names or INFO keys used assess if rows are duplicates")
    parser.add_argument("--info-fields", required=False, nargs="+", default=None,
        help="List of key names in the VCF INFO field to consider when merging."
             "Default is to ignore all INFO keys")
    parser.add_argument("--ignore-fields", required=False, nargs="+", default=[],
        help="List of column names to ignore Default is to ignore no columns")
    
    parser.add_argument("--silent", required=False, default=False)
    args = parser.parse_args()
    
    vcf_a, _, vcf_columns_a = read_vcf(args.vcf_a)
    vcf_b, _, vcf_columns_b = read_vcf(args.vcf_b)
    vcf_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'DATA']
    for col in vcf_columns:
        if col not in vcf_columns_b:
            vcf_b[col] = "."
        if col not in vcf_columns_a:
            vcf_a[col] = "."

    # sort each VCF
    vcf_a = vcf_a.sort(["CHROM", "POS"])
    vcf_b = vcf_b.sort(["CHROM", "POS"])

    # append a-->b
    vcf = vcf_a.append(vcf_b)
    vcf = vcf.reset_index(drop=True)
    vcf_columns = vcf.columns
    # extract INFO field into parsable dictionary
    info_fields = args.info_fields
    INFO_data = map(lambda x: extract_info_field(x), vcf.INFO.values)
    # and now into columns
    for i in info_fields:
        vcf[i] = map(lambda x: x.get(i, None), INFO_data)

    # drop all exact duplicates
    vcf = vcf.drop_duplicates()

    # find duplicated rows based on Chrom, Pos, Alt and Sample
    # this could be made more or less strict 
    keys = args.merge_keys

    # Ignore the INFO column proper (since relevant keys have already been copied into new columns)
    ignore_cols = ["INFO"]
    ignore_cols.extend(args.ignore_fields)
    cols = filter(lambda x: x not in ignore_cols, vcf.columns)
    print cols
    out = []
    to_resolve = []
    unresolved_line_count = 0
    for cnt, g in vcf.groupby(keys):
        if len(g) == 1:
            out.append(g)
        elif len(g) == 2:
            unresolved_line_count += 2
            resolve_rows = ""
            if not args.silent:
                a = g[cols].values[0]
                b = g[cols].values[1]
                matcher = difflib.SequenceMatcher(lambda x: x in cols, a, b)
                for tag, i1, i2, j1, j2 in matcher.get_opcodes(): 
                    if tag != "equal":
                        resolve_rows += "%s\tRow 1:%s\n\tRow 2:%s\n" % (cols[i1], a[i1:i2][0],b[j1:j2][0])
            if resolve_rows != "":
                print "=========" * 10
                print "\t" + "\t".join(cols)
                print "Row 1:\t" + "\t".join([str(x) for x in a])
                print "Row 2:\t" + "\t".join([str(x) for x in b])
                print resolve_rows
                to_resolve.append(g)
            else:
                out.append(g[cols])
        else:
            unresolved_line_count += len(g)
            if not args.silent:
                print "=========" * 10 
                print "More than 2 duplicates found!"
                for i in range(len(g)):
                    print "Row %d:" % i + "\t".join([str(x) for x in a])
            to_resolve.append(g)
    
    if len(to_resolve) > 0:
        pd.concat(to_resolve)[vcf_columns_a].to_csv(args.outfile + ".unresolved", sep="\t", index=False)
        print "Wrote %d unresolved records (%d lines) to %s" % (len(to_resolve), unresolved_line_count, args.outfile + ".unresolved")
    if len(out) > 0:
        pd.concat(out)[vcf_columns_a].to_csv(args.outfile, sep="\t", index=False)
        print "Wrote %d merged records to %s" % (len(out), args.outfile)
