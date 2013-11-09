import pandas as pd
import argparse
from cStringIO import StringIO
import difflib
import yaml

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

def parse_config(config, args):
    args.merge_keys = config.get("merge").get("keys") or args.merge_keys
    args.ignore_fields = config.get("merge").get("ignore") or args.ignore_fields
    args.info_fields = [] #map(lambda x: x.lstrip("INFO."), filter(lambda x: x.startswith("INFO"), config.get("rules").keys()))
    for rule in config.get("rules"):
        if rule.keys()[0].startswith("INFO."):
            args.info_fields.append(rule.keys()[0].lstrip("INFO."))
    rules = config.get("rules", None)
    return args, rules

def fuzzy_pos(positions, fuzzy=2):
    last_pos = 0
    key = 0
    out = []
    for pos in positions:
        if pos - last_pos > fuzzy:
            key += 1
        out.append(key)
        last_pos = pos
    return out

def compare_rule_kw(kw, kw_list):
    if kw_list == "*":
        return 0
    else:
        try:
            return kw_list.index(kw)
        except ValueError:
            return None

def apply_rule(a, b, rule):
    rule_name = rule.keys()[0]
    rule = rule[rule_name]
    short_rule_name = rule_name.lstrip("INFO.")
    for r in rule["order"]:
        a_ix = compare_rule_kw(a[short_rule_name], rule[r])
        b_ix = compare_rule_kw(b[short_rule_name], rule[r])
        if r == "accept":
            if a_ix is None and b_ix is None:
                continue    
            elif a_ix >= b_ix:
                return [a], True
            elif b_ix > a_ix:
                return [b], True
            else:
                # accept kw not found, continue on
                continue
        elif r == "reject":
            if (a_ix is not None) or (b_ix is not None):
                # reject kw triggered in either row
                # do not merge and return both
                return [a, b], True
            else:
                # reject kw not found, continue on
                continue
        else:
            Exception("unknown rule!")
    # no rules applied, return both
    return [a,b], False

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_a")
    parser.add_argument("vcf_b")
    parser.add_argument("outfile")
    parser.add_argument("--config", required=False, default=None)
    parser.add_argument("--merge-keys", required=False, nargs="?", 
        default=["CHROM","POS","SAMPLES"],
        help="List of VCF column names or INFO keys used assess if rows are duplicates")
    parser.add_argument("--info-fields", required=False, nargs="+", default=None,
        help="List of key names in the VCF INFO field to consider when merging."
             "Default is to ignore all INFO keys")
    parser.add_argument("--ignore-fields", required=False, nargs="+", default=[],
        help="List of column names to ignore Default is to ignore no columns")
    parser.add_argument("--fuzzy", required=False, default=0, type=int,
        help="Consider POS values +/- this value to be equivalent")
    parser.add_argument("--silent", required=False, default=False)

    # Parse args and optional YAML config file
    args = parser.parse_args()
    rules = None
    if args.config:
        args, rules = parse_config(yaml.load(open(args.config,'r')), args)

    rule_list = map(lambda x: x.keys()[0], rules)
    
    # Load each vcf
    vcf_a, _, vcf_columns_a = read_vcf(args.vcf_a)
    vcf_b, _, vcf_columns_b = read_vcf(args.vcf_b)

    # normalize VCFs if missing columns
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

    # extract INFO field into columns
    INFO_data = map(lambda x: extract_info_field(x), vcf.INFO.values)
    for i in args.info_fields:
        vcf[i] = map(lambda x: x.get(i, None), INFO_data)

    # drop all exact duplicates
    vcf = vcf.drop_duplicates()
    
    # find duplicated rows based on Chrom, Pos, Alt and Sample
    # this could be made more or less strict
    if args.fuzzy > 0:
        args.merge_keys = [x if x != "POS" else fuzzy_pos(vcf.POS, args.fuzzy) for x in args.merge_keys]
        grouped_vcf = vcf.groupby(args.merge_keys)
    else:
        grouped_vcf = vcf.groupby(args.merge_keys)

    # Ignore the INFO column proper (since relevant keys have already been copied into new columns)
    ignore_cols = ["INFO"]
    ignore_cols.extend(args.ignore_fields)
    cols = filter(lambda x: x not in ignore_cols, vcf.columns)
    print cols
    out = []
    to_resolve = []
    for cnt, g in grouped_vcf:
        if len(g) == 1:
            # Simple case-- nothing to resolve (unique record)
            out.append(g.ix[g.index[0]])
        elif len(g) == 2:
            a = g[cols].values[0]
            b = g[cols].values[1]
            matcher = difflib.SequenceMatcher(lambda x: x in cols, a, b)
            resolve_rows = ""
            resolve_fields = []
            for tag, i1, i2, j1, j2 in matcher.get_opcodes(): 
                if tag != "equal":
                    resolve_rows += "%s\tRow 1:%s\n\tRow 2:%s\n" % (cols[i1], a[i1:i2][0],b[j1:j2][0])
                    if cols[i1] in args.info_fields:
                        resolve_fields.append("INFO." + cols[i1])
                    else:
                        resolve_fields.append(cols[i1])
            rule_success = False
            if resolve_rows != "":
                # Apply rules here if available to resolve if possible
                for field in resolve_fields:
                    if field in rule_list:
                        rule = rules[rule_list.index(field)]
                        c, s = apply_rule(g.ix[g.index[0]][cols], g.ix[g.index[1]][cols], rule=rule)
                        if s == False:
                            continue
                        else:
                            if len(c) == 1:
                                # successfully resolved fields
                                out.append(g.ix[c[0].name])
                                rule_success = True
                                resolved = True
                                break
                            elif len(c) == 2:
                                to_resolve.append(g)
                                #to_resolve.append(c[1])
                                rule_success = True
                                resolved = False
                                break
                            else:
                                raise Exception("Rule ERROR")
                if not rule_success:
                    to_resolve.append(g)
                    resolved=False
                if not resolved and not args.silent:
                    print "=========" * 10
                    print "\t" + "\t".join(cols)
                    print "Row 1:\t" + "\t".join([str(x) for x in a])
                    print "Row 2:\t" + "\t".join([str(x) for x in b])
                    print resolve_rows
            else:
                # no differences found in relevant fields, append the first row by default
                out.append(g.ix[g.index[0]])
        else:
            if not args.silent:
                print "=========" * 10 
                print "More than 2 duplicates found!"
                for i in range(len(g)):
                    print "Row %d:" % i + "\t".join([str(x) for x in g.ix[g.index[i]][cols]])
            to_resolve.append(g)
    
    if len(to_resolve) > 0:
        pd.concat(to_resolve)[vcf_columns_a].to_csv(args.outfile + ".unresolved", sep="\t", index=False)
        print "Wrote %d unresolved records to %s" % (len(to_resolve), args.outfile + ".unresolved")
    if len(out) > 0:
        pd.DataFrame(out)[vcf_columns_a].to_csv(args.outfile, sep="\t", index=False)
        print "Wrote %d merged records to %s" % (len(out), args.outfile)
