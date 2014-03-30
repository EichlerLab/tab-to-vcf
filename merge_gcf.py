import pandas as pd
import argparse
from cStringIO import StringIO
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
    try:
        info_field = info_field.rstrip(";").split(";") or []
    except:
        return {}
    out = {}
    for i in info_field:
        k,v = i.split("=")
	out[k] = v
    return out

def parse_config(config, args):
    args.merge_keys = config.get("merge").get("keys") or args.merge_keys
    args.ignore_fields = config.get("merge").get("ignore") or args.ignore_fields
    args.fuzzy = config.get("merge").get("fuzzy") or args.fuzzy
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

def apply_rule_multi(rows, rule):
    rule_name = rule.keys()[0]
    rule = rule[rule_name]
    short_rule_name = rule_name.lstrip("INFO.")
    for r in rule["order"]:
        rows["ix"] = map(lambda x: compare_rule_kw(x[1][short_rule_name], rule[r]), rows.iterrows())
        ixs = rows["ix"].values
        if r == "accept":
            if set(ixs) == set([None]):
                # accept kw not found in any rows, continue on
                continue
            else:
                # some rows contained the accept keyword, return that row
                return rows.sort("ix").head(1), True
        elif r == "reject":
            if len(set(ixs) - set([None])) > 0:
                # reject kw triggered in any row
                # do not merge and return all
                return rows, True
            else:
                # reject kw not found, continue on
                continue
        else:
            Exception("unknown rule!")
    # no rules applied, return all
    return rows, False

def diff(cols, rows):
    res = []
    for c in cols:
        if len(set(rows[c].values)) > 1:
            res.append(c)
    return res

DEFAULT_VCF_COLS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'DATA']

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge GCF files")
    parser.add_argument('vcf1', metavar='input', help="first GCF file")
    parser.add_argument('vcf2', nargs='+', metavar='input', help="additional GCF files to merge")
    parser.add_argument("outfile")
    parser.add_argument("--config", required=False, default=None)
    parser.add_argument("--silent", required=False, default=False)

    # Parse args and optional YAML config file
    args = parser.parse_args()
    rules = None
    if args.config:
        args, rules = parse_config(yaml.load(open(args.config,'r')), args)

    rule_list = map(lambda x: x.keys()[0], rules)
    print args
    # Load each vcf
    all_inputs = [args.vcf1] + args.vcf2
    print all_inputs
    all_vcf = []
    all_columns = []

    vcf = pd.DataFrame()
    
    for i in all_inputs:
        # read file
        vcf_in, _, cols = read_vcf(i)
        # normalize to standard columns
        for c in DEFAULT_VCF_COLS:
            if c not in cols:
                vcf_in[c] = "."
        # save to list
        vcf_in = vcf_in.sort(["CHROM", "POS"])
        all_vcf.append(vcf_in)
        all_columns.append(cols)
        vcf = vcf.append(vcf_in)
    
    # resort and re-create index

    vcf = vcf.sort(["CHROM", "POS"])
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
        merge_keys = [x if x != "POS" else fuzzy_pos(vcf.POS, args.fuzzy) for x in args.merge_keys]
        grouped_vcf = vcf.groupby(merge_keys)
    else:
        grouped_vcf = vcf.groupby(args.merge_keys)

    # Ignore the INFO column proper (since relevant keys have already been copied into new columns)
    ignore_cols = ["INFO"]
    ignore_cols.extend(args.ignore_fields)
    cols = filter(lambda x: x not in ignore_cols, vcf.columns)
    out = []
    to_resolve = []
    for cnt, g in grouped_vcf:
        if len(g) == 1:
            # Simple case-- nothing to resolve (unique record)
            out.append(g.ix[g.index[0]])
        else:
            # get list of differing columns
            diff_cols = diff(cols, g[cols])
            # convert the column names to INFO.<col> where needed
            resolve_fields = [tag if tag not in args.info_fields else "INFO.%s" % tag for tag in diff_cols]
            rule_success = False
            
            if len(resolve_fields) > 0:
                # Apply rules here if available to resolve if possible
                for field in resolve_fields:
                    if field in rule_list:
                        # a rule exists for this field
                        rule = rules[rule_list.index(field)]
                        # apply the rule
                        c, s = apply_rule_multi(g, rule=rule)
                        if s == False:
                            continue
                        else:
                            if len(c) == 1:
                                # successfully resolved fields
                                out.append(c.ix[c.index[0]])
                                rule_success = True
                                resolved = True
                                break
                            elif len(c) >= 2:
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
                    print g
                    print ", ".join(resolve_fields)
            else:
                # no differences found in relevant fields, append the first row by default
                out.append(g.ix[g.index[0]])
    #import IPython; IPython.embed()
    if len(to_resolve) > 0:
        pd.concat(to_resolve)[DEFAULT_VCF_COLS] \
            .sort(["CHROM", "POS"]) \
            .to_csv(args.outfile + ".unresolved", sep="\t", index=False)
        print "Wrote %d unresolved records to %s" % (len(to_resolve), args.outfile + ".unresolved")
    if len(out) > 0:
        out = pd.DataFrame(out)
        out.POS = out.POS.astype(int)
        out[DEFAULT_VCF_COLS] \
            .sort(["CHROM", "POS"]) \
            .to_csv(args.outfile, sep="\t", index=False)
        print "Wrote %d merged records to %s" % (len(out), args.outfile)
