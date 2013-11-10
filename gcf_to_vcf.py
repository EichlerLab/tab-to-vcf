import pandas as pd
import argparse
from cStringIO import StringIO

def read_vcf(vcf_filename, columns=None):
    columns = None
    s = StringIO()
    vcf_header_lines = ""
    with open(vcf_filename) as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    columns = line.lstrip("#").split()
                else:
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

def generate_info_field(vcf, info_keys):
    info_out = []
    for ix, row in vcf.iterrows():
        info_out.append(";".join(["%s=%s" % (k,v) for k, v in zip(info_keys, row[info_keys].values)]))
    return info_out

def group_concat(l):
    return [item for sublist in [l] for item in sublist]

def group_concat_join(l):
    return ",".join([str(item) for sublist in [l] for item in sublist])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("gcf_in", help="Input gcf file to convert")
    parser.add_argument("vcf_out", help="Path for output vcf file")
    # Not implemented
    # parser.add_argument("--create-data-columns", default=False, action="store_true")
    args = parser.parse_args()

    gcf, gcf_header, gcf_cols = read_vcf(args.gcf_in)
    print gcf.head()
    INFO_data = map(lambda x: extract_info_field(x), gcf.INFO.values)
    INFO_keys = set()
    for i in INFO_data:
        INFO_keys.update(i.keys())
    INFO_keys = list(INFO_keys)
    for k in INFO_keys:
        gcf[k] = map(lambda x: x.get(k, None), INFO_data)

    agg_dict = {k: group_concat_join for k in INFO_keys}

    vcf = gcf.groupby(["CHROM","POS","ALT","REF"], as_index=False)\
             .agg(agg_dict)
    
    vcf["INFO"] = generate_info_field(vcf, INFO_keys)
    vcf["FORMAT"] = "."
    vcf["DATA"] = "."
    vcf["ID"] = "."
    vcf["QUAL"] = "."
    vcf["FILTER"] = "."
    vcf = vcf.rename(columns={"CHROM":"#CHROM"})
    vcf = vcf[["#CHROM","POS","ID","REF","ALT","QUAL", "FILTER","INFO","FORMAT","DATA"]]
    vcf = vcf.sort(["#CHROM","POS"])
    vcf.to_csv(args.vcf_out, sep="\t", index=False)
