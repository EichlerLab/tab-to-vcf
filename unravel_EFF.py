import pandas as pd
from cStringIO import StringIO
import numpy as np
import argparse
import copy
import operator
from collections import defaultdict
EFF_LEVELS = {"SPLICE_SITE_ACCEPTOR": 4, 
            "SPLICE_SITE_DONOR": 4, 
            "START_LOST": 4, 
            "EXON_DELETED": 4, 
            "FRAME_SHIFT": 4, 
            "STOP_GAINED": 4, 
            "STOP_LOST": 4, 
            "RARE_AMINO_ACID": 4, 
            "NON_SYNONYMOUS_CODING": 3, 
            "CODON_CHANGE": 3, 
            "CODON_INSERTION": 3, 
            "CODON_CHANGE_PLUS_CODON_INSERTION": 3, 
            "CODON_DELETION": 3, 
            "CODON_CHANGE_PLUS_CODON_DELETION": 3, 
            "UTR_5_DELETED": 3, 
            "UTR_3_DELETED": 3, 
            "SYNONYMOUS_START": 2, 
            "NON_SYNONYMOUS_START": 2, 
            "START_GAINED": 2, 
            "SYNONYMOUS_CODING": 2, 
            "SYNONYMOUS_STOP": 2, 
            "UTR_5_PRIME": 1, 
            "UTR_3_PRIME": 1, 
            "REGULATION": 1, 
            "UPSTREAM": 1, 
            "DOWNSTREAM": 1, 
            "GENE": 1, 
            "TRANSCRIPT": 1, 
            "EXON": 1, 
            "INTRON_CONSERVED": 1, 
            "INTRON": 1, 
            "INTRAGENIC": 1, 
            "INTERGENIC": 1, 
            "INTERGENIC_CONSERVED": 1, 
            "NONE": 1, 
            "CHROMOSOME": 1, 
            "CUSTOM": 1, 
            "CDS": 1}      


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

def parse_INFO(info):
    field_list = info.split(";")
    out = {}
    for field in field_list:
        try:
            key, value = field.split("=")
            flag = False
        except:
            #This is a VCF flag
            key, value = field, True
        if key.startswith("EFF"):
            out["EFF"] = parse_EFF(value, 20)
        else:
            out[key] = value
    return out

def parse_EFF(value, max_num_effects=1):
    # this is the SNPEFF field, parse it appropriately
    #NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Gtt/Att|V5I|293|HNRNPCL1||CODING|NM_001013631.1|2|1),
    #MODERATE|MISSENSE|cGc/cCc|R1113P|1159|INPP5D||CODING|NM_005541.3|25|1|WARNING_TRANSCRIPT_INCOMPLETE
    EFF_LIST = []
    for effect in value.split(","):
        EFF = {}    
        EFF["e"], t = effect.split("(",1)
        try:
            # no optional warning field
            # HIGH||-/-|-177|305|OR4F5||CODING|NM_001005484.1|1|1)
            EFF["class"], EFF["group"], EFF["cc"], EFF["aa"], _, EFF["g"], _, _, EFF["tx"], EFF["r"], _ = t.split("|", 10)
        except ValueError:
            print value
            raise
            # clear out any empty fields!
        EFF_LIST.append({k:v for k,v in EFF.iteritems() if v is not ''})
    
    eff_levels = np.array([EFF_LEVELS[eff["e"]] for eff in EFF_LIST])
    eff_argsort = np.argsort(eff_levels)[::-1]
    eff_sorted = np.array(EFF_LIST)[eff_argsort][0:max_num_effects]
    already_added_changes = []
    i = 1
    out = []
    for e_ix, e in enumerate(eff_sorted[0:max_num_effects]):
        if ("aa" not in e) or (e["aa"] not in already_added_changes):
            d = {}
            d["i"] = i
            d["gene"] = e.get("g", "")
            d["class"] = e.get("class", "")
            d["codon"] = e.get("cc", "")
            d["effect"] = e.get("e","")
            d["group"] = e.get("group","")
            d["exon"] = e.get("r", "")
            d["AA"] = e.get("aa","")
            d["transcript"] = e.get("tx", "")
            already_added_changes.append(e.get("aa",""))
            i += 1
            out.append(d)
    return out

def format_info(info_d, order="alphabetical"):
    if order == "alphabetical":
        _order = sorted(info_d.keys())
    elif type(order) == list:
        _order = order
    else:
        _order = info_d.keys()
    out = []
    for k in _order:
        v = info_d[k]
        if k == "EFF":
            # format EFF field based on dictionary
            out += ["EFF={effect}({class}|{group}|{codon}|{AA}||{gene}|||{transcript}|{exon}||)".format(**v)]
        elif type(v) == bool:
            # is a flag
            out += ["%s" % k]
        else:
            # is a k=v pair
            out += ["%s=%s" % (k,v)]
    return ";".join(out)

def write_vcf(filename, df, header_lines, columns):
    with open(filename, 'w') as out_f:
        out_f.writelines(header_lines)
        df[columns].to_csv(out_f, sep="\t", index=False)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf")
    parser.add_argument("outvcf")
    args = parser.parse_args()

    vcf, header, cols = read_vcf(args.vcf)
    vcf["INFO_d"] = map(parse_INFO, vcf.INFO.values)

    vcf["EFF"] = map(lambda x: x.get("EFF"), vcf.INFO_d.values)
    cols_to_copy = ["CHROM","POS", "ID","REF","ALT","QUAL","FILTER","INFO_d"]
    out = []
    for ix, row in vcf.iterrows():
        #isoforms = [x.get("transcript") for x in row.EFF]
        for i in row.EFF:
            d = copy.deepcopy(dict(row[cols_to_copy]))
            d["INFO_d"]["EFF"] = i
            out.append(d)

    out = pd.DataFrame(out)

    out["INFO"] = map(format_info, out.INFO_d.values)
    del out["INFO_d"]

    out = out.rename(columns={"CHROM": "#CHROM"})
    out_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    header += '##unravel_EFF.py=<Nik Krumm 2014 nkrumm@uw.edu>\n'

    write_vcf(args.outvcf, out, header, out_cols)




