import pandas as pd
import argparse
from cStringIO import StringIO
import numpy as np
import yaml

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

def parse_EFF(value):
    # this is the SNPEFF field, parse it appropriately
    #NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Gtt/Att|V5I|293|HNRNPCL1||CODING|NM_001013631.1|2|1),
    #MODERATE|MISSENSE|cGc/cCc|R1113P|1159|INPP5D||CODING|NM_005541.3|25|1|WARNING_TRANSCRIPT_INCOMPLETE
    EFF_LIST = []
    for effect in value.split(","):
        EFF = {}    
        EFF["e"], t = effect.split("(",1)
        try:
            # no optional warning field
            _, EFF["f"], EFF["cc"], EFF["aa"], _, EFF["g"], _, _, EFF["tx"], EFF["r"], _ = t.split("|")
        except:
            _, EFF["f"], EFF["cc"], EFF["aa"], _, EFF["g"], _, _, EFF["tx"], EFF["r"], _, EFF["err"] = t[:-1].split("|") #-1 removes trailing ")"
            # clear out any empty fields!
        EFF_LIST.append({k:v for k,v in EFF.iteritems() if v is not ''})
    return EFF_LIST


def parse_annotations(info_field, config_df):
    field_list = info_field.split(";")
    out = {}
    for field in field_list:
        try:
            key,value = field.split("=")
        except:
            #print "COULD NOT PARSE INFO FIELD", field
            continue
        for ix, key_config in config_df[config_df["vcf-name"] == "INFO.%s" % key].iterrows():
            out_column_name = key_config["col"]
            if out_column_name in ["EFF"]:
                out_value = parse_EFF(value)
            elif type(key_config["formatter"]) == str:
                out_value = eval(key_config["formatter"])(value)
            else:
                out_value = key_config["formatter"](value)
            out[out_column_name] = out_value
    return out

def load_config(config_file):
    y = yaml.load(open(args.config))
    for ix, col in enumerate(y["output"]):
        if type(col) != dict:
            y["output"][ix] = {"col":col, "formatter":str,"vcf-name":col}
        else:
            col_name = col.keys()[0]
            d = {"col": col_name}
            if col_name in FIELDSETS:
                d.update(col[col_name])
            else:
                d["vcf-name"] = y["output"][ix][col_name].get("vcf-name",col_name)
                if "formatter" not in y["output"][ix][col_name]:
                    d["formatter"] = str
                elif y["output"][ix][col_name]["formatter"].startswith("lambda"):
                    d["formatter"] = eval(y["output"][ix][col_name]["formatter"])
                else:
                    d["formatter"] = y["output"][ix][col_name]["formatter"]
            y["output"][ix] = d
    return y

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("in_vcf", help="Input VCF file to process")
    parser.add_argument("out_tab", help="Output tab-delimited file")
    parser.add_argument("--config", required=False, default=None)
    parser.add_argument("--max-num-effects", required=False, default=3, type=int, help="Maximum number of effects to output per variant")
    #parser.add_argument("--panda-df", required=False, help="Optional pickled ouput panda dataframe")
    args = parser.parse_args()

    config = load_config(args.config)

    vcf,_,_ = read_vcf(args.in_vcf)

    MAX_EFF = args.max_num_effects
    FIELDSETS = ["EFF"]

    out = []
    for ix, row in vcf.iterrows():
        info = parse_annotations(row["INFO"], config)
        # TODO, start here with refactor
        for s_ix, s in enumerate(info["SAMPLE"]):
            d = dict(row).copy()
            d["sample_id"] = s
            d.update(info)
            d["study"] = info["STUDY"][s_ix]
            if "VALIDATION" in info:
                d["VALIDATION"] = info["VALIDATION"][s_ix]
            del d["SAMPLE"]
            del d["INFO"]
            if "EFF" in d:
                eff_levels = np.array([EFF_LEVELS[eff["e"]] for eff in d["EFF"]])
                eff_argsort = np.argsort(eff_levels)[::-1]
                eff_sorted = np.array(d["EFF"])[eff_argsort][0:MAX_EFF]
                already_added_changes = []
                i = 1
                for e_ix, e in enumerate(eff_sorted):
                    if ("aa" not in e) or (e["aa"] not in already_added_changes):
                        d["EFF_%d_gene" % i] = e["g"]
                        d["EFF_%d_effect" % i] = e["e"]
                        d["EFF_%d_group" % i] = e.get("f",None)
                        d["EFF_%d_exon" % i] = e.get("r", None)
                        d["EFF_%d_AA" % i] = e.get("aa",None)
                        d["EFF_%d_transcript" % i] = e.get("tx", None)
                        already_added_changes.append(e.get("aa",None))
                        i += 1
                del d["EFF"]
            del d["FORMAT"]
            del d["QUAL"]
            del d["DATA"]
            del d["FILTER"]        
            out.append(d)

    eff_cols = ["EFF_%d_gene","EFF_%d_effect","EFF_%d_group","EFF_%d_AA","EFF_%d_exon","EFF_%d_transcript"]
    eff_out_cols = []
    for i in range(1, MAX_EFF+1):
        eff_out_cols.extend([s % i for s in eff_cols])

    out_cols =  ["1000Genomes_Total_Count", "ESP_AA_Allele_Fraction", "ESP_EA_Allele_Fraction", 'SIFT', 'PolyPhen2_hdiv', 'PolyPhen2_hvar', 'Mutation_Assessor', "LRT_score", 'Mutation_Taster', 'GERPrs', 'GERPnr', 'PhyloP_Score', "SiPhy", 'FATHMM_score', "Ancestral_Allele", "Uniprot_id"]
    print_cols = ["CHROM","POS","REF","ALT","sample_id","study","VALIDATION", "ID","dbSNPBuildID"] + eff_out_cols + out_cols
    out = pd.DataFrame(out)[print_cols]
    out.to_csv(args.out_tab, sep="\t", index=False)
