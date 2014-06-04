import pandas as pd
import argparse
from cStringIO import StringIO
import numpy as np
import yaml
import re

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

class FormatterManager(object):
    def __init__(self):
        super(FormatterManager, self).__init__()
        prefix = "formatters."
        self.formatters = {
            "%sEFF" % prefix: self.parse_EFF
        }
        self.columns = {
            "%sEFF" % prefix: self.cols_EFF
        }
        self.options = {}

    def get_formatter(self, name):
        if name in self.formatters:
            return self.formatters[name]
        else:
            return None
    
    def get_columns(self, name):
        if name in self.columns:
            return self.columns[name]
        else:
            return None
    
    def parse_EFF(self, value, max_num_effects=1):
        # this is the SNPEFF field, parse it appropriately
        #NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Gtt/Att|V5I|293|HNRNPCL1||CODING|NM_001013631.1|2|1),
        #MODERATE|MISSENSE|cGc/cCc|R1113P|1159|INPP5D||CODING|NM_005541.3|25|1|WARNING_TRANSCRIPT_INCOMPLETE
        #Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] 
        EFF_LIST = []
        for effect in value.split(","):
            EFF = {}    
            EFF["e"], t = effect.split("(",1)
            try:
                # no optional warning field
                _, EFF["f"], EFF["cc"], EFF["aa"], EFF["aa_len"], EFF["g"], _, _, EFF["tx"], EFF["r"], _ = t.split("|")
            except:
                _, EFF["f"], EFF["cc"], EFF["aa"], EFF["aa_len"], EFF["g"], _, _, EFF["tx"], EFF["r"], _, EFF["err"] = t[:-1].split("|") #-1 removes trailing ")"
                # clear out any empty fields!
            EFF_LIST.append({k:v for k,v in EFF.iteritems() if v is not ''})
        
        if max_num_effects == "all":
            max_num_effects = len(EFF_LIST)
        eff_levels = np.array([EFF_LEVELS[eff["e"]] for eff in EFF_LIST])
        eff_argsort = np.argsort(eff_levels)[::-1]
        eff_sorted = np.array(EFF_LIST)[eff_argsort][0:max_num_effects]
        already_added_changes = []
        i = 1
        d = {}
        for e_ix, e in enumerate(eff_sorted[0:max_num_effects]):
            if ("aa" not in e) or (e["aa"] not in already_added_changes):
                d["EFF_%d_gene" % i] = e["g"]
                d["EFF_%d_effect" % i] = e["e"]
                d["EFF_%d_group" % i] = e.get("f",None)
                d["EFF_%d_exon" % i] = e.get("r", None)
                d["EFF_%d_AA" % i] = e.get("aa",None)
                d["EFF_%d_AA_len" % i] = e.get("aa_len",None)
                d["EFF_%d_transcript" % i] = e.get("tx", None)
                #already_added_changes.append(e.get("aa",None))
                i += 1
        if (i-1) > self.options.get("max_num_effects",0):
            self.options["max_num_effects"] = (i-1)
        return d

    def cols_EFF(self, max_num_effects=1):
        if max_num_effects == "all":
            mne = self.options.get("max_num_effects", max_num_effects)
        else:
            mne = max_num_effects
        eff_cols = ["EFF_%d_gene","EFF_%d_effect","EFF_%d_group","EFF_%d_AA","EFF_%d_AA_len","EFF_%d_exon","EFF_%d_transcript"]
        eff_out_cols = []
        for i in range(1, mne+1):
            eff_out_cols.extend([s % i for s in eff_cols])
        return eff_out_cols

def parse_annotations(info_field, config_df, formatter_manager):
    field_list = info_field.split(";")
    out = {}
    for field in field_list:
        try:
            key, value = field.split("=")
            flag = False
        except:
            #This is a VCF flag
            key, value = field, None
            flag = True
        for ix, key_config in config_df[config_df["vcf-name"] == "INFO.%s" % key].iterrows():
            out_column_name = key_config["col"]
            if key_config["formatter"] in formatter_manager.formatters:
                kwargs = key_config.get("options", None)
                out_value = formatter_manager.get_formatter(key_config["formatter"])(value, **kwargs)
            elif key_config["ungroup"]:
                out_value = key_config["ungroup"](value)
            elif flag:
                out_value = key_config.get("flag", "True")
            elif type(key_config["formatter"]) == str:
                try:
                    out_value = eval(key_config["formatter"])(value)
                except:
                    print "Could not convert value", value, " in column: ", key_config["col"]
                    out_value = ""
            else:
                try:
                    out_value = key_config["formatter"](value)
                except:
                    print "Could not convert value", value, " in column: ", key_config["col"]
                    out_value = ""
            if type(out_value) == dict:
                out.update(out_value)
            else:
                out[out_column_name] = out_value
    return out

def load_config(config_file, defaults):
    y = yaml.load(open(config_file))
    for ix, col in enumerate(y["output"]):
        if type(col) != dict:
            y["output"][ix] = {"col":col, 
                               "formatter":str,
                               "vcf-name":col,
                               "ungroup": False}
        else:
            col_name = col.keys()[0]
            d = {"col": col_name}
            d["ungroup"] = y["output"][ix][col_name].get("ungroup", False)
            if type(d["ungroup"]) == str:
                d["ungroup"] = eval(d["ungroup"])
            if y["output"][ix][col_name].get("arg-name", False):
                default_val = defaults.get(y["output"][ix][col_name]["arg-name"], False)
                if default_val:
                    d.update({"default-value": default_val})
            elif col_name in ["EFF"]:
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
    return pd.DataFrame(y["output"])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("in_vcf", help="Input VCF file to process")
    parser.add_argument("out_tab", help="Output tab-delimited file")
    parser.add_argument("config", help="config.yaml file")
    parser.add_argument("--defaults", nargs="+", help="Default values for arg-name keys to add as columns in output", default=None)
    args = parser.parse_args()

    if args.defaults:
        defaults = {d.split("=")[0]: d.split("=")[1] for d in args.defaults}
    else:
        defaults = {}
    config_df = load_config(args.config, defaults)
    formatter_mgr = FormatterManager()

    vcf,_,_ = read_vcf(args.in_vcf)

    UNGROUP_KEYS = list(config_df[config_df["ungroup"] != False]["col"].values)
    if len(UNGROUP_KEYS) > 0:
        ungroup = True
    else:
        ungroup = False

    out = []
    for ix, row in vcf.iterrows():
        info = parse_annotations(row["INFO"], config_df, formatter_mgr)
        if ungroup:
            for keys in zip(*[info[k] for k in UNGROUP_KEYS]):
                d = dict(row).copy()
                d.update(info)
                d.update(zip(UNGROUP_KEYS, keys))
                del d["INFO"]
                out.append(d)
        else:
            d = dict(row).copy()
            d.update(info)
            del d["INFO"]
            out.append(d)

    out = pd.DataFrame(out)
    print_cols = []
    for ix, col in config_df.iterrows():
        if col.get("default-value", None) not in [None, np.nan]:
            out[col["col"]] = col["default-value"]
            print_cols.append(col["col"])
        elif col["formatter"] in formatter_mgr.formatters:
            kwargs = col.get("options", {})
            columns = formatter_mgr.get_columns(col["formatter"])(**kwargs)
            for c in columns:
                if c not in out:
                    out[c] = None
            print_cols.extend(columns)
        elif col["col"] not in out:
            out[col["col"]] = None
            print_cols.append(col["col"])
        else:
            print_cols.append(col["col"])

    
    out[print_cols].to_csv(args.out_tab, sep="\t", index=False)
