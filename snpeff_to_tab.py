import pandas as pd
import argparse
from cStringIO import StringIO
import numpy as np

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

formatters = {'AB': float,
              'AC': float,
              'AF': float,
              'AN': float,
              'Alignability': float,
              'BaseQRankSum': float,
              'DB': float,
              'DP': float,
              'Dels': float,
              'EFF': None,
              'HRun': float,
              'LowMQ': lambda x: [float(y) for y in x.split(",")],
              'MQ': float,
              'MQ0': float,
              'MQRankSum': float,
              'QD': float,
              'SB': float,
              'dbNSFP_1000Gp1_AC': float,
              'dbNSFP_1000Gp1_AFR_AC': float,
              'dbNSFP_1000Gp1_AMR_AC': float,
              'dbNSFP_1000Gp1_ASN_AC': float,
              'dbNSFP_1000Gp1_EUR_AC': float,
              'dbNSFP_29way_logOdds': float,
              'dbNSFP_29way_pi': lambda x: [float(y) for y in x.split(":")],
              'dbNSFP_Ancestral_allele': str,
              'dbNSFP_ESP6500_AA_AF': float,
              'dbNSFP_ESP6500_EA_AF': float,
              'dbNSFP_FATHMM_score': float,
              'dbNSFP_GERP++_NR': float,
              'dbNSFP_GERP++_RS': float,
              'dbNSFP_GERP++_RS': float,
              'dbNSFP_GERP++_NR': float,
              'dbNSFP_LRT_Omega': float,
              'dbNSFP_LRT_score': float,
              'dbNSFP_MutationAssessor_score': float,
              'dbNSFP_MutationTaster_score': float,
              #'dbNSFP_Polyphen2_HDIV_score': lambda x: [float(y) for y in x.split(",")],
              #'dbNSFP_Polyphen2_HVAR_score': lambda x: [float(y) for y in x.split(",")],
              'dbNSFP_Polyphen2_HDIV_score': lambda x: max([float(y) for y in x.split(",")]),
              'dbNSFP_Polyphen2_HVAR_score': lambda x: max([float(y) for y in x.split(",")]),
              'dbNSFP_SIFT_score': float,
              'dbNSFP_SLR_test_statistic': float,
              'dbNSFP_Uniprot_id': lambda x: ", ".join(filter(lambda x: x!=".", x.split(","))),
              'dbNSFP_Ensembl_transcriptid': lambda x: x.split(","),
              'dbNSFP_phyloP': float,
              'dbNSFP_UniSNP_ids': str,
              'dbSNPBuildID':int,
              'VALIDATION': lambda x: x.split(","),
              'STUDY': lambda x: x.split(",")}

short_key_names = {'dbNSFP_1000Gp1_AC': '1000Genomes_Total_Count',
                   'dbNSFP_1000Gp1_AFR_AC': '1000Genomes_AFR_Count',
                   'dbNSFP_1000Gp1_AMR_AC': '1000Genomes_AMR_Count',
                   'dbNSFP_1000Gp1_ASN_AC': '1000Genomes_ASN_Count',
                   'dbNSFP_1000Gp1_EUR_AC': '1000Genomes_EUR_Count',
                   'dbNSFP_Ancestral_allele': 'Ancestral_Allele',
                   'dbNSFP_ESP6500_AA_AF': 'ESP_AA_Allele_Fraction',
                   'dbNSFP_ESP6500_EA_AF': 'ESP_EA_Allele_Fraction',
                   'dbNSFP_FATHMM_score': 'FATHMM_score',
                   'dbNSFP_GERP++_RS': 'GERPrs',
                   'dbNSFP_GERP++_NR': 'GERPnr',
                   'dbNSFP_LRT_score': 'LRT_score',
                   'dbNSFP_MutationAssessor_score': 'Mutation_Assessor',
                   'dbNSFP_MutationTaster_score': 'Mutation_Taster',
                   'dbNSFP_Polyphen2_HDIV_score': 'PolyPhen2_hdiv',
                   'dbNSFP_Polyphen2_HVAR_score': 'PolyPhen2_hvar',
                   'dbNSFP_SIFT_score': 'SIFT',
                   'dbNSFP_SLR_test_statistic': 'SLR',
                   'dbNSFP_29way_logOdds': 'SiPhy',
                   'dbNSFP_phyloP': 'PhyloP_Score',
                   'dbNSFP_Uniprot_id':'Uniprot_id'}

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

def parse_annotations(info_field):
    field_list = info_field.split(";")
    out = {}
    for field in field_list:
        try:
            key,value = field.split("=")
        except:
            #print "COULD NOT PARSE INFO FIELD", field
            continue
        if key == "SAMPLE":
            out["SAMPLE"] = value.split(",")
        elif key == "EFF":
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
            out["EFF"] = EFF_LIST
        elif key[0:6] == "dbnsfp":
            # These are tne dbNSFP fields, process them appropriately!
            if "dbNSFP" not in out:
                out["dbNSFP"] = {}
            try:
                k = short_key_names.get(key,key)
                out["dbNSFP"][k] = formatters[key](value)
            except ValueError:
                print key, value
        elif key == "LOF":
            # this is the output from the -lof function
            # value will be (Gene | ID | num_transcripts | percent_affected)
            value = value.split(",")[0]
            _, _, total_tx, lof_percent = value.rstrip(")").split("|")
            out["LOF_tx"] = int(float(lof_percent) * int(total_tx))
            out["Num_tx"] = int(total_tx)
        else:
            if key in formatters:
                k = short_key_names.get(key,key)
                if k is not None:
                    out[k] = formatters[key](value)
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("in_vcf", help="Input VCF file to process")
    parser.add_argument("out_tab", help="Output tab-delimited file")
    parser.add_argument("--max-num-effects", required=False, default=3, type=int, help="Maximum number of effects to output per variant")
    #parser.add_argument("--panda-df", required=False, help="Optional pickled ouput panda dataframe")
    args = parser.parse_args()


    vcf,_,_ = read_vcf(args.in_vcf)

    MAX_EFF = args.max_num_effects

    out = []
    for ix, row in vcf.iterrows():
        info = parse_annotations(row["INFO"])
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
