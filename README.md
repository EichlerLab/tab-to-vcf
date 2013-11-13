vcf-flow
=======
### a yaml-based toolkit for managing VCF files


Vcf-flow is a set of tools designed for working with [Variant Call Files](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41). Unique to vcf-flow is the ability to configure each tool using a [YAML](http://en.wikipedia.org/wiki/YAML) configuration file, which makes repeatable workflows or pipelines possible.  Currently there are several available tools:

There are three basic file types involved in vcf-flow:

 *  .tab/.txt: Tab-delimited files. Often the first input or final output format of a workflow.
 *  .VCF files: Each line corresponds to a unique Chromosome, Position and alternate allele. Additional info is represented in the INFO column.
 *  .GCF files: Similar to a VCF file, but instead of "Variants" on each line, a GCF represents a single "Genotype" on each line-- Making each line unique in Chrom, Pos, Alternate *and* Sample

Several programs are available for vcf-flow:

 *  *tab_to_gcf.py* Converts a tab-delimited file into a standardized GCF (see below) formatted file, and convert GATK-style reference/alternate alleles to VCF format
 *  *convert_iupac.py* Convert IUPAC nucleotide names in a gcf/vcf to standard A/T/C/Gs only (biallelic genotypes only)
 *  *merge_gcf.py* Uses a YAML configuration file to define a set of rules to merge two GCF files. Reports errors to a separate file.
 *  *gcf_to_vcf.py* Convert a GCF file into a VCF file, grouping relevant attributes from each GCF record.
 *  *vcf_to_tab.py* Converts a VCF file into a tab-delimited file based on a YAML config file

### Install

    git clone https://github.com/EichlerLab/vcf-flow.git
    cd vcf-flow
    virtualenv vcf
    source vcf/bin/activate
    pip install -r requirements.txt
    pip install -e git+https://github.com/brentp/fastahack-python.git#egg=fastahack-python


## Sample Workflow(s)
### Convert tab-delimited file to GCF:

1. Download reference genome

        wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta*
        gunzip human_1kg_v37.fasta.gz

2.  Run tab_to_gcf.py

    ```shell
    python tab_to_gcf.py variants.tab variants.gcf human_1kg_v37.fasta \
         --convert-iupac
        --info-fields Sample:SAMPLE,Study:STUDY,Validation_Status:VALIDATION
    ```

    the `--info-fields` specifies which tab-delimited fields int the input variants.tab should be converted to  `KEY=VALUE` pairs in the INFO field of the VCF. The `--info-fields` are given as a comma-separated mapping of `input:output`, where `input` refers to the column name in the input file and `output` referes to the `KEY` of the INFO field.

### Merge two GCF files (with optional fuzziness parameter)

1. Set up merge.yaml configure file (See comments below for meaning of options)

    ```yaml
    merge:                                                # This first section describes which parts of the input GCFs to check for merge conflicts
        keys: ["CHROM","POS","SAMPLE"]                    # Merge the GCF on these fields (default)
        ignore: ["QUAL", "DATA", "FORMAT"]                # Do not consider these fields when merging

    rules:                                                # the "rules" section specifies which GCF record to put into final VCF file
                                                          # rules can specifiy a "reject" pattern and an "accept" pattern
                                                          # depending on the "order" of and the order of keywords, GCF records with matching
                                                          # merge-keys are either merged into a single VCF record or to placed into the .unresolved file
        - REF:                                            # First rule (Rules are given in order) operates on the REF field
            order: [reject]                               # Always reject...
            reject: "*"                                   # Any differences in the REF field
        - INFO.SAMPLE:                                    # Second rule operates on the SAMPLE key of the INFO field
            order: [reject]                               # Always reject...
            reject: "*"                                   # all differences.
        - INFO.VALIDATION:                                # Third rule (if REF and SAMPLE are not differing)
            order: [reject, accept]                       # Apply "reject" rules first, then "accept"
            reject: [invalid, failed, fail]               # Do any of these keywords match the value of the VALIDATION key? If so, reject merge
            accept: [valid, "denovo"]                     # Otherwise, does VALIDATION match any of these keys? If so, accept that GCF record
    ```

    Some other notes:
     * For each set of merged (i.e., duplicate) GCF records, the parser will attempt to apply rules in the order of the config file. Once either a "reject" or "accept" rule is triggered, the matching GCF record is merged into the VCF file (for accept) or the `.unresolved` file (for reject)
     * If no suitable rules are found, the GCF records are appended to the unmerged file
     * If more than two records are found for a set of merge-keys, no rules are evaluated and all records are appended to the unmerged file

2. Run merge_gcf.py

        python merge_gcf.py \
            --config merge.yaml \
            [--fuzzy 5 ] \
            input_A.gcf \
            input_B.gcf \
            merged.gcf > merge.log
    
    Here we are taking to input GCF files (input_A.gcf, input_B.gcf) which may have overlapping records based on CHROM, POS and SAMPLE. An optional `--fuzzy` parameter specifies if we should tolerate up to `n` bases in the POS field (useful for finding off-by-one or poorly mapped duplicated indel calls or complex variants). 

    Output files:

    * merged.gcf: Merged output. Contains only a single record for each CHROM/POS/SAMPLE. The 
    * merged.gcf.unresolved: Unresolved merges (GCF records which triggered a reject rule or had more than 2 GCF records for a (fuzzy) position)
    * merge.log: Human readable list of merge conflicts

### Convert a GCF file to a VCF file:

    python tab-to-vcf/gcf_to_vcf.py merged.gcf merged.vcf


### Create a tab-delimited file from a VCF file:

1. Set up an `output.yaml` configuration file, which specifies how to map VCF fields into a tab-delimited format

    ```yaml
    output:
        - CHROM                                  # each entry corresponds to the name of the OUTPUT column
        - POS                                    # The default is that VCF and OUTPUT column names are identical
        - REF                                    # Here we simply are copying the initial VCF fields into our output
        - ALT
        - SAMPLE:                                # If you wish to export a KEY=VALUE field from the INFO field, 
            vcf-name: INFO.SAMPLE                # use the vcf-name option and the "INFO.<KEY>" scheme
            ungroup: 'lambda x: x.split(",")'    # This will "ungroup" or "unravel" the VCF based on this key 
                                                 # and will split the comma-separated list of SAMPLES
        - STUDY:
            vcf-name: INFO.STUDY
            ungroup: 'lambda x: x.split(",")'    # Again, split the comma-separated list in INFO.STUDY into new output lines
        - VALIDATION:                   
            vcf-name: INFO.VALIDATION
            ungroup: 'lambda x: x.split(",")'
        - ID
        - dbSNPBuildID:                          # Here we are renaming the "dbSNPBuildID" key in the VCF INFO field to
            vcf-name: INFO.dbSNPBuildID          # simply "dbSNPBuildID" as an output column
        - EFF:
            vcf-name: INFO.EFF
            formatter: formatters.EFF            # Make use of the built-in formatter for SnpEff fields
            options:                             # these formatters are in the FormatterManager class
                max_num_effects: 1               # output only the most deleterious effect
        - LOF:
            vcf-name: INFO.LOF
            formatter: 'lambda x: x.split("|")[-1].rstrip(")")'  # Use a custom formatter to process the LOF key in the INFO field
        - 1000Genomes_Total_Count:
            vcf-name: INFO.dbNSFP_1000Gp1_AC
            formatter: int
        - ESP_AA_Allele_Fraction:
            vcf-name: INFO.dbNSFP_ESP6500_AA_AF
            formatter: float
    ```

2. Run the vcf_to_tab.py program:
    
        python vcf_to_tab.py --config output.yaml merged.vcf output.tab
    
