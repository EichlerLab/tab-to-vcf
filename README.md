tab-to-vcf
==========

Convert tab-delimited file into a VCF

Install
-------

    git clone https://github.com/huddlej/tab-to-vcf.git
    cd tab-to-vcf
    virtualenv vcf
    source vcf/bin/activate
    pip install -r requirements.txt
    pip install -e git+https://github.com/brentp/fastahack-python.git#egg=fastahack-python

Run
---

    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta*
    gunzip human_1kg_v37.fasta.gz
    python tab_to_vcf.py variants.tab variants.vcf human_1kg_v37.fasta
