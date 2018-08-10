# Vcf_coding_challenge
Coding challenge to extract info and re-annotate a vcf file

## Challenge Text

We frequently encounter situations where we need a custom solution for a particular
investigator’s research project. Either the tool doesn’t exist or an existing tool doesn’t
quite do what we need. Thus, we’d like to see how you approach a coding problem.

**Problem:** An investigator has a VCF that they would like annotated in a specific way.

**Coding challenge:** We will provide you with the VCF file from the investigator, and we
would like you to create a tool to output a table annotating each variant in the file. Every
variant (if multiallelic, decompose into individual mutations) should be annotated with
the following information derived from the VCF and querying the ExAC database via the
API (documentation can be found at http://exac.hms.harvard.edu/):

1. Variant type (e.g. insertion, deletion, etc.). 

2. Variant effect (e.g. missense,synonymous, etc.).

*Note: If multiple variant types exist in the ExAC database, annotate with the most
deleterious possibility.*

3. Read depth at the site of variation.

4. Number of reads supporting the variant.

5. Percentage of reads supporting the variant versus those supporting reference reads.

6. Allele frequency of variant

7. (Optional) Any other information from ExAC that you feel might be relevant.

## 7-30-2018_coding_challenge_vcf_filter_v1.sh Overview

This script will take the input vcf file coding_challenge_final.vcf and extract the following info from
both the original vcf file and from a ExAC API database query (http://exac.hms.harvard.edu/):

1.	Variant type (e.g. insertion, deletion, etc.).
2.	Variant effect (e.g. missense, synonymous, etc.).
3. 	Read depth at the site of variation.
4. 	Number of reads supporting the variant.
5.	Percentage of reads supporting the variant versus those supporting reference reads.
6. 	Allele frequency of variant
7. 	Gene ID
8.	Gene Symbol

**These details are outputed as a tsv and re-annotated vcf file**

## Required Files

The two files coding_challenge_final_additional_INFO_part.vcf and Allele_consequence_severity_rank.csv are
required for full function of the script.  **Both files must be copied to the working directory**

The file coding_challenge_final_additional_INFO_part.vcf includes vcf INFOlines to update the output vcf file.

The file Allele_consequence_severity_rank.csv is a list of allele consequences ranked according to results
found in Kircher, Martin, et al. "A general framework for estimating the relative pathogenicity of human genetic variants." Nature genetics 46.3 (2014): 310.

## Command Line Inputs

The user will supply four arguments when running the 7-30-2018_coding_challenge_vcf_filter_v1.sh:

The first read in $1 is the working folder

The second read in $2 is the starting vcf file

The third read in $3 is the output tsv file

The forth read in $4 is the output re-annotated vcf file

## Output Files

The 7-30-2018_coding_challenge_vcf_filter_v1.sh script will output both a tsv and vcf file

### TSV Output File

Each variant from the original vcf file is represented by a single row in the table, with the exception of multiallelic loci that are split into one row for each allele.

The column names and descriptions are shown below:

**chromosome:** Human chromosome number

**location:** Base number location in chromosome

**reference_allele:** Human genome reference allele

**alternate_allele:** Sample alternate allele

**ExAC_name:** Name of allele in ExAC database query format

**variant_type:** Variant type (e.g. insertion, deletion, etc.)

**ExAC_variant_effect_most_Severe_Consequence:** Variant effect (e.g. missense, synonymous, etc.) derived from ExAC database query

**read_depth:** Read depth at the site of variation

**variant_read_support:** Number of reads supporting the variant.  Note the pair of numbers separated by and underscore
represent the reference allele count first and alternate allele second.

**percentage_read_variant_support:** Percentage of reads supporting the variant versus those supporting reference reads

**ExAC_allele_frequency:** Allele frequency of variant derived from ExAC database query

**gene_ID:** Ensembl Gene ID associated with allele genomic location

**gene_Symbol:** Official HUGO Gene Nomenclature Committee Symbol

### VCF Output File

The vcf output file includes additional INFO features

**VT:** Variant type (e.g. insertion, deletion, etc.)

**VC:** Variant consequence (e.g. missense, synonymous, etc.) derived from ExAC database query

**DPN:** Read depth at the site of variation

**VS:** Number of reads supporting the variant.  Note the pair of numbers separated by and underscore
represent the reference allele count first and alternate allele second.

**PS:** Percentage of reads supporting the variant versus those supporting reference reads

**EAF:** Allele frequency of variant derived from ExAC database query

**ID:** Ensembl Gene ID associated with allele genomic location

**GS:** Official HUGO Gene Nomenclature Committee Symbol
