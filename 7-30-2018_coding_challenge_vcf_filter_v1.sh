#!/bin/bash

#	This script will take the input vcf file coding_challenge_final.vcf and extract the following info from
#	both the original vcf file and from a ExAC API database query (http://exac.hms.harvard.edu/):

#	1.	Variant type (e.g. insertion, deletion, etc.).
#	2.	Variant effect (e.g. missense, synonymous, etc.).
#	3. 	Read depth at the site of variation.
#	4. 	Number of reads supporting the variant.
#	5.	Percentage of reads supporting the variant versus those supporting reference reads.
#	6. 	Allele frequency of variant
#	7. 	Gene ID
#	8.	Gene Symbol

#	These details will be outputed in a tsv and re-annotated vcf file

#	NOTE that this script dependes on two files (coding_challenge_final_additional_INFO_part.vcf and Allele_consequence_severity_rank.csv)
#	installed to the  working directory to fully function


#	The following are user supplied read in details

#	The first read in $1 is the working folder
#	The second read in $2 is the starting vcf file
#	The third read in $3 is the output tsv file
#	The forth read in $4 is the output re-annotated vcf file


#	Changing directories to user defined working directory 

cd  "$1"

#	Making folder to deposit ExAC database output JSON files

mkdir ExAC_JSON_output_files

#	Making header file specific to coding_challenge_final.vcf
#	Make split header files to add new INFO lines

cat coding_challenge_final.vcf | head -n 56 > coding_challenge_final_first_header_portion.vcf
cat coding_challenge_final.vcf | tail -n 7098 | head -n 98 > coding_challenge_final_last_header_portion.vcf

#	Add in new INFO lines and preserve original coding_challenge_final_first_header_portion.vcf

cp coding_challenge_final_first_header_portion.vcf coding_challenge_final_first_header_portion_copy.vcf
cat coding_challenge_final_additional_INFO_part.vcf >> coding_challenge_final_first_header_portion.vcf
cat coding_challenge_final_last_header_portion.vcf >> coding_challenge_final_first_header_portion.vcf
mv coding_challenge_final_first_header_portion.vcf coding_challenge_final_header.vcf
mv coding_challenge_final_first_header_portion_copy.vcf coding_challenge_final_first_header_portion.vcf

#	Making temporary vcf file with no header for use in making output tsv and re-annotated vcf


cat "$2" | grep -v "##" |grep -v "#" > temp.vcf

#	Creating temp and output files

touch temp_output.vcf
touch "$3"


#	The following is a list of variables and their purpose for this program

#	chrom: 					Chromosome number for variant entry
#	location:				Exact numeric chromosome position for variant entry
#	refAllele:				Reference allele for variant entry
#	altAllele:				Alternative allele for variant entry
#	varType:				Variant type (e.g. insertion, deletion, etc.)
#	readDepth:				Read depth at the site of variation
#	exacEntry:				ExAC database format query 
#	multiAlleleCheck:		Check variable to determine if multiple alleles are assocated with a variant entry
#	exacEntryFirstAlt:		ExAC database format query from first of two alleles for a varient entry
#	exacEntrySecondAlt:		ExAC database format query from second of two alleles for a varient entry
#	nullTest:				Check variable to determine if the ExAC database output is associated with a known allele
#	varEffect:				All listed consequences listed in ExAC database output (null if no ExAC datbase match)
#	varEffectCount:			Check varible to determine if multiple consequences are listed in the ExACX database output
#	effectList:				List of consequences listed in the ExACX database output
#	mostSevereConsequence:	Most severe consequence of variant allele based on scale from Kircher et al. 2014
#	alleleFreq:				Allele frequence listed in ExAC database output (0 if no ExAC database match)
#	geneID:					Ensembl gene ID associated with variant location (null if no ExAC datbase match)
#	geneSymbol:				Human gene symbol (null if no ExAC datbase match)


#	Reading in vcf lines one at a time to check for multiple alleles (max of two) and multiple
#	variant types present in the ExAC database

filename="temp.vcf"
filelines=`cat $filename`


IFS=$'\n'
for line in $filelines ; do
	chrom=$(echo -e "$line" | awk '{print $1}' | sed 's/chr//g')
	location=$(echo -e "$line" | awk '{print $2}')
	refAllele=$(echo -e "$line" | awk '{print $4}')
	altAllele=$(echo -e "$line" | awk '{print $5}')
	varType=$(echo -e "$line" | awk 'BEGIN{FS=";"}{print $41}' | awk 'BEGIN{FS="="}{print $2}')
	readDepth=$(echo -e "$line" | awk 'BEGIN{FS=";"}{print $8}' | awk 'BEGIN{FS="="}{print $2}')
	exacEntry=$(echo -e "$line" | awk -v var1="$chrom" 'BEGIN{OFS="-"}{print var1,$2,$4,$5}')

#	Test for multiple alleles

	multiAlleleCheck=$(echo -e "$exacEntry" | awk '{print $1}' | grep -c ",")
	echo "$multiAlleleCheck"
	if [[ "$multiAlleleCheck" -gt "0" ]]; then

#	First variant allele variables if multiple alleles present
		altAllele=$(echo -e "$altAllele" | awk -F, '{print $1}')
		exacEntryFirstAlt=$(echo -e "$exacEntry" | awk '{print $1}' | awk -F, '{print $1}')
		wget http://exac.hms.harvard.edu/rest/variant/$exacEntryFirstAlt
		wait
		nullTest=$(cat $exacEntryFirstAlt | head -c 40 | grep -c "null")
		if [[ "$nullTest" -lt "1" ]]; then
			varEffect=$(cat $exacEntryFirstAlt | grep -o -P -m 1 '", "Consequence":.{0,300}' | head -n 1 | cut -d '"' -f5)

#	Filtering entries with multiple variant types by severity based on scale from Kircher et al. 2014

			varEffectCount=$(echo "$varEffect" | grep -c "&")
			if [[ "$varEffectCount" -gt "0" ]]; then
				effectList=$(echo -e "$varEffect" | sed 's/&/\n/g')

#	Comparing multiple variants from most the least severe
				touch effectCompare_temp.txt
				for line2 in $effectList ; do
					cat Allele_consequence_severity_rank.csv | grep "$line2" >> effectCompare_temp.txt
				done
				mostSevereConsequence=$(cat effectCompare_temp.txt | sort -k1,1n | head -n 1 | awk -F, '{print $2}')
				rm effectCompare_temp.txt
			else
				mostSevereConsequence=$(echo -e "$varEffect")
			fi	

			alleleFreq=$(cat $exacEntryFirstAlt | grep -o -P 'allele_freq.{0,30}' | cut -d '"' -f2 | sed 's/^:\ //g' | sed 's/,//g' | sed 's/\ //g')
			geneID=$(cat $exacEntryFirstAlt | head -c 100 | cut -d '"' -f6)
			geneSymbol=$(cat $exacEntryFirstAlt | head -c 100 | cut -d '"' -f10)
		else
			varEffect="null"
			mostSevereConsequence=$(echo -e "$varEffect")
			alleleFreq="0"
			geneID="null"
			geneSymbol="null"
		fi
		readSupport_ref_and_var=$(echo -e "$line" | awk 'BEGIN{FS=":"}{print $12}' | awk -F, 'BEGIN{OFS="_"}{print $1,$2}')
		variantFreq=$(echo -e "$line" | awk 'BEGIN{FS=":"}{print $12}' | awk -F, 'BEGIN{OFS=","}{sampleDepth= $1 + $2 + $3} {print ($2/sampleDepth)*100} ')

#	Writing first variant allele (if present) to tsv and vcf output files

		echo -e "$chrom\t$location\t$refAllele\t$altAllele\t$exacEntryFirstAlt\t$varType\t$mostSevereConsequence\t$readDepth\t$readSupport_ref_and_var\t$variantFreq\t$alleleFreq\t$geneID\t$geneSymbol" >> "$3"
		vcfLine=$(echo -e "$line" | awk -v var1="$varType" -v var2="$mostSevereConsequence" -v var3="$readDepth" -v var4="$readSupport_ref_and_var" -v var5="$variantFreq" -v var6="$alleleFreq" -v var7="$geneID" -v var8="$geneSymbol"  'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8";VT="var1";VC="var2";DPN="var3";VS="var4";PS="var5";EAF="var6";ID="var7";GS="var8,$9,$10}')
		echo -e "$vcfLine" >> temp_output.vcf
		mv $exacEntryFirstAlt ./ExAC_JSON_output_files

#	Second variant allele variables if multiple alleles present

		altAllele=$(echo -e "$line" | awk '{print $5}' | awk -F, '{print $2}')
		exacEntrySecondAlt=$(echo -e "$exacEntry" | sed 's/-*,/-/g' | awk -F- 'BEGIN{OFS="-"}{print $1,$2,$3,$5}')
		wget http://exac.hms.harvard.edu/rest/variant/$exacEntrySecondAlt
		wait
		nullTest=$(cat $exacEntrySecondAlt | head -c 40 | grep -c "null")
		if [[ "$nullTest" -lt "1" ]]; then
			varEffect=$(cat $exacEntrySecondAlt | grep -o -P -m 1 '", "Consequence":.{0,300}' | head -n 1 | cut -d '"' -f5)

#	Filtering entries with multiple variant types by severity based on scale from Kircher et al. 2014

			varEffectCount=$(echo "$varEffect" | grep -c "&")
			if [[ "$varEffectCount" -gt "0" ]]; then
				effectList=$(echo -e "$varEffect" | sed 's/&/\n/g')

#	Comparing multiple variants from most the least severe
				touch effectCompare_temp.txt
				for line2 in $effectList ; do
					cat Allele_consequence_severity_rank.csv | grep "$line2" >> effectCompare_temp.txt
				done
				mostSevereConsequence=$(cat effectCompare_temp.txt | sort -k1,1n | head -n 1 | awk -F, '{print $2}')
				rm effectCompare_temp.txt
			else
				mostSevereConsequence=$(echo -e "$varEffect")
			fi	

			alleleFreq=$(cat $exacEntrySecondAlt | grep -o -P 'allele_freq.{0,30}' | cut -d '"' -f2 | sed 's/^:\ //g' | sed 's/,//g' | sed 's/\ //g')
			geneID=$(cat $exacEntrySecondAlt | head -c 100 | cut -d '"' -f6)
			geneSymbol=$(cat $exacEntrySecondAlt | head -c 100 | cut -d '"' -f10)
		else
			varEffect="null"
			mostSevereConsequence=$(echo -e "$varEffect")
			alleleFreq="0"
			geneID="null"
			geneSymbol="null"
		fi
		readSupport_ref_and_var=$(echo -e "$line" | awk 'BEGIN{FS=":"}{print $12}' | awk -F, 'BEGIN{OFS="_"}{print $1,$3}')
		variantFreq=$(echo -e "$line" | awk 'BEGIN{FS=":"}{print $12}' | awk -F, 'BEGIN{OFS=","}{sampleDepth= $1 + $2 + $3} {print ($3/sampleDepth)*100} ')

#	Writing second variant allele (if present) to tsv and vcf output files

		echo -e "$chrom\t$location\t$refAllele\t$altAllele\t$exacEntrySecondAlt\t$varType\t$mostSevereConsequence\t$readDepth\t$readSupport_ref_and_var\t$variantFreq\t$alleleFreq\t$geneID\t$geneSymbol" >> "$3"
		vcfLine=$(echo -e "$line" | awk -v var1="$varType" -v var2="$mostSevereConsequence" -v var3="$readDepth" -v var4="$readSupport_ref_and_var" -v var5="$variantFreq" -v var6="$alleleFreq" -v var7="$geneID" -v var8="$geneSymbol"  'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8";VT="var1";VC="var2";DPN="var3";VS="var4";PS="var5";EAF="var6";ID="var7";GS="var8,$9,$10}')
		echo -e "$vcfLine" >> temp_output.vcf
		mv $exacEntrySecondAlt ./ExAC_JSON_output_files

	else

#	Single allele process
		
		wget http://exac.hms.harvard.edu/rest/variant/$exacEntry
		wait
		nullTest=$(cat $exacEntry | head -c 40 | grep -c "null")
		if [[ "$nullTest" -lt "1" ]]; then
			varEffect=$(cat $exacEntry | grep -o -P -m 1 '", "Consequence":.{0,300}' | head -n 1 | cut -d '"' -f5)

#	Filtering entries with multiple variant types by severity based on scale from Kircher et al. 2014

			varEffectCount=$(echo "$varEffect" | grep -c "&")
			if [[ "$varEffectCount" -gt "0" ]]; then
				effectList=$(echo -e "$varEffect" | sed 's/&/\n/g')

#	Comparing multiple variants from most the least severe
				touch effectCompare_temp.txt
				for line2 in $effectList ; do
					cat Allele_consequence_severity_rank.csv | grep "$line2" >> effectCompare_temp.txt
				done
				mostSevereConsequence=$(cat effectCompare_temp.txt | sort -k1,1n | head -n 1 | awk -F, '{print $2}')
				rm effectCompare_temp.txt
			else
				mostSevereConsequence=$(echo -e "$varEffect")
			fi	

			alleleFreq=$(cat $exacEntry | grep -o -P 'allele_freq.{0,30}' | cut -d '"' -f2 | sed 's/^:\ //g' | sed 's/,//g' | sed 's/\ //g')
			geneID=$(cat $exacEntry | head -c 100 | cut -d '"' -f6)
			geneSymbol=$(cat $exacEntry | head -c 100 | cut -d '"' -f10)
		else
			varEffect="null"
			mostSevereConsequence=$(echo -e "$varEffect")
			alleleFreq="0"
			geneID="null"
			geneSymbol="null"
		fi
		readSupport_ref_and_var=$(echo -e "$line" | awk 'BEGIN{FS=":"}{print $12}' | awk -F, 'BEGIN{OFS="_"}{print $1,$2}')
		variantFreq=$(echo -e "$line" | awk 'BEGIN{FS=":"}{print $12}' | awk -F, 'BEGIN{OFS=","}{sampleDepth= $1 + $2 + $3} {if (length($3) == 0) {print ($2/sampleDepth)*100} }')

#	Writing (single) variant allele to tsv and vcf output files

		echo -e "$chrom\t$location\t$refAllele\t$altAllele\t$exacEntry\t$varType\t$mostSevereConsequence\t$readDepth\t$readSupport_ref_and_var\t$variantFreq\t$alleleFreq\t$geneID\t$geneSymbol" >> "$3"
		vcfLine=$(echo -e "$line" | awk -v var1="$varType" -v var2="$mostSevereConsequence" -v var3="$readDepth" -v var4="$readSupport_ref_and_var" -v var5="$variantFreq" -v var6="$alleleFreq" -v var7="$geneID" -v var8="$geneSymbol"  'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8";VT="var1";VC="var2";DPN="var3";VS="var4";PS="var5";EAF="var6";ID="var7";GS="var8,$9,$10}')
		echo -e "$vcfLine" >> temp_output.vcf
		mv $exacEntry ./ExAC_JSON_output_files
	fi
done
unset IFS

#	Adding header to tsv file

sed -i '1 i\chromosome\tlocation\treference_allele\talternate_allele\tExAC_name\tvariant_type\tExAC_variant_effect_most_Severe_Consequence\tread_depth\tvariant_read_support\tpercentage_read_variant_support\tExAC_allele_frequency\tgene_ID\tgene_Symbol' "$3"


#	Putting header file (coding_challenge_final_header.vcf) and re-annotated variants file (temp_output.vcf) together to create final output vcf
#	and removing temporary files


cat temp_output.vcf >> coding_challenge_final_header.vcf
mv coding_challenge_final_header.vcf "$4"
#rm temp.vcf 
rm coding_challenge_final_first_header_portion.vcf
rm coding_challenge_final_last_header_portion.vcf
