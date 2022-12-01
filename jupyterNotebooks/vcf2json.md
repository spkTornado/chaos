```python
'''
Description: Convert any annotated (by snpEff and snpShift)
into.json format. Will be very helpful for the integration
of vcf into other workflows.
'''
```




    '\nDescription: Convert any annotated (by snpEff and snpShift)\ninto.json format. Will be very helpful for the integration\nof vcf into other workflows.\n'




```python
import os
import sys
import numpy as np
import json
```


```python
#iVCF = '/Users/shashank.katiyar/projects/vcfAnnotation/antiADARhits_ENST00000226574.snpeff.snpsift.dbsnp.vcf'
#iVCF='/Users/shashank.katiyar/projects/vcfAnnotation/ADARhits_ENST00000226574.snpeff.snpsift.dbsnp.clinvar.vcf'
iVCF='/Users/shashank.katiyar/projects/vcfAnnotation/antiADARhits_ENST00000226574.snpeff.snpsift.dbsnp.clinvar.vcf'
odir = '/Users/shashank.katiyar/projects/vcfAnnotation/'
```


```python
def getID(idData):
    '''
    There can be multiple varrIDs after annotating with dbSNP;
    One from transcript2vcf and another from dbSNP
    This function returns the transcript2vcf generated id.
    This can be helpful later to integrate with the Benchling/Intake sheet
    '''
    return idData.split(';')[0]

def getHighestPathogenicityPrediction(inputPredictionsStates):
    '''
    Given the multiple predictions from a pathgenicity predictor,
    return the highest. Eg. 'D,D,T,T,D' --> 'D'
    Order is : 'A'|'D' > 'P'|'B'|'T'|'N' > 'U'
    '''
    
    #damagin_chars_dict = {'A','D'}
    #tolerated_chars_dict = {'P','B','T','N'}
    
    if 'D' in inputPredictionsStates or 'A' in inputPredictionsStates:
        return 'D'
    elif 'T' in inputPredictionsStates or\
         'P' in inputPredictionsStates or\
         'B' in inputPredictionsStates or\
         'N' in inputPredictionsStates:
        return 'T'
    else:
        return 'U'

def verifyVarID(iVarID, varID_alt):
    if not iVarID:    # This should not happen
        print("\nNo varID identified for the variant: ", varID_alt)
        sys.exit()
        
varDict = {}
concise_varDict = {}


with open(iVCF) as ivcf:
    for line in ivcf:
        if line.startswith('#'):
            continue
        lineData = line.replace(" ","").replace("\n","").replace("\r","").split("\t")
        
        chrom, pos = lineData[:2]
        ref, alt = lineData[3:5]
        
        varID_alt = getID(lineData[2])    #we will use c. to use an id
        
        varID = False
        
        infoData_list = lineData[7].split(";")
        
        for infoField in infoData_list:
            if infoField.startswith('ANN'):    # Info annotation by snpEff
                annoFields_list = infoField.split(',')
                
                # Because there can be multiple annotations for the same variant
                effect = []
                impact = []
                gene = []
                geneID = []
                pdot = []
                cdot = []
                
                
                for multiAnnoAnnotation in annoFields_list:
                    annoField_list = multiAnnoAnnotation.split('|')

                    effect.append(annoField_list[1])
                    impact.append(annoField_list[2])
                    gene.append(annoField_list[3])
                    geneID.append(annoField_list[4])
                    cdot.append(annoField_list[9])
                    pdot.append(annoField_list[10])
                    
                
                effect = ";".join(effect)
                impact = ";".join(impact)
                gene = ";".join(gene)
                geneID = ";".join(geneID)
                pdot = ";".join(pdot)
                cdot = ";".join(cdot)
                
                varID = cdot
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                # update the dicts with coordinates
                varDict[varID] = {'chrom': chrom,
                              'tpos': varID_alt.split("-")[-1],
                              'gpos': pos,
                              'ref': ref,
                              'alt': alt}

                concise_varDict[varID] = {'chrom': chrom,
                              'tpos': varID_alt.split("-")[-1],
                              'gpos': pos,
                              'ref': ref,
                              'alt': alt} 
                
                # initialize main categories
                varDict[varID]['pathogenicity'] = []
                varDict[varID]['conserved'] = []
                varDict[varID]['popAF'] = []

                # add to the dictionary
                varDict[varID]['effect'] =  effect
                varDict[varID]['impact'] = impact
                varDict[varID]['gene'] = gene                
                varDict[varID]['geneID'] = geneID               
                varDict[varID]['pdot'] = pdot
                varDict[varID]['cdot'] = cdot
                
                
                concise_varDict[varID]['effect'] = effect
                concise_varDict[varID]['pdot'] = pdot
                concise_varDict[varID]['cdot'] = cdot
                concise_varDict[varID]['impact'] = impact 
                               
                
            elif infoField.startswith('RS='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                varDict[varID]['rsid'] = infoField.replace("RS=",'')
            
            # Note about annotations:
            # From here onwards, 'https://varianttools.sourceforge.net/Annotation/DbNSFP' to see
            # the explanations of the annotations
            
            elif infoField.startswith('dbNSFP_GERP___NR='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                varDict[varID]['GERP_NR'] = infoField.replace("dbNSFP_GERP___NR=",'')
            
            
            elif infoField.startswith('dbNSFP_GERP___RS='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                gerpRS_raw = np.mean([float(x) for x in infoField.replace("dbNSFP_GERP___RS=",'').split(',')])
                # range is: -12.3 - 6.17
                # let's normalize it between 0-1; minMax normalization
                gerpRS = (gerpRS_raw - (-12.3)) / (6.17 - (-12.3))
                
                varDict[varID]['GERP_RS'] = gerpRS
                varDict[varID]['conserved'].append(gerpRS)
                
            elif infoField.startswith('dbNSFP_MutationTaster_pred='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                mutationTester_pred = getHighestPathogenicityPrediction(infoField.replace("dbNSFP_MutationTaster_pred=",''))
                varDict[varID]['MutationTaster_pred'] = mutationTester_pred
                varDict[varID]['pathogenicity'].append(mutationTester_pred)
                
            elif infoField.startswith('dbNSFP_PROVEAN_pred='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                priveanPred = getHighestPathogenicityPrediction(infoField.replace("dbNSFP_PROVEAN_pred=",''))
                varDict[varID]['PROVEAN_pred'] = priveanPred
                varDict[varID]['pathogenicity'].append(priveanPred)
                
            elif infoField.startswith('dbNSFP_SIFT_pred='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                siftPred = getHighestPathogenicityPrediction(infoField.replace("dbNSFP_SIFT_pred=",''))
                varDict[varID]['SIFT_pred'] = siftPred
                varDict[varID]['pathogenicity'].append(siftPred)
            
            elif infoField.startswith('dbNSFP_MetaSVM_pred='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                metaSVM_pred = getHighestPathogenicityPrediction(infoField.replace("dbNSFP_MetaSVM_pred=",''))
                varDict[varID]['MetaSVM_pred'] = metaSVM_pred
                varDict[varID]['pathogenicity'].append(metaSVM_pred)
                
            elif infoField.startswith('dbNSFP_Uniprot_acc='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                varDict[varID]['Uniprot_ids'] = infoField.replace("dbNSFP_Uniprot_acc=",'')
                concise_varDict['Uniprot_ids'] = infoField.replace("dbNSFP_Uniprot_acc=",'')
                
            elif infoField.startswith('dbNSFP_phastCons100way_vertebrate='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                # Source: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=compGeno&hgta_track=cons100way&hgta_table=phastCons100way&hgta_doSchema=describe+table+schema
                # Source: http://compgen.cshl.edu/phast/phastCons-HOWTO.html
                # Treat is as the conservation score of the nucleotide base in vertibrates
                # Higher the score, more conserved the nucleotide base is
                phastCons100 = np.mean([float(x) for x in infoField.replace("dbNSFP_phastCons100way_vertebrate=",'').split(',')])
                varDict[varID]['phastCons100way_conservationScore'] = phastCons100
                varDict[varID]['conserved'].append(phastCons100)
                
            elif infoField.startswith('dbNSFP_LRT_pred='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                lrtPred = infoField.replace("dbNSFP_LRT_pred=",'')
                varDict[varID]['LRT_pred'] = lrtPred
                varDict[varID]['pathogenicity'].append(lrtPred)
                
            elif infoField.startswith('dbNSFP_MutationAssessor_pred='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                mutationAssessorPred = infoField.replace("dbNSFP_MutationAssessor_pred=",'')
                varDict[varID]['MutationAssessor_pred'] = mutationAssessorPred
                varDict[varID]['pathogenicity'].append(mutationAssessorPred)
                
            elif infoField.startswith('dbNSFP_FATHMM_pred='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                fathmmPred = infoField.replace("dbNSFP_FATHMM_pred=",'')
                varDict[varID]['FATHMM_pred'] = fathmmPred
                varDict[varID]['pathogenicity'].append(fathmmPred)
                
            elif infoField.startswith('dbNSFP_Polyphen2_HVAR_pred='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                polyPhen2_hvarPred = infoField.replace("dbNSFP_Polyphen2_HVAR_pred=",'')
                varDict[varID]['Polyphen2_HVAR_pred'] = polyPhen2_hvarPred
                varDict[varID]['pathogenicity'].append(polyPhen2_hvarPred)
                
            elif infoField.startswith('dbNSFP_Polyphen2_HDIV_pred='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                polyPhen2_hdivPred = infoField.replace("dbNSFP_Polyphen2_HDIV_pred=",'')
                varDict[varID]['Polyphen2_HDIV_pred'] = polyPhen2_hdivPred
                varDict[varID]['pathogenicity'].append(polyPhen2_hdivPred)
                
            # Population AF
            elif infoField.startswith('dbNSFP_ESP6500_AA_AF='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                pop_ESP6500_AA_AF = float(infoField.replace("dbNSFP_ESP6500_AA_AF=",''))
                varDict[varID]['ESP6500_AA_AF'] = pop_ESP6500_AA_AF
                varDict[varID]['popAF'].append(pop_ESP6500_AA_AF)
                
            elif infoField.startswith('dbNSFP_ESP6500_EA_AF='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                pop_ESP6500_EA_AF = float(infoField.replace("dbNSFP_ESP6500_EA_AF=",''))
                varDict[varID]['ESP6500_EA_AF'] = pop_ESP6500_EA_AF
                varDict[varID]['popAF'].append(pop_ESP6500_EA_AF)
                
            elif infoField.startswith('dbNSFP_ExAC_AF='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                pop_ExAC_AF = float(infoField.replace("dbNSFP_ExAC_AF=",''))
                varDict[varID]['ExAC_AF'] = pop_ExAC_AF
                varDict[varID]['popAF'].append(pop_ExAC_AF)
                
            elif infoField.startswith('dbNSFP_1000Gp3_AF='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                pop_1000Gp3_AF = float(infoField.replace("dbNSFP_1000Gp3_AF=",''))
                varDict[varID]['1000Gp3_AF'] = pop_1000Gp3_AF
                try:
                    varDict[varID]['popAF'].append(pop_1000Gp3_AF)
                except KeyError:
                    print(varID_alt)
                    print(varDict[varID])
                    raise
                          
            # Domain effects
            elif infoField.startswith('dbNSFP_Interpro_domain='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                domainInfo = infoField.replace("dbNSFP_Interpro_domain=",'')
                if 'DNA-binding' in domainInfo:
                    varDict[varID]['Interpro_domain'] = 'DNA-binding'
                    concise_varDict[varID]['Interpro_domain'] = 'DNA-binding'
                elif  'protein-binding' in domainInfo.lower():
                    varDict[varID]['Interpro_domain'] = 'Protein-binding'
                    concise_varDict[varID]['Interpro_domain'] = 'Protein-binding'
                else:
                    varDict[varID]['Interpro_domain'] = infoField
                    concise_varDict[varID]['Interpro_domain'] = infoField
            
            # Clinvar annotations
            elif infoField.startswith('ALLELEID='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                clinvarID = int(infoField.replace("ALLELEID=",''))
                varDict[varID]['clinvarID'] = clinvarID
            
            elif infoField.startswith('CLNDISDB='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                varDict[varID]['clinvarSources'] = infoField.replace("CLNDISDB=",'')
            
            elif infoField.startswith('CLNDN='):
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                clinvar_disease = infoField.replace("CLNDN=",'')
                if 'not_' in clinvar_disease:
                    varDict[varID]['clinvarDisease'] = 'NA'
                else:
                    varDict[varID]['clinvarDisease'] = clinvar_disease
            elif infoField.startswith('CLNSIG='): 
                # we can not move forward without identifying the varid
                verifyVarID(varID, varID_alt)
                
                varDict[varID]['clinvarSignificance'] = infoField.replace("CLNSIG=",'')


            
            '''    # Disabling as I could not find literature on this
            elif infoField.startswith('dbNSFP_CADD_phred='):
                varDict[varID]['CADD_phred'] = infoField.replace("dbNSFP_CADD_phred=",'')
            ''' 
  
            
            '''    # Reserved for the future updates
            elif infoField.startswith('='): 
                varDict[varID][''] = infoField.replace("=",'')
            elif infoField.startswith('='): 
                varDict[varID][''] = infoField.replace("=",'')
            elif infoField.startswith('='): 
                varDict[varID][''] = infoField.replace("=",'')
            '''

```


```python
# Update the concise_varDict by using the varDict

#pathogenicity_keys = ['MutationTaster_pred', 'Polyphen2_HVAR_pred', 'Polyphen2_HDIV_pred', 'dbNSFP_PROVEAN_pred', 'dbNSFP_SIFT_pred', 'dbNSFP_MetaSVM_pred', 'dbNSFP_LRT_pred']
#conservedness_keys = ['GERP_RS', 'phastCons100way_conservationScore']
#popAF_keys = ['dbNSFP_ESP6500_AA_AF', ]
#print(concise_varDict)


highImpactVar_list = []
pathogenicVar_list = []
pathogenicityScore_dict = {}

for varID in varDict:
    
    # set the defaults
    concise_varDict[varID]['isSNP'] = 'False'
    
    # dnsnp
    try:
        concise_varDict[varID]['rsid'] = varDict[varID]['rsid']
        concise_varDict[varID]['isSNP'] = 'True'
    except KeyError:
        concise_varDict[varID]['rsid'] = '-'
    
    # conservedness (this is a even a word?)
    varConservedness = round(np.mean(varDict[varID]['conserved']), 2)

    if str(varConservedness) == 'nan':
        varConservedness = 'NA'
        
    # population AF
    try:
        varPopAF = max(varDict[varID]['popAF'])
    except ValueError:
        varPopAF = 0
        
    # Pathogenicity
    varPathogenicity_list = varDict[varID]['pathogenicity']
        
    dCount = 0    # Count of 'D' as Damagin variant
    tCount = 0    # Count of 'T' as Tolerated variant
    uCount = 0    # Count of 'U' as Unknown variant pathogenicity
        
    for eachEffect in varPathogenicity_list:
        if 'D' in eachEffect:
            dCount += 1
        elif 'T' in eachEffect:
            tCount += 1
        else:
            uCount += 1
    pathogenicity_finalTerm = "D%sT%sU%s" %(dCount,
                                            tCount,
                                            uCount)
    if dCount+uCount+tCount == 0:
        pathgenicity_finalScore = 0
    else:
        pathgenicity_finalScore = round((dCount+1)/float(tCount + 1), 4)
        
    # update concise_varDict
    concise_varDict[varID]['conserved'] = varConservedness
    concise_varDict[varID]['popAF'] = round(varPopAF, 4)

    if varPopAF >= .02:
        concise_varDict[varID]['isSNP'] = 'True'
    concise_varDict[varID]['pathogenicity_occurances'] = pathogenicity_finalTerm
    concise_varDict[varID]['pathogenicity_score'] = pathgenicity_finalScore
    
    # clinvar stuff; this is very important
    try:
        concise_varDict[varID]['clinvarID'] = varDict[varID]['clinvarID']
        concise_varDict[varID]['clinvarDisease'] = varDict[varID]['clinvarDisease']
        concise_varDict[varID]['clinvarSignificance'] = varDict[varID]['clinvarSignificance']
        if 'pathogenic' in varDict[varID]['clinvarSignificance']:
            pathogenicVar_list.append(varID)
    except KeyError:
        concise_varDict[varID]['clinvarID'] = 'NA'
        concise_varDict[varID]['clinvarDisease'] = 'NA'
        concise_varDict[varID]['clinvarSignificance'] = 'NA'
    
    # Update the significant variant data
    if 'high' in varDict[varID]['impact'].lower():
        highImpactVar_list.append(varID)
    
    try:
        pathogenicityScore_dict[pathgenicity_finalScore].append(varID)
    except KeyError:
        pathogenicityScore_dict[pathgenicity_finalScore] = [varID]
```


```python
# let's dump to the json
#varDict = {}
#concise_varDict = {}

varDict_opath = "%s/%s" %(odir, os.path.basename(iVCF).replace('.vcf','.full.json'))
conciseVarDict_opath = "%s/%s" %(odir, os.path.basename(iVCF).replace('.vcf','.concise.json'))

with open(varDict_opath, 'w+') as vdf:
    json.dump(varDict, vdf, indent = 4)
    
with open(conciseVarDict_opath, 'w+') as cvdf:
    json.dump(concise_varDict, cvdf, indent = 4)
```


```python
# sort the variants by their pathoginicity

pathogenicVariants = ['#CHR\tPOS\tcdot\taminoAcidAlteration\tpredictedPathogenicity\tImpact\tconservedness\tClinvarSignificance\tClinvarDisease\tdomainInfo']
for pathogenicityScore in sorted(pathogenicityScore_dict.keys(), reverse=True):
    #print("==", pathogenicityScore)
    #print(pathogenicityScore_dict[pathogenicityScore])
    
    
    variantList = pathogenicityScore_dict[pathogenicityScore]
    
    for eachVar in variantList:
        
        # DO NOT DISCARD CLINVAR REPORTED VARIANTS
        isPathogenic = False
        if 'Pathogenic' in concise_varDict[eachVar]["clinvarSignificance"]:
            isPathogenic = True

        if 'True' in concise_varDict[eachVar]["isSNP"] and (not isPathogenic):
            continue

        if pathogenicityScore < 1.0 and (not isPathogenic):
            continue

        try:
            domainInfo = concise_varDict[eachVar]["Interpro_domain"]
        except KeyError:
            domainInfo = '-'
        
        varImpact = concise_varDict[eachVar]["impact"]
        
        # filter variants of less importance
        if (not varImpact == 'HIGH') and (pathogenicityScore < 3) and (domainInfo == '-') and (not isPathogenic):
            continue
          
        pathogenicVariants.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(concise_varDict[eachVar]["chrom"],
                                 concise_varDict[eachVar]["gpos"],
                                 eachVar,
                                 concise_varDict[eachVar]["pdot"].replace("p.",""),
                                 pathogenicityScore,
                                 varImpact,
                                 concise_varDict[eachVar]["conserved"],
                                 concise_varDict[eachVar]["clinvarSignificance"],
                                 concise_varDict[eachVar]["clinvarDisease"],
                                 domainInfo))
    
```


```python
pathogenicVars_opath = "%s/%s" %(odir, os.path.basename(iVCF).replace('.vcf','.pathogenicVars.tsv'))
writePathogenicVars_FH = open(pathogenicVars_opath, 'w+') 
writePathogenicVars_FH.write("\n".join(pathogenicVariants))
writePathogenicVars_FH.close()
```


```python

```
