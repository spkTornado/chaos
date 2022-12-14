```python
'''Description: Accept any transcript and generate 'A/G' and 'G/A'
mutation VCF. Usefull to explore the ADAR targets.
'''
```




    "Description: Accept any transcript and generate 'A/G' and 'G/A'\nmutation VCF. Usefull to explore the ADAR targets.\n"




```python
import os
import sys
import gget
import pandas as pd
import numpy as np
from datetime import datetime, timezone, timedelta
```


```python
def getTimeStamp():
    '''
    Return a timestap as: MM/DD/YY, HH:MM
    '''
    timezone_offset = -8.0  # Pacific Standard Time (UTC−08:00)
    tzinfo = timezone(timedelta(hours=timezone_offset))
    return datetime.now(tzinfo).strftime("Date: %m/%d/%Y, time: %H:%M")

print(getTimeStamp())
```

    Date: 12/01/2022, time: 10:42



```python

```


```python
def getHumanReference(refGenome):
    '''
    Set the reference genome
    '''
    gget.ref(refGenome)
```


```python
def ensembleID2seq(ensembleID):
    '''
    get the sequence from ensembleID
    default human genome assemble is: GRCh38
        Parameters:
            'ensembleID'(str): Transcript ensemble id
        Returns:
            'assembly' (str): Assembly of the reference genome
            'genomic_coordinates' (list): Chromosome, start, and 'end'
            'sequence' (str): Transcript sequence
    '''

    # Set the reference
    keyword2refName_dict = {'human': 'homo_sapiens'}
    refGenome = 'human'
    #getHumanReference(keyword2refName_dict[refGenome])
    gget.ref('homo_sapiens')
    
    # Note: gget.* returns a PandasDataframe
    seqMetaData, sequence = gget.seq([ensembleID])
    assembly, chrom, start, end = seqMetaData.split(':')[1:-1]
    
    return assembly, [chrom, int(start), int(end)], sequence
```


```python
def hitByADAR(iBase):
    '''
    Check and mutate A/G
        Paramter:
            'iBase' (str): input nucleotide base
        Returns:
            '(str) or 'False': Either returns 'G' or 'False' based on the input
    '''
    fitEffect_dict = {'A':'G',
                     'a': 'g'}
    try:
        return fitEffect_dict[iBase]
    except KeyError:
        return False
```


```python
ensembleID = 'ENST00000226574.9'
#print(gget.info([ensembleID]))
```


```python
print(gget.ref('homo_sapiens'))
```

    Thu Dec  1 10:42:32 2022 INFO Fetching reference information for homo_sapiens from Ensembl release: 108.


    {'homo_sapiens': {'transcriptome_cdna': {'ftp': 'http://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz', 'ensembl_release': 108, 'release_date': '2022-10-04', 'release_time': '20:17', 'bytes': '74M'}, 'genome_dna': {'ftp': 'http://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz', 'ensembl_release': 108, 'release_date': '2022-10-04', 'release_time': '18:37', 'bytes': '840M'}, 'annotation_gtf': {'ftp': 'http://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz', 'ensembl_release': 108, 'release_date': '2022-10-04', 'release_time': '21:10', 'bytes': '52M'}, 'coding_seq_cds': {'ftp': 'http://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz', 'ensembl_release': 108, 'release_date': '2022-10-04', 'release_time': '20:18', 'bytes': '21M'}, 'non-coding_seq_ncRNA': {'ftp': 'http://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz', 'ensembl_release': 108, 'release_date': '2022-10-04', 'release_time': '21:05', 'bytes': '18M'}, 'protein_translation_pep': {'ftp': 'http://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz', 'ensembl_release': 108, 'release_date': '2022-10-04', 'release_time': '20:18', 'bytes': '14M'}}}



```python

assembly, ensembleID_coordinates, sequence = ensembleID2seq(ensembleID)
print(assembly)
print(ensembleID_coordinates)

```

    Thu Dec  1 10:42:39 2022 INFO Fetching reference information for homo_sapiens from Ensembl release: 108.
    Thu Dec  1 10:42:39 2022 INFO We noticed that you may have passed a version number with your Ensembl ID.
    Please note that gget seq will return information linked to the latest Ensembl ID version.



```python
#print(gget.seq([ensembleID]))
```


```python
# gget() outputs GRCh37 sequence, is it correct? let's verify using ensemble rest api.
# let's try:
# source: https://rest.ensembl.org/documentation/info/sequence_id
import requests, sys

ensembleID = 'ENST00000226574'
server = "https://rest.ensembl.org"
ext = "/sequence/id/%s?" %ensembleID
 
r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
 
#print(r.text)
```


```python
# We got the sequence, let's generate a vcf file

vcfMetaData = '''##fileformat=VCFv4.3
##fileDate=%s
##source=transcript2vcf.py
##reference=Homo sapiens
##contig=<ID=-,length=%s,assembly=%s,md5=-,species="Homo sapiens",taxonomy=x>
##phasing=partial
#ID: mut<num>-<num>,Type=String,Description="Informs about the mutation number and postion on the transcript."
##INFO=<ID=DP,Number=1,Type=Integer,Description="Artificial Depth, will be the same for all; 500">
##INFO=<ID=AF,Number=A,Type=Float,Description="Artificial AF, will be the same for all; 1.0">
##FILTER=<ID=PASS,Artificial mutation data, ALl PASS>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Artificial Read Depth, will be the same for all; 500">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Artificial data quality, will be the same for all; 60">''' %(getTimeStamp(),
                                                                            len(sequence),
                                                                           assembly,)


```


```python

```


```python
vcfData_list = ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' %ensembleID]

chrom, start, end = ensembleID_coordinates
currentCoordinates = start - 1

hitFound = 0
hitPosOnTranscript = 0
for eachPos in sequence:
    
    # counter's update
    currentCoordinates += 1
    hitPosOnTranscript += 1
    
    # check if ADAR will attack here
    if hitByADAR(eachPos):
        hitFound += 1
        
        # add it to the VCF file
        hitID = 'mut%s-%s' %(hitFound,
                                    hitPosOnTranscript)
        vcfLine = "%s\t%s\t%s\tA\tG\t60\tPASS\tDP=500;AF=1.0\tGT:GQ:DP:HQ\t0|1:60:500:60,60" %(chrom,
                        currentCoordinates,
                         hitID)
        vcfData_list.append(vcfLine)
```


```python
# save the vcf file
odir = '%s/data' %(os.getcwd())
vcfFileName = "%s/ADARhits_%s.vcf" %(odir,
                                    ensembleID.split(".")[0])
wrtieVCF_FH = open(vcfFileName, 'w+')
wrtieVCF_FH.write("%s\n" %vcfMetaData)
wrtieVCF_FH.write("\n".join(vcfData_list))
wrtieVCF_FH.close()
```


```python

```
