<!DOCTYPE html>
<html lang="en">
<head> 
    <meta charset="UTF-8"> 
    <title>Title</title> 
</head> 
<body>

</body>
</html>

# VCF VarSelect
This is a python package to select whole-exome sequencing variants from VCF files. The parse of VCF file is based on [VCF Parser](https://github.com/moonso/vcf_parser). Paper is coming soon. 

## Installation
```shell script
git clone http://github.com/DanyangLi107/vcf_varselect.git
cd vcf_varselect
pip3 install .
```

## File preparation
vcf file: a text file with extension .vcf or .vcf.gz. A tiny example can be found in: /example/example.vcf

gene_file: disorder related gene list, an example is in: /example/gene_list.csv

gender_file: male samples list, an example is in: /example/male_list.txt

innerfreq_file: file of variant inner-freq from all samples, an example is in: /example/inner_freq.json

## Basic function
### Return a dictionary of selected variants from one sample VCF file
Total variants in VCF:
```python
from vcf_varselect import VariantSelection
vcf = VariantSelection(infile='file.vcf')
for var in vcf:
    print (vcf.variant)
```
VCF file information:
```python
vcf.header         # header information in vcf
vcf.id_dict        # information of INFO, FORMAT and FILTER in vcf
vcf.vep_columns    # information VEP annotation in vcf
vcf.sample         # sample ID in vcf
```
Select variants with good quality:
```python
quality = vcf.quality_selection(FILTER='PASS', DP=10.0, QD=2.0, MQ=40.0)
for var in vcf:
    print(quality)
``` 
Select rare variants:
```python
freq = vcf.freq_selection(KG=0.001, EXAC=0.001, GNOMAD=0.001, SWEGEN=0.001, innerfreqfile=innerfreq_file)
for var in vcf:
    print(freq)
```
Select damaging variants:
```python
damaging, lof, mis_damage = vcf.damaging_selection(criteria=['SIFT', 'POLYPHEN', 'MPC', 'CADD', 'SPIDEX', 'PHYLOP'])
for var in vcf:
    print(damaging)
    print(lof)
    print(mis_damage)
```
Select rare damaging, rare loss-of-function and rare missense variants:
```python
damaging_var, lof_var, mis_var = vcf.comb_selection(FILTER='PASS', DP=10.0, QD=2.0, MQ=40.0,
                                          KG=0.001, EXAC=0.001, GNOMAD=0.001, SWEGEN=0.001, innerfreqfile=innerfreq_file,
                                          criteria=['SIFT', 'POLYPHEN', 'MPC', 'CADD', 'SPIDEX', 'PHYLOP'])
for var in vcf:
    print(damaging_var)
    print(lof_var)
    print(mis_var)
```
Select variants of disorder related genes:
```python
from vcf_varselect import match_gene
gene_var = match_gene(damaging_var, gene_file, gender_file)
```
### Return a dataframe of disorder related rare damaging variants from multiple samples
```python
from vcf_varselect import sample_combine
df = sample_combine(dir, innerfreq_file, gene_file, gender_file,
                    FILTER='PASS', DP=10.0, QD=2.0, MQ=40.0,
                    KG=0.001, EXAC=0.001, GNOMAD=0.001, SWEGEN=0.001,
                    criteria=['SIFT', 'POLYPHEN', 'MPC', 'CADD', 'SPIDEX', 'PHYLOP'])
```

