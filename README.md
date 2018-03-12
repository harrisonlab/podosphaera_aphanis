# podosphaera_aphanis
Commands used in the assembly and annotation of the podosphaera aphanis genome



#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```bash
ProjectDir=/home/groups/harrisonlab/project_files/podosphaera
mkdir -p $ProjectDir
```

```bash
  ProjectDir=/home/groups/harrisonlab/project_files/podosphaera
  RawDat=/data/seq_data/miseq/2018/ANALYSIS/180227_M04465_0070_000000000-B86MV/Data/Intensities/BaseCalls
  OutDir=raw_dna/paired/P.aphanis/C1/F
  mkdir -p $OutDir
  cd $OutDir
  cp -s $RawDat/MildewC1_S1_L001_R1_001.fastq.gz .
  cd $ProjectDir
  OutDir=raw_dna/paired/P.aphanis/C1/R
  mkdir -p $OutDir
  cd $OutDir
  cp -s $RawDat/MildewC1_S1_L001_R2_001.fastq.gz .
  cd $ProjectDir
```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:


```bash
  for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```


Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf

```bash
for StrainPath in $(ls -d raw_dna/paired/P.*/*); do
echo $StrainPath
Read_F=$(ls $StrainPath/F/*.fastq.gz)
Read_R=$(ls $StrainPath/R/*.fastq.gz)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```


```bash
  for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```


#Assembly

Assembly was performed with:
* Spades

## Spades Assembly

Assembly was submitted for genomes with a single run of data

```bash
for StrainPath in $(ls -d qc_dna/paired/P.*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
OutDir=assembly/spades/$Organism/$Strain
echo $F_Read
echo $R_Read
qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct
done
```
