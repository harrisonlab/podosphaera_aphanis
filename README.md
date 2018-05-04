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

Reads were aligned to the vesca assembly to assess presence of strawberry DNA in
the sample:

```bash
Reference=$(ls ../../../../../data/scratch/armita/Fvesca-genome.v4.0.a1/assembly/Fragaria_vesca_v4.0.a1.fasta)
for StrainPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_vesca
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
done
```

```bash
mkdir -p assembly/external_group/F.annanassa/redgauntlet_2016-11-25/contigs
scp -r ../../../../../home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq/assemblies/redgauntlet_contigs_2016-11-25_ei_version.fasta.gz assembly/external_group/F.annanassa/redgauntlet_2016-11-25/contigs/.
gunzip assembly/external_group/F.annanassa/redgauntlet_2016-11-25/contigs/redgauntlet_contigs_2016-11-25_ei_version.fasta.gz

Reference=$(ls assembly/external_group/F.annanassa/redgauntlet_2016-11-25/contigs/redgauntlet_contigs_2016-11-25_ei_version.fasta)
for StrainPath in $(ls -d qc_dna/paired/*/*); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_F.ananassa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie_unaligned.sh $Reference $F_Read $R_Read $OutDir
done
```

Reads were also aligned to the podosphaera xanthii transcriptome to assess
the relative level of podosphaera reads in the dataset:

```bash
DownloadDir=assembly/external_group/P.xanthii/SF2086/transcriptome
mkdir -p $DownloadDir
cd $DownloadDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/GE/UO/GEUO01/GEUO01.1.fsa_nt.gz
gunzip *.gz
cp -s GEUO01.1.fsa_nt P.xanthii_SF2086_transcriptome.fa
cd /home/groups/harrisonlab/project_files/podosphaera

Reference=$(ls $DownloadDir/P.xanthii_SF2086_transcriptome.fa)
for StrainPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_P.xanthii
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
done
```

Reads were also aligned to the grape mildew (Erysiphe necator) genome to assess
the relative level of podosphaera reads in the dataset:
https://genome.jgi.doe.gov/Erynec1/Erynec1.home.html

```bash
DownloadDir=assembly/external_group/E.necator/strain_c/assembly
mkdir -p $DownloadDir
cd $DownloadDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/JN/VN/JNVN01/JNVN01.1.fsa_nt.gz
gunzip *.gz
cp -s JNVN01.1.fsa_nt E.necator_strain_c_assembly.fa
cd /home/groups/harrisonlab/project_files/podosphaera

Reference=$(ls $DownloadDir/E.necator_strain_c_assembly.fa)
for StrainPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_E.necator
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
done
```

Reads were also aligned to the wheat mildew (Blumeria graminis f. sp. tritici)
genome to assess the relative level of podosphaera reads in the dataset:
https://genome.jgi.doe.gov/Erynec1/Erynec1.home.html

```bash
DownloadDir=assembly/external_group/B.graminis_f.sp._tritici/96224/assembly
mkdir -p $DownloadDir
cd $DownloadDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AN/ZE/ANZE01/ANZE01.1.fsa_nt.gz
gunzip *.gz
cp -s ANZE01.1.fsa_nt B.graminis_f.sp._tritici_96224_assembly.fa
cd /home/groups/harrisonlab/project_files/podosphaera

Reference=$(ls $DownloadDir/B.graminis_f.sp._tritici_96224_assembly.fa)
for StrainPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_B.graminis
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```bash
  for TrimPath in $(ls -d qc_dna/paired/P.*/*); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    TrimF=$(ls $TrimPath/F/*.fq.gz)
    TrimR=$(ls $TrimPath/R/*.fq.gz)
    echo $TrimF
    echo $TrimR
    qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
  done
```

Reads were also aligned to PhiX

```bash
Reference=$(ls ../../../../../home/harrir/git_master/seq_tools/phix.fa)
for StrainPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_phiX
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
done
```


#Assembly

Assembly was performed with:
* Spades

## Spades Assembly

Assembly was submitted for genomes with a single run of data
<!--
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
``` -->

```bash
  for RawData in $(ls analysis/genome_alignment/bowtie/P.aphanis/C1/vs_F.ananassa/*.fq.*.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```


```bash
for StrainPath in $(ls -d analysis/genome_alignment/bowtie/P.aphanis/C1/vs_F.ananassa); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f3 -d '/' | rev)
F_Read=$(ls $StrainPath/*.fq.1.gz)
R_Read=$(ls $StrainPath/*.fq.2.gz)
OutDir=assembly/spades/$Organism/${Strain}_no_strawberry
# OutDir=assembly/spades/$Organism/${Strain}_no_strawberry_meta
echo $F_Read
echo $R_Read
qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct
# qsub $ProgDir/submit_metaSPAdes.sh $F_Read $R_Read $OutDir
done
```


Quast and busco were run to assess assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep '_no_strawberry'); do
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
  OutDir=gene_pred/busco/$Organism/$Strain/assembly
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


Potential contmainant contigs were identified using deconseq


```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep '_no_strawberry'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    # for Exclude_db in "bacillus" "delftia" "paenibacillus" "stenotrophomonas" "vesca" "ananassa" "rickettsia"; do
    for Exclude_db in "rickettsia"; do
      Good_db="blumeria"
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
    done
  done
```

```bash
for Exclude_db in "bacillus" "delftia" "paenibacillus" "stenotrophomonas" "vesca" "ananassa"; do
echo $Exclude_db
for File in $(ls assembly/spades/P.*/*/*/log.txt | grep "_${Exclude_db}"); do
Name=$(echo $File | rev | cut -f3 -d '/' | rev);
Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
Both=$(cat $File |cut -f2 | head -n2 | tail -n1);
Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
printf "$Name\t$Good\t$Both\t$Bad\n";
done
done
```
```
bacillus
C1_no_strawberry	150998	77	56
delftia
C1_no_strawberry	151117	1	13
paenibacillus
C1_no_strawberry	151113	1	17
stenotrophomonas
C1_no_strawberry	151109	4	18
vesca
C1_no_strawberry	21145	33	129953
ananassa
C1_no_strawberry	17789	17	133325
```

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep '_no_strawberry'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
StrainDir=$(ls -d assembly/spades/$Organism/$Strain)
mkdir -p $StrainDir/deconseq_appended
for File in $(ls $StrainDir/deconseq_*/*cont.fa | grep -e "bacillus" -e "delftia" -e "paenibacillus" -e "stenotrophomonas" -e "vesca" -e "ananassa" ); do
cat $File | grep '>'
done | sort | uniq | tr -d '>' > $StrainDir/deconseq_appended/exclude_list.txt
Instructions=$StrainDir/deconseq_appended/exclude_instructions.txt
printf "Exclude:\nSequence name, length, apparent source\n" > $Instructions
cat $StrainDir/deconseq_appended/exclude_list.txt | sed -r 's/$/\t.\tcontaminant/g' >> $Instructions
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $StrainDir/deconseq_appended/contigs_min_500bp_renamed.fasta --coord_file $Instructions
done
```

Quast and busco were run to assess assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/P.aphanis/C1_no_strawberry/deconseq_appended/contigs_min_500bp_renamed.fasta); do
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
  # OutDir=gene_pred/busco/$Organism/$Strain/assembly
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```
Assembly                   contigs_min_500bp_renamed
# contigs (>= 0 bp)        17365                    
# contigs (>= 1000 bp)     8764                     
Total length (>= 0 bp)     39188020                 
Total length (>= 1000 bp)  33700446                 
# contigs                  17365                    
Largest contig             30706                    
Total length               39188020                 
GC (%)                     42.87                    
N50                        4533                     
N75                        1971                     
L50                        2539                     
L75                        5778                     
# N's per 100 kbp          3.88     
```

```
1117	Complete BUSCOs (C)
1115	Complete and single-copy BUSCOs (S)
2	Complete and duplicated BUSCOs (D)
101	Fragmented BUSCOs (F)
97	Missing BUSCOs (M)
1315	Total BUSCO groups searched
```

Next step, I will run busco vs Blumeria and Vesca and extract both of their busco hits. I will align the busco hits from the 3 organisms and see wether the hits to the remiaing contigs look more fungal or plant-al


```bash
# P. aphanis
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/P.aphanis/C1_no_strawberry/deconseq_appended/contigs_min_500bp_renamed.fasta); do
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=analysis/busco_eukaryote/$Organism/$Strain/assembly
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
# Mildew genomes
  for Assembly in $(ls assembly/external_group/*/*/assembly/*.fa); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/busco_eukaryote/$Organism/$Strain/assembly
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
# Strawberry genomes
for Assembly in $(ls /data/scratch/armita/Fvesca-genome.v4.0.a1/assembly/Fragaria_vesca_v4.0.a1.fasta); do
  Organism="F.vesca"
  Strain="v4.0"
  echo "$Organism - $Strain"
  OutDir=analysis/busco_eukaryote/$Organism/$Strain/assembly
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
for Assembly in $(ls assembly/external_group/F.annanassa/redgauntlet_2016-11-25/contigs/redgauntlet_contigs_2016-11-25_ei_version.fasta); do
  Organism="F.ananassa"
  Strain="redgauntlet"
  echo "$Organism - $Strain"
  OutDir=analysis/busco_eukaryote/$Organism/$Strain/assembly
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls analysis/busco_eukaryote/*/*/assembly/*/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

```
B.graminis_f.sp._tritici	96224	290	4	9	303
E.necator	strain_c	293	4	6	303
F.ananassa	redgauntlet	221	37	45	303
F.vesca	v4.0	245	9	49	303
P.aphanis	C1_no_strawberry	282	8	13	303
```


Find single copy busco genes in assemblies

Create a list of all BUSCO IDs

```bash
cd /home/groups/harrisonlab/project_files/podosphaera

OutDir=analysis/popgen/busco_phylogeny
mkdir -p $OutDir
BuscoDb="eukaryota_odb9"
ls -1 /home/groups/harrisonlab/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

For each busco gene create a folder and move all single copy busco hits from
each assembly to the folder.
Then create a fasta file containing all the aligned reads for each busco gene for
alignment later.

```bash
printf "" > analysis/popgen/busco_phylogeny/single_hits.txt
for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
OutDir=analysis/popgen/busco_phylogeny/$Busco
mkdir -p $OutDir
for Fasta in $(ls analysis/busco_eukaryote/*/*/assembly/*/single_copy_busco_sequences/$Busco*.fna); do
Strain=$(echo $Fasta | rev | cut -f5 -d '/' | rev)
Organism=$(echo $Fasta | rev | cut -f6 -d '/' | rev)
FileName=$(basename $Fasta)
cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
done
cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny/single_hits.txt
done
```

If all isolates have a single copy of a busco gene, move the appended fasta to
a new folder

```bash
  OutDir=analysis/popgen/busco_phylogeny/alignments
  mkdir -p $OutDir
  OrganismNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
  for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
  echo $Busco
  HitNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
  if [ $HitNum == $OrganismNum ]; then
    cp analysis/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
  fi
  done
```

Submit alignment for single copy busco genes with a hit in each organism


```bash
  AlignDir=analysis/popgen/busco_phylogeny/alignments
  CurDir=$PWD
  cd $AlignDir
  ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
  qsub $ProgDir/sub_mafft_alignment.sh
  cd $CurDir
```


Results confirmed that busco genes in the assembly were not strawberry. Blast
searches of publically available P. aphanis ITS sequence vs the genome and
reciprocal blast search of the hit vs ncbi confirmed that the genome was a
podosphaera spp.. Interestingly there werea few SNPs that were distinct from the
reference sequence.

Coverage of the assembled genome was determined through alignment of trimmed
reads vs the assembly.

```bash

cd /home/groups/harrisonlab/project_files/podosphaera

Reference=$(ls assembly/spades/P.aphanis/C1_no_strawberry/deconseq_appended/contigs_min_500bp_renamed.fasta)
for StrainPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_P.aphanis_no_strawberry
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
done
```

```bash
cat analysis/genome_alignment/bowtie/P.aphanis/C1/vs_P.aphanis_no_strawberry/bowtie_log.txt
```

```
23687450 reads; of these:
  23687450 (100.00%) were paired; of these:
    20680676 (87.31%) aligned concordantly 0 times
    2682504 (11.32%) aligned concordantly exactly 1 time
    324270 (1.37%) aligned concordantly >1 times
    ----
    20680676 pairs aligned concordantly 0 times; of these:
      1109601 (5.37%) aligned discordantly 1 time
17.38% overall alignment rate
```


### Read coverage

Identify read coverage over each bp

```bash
  for Bam in $(ls analysis/genome_alignment/bowtie/*/*/vs_P.aphanis_no_strawberry/*_aligned_sorted.bam); do
    Target=$(echo $Bam | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain - $Target"
    OutDir=$(dirname $Bam)
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_${Target}_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_${Target}_depth.tsv > $OutDir/${Organism}_${Strain}_${Target}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_${Target}_depth_10kb.tsv
  done
  for Target in "vs_P.aphanis_no_strawberry"; do
    OutDir=analysis/genome_alignment/bowtie/grouped_${Target}
    mkdir -p $OutDir
    cat analysis/genome_alignment/bowtie/*/*/*/*_*_${Target}_depth_10kb.tsv > $OutDir/${Target}_grouped_depth.tsv
  done
```
