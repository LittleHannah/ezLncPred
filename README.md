# ***LncPac***: An integrated python package and web server for LncRNA identification



*LncPac* is an comprehensive python package for LncRNA identification which integrate 9 state-of-the-art lncRNA prediction models. *LncPac* python package provides a convenient command line method for researchers who intends to identify LncRNA. The *LncPac* web server tool provides a straight-to-the-point answer for input `fasta` file, or an email for time-consuming identification task.

##Integration



*LncPac* currently provides 9 LncRNA prediction models, which are listed as follows. 

 - CNCI
 - CPC2
 - lgc
 - PLEK
 - CPAT
 - CPPred
 - longdist
 - PredLnc-GFStack
 - LncADeep

##Web server
------------


##Python package installation



 - Prerequisite
    - python 3.0 version (or above)
    - linux operating system
    - Biopython
    - C/C++ compiler(for PLEK)
 - Download *LncPac* by

```bash
pip install LncPac
```

##Help
--------------------

For detailed message of *LncPac*, run

```bash
LncPac --manual
```

For detailed message of each model and their parameters, run

```bash
LncPac  --manual -m [model]
```

##Usage
----------------------

*LncPac* offers a total of 9 LncRNA prediction models, each with a different variety of parameter choices, users can refer to the list below to customarize your prediction procedure.
First, *LncPac* **must** receive at least three parameters to specify the `input file` `output directory` and `prediction model`

 - `-i` or `--input`  <a href= >fasta</a> format input files
 - `-o` or `--output` the output directory to store the identification results
 - `-m` or `--model` one out of nine LncRNA identification models to choose, careful not to make case error or spelling error

Individual model parameters

>CNCI

 - `-g` or `--gtf`  if your input file is gtf format please use this parameter
 - `-d` or `--directory` specify the path of your reference genome if your input file is gtf format
 - `--parallel` assign the running CUP numbers
<br>
example
```bash
LncPac -m CNCI -i example.fa -o results
```

> CPC2

example
```bash
LncPac -m CPC2 -i example.fa -o results/CPC2_outfile
```
> lgc

example
```bash
LncPac -m lgc -i example.fa -o results/lgc_outfile
```
> PLEK   

 - `--thread` the number of threads to run the PLEK programme
 - `-z` or `--size` minimum sequence size to consider, default is 200
 - `--isoutmsg` output message to stdout or not, the existence of this parameter means that PLEK will be run quietly
 - `isrmtempfile` remove temporary files or not, the existence of this parameter means that PLEK programme will remove all temporary files
<br>
example
```bash
LncPac -m PLEK -i example.fa -o results/PLEK_outfile
```
>CPAT
  
 - `-p` or `--species` specify the species of the LncRNAs choose from `Human` `Mouse` `Fly` `Zebrafish` (note that the first character is upper case)
 - `-s` or `--start` Start codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF), default is ATG
 - `-t` or `--stop` Stop codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). Multiple stop codons should be separated by ',' default is TAG,TAA,TGA
<br>
example
```bash
LncPac -i example.fa -m CPAT -o results/CPAT_outfile -p Human
```
>CPPred

 - 
