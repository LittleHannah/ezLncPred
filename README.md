
# ***ezLncPred***: A Python package with integrated lncRNA prediction methods



*ezLncPred* is an comprehensive python package for LncRNA identification which integrates *9* state-of-the-art lncRNA prediction models. *ezLncPred* python package provides a convenient command line method for researchers who intends to identify lncRNAs. The *ezLncPred* web server tool provides a straight-to-the-point answer for input fasta file, or an email for time-consuming identification task.

Integration
---------------------------------------------


*ezLncPred* currently provides 9 LncRNA prediction models, which are listed as follows. 

 - CNCI
 - CPC2
 - lgc
 - PLEK
 - CPAT
 - CPPred
 - longdist
 - PredLnc-GFStack
 - LncADeep

Web server
---------------------------------------------


Python package installation
---------------------------------------------


 - Prerequisite
    - python 3.0 version (or above)
    - linux operating system
    - C/C++ compiler(for PLEK)
 - Download *ezLncPred* by

```bash
pip install ezLncPred
```

Help
---------------------------------------------

For detailed message of *ezLncPred*, run

```bash
ezLncPred --manual
```

For detailed message of each model and their parameters, run

```bash
ezLncPred  --manual [model]{CNCI,CPC2,CPAT,lgc,CPPred,GFStack,longdist,PLEK,LncADeep}
```

Usage
---------------------------------------------

*ezLncPred* offers a total of 9 LncRNA prediction models, each with a different variety of parameter choices, users can refer to the list below to customarize your prediction procedure.
First, *ezLncPred* **must** receive at least three parameters to specify the `input file` `output directory` and `prediction model`

	-i --input		fasta format input files
	
	-o --output		the output directory to store the identification results
	
	-m --manual		show manuals
	
	-v --version		show program's version number and exit

For example, run
```bash
ezLncPred -i your_fasta_file -o output_directory model [parameters]
```


Individual model parameters
---------------------------------------------

>CNCI

	-h --help		show this help message and exit
	
	--parallel		assign the running CUP numbers
	
	-p {ve,pl} --species {ve,pl}
				assign the classification models ("ve" for vertebrate species, "pl" for plat species)
 
**example**
```bash
	ezLncPred CNCI -h
	ezLncPred -i example.fa -o results CNCI
	ezLncPred -i example.fa -o results CNCI --parallel
	ezLncPred -i example.fa -o results CNCI -p ve
```

>CPC2

	-h --help		show this help message and exit
	
	-r REVERSE --reverse
				REVERSE also check the reverse strand

**example**
```bash
	ezLncPred CPC2 -h
	ezLncPred -i example.fa -o results CPC2
	ezLncPred -i example.fa -o results CPC2 -r REVERSE
```

>lgc

	-h --help		show this help message and exit

**example**
```bash
	ezLncPred lgc -h
	ezLncPred -i example.fa -o results lgc
```


>PLEK   

	-h --help		show this help message and exit
	
	--thread		the number of threads to run the PLEK programme
	
	--isoutmsg		Output messages to stdout(screen) or not. "0" means 
				that PLEK be run quietly. "1" means that PLEK outputs
				the details of processing. Default value: 0
					
	--isrmtempfile		Remove temporary files or not. "0" means that PLEK 
				retains temporary files. "1" means that PLEK remove 
				temporary files. Default value: 1

**example**
```bash
	ezLncPred PLEK -h
	ezLncPred -i example.fa -o results PLEK
	ezLncPred -i example.fa -o results PLEK --thread 4
	ezLncPred -i example.fa -o results PLEK --isoutmsg 1
	ezLncPred -i example.fa -o results PLEK --isrmtempfile 0
```


>CPAT

	-h --help		show this help message and exit
	
	-p --species    	{Human,Mouse,Fly,Zebrafish}
				specify the species of the LncRNAs choose from Human 
				Mouse Fly Zebrafish (note that the first character 
				is upper case)
					
	-s --start		Start codon (DNA sequence, so use 'T' instead of 'U')
				used to define open reading frame (ORF), default is ATG
					
	-t --stop		Stop codon (DNA sequence, so use 'T' instead of 'U')
				used to define open reading frame (ORF). Multiple stop
				codons should be separated by ',' default is TAG,TAA,TGA

**example**
```bash
	ezLncPred CPAT -h
	ezLncPred -i example.fa -o results CPAT
	ezLncPred -i example.fa -o results CPAT -p Human
	ezLncPred -i example.fa -o results CPAT -s TAG
	ezLncPred -i example.fa -o results CPAT -t ATG,TGA,TTA
```


>CPPred

	-h --help		show this help message and exit
	
	-p --species		{Human,Integrated}
				the model of the species to choose (Human,Integrated).
	
**example**
```bash
	ezLncPred CPPred -h
	ezLncPred -i example.fa -o results CPPred
	ezLncPred -i example.fa -o results CPPred -p Integrated
```


>longdist

	-h --help		show this help message and exit
	
	-z <200>, --size <200>
				Minimun sequence size to consider. Default is 200.

	-p --species	{Human,Mouse}
				the model of the species to choose (human,mouse).

**example**
```bash
	ezLncPred longdist -h
	ezLncPred -i example.fa -o results longdist
	ezLncPred -i example.fa -o results longdist -z 150
	ezLncPred -i example.fa -o results longdist -p Human
```
	

>PredLnc-GFStack

	-h --help		show this help message and exit
	
	-p --species		{human,mouse}
				choose a species type from Human and Mouse

**example**
```bash
	ezLncPred GFStack -h
	ezLncPred -i example.fa -o results GFStack
	ezLncPred -i example.fa -o results GFStack -p human
```


>LncADeep

	-h --help		show this help message and exit

	-mt --modeltype	{full,partial}
				the model used for lncRNA identification,
				choose from partial full default is partial
				default is "partial".
					
	-HMM --HMMthread
				the thread number of using HMMER, default is 8
					
	-p --species	{human,mouse}
				the model of the species to choose (human,mouse).
				default is "human".

	-t --thread	THREAD
                    		The number of threads for running the LncADeep program.default is 8.

**example**
```bash
	ezLncPred LncADeep -h
	ezLncPred -i example.fa -o results LncADeep
	ezLncPred -i example.fa -o results LncADeep -mt full
	ezLncPred -i example.fa -o results LncADeep -HMM 4
	ezLncPred -i example.fa -o results LncADeep -p human
	ezLncPred -i example.fa -o results LncADeep -t 4
```

Test case
---------------------------------------------
We used lncRNA fasta file (1.0M) to test the cost time of prediction by each model. 
We run the default command without any additional parameters on *Intel(R) Xeon(R) Gold 6146 CPU @ 3.20GHz*.
The test results are for reference only, and the specific values ​​depend on the specific situation.

| Model | Real Time | User Time | Sys Time |
|:-----:|:---------:|:---------:|:--------:|
|CNCI|2m14.578s|2m13.992s|0m0.335s|
|CPC2|0m1.928s|0m1.756s|0m0.163s|
|CPAT|0m2.076s|0m1.896s|0m0.182s|
|lgc|0m1.640s|0m1.476s|0m0.162s|
|CPPred|0m13.973s|0m13.767s|0m0.200s|
|GFStack|0m44.763s|4m35.646s|0m14.584s|
|longdist|0m2.689s|0m2.569s|0m0.119s|
|PLEK|0m21.603s|1m30.905s|0m1.431s|
|LncADeep|4m1.667s|6m12.511s|0m23.426s|

Papers
---------------------------------------------
 - CNCI : *“Utilizing sequence intrinsic composition to classify protein-coding and long non-coding transcripts”*, Sun et al. (2013).
 - CPC2 : *“CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features”*, Kang, J., et al. (2017).
 - CPAT : *“CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model”*, Wang et al. (2013).
 - lgc : *“Characterization and identification of long noncoding RNAs based on feature relationship”*, Wang et al. (2019).
 - CPPred : *“CPPred: coding potential prediction based on the global description of RNA sequence”*, Tong et al. (2019).
 - GFStack : *“PredLnc-GFStack: A Global Sequence Feature Based on a Stacked Ensemble Learning Method for Predicting lncRNAs from Transcripts”*, Liu et al. (2019).
 - longdist : *“A Support Vector Machine based method to distinguish long non-coding RNAs from protein coding transcripts”*, W.Schneider, H., et al. (2017).
 - PLEK : *“PLEK: a tool for predicting long non-coding RNAs and messenger RNAs based on an improved k-mer scheme”*, Li, A., et al. (2014).
 - LncADeep : *“LncADeep: An ab initio lncRNA identification and annotation tool based on deep learning”*, Yang et al. (2017).


















