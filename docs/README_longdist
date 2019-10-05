# longdist

## Requirements

Python 3, Tkinter and pip are required for running this tool. Also, the
following packages will be installed during the installation process:
- NumPy - ``pip install numpy``
- Biopython - ``pip instal biopython``
- scipy - ``pip install scipy``
- scikit-learn - ``pip install scikit-learn``
- configparser - ``pip install configparser``
- matplotlib - ``pip install matplotlib``

These packages will be automatically installed during the installation process

## Installation

```
wget https://github.com/hugowschneider/longdist.py/archive/v1.0.3.tar.gz
tar zxvf v1.0.3.tar.gz
cd longdist.py-1.0.3
pip3 install .
```

If you install these package in a local repository, please remember to set the ``PYTHONPATH``
and ``PATH`` environment variables correctly.

## Usage
```
usage: longdist.py [-h] [--citation] [--version]
                   [--longs <longs1.fa longs2.fa ...> [<longs1.fa longs2.fa ...> ...]]
                   [--pcts <pcts1.fa pcts2.fa ...> [<pcts1.fa pcts2.fa ...> ...]]
                   [--input <input.fa>] [--kmers <50>] [--orf <1>]
                   [--ratio <0.75>] [--size <200>] [--cv <10>]
                   [--log2c <-5,15,2>] [--log2g <3,-15,-2>] [--processes <1>]
                   [--out_roc <"lncRNA file"x"PCT file"x"kmers"_roc.eps>]
                   [--out_csv <"lncRNA file"x"PCT file"x"kmers".csv>]
                   [--out_model <"lncRNA file"x"PCT file"x"kmers".plk>]
                   [--predict]
                   [--model_config <"lncRNA file"x"PCT file"x"kmers".plk>]
                   [--out <"Input File".csv>] [--purge]

longdist: Method implementation for long ncRNAs and PCT distinction. This
application can create and use models base on the method by Schneider et al
(2017).

optional arguments:
  -h, --help            show this help message and exit
  --citation            Prints bibtex citation.
  --version             Prints version number.

Method Paramenters:
  --longs <longs1.fa longs2.fa ...> [<longs1.fa longs2.fa ...> ...]
                        List of fasta files containing long non-coding RNAs.
                        The files should the separated by spaces and should be
                        ordered by species in same order as the pct file list.
                        This argument is required.
  --pcts <pcts1.fa pcts2.fa ...> [<pcts1.fa pcts2.fa ...> ...]
                        List of fasta files containing protein coding
                        transcripts. The files should the separated by spaces
                        and should be ordered by species in same order as the
                        lcnRNA file list. This argument is required.
  --input <input.fa>    Fasta file containing transcripts to predict with the
                        model
  --kmers <50>          Number of nucleotide pattern frequencies to consider
                        in the model. Default is 50.
  --orf <1>             The orf feature to be used by the model. Default is 1.
                        Possible values are: 0 - No orf feature; 1 - First ORF
                        relative length; 2 - Longest ORF relative length
  --ratio <0.75>        The ratio of whole dataset that should be used for
                        training. Default is 0.75.
  --size <200>          Minimum sequence size to consider. Default is 200.
  --cv <10>             Number of folds in cross-validation. Default is 10.
  --log2c <-5,15,2>     Set the range of c to
                        2^{begin,...,begin+k*step,...,end}. Default is
                        -5,15,2.
  --log2g <3,-15,-2>    Set the range of g to
                        2^{begin,...,begin+k*step,...,end}. Default is
                        3,-15,-2.
  --processes <1>       Number of parallel processes for parameters search.
                        Default is 1.
  --out_roc <"lncRNA file"x"PCT file"x"kmers"_roc.eps>
                        Name of the output file for the roc Curve. Default is
                        roc.eps.
  --out_csv <"lncRNA file"x"PCT file"x"kmers".csv>
                        Name of the output CSV file containing the results.
                        Default is a name built from the names of both fasta
                        files.
  --out_model <"lncRNA file"x"PCT file"x"kmers".plk>
                        Name of the output file containing the SVM Model.
                        Default is a name built from the names of both fasta
                        files.
  --predict             Just use a predefined model to distinguish long ncRNAs
                        and PCTs in the input fasta file
  --model_config <"lncRNA file"x"PCT file"x"kmers".plk>
                        The file name containing the model configuration
                        properties for prediction.
  --out <"Input File".csv>
                        The output CSV file for the result of the distinction
                        made in the input file.
  --purge               Purge all intermediate files. Intermediate files have
                        ".longdist." in their names. All intermediate files
                        are used to accelerate consecutive runs of the method.
                        Don't purge them if you want to run this method a
                        second time with the same data
```
### Building Models

The following command builds a model with data from two fasta files, one containing
long noncoding RNAs (``test/GRCm38.lncRNA.fa``) and other containing protein coding
(``test/GRCm38.pct.fa``):

```
longdist --longs test/GRCm38.lncRNA.fa --pcts test/GRCm38.pct.fa
```

This command line will build a model with 50 nucleotide patterns frequencies (kmers)
and the first ORF relative length. To change the number of kmers, should be included
the parameter ``--kmers``, for example:

```
longdist --longs test/GRCm38.lncRNA.fa --pcts test/GRCm38.pct.fa --kmers 10
```

It is also possible to change the training data ratio for building the model with
the parameter ``--ratio``, for example:

```
longdist --longs test/GRCm38.lncRNA.fa --pcts test/GRCm38.pct.fa --ratio 0.5
```

The paramenters ``--log2c`` and ``--log2g`` changes the search space for the C and
gamma parameter, for example:

```
longdist --longs test/GRCm38.lncRNA.fa --pcts test/GRCm38.pct.fa --log2c 1,15,2 --log2g 3,-1,-1
```

The model build process creates intermediate files to accelerate the build of the
following models. This intermediate files can be deleted with the parameter ``--purge``.
Also, some output files will be created to evaluated model performances. All file names
are based on the input files, for example ``GRCm38.lncRNA_x_GRCm38.pct_50.csv``.

The process outputs the following output files:
- A CSV file with the prediction results;
- An EPS file with the ROC Curve;
- A PLK and a PLK.CONF files. PLK is the SVM Model and the PLK.CONF is the
configuration to use this model.

### Distinguishing lncRNAs and PCTs

To distinguish PCTs and lncRNAs from a input fasta file (``input.fa``) using a pre-built
model (``model.plk.conf``), the following command should be used:

```
longdist --predict --input input.fa --model_config model.plk.conf
```

This command will create a csv file named ``input.fa.csv`` with the prediction results.

### Models

Models are PLK files built with scikit-learn package and they are only compatible with
the package's SVM implementation. To use the model the PLK.CONF file should be used together
with the model. To move the model to another folder, the PLK.CONF should be moved to.

The selected attributes for the model are listed in the PLK.CONF file.

## Copyright
The work herein is Copyright 2013--2017 Hugo Wruck Schneider and Universidade de Brasília (UnB). **No rights are given to reproduce or modify this work**.
