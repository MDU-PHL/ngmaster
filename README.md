[![Build Status](https://travis-ci.org/MDU-PHL/ngmaster.svg?branch=master)](https://travis-ci.org/MDU-PHL/ngmaster)
[![License: GPLv2](https://img.shields.io/badge/License-GPL_2.0-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
![Python 3.6](https://img.shields.io/badge/Language-Python_3.6-steelblue.svg)
[![DOI:10.1099/mgen.0.000076](https://zenodo.org/badge/DOI/10.1099/mgen.0.000076.svg)](https://doi.org/10.1099/mgen.0.000076)

# ngmaster

*In silico* multi-antigen sequence typing for *Neisseria gonorrhoeae* (NG-MAST).  

## Synopsis
```
% ngmaster gono.fa
ID         NG-MAST    POR    TBPB
gono.fa    10699      6277   4
```

## Dependencies

* [Python >= 3.6](https://www.python.org/)
* [BioPython](http://biopython.org/)
* [isPcr >= v33x2](http://hgwdev.cse.ucsc.edu/~kent/src/) by Jim Kent

## Installation

#### PiPy
```
pip3 install ngmaster
```
#### Brew
```
brew install brewsci/bio/ngmaster
```
#### Conda
```
conda install -c conda-forge -c bioconda -c defaults ngmaster  # COMING SOON
```

## Test

Once installed, you can run the following to ensure `ngmaster` is successfully working:

    $ ngmaster --test

If everything works, you will see the following:

```
Running ngmaster.py on test example (NG-MAST 10699) ...
$ ngmaster.py test/test.fa
ID    NG-MAST    POR    TBPB
test.fa    10699    6277    4
... Test successful.
```

## Usage

    $ ngmaster -h

    usage:
      ngmaster [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>

    In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)

    Please cite as:
      Kwong JC, Goncalves da Silva A, Howden BP and Seemann T.
      NGMASTER: in silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)
      GitHub: https://github.com/MDU-PHL/ngmaster

    positional arguments:
      FASTA            input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN

    optional arguments:
      -h, --help       show this help message and exit
      --db DB          specify custom directory containing allele databases
                       directory must contain database files "POR.tfa", "TBPB.tfa", and "ng_mast.txt"
      --csv               output comma-separated format (CSV) rather than tab-separated
      --printseq FILE  specify filename to save allele sequences to (default=off)
       --updatedb       update allele database from <www.ng-mast.net>
      --test           run test example
      --version        show program's version number and exit


## Quick start

**To perform *in silico* NG-MAST on FASTA files:**

`$ ngmaster <fasta1> <fasta2> <fasta3> ... <fastaN>`

The NG-MAST result and allele numbers are printed in tab-separated format to `stdout`.
* If an allele is not found (ie. unable to located with primers), the allele result is "`–`".
* If an allele is found (ie. located with primers), but the conserved region containing the starting key motif required for sequence trimming cannot be located, the allele result is "`no_key`".
* If an allele is found (ie. located with primers), but the trimmed sequence is novel, and not in the current database, the allele result is "`new`".

**To save results to a tab-separated text file, redirect `stdout`:**

`$ ngmaster <fasta1> <fasta2> <fasta3> ... <fastaN>  > results.txt`

**To display results in comma-separated format, use the `--csv` option:**

`$ ngmaster --csv <fasta1> <fasta2> <fasta3> ... <fastaN>`

**To save sequences of the alleles to a file (eg. for uploading "new" sequences to [http://www.ng-mast.net](http://www.ng-mast.net/)):**

`$ ngmaster --printseq [filename] <fasta1> <fasta2> <fasta3> ... <fastaN>`

## Updating the allele databases

**To update the allele databases from http://www.ng-mast.net :**  
*Warning: This will overwrite the existing databases so ensure you back them up if you wish to keep them.*

    $ ngmaster.py --updatedb

A copy of the old database is saved just in case, but is overwritten with each subsequent   ```--updatedb```.

**To update the allele databases into a different folder (ie. not the /db folder in the ngmaster directory):**

    $ ngmaster.py --updatedb --db path/to/folder

This will download the database files into the folder ```path/to/folder```.
This can then be specified when running ngmaster using the ```--db  path/to/folder``` option.

## Creating a custom allele database

1. Create custom database files: `POR.tfa`, `TBPB.tfa`, `ng_mast.txt`  
   See default `db` directory for examples.  
   `POR.tfa` and `TBPB.tfa` contain the respective allele sequences in FASTA format.  
   `ng_mast.txt` contains a list of NG-MAST types and the corresponding allele types.

2. Place the custom database files in a folder.

3. Specify the path to that custom database folder:  
   `$ ngmaster --db [/path/to/custom/folder/] <fasta1> <fasta2> <fasta3> ... <fastaN>`

## Citation

Kwong JC, Gonçalves da Silva A, Dyet K, Williamson DA, Stinear TP, Howden BP and Seemann T.  
*NGMASTER: in silico multi-antigen sequence typing for Neisseria gonorrhoeae*
**Microbial Genomics**
2016 Aug 25;2(8):e000076.
PMID:[28348871](https://www.ncbi.nlm.nih.gov/pubmed/28348871)
DOI:[10.1099/mgen.0.000076](https://doi.org/10.1099/mgen.0.000076)

## Bugs

### Software
Please submit via the [GitHub issues page](https://github.com/MDU-PHL/ngmaster/issues).  

### Database
Note that the NG-MAST databases and website are curated and hosted at the
Department of Infectious Disease Epidemiology, Imperial College London.  For
issues with the NG-MAST databases, please contact the [NG-MAST
curator](mailto:d.aanensen@imperial.ac.uk).

## Software Licence

[GPLv2](https://github.com/MDU-PHL/ngmaster/blob/master/LICENSE)

## References

* Martin et al. J Infect Dis, 2004 Apr 15; 189(8): 1497-1505.  
* See also [http://www.ng-mast.net](http://www.ng-mast.net/).

## Authors

* Jason Kwong (@kwongjc)
* Anders Gonçalves da Silva (@drandersgs)
* Mark Schultz (@schultzm)
* Torsten Seemann (@torstenseemann)

## Development

When incrementing version (i.e., minor patch), run the following:

```
bumpversion --verbose --dry-run --new-version <major.minor.patch> patch
bumpversion --new-version <new.version.number> patch
```

The same can be done for `minor` and `major` numbers.

This will automatically commit and tag the commit with the new version number.
It will also update the necessary location in the file.

## Pushing to pypi

**Must be uploaded to maintainer's account.**

```
bumpversion --new-version <new.version.number> <patch|minor|major>
git push
# create distribution
python3 setup.py sdist bdist_wheel
# run some checkes
twine check dist/*
# upload to test pypi to see if everything works
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
# upload to pypi
twine upload dist/*
```
