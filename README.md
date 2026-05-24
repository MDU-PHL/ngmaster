[![Tests](https://github.com/MDU-PHL/ngmaster/actions/workflows/test.yml/badge.svg?branch=master)](https://github.com/MDU-PHL/ngmaster/actions/workflows/test.yml)
[![Tests](https://github.com/MDU-PHL/ngmaster/actions/workflows/ci.yml/badge.svg)](https://github.com/MDU-PHL/ngmaster/actions/workflows/ci.yml)
[![GitHub release](https://img.shields.io/github/v/release/MDU-PHL/ngmaster)](https://github.com/MDU-PHL/ngmaster/releases)
![PyPI - License](https://img.shields.io/pypi/l/ngmaster)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/ngmaster)](https://anaconda.org/bioconda/ngmaster)
[![DOI:10.1099/mgen.0.000076](https://zenodo.org/badge/DOI/10.1099/mgen.0.000076.svg)](https://doi.org/10.1099/mgen.0.000076)

# ngmaster

*In silico* **m**ulti-**a**ntigen **s**equence **t**yping for ***N**eisseria **g**onorrhoeae* (NG-MAST) and
_**N**eisseria **g**onorrhoeae_ **s**equence **t**yping for **a**ntimicrobial **r**esistance (NG-STAR).  

## Synopsis
```
ngmaster gono.fa
FILE     SCHEME      NG-MAST/NG-STAR  porB_NG-MAST  tbpB  penA  mtrR  porB_NG-STAR  ponA  gyrA  parC  23S  CC
gono.fa  ngmaSTar    4186/231         2569           241   23    42    100           100   10    2     100  231
```

## What's New

- **Version and database version info**: `ngmaster --version` now reports both the tool version and the bundled database version. Note: `unauthenticated` database automatically means database pre-dating 2025.
- **Clonal complex (CC)**: output now includes a `CC` column from the updated NG-STAR database.
- **Duplicate 23S alleles**: multiple identical 23S alleles detected in a sample are now collapsed into a single call.
- **Large database fix**: requires `mlst >= 2.25.0`, which fixes alleles going undetected when the BLAST database is large (see [mlst v2.25.0](https://github.com/tseemann/mlst/releases/tag/v2.25.0)).
- **Database licensing notice**: the bundled database reflects the last freely redistributable snapshot of PubMLST. The PubMLST database after 2024-12-31 is subject to [PubMLST terms and conditions](https://pubmlst.org/terms-conditions) and cannot be redistributed without appropriate licensing. Use `--updatedb` (with authenticated access) to get the latest data. For authenticated access, ensure you are registered with PubMLST and have an API key to access the database update endpoint: `mlstdb connect -d pubmlst`.


## Dependencies

* [Python >= 3.7](https://www.python.org/)
* [BioPython](http://biopython.org/)
* [mlst >= 2.25.0](https://github.com/tseemann/mlst)
* [mlstdb](https://github.com/MDU-PHL/mlstdb) — for authenticated database updates

## Installation

#### PiPy

```
# mlst must be installed separately via conda
conda install -c bioconda mlst
pip install ngmaster
```

#### pixi (fast reproducible environment)

```bash
# Install pixi if not already installed
curl -fsSL https://pixi.sh/install.sh | bash

# Clone the repo and set up the environment
git clone https://github.com/MDU-PHL/ngmaster
cd ngmaster
pixi install
pixi run pip install -e .
```

#### Brew

```
# TODO how to integrate mlst dependency
brew install brewsci/bio/ngmaster
```

#### Conda/Mamba

```
conda install -c bioconda ngmaster
```

## Test

Once installed, you can run the following to ensure `ngmaster` is successfully working:

    ngmaster --test

If everything works, you will see the following:

```
Running ngmaster.py on test example (NG-MAST 4186 / NG-STAR 231) ...
FILE    SCHEME  NG-MAST/NG-STAR porB_NG-MAST    tbpB    penA    mtrR    porB_NG-STAR    ponA    gyrA    parC 23S      CC
test.fa     ngmaSTar        4186/231     2569     241     23      42      100     100     10      2       100     231
... Test successful.
```

## Usage

    ngmaster -h

    usage:
      ngmaster [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>

    In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)
    and Neisseria gonorrhoeae Sequence Typing for Antimicrobial Resistance (NG-STAR)

    Please cite as:
      Kwong JC, Goncalves da Silva A, Howden BP and Seemann T.
      NGMASTER: in silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)
      GitHub: https://github.com/MDU-PHL/ngmaster

    positional arguments:
      FASTA            input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN

    optional arguments:
      -h, --help       show this help message and exit
      --db DB          specify custom directory containing allele databases
                       directory must contain database sequence files (.tfa) and allele profile files (ngmast.txt / ngstar.txt)
                       in mlst format (see <https://github.com/tseemann/mlst#adding-a-new-scheme>)
                       default: <path/to/ngmaster/package/db>
      --csv            output comma-separated format (CSV) rather than tab-separated
      --printseq FILE  specify filename to save allele sequences to
      --minid MINID    DNA percent identity of full allele to consider 'similar' [~]
      --mincov MINCOV  DNA percent coverage to report partial allele at [?]
      --updatedb       update NG-MAST and NG-STAR allele databases from <https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef>
      --assumeyes      assume you are certain you wish to update db
      --test           run test example
      --comments       include NG-STAR comments for each allele in output
      --version        show program's version number and exit

## Quick start

**To perform *in silico* NG-MAST and NG-STAR typing on FASTA files:**

`$ ngmaster <fasta1> <fasta2> <fasta3> ... <fastaN>`

The NG-MAST and NG-STAR results and allele numbers are printed in tab-separated format to `stdout`.

* `ngmaster` reports alleles according to the same rules that are implemented in `mlst`.
* `mlst`'s arguments `--minid` and `--mincov` are available directly in `ngmaster` 
* For each allele n:

Symbol | Meaning | Length | Identity
---   | --- | --- | ---
`n`   | exact intact allele                   | 100%            | 100%
`~n`  | novel full length allele similar to n | 100%            | &ge; `--minid`
`n?`  | partial match to known allele         | &ge; `--mincov` | &ge; `--minid`
`-`   | allele missing                        | &lt; `--mincov` | &lt; `--minid`
`n,m` | multiple alleles                      | &nbsp;          | &nbsp;

**To save results to a tab-separated text file, redirect `stdout`:**

`$ ngmaster <fasta1> <fasta2> <fasta3> ... <fastaN>  > results.txt`

**To display results in comma-separated format, use the `--csv` option:**

`$ ngmaster --csv <fasta1> <fasta2> <fasta3> ... <fastaN>`

**To save sequences of the alleles to a file (eg. for uploading "new" sequences to [PubMLST](https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/)):**

`$ ngmaster --printseq [filename] <fasta1> <fasta2> <fasta3> ... <fastaN>`

This will create two files:

1. `NGMAST__filename`
2. `NGSTAR__filename`

## Updating the allele databases

**To update the allele databases from [PubMLST](https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/):**  
*Warning: This will overwrite the existing databases so ensure you back them up if you wish to keep them.*

    $ ngmaster.py --updatedb

A copy of the old database is saved just in case, but is overwritten with each subsequent   ```--updatedb```.

**To update the allele databases into a different folder (ie. not the /db folder in the ngmaster directory):**

    $ ngmaster.py --updatedb --db path/to/folder

This will download the database files into the folder ```path/to/folder```.
This can then be specified when running ngmaster using the ```--db  path/to/folder``` option.

**Check version and database version:**

```bash
ngmaster --version
```

## Updating the allele databases

The bundled database is the last freely redistributable snapshot of the PubMLST NG-MAST and NG-STAR schemes. The [PubMLST terms and conditions](https://pubmlst.org/terms-conditions) do not permit redistribution of database content after 2024-12-31 without appropriate licensing.

To update to the latest alleles, you need an authenticated PubMLST account and API credentials via `mlstdb`:

```bash
# Register and connect to PubMLST (one-time setup)
mlstdb connect -d pubmlst

# Update the bundled database
ngmaster --updatedb

# Or update into a custom directory
ngmaster --updatedb --db path/to/custom/db
```

> **Warning:** `--updatedb` overwrites existing database files. Back up your current `db/` directory first if needed.


## Creating a custom allele database

To create a custom allele database please follow the instructions for creating a custom ```mlst``` database
described [here](https://github.com/tseemann/mlst#adding-a-new-scheme).
Usually, this should not be necessary, simply run `ngmaster --update` to update to the latest NG-MAST and NG-STAR schemes from PubMLST.

## Citation

Kwong JC, Gonçalves da Silva A, Dyet K, Williamson DA, Stinear TP, Howden BP and Seemann T.  
*NGMASTER: in silico multi-antigen sequence typing for Neisseria gonorrhoeae*
**Microbial Genomics**
2016 Aug 25;2(8):e000076.
PMID:[28348871](https://www.ncbi.nlm.nih.gov/pubmed/28348871)
DOI:[10.1099/mgen.0.000076](https://doi.org/10.1099/mgen.0.000076)

## Bugs

Please submit via the [GitHub issues page](https://github.com/MDU-PHL/ngmaster/issues).  

## Software Licence

[GPLv2](https://github.com/MDU-PHL/ngmaster/blob/master/LICENSE)

## References

* [Martin et al. J Infect Dis, 2004 Apr 15; 189(8): 1497-1505](https://doi.org/10.1086/383047).
* [Demczuk et al. J Clin Microbiol, 2017 May; 55(5): 1454-1468](https://doi.org/10.1128/jcm.00100-17)
* See also [PubMLST](https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/).

## Authors

* Jason Kwong (@kwongjc)
* Anders Gonçalves da Silva (@drandersgs)
* Mark Schultz (@schultzm)
* Torsten Seemann (@torstenseemann)
* Andreas Stroehlein (@stroehleina)

## Development

When incrementing version (i.e., minor patch), run the following:

```
bumpversion --verbose --dry-run --new-version <major.minor.patch> patch
bumpversion --new-version <new.version.number> patch
```

The same can be done for `minor` and `major` numbers.

This will automatically commit and tag the commit with the new version number.
It will also update the necessary location in the file.

