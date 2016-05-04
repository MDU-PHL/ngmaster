#ngmaster

*In silico* multi-antigen sequence typing for *Neisseria gonorrhoeae* (NG-MAST).  

##Authors

* Jason Kwong (@kwongjc)
* Torsten Seemann (@torstenseemann)

##Dependencies

* BioPython
* isPcr

##Installation

The easiest way of installing `ngmaster` is using `pip`:

    pip install --user git+https://github.com/MDU-PHL/ngmaster.git
    
The `--user` option will install the package locally, rather than in the global `python` directory. 

Thus, by default, this will install the package in `$HOME/.local/`, and the executable in `$HOME/.local/bin/`. To install the executable in a custom location (e.g., `$HOME/bin`), use the following:

    pip install --install-option="--install-scripts=$HOME/bin" --user git+https://github.com/MDU-PHL/ngmaster.git

To upgrade to a newer version: 

    pip install --upgrade --install-option="--install-scripts=$HOME/bin" --user git+https://github.com/MDU-PHL/ngmaster.git

The simplest way to install dependencies is to use the Brew (Mac OS X) or LinuxBrew (Linux) system.
```
brew tap homebrew/science
brew tap chapmanb/cbl
brew tap tseemann/homebrew-bioinformatics-linux
```

### To test installation

Once installed, you can run the following to ensure `ngmaster` is successfully working:

    ngmaster --test

If everything works, you will see the following:

```
Running ngmaster.py on test example (NG-MAST 10699) ...
$ ngmaster.py test/test.fa
ID	NG-MAST	POR	TBPB
test.fa	10699	6277	4
... Test successful.
```

##Usage

	$ ngmaster -h
        
	usage: 
	  ngmaster [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>
	
	In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)
	
	Ref: Martin et al. J Infect Dis, 2004 Apr 15;189(8):1497-1505.
	See also http://www.ng-mast.net/
	
	positional arguments:
	  FASTA            input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN
	
	optional arguments:
	  -h, --help       show this help message and exit
	  --db DB          specify custom directory containing allele databases
	                   directory must contain database files "POR.tfa", "TBPB.tfa", and "ng_mast.txt"
	  --csv			   output comma-separated format (CSV) rather than tab-separated
	  --printseq FILE  specify filename to save allele sequences to (default=off)
 	  --updatedb       update allele database from <www.ng-mast.net>
	  --test		   run test example
	  --version        show program's version number and exit


##Quick start

**To perform *in silico* NG-MAST on FASTA files:**

`$ ngmaster <fasta1> <fasta2> <fasta3> ... <fastaN>`

##Updating the allele databases

**To update the allele databases from http://www.ng-mast.net :**  
*Warning: This will overwrite the existing databases so ensure you back them up if you wish to keep them.*

	$ ngmaster.py --updatedb

A copy of the old database is saved just in case, but is overwritten with each subsequent   ```--updatedb```.

**To update the allele databases into a different folder (ie. not the /db folder in the ngmaster directory):**

	$ ngmaster.py --updatedb --db path/to/folder

This will download the database files into the folder ```path/to/folder```.
This can then be specified when running ngmaster using the ```--db  path/to/folder``` option.

##Creating a custom allele database

1. Create custom database files: `POR.tfa`, `TBPB.tfa`, `ng_mast.txt`  
   See default `db` directory for examples.  
   `POR.tfa` and `TBPB.tfa` contain the respective allele sequences in FASTA format.  
   `ng_mast.txt` contains a list of NG-MAST types and the corresponding allele types.

2. Place the custom database files in a folder.

3. Specify the path to that custom database folder:  
   `$ ngmaster --db [/path/to/custom/folder/] <fasta1> <fasta2> <fasta3> ... <fastaN>`

**To save sequences of the alleles to a file (eg. for uploading to [http://www.ng-mast.net](http://www.ng-mast.net/)):**

`$ ngmaster --printseq [filename] <fasta1> <fasta2> <fasta3> ... <fastaN>`

##Citation

Please cite as:

Kwong JC, Goncalves da Silva A, Howden BP and Seemann T.  
*NGMASTER: in silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)*  
GitHub: https://github.com/MDU-PHL/ngmaster  

##Bugs

Please submit via the [GitHub issues page](https://github.com/MDU-PHL/ngmaster/issues)  

Note that the NG-MAST databases and website are curated and hosted at the Department of Infectious Diseases Epidemiology, Imperial College London. For issues with the NG-MAST databases, please contact the [NG-MAST curator](mailto:d.aanensen@imperial.ac.uk?subject=N. gonorrhoeae homepage).

##Software Licence

[GPLv2](https://github.com/MDU-PHL/ngmaster/blob/master/LICENSE)

##References

* Martin et al. J Infect Dis, 2004 Apr 15; 189(8): 1497-1505.  
* See also [http://www.ng-mast.net](http://www.ng-mast.net/).
