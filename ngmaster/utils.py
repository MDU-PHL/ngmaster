# Modules and Functions
import sys
import os
import os.path
import shutil
import re
from Bio import SeqIO
import subprocess
import requests
from pkg_resources import resource_filename

class MlstRecord:
    '''A class that defines an NG-MAST or NG-STAR
    output record for a single FASTA file
    analysed by the mlst tool'''
    def __init__(self, fname, scheme, st, alleles):
        self.fname = fname
        self.scheme = scheme
        self.st = st
        self.alleles = alleles
        # self.porb = ["-"]
        # msg("self.alleles is of type ", str(type(self.alleles))) # a list
        if self.scheme == 'ngstar':
            self.porb = self.alleles[2]
        else:
            self.porb = self.alleles[0]

        # msg("self.porb ", self.porb)


        # self.simil = "" # similar "~"
        # self.part = "" # partial "?"

        # FIXME how to handle multi allele porB?
        # INFO if there is a multi hit then the approx hits won't be present
        # INFO if there are approximate hits then no multi hits will be present
        # TODO make self.porb a list (so we can do multi match to ngmast via ttable)
        # INFO for all other alleles this is not needed

        if re.match('~', self.porb):
            self.porb = self.porb[1:]
            self.simil = "~"
        if re.search(r'\?$', self.porb):
            self.porb = self.porb[:-1]
            self.part = "?"
        if re.search(',', self.porb):
            self.porb = self.porb.split(',')

        # TODO account for "^~", "?$" and multi alleles (comma separated list)
        # n       exact intact allele                     100%        100%
        # ~n      novel full length allele similar to n   100%        ≥ --minid
        # n?      partial match to known allele           ≥ --mincov  ≥ --minid
        # -       allele missing                          < --mincov  < --minid
        # n,m     multiple alleles         

    def get_record(self, sep):

        record = sep.join([self.fname, self.scheme, self.st] + [a for a in self.alleles])
        return record


# Log a message to stderr
def msg(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
    msg(*args, **kwargs)
    sys.exit(1)

# Update DB function
def update_db(db_folder, db):
    '''
    A function to update the database
    Has four inputs:
     db_folder: the base path for all DBs
     db: a dict with keys 'url' (the url to obtain the new db from) and 'db' (the filename to save to)
    '''

    #check if db folder exists and try to create it if not
    for dir in ["/blast","/pubmlst","/pubmlst/ngmast","/pubmlst/ngmast"]:
        try:
            if not os.path.exists(db_folder):
                os.makedirs(db_folder)
        except:
            err("Could not find/create db folder:'{}'".format( db_folder))

    #download and process db information
    try:
        if os.path.isfile(db['db']):
            shutil.copy(db['db'], db['db'] + '.old')

        try:
            pubmlst = requests.get(db['url']).text
        except requests.exceptions.RequestException as e:
            raise SystemExit(e)
        
        # Clean up names from PubMLST so they work well with mlst's mlst-make_blast_db
        pubmlst = pubmlst.replace('NG-MAST_','')
        pubmlst = pubmlst.replace('NG_','')
        pubmlst = pubmlst.replace('NEIS1753','penA')
        pubmlst = pubmlst.replace('\'mtrR','mtrR')
        new_db = pubmlst


    except:
        err("Unable to download/process URL:'{}'".format(db['url']))

    # save new db to db['db']
    with open(db['db'],'w') as f:
        f.write(new_db)
        
    msg(db['db'] + ' ... Done.')

# Create mlst database
def make_mlst_db(DBpath, mkblastdbpath):
    '''
    A function to run mlst's mlst-make_blast_db script to create custom ngstar and ngmast databases
    in the db folder (copy to scripts folder)
    Has two inputs:
    DBpath: base directory with blast and pubmlst/ngmast and pubmlst/ngstar folders
    mkblastdbpath: path to mlst-make_blast_db script
    '''

    if os.path.exists(mkblastdbpath):
        # msg('mlst installed (' + mlstpath.decode('utf-8').strip() + ') and mlst-make_blast_db (' + mkblastdbpath + ') found.')
        #copied make-blast-db == cpmbdb
        # cpmbdb = resource_filename(__name__, 'scripts/') + 'mlst-make_blast_db'
        cpmbdb = resource_filename(__name__, 'scripts/mlst-make_blast_db')
        #msg(cpmbdb)
        try:
            shutil.copy(mkblastdbpath, cpmbdb)
        except:
            err("Could not copy mlst-make_blast_db script to :'{}'".format(cpmbdb))

        try:
            subprocess.run(cpmbdb, check=True)
        except subprocess.CalledProcessError as e:
            err("Could not run mlst-make_blast_db script :'{}'".format(e))
            
    else:
        err('ERROR: Could not find mlst-make_blast_db script in ' + mkblastdbpath + '. Check mlst (https://github.com/tseemann/mlst) is installed correctly and in $PATH.')

# FIXME download comments for each allele
# #!/usr/bin/python3

# TODO you have never gone through the sequence files one by one?

# TODO match these to 'clean' names
# NG_porB 
# NG_ponA
# NG_parC
# NG_gyrA
# NG_23S
# NEIS1753
# 'mtrR

# porB 
# ponA
# parC
# gyrA
# 23S
# penA
# mtrR

#  # ls *.fa | 
#  # while read line; do head $line | 
#  # fa2tab | 
#  # awk -F"\t" '{split($1,a,"_"); printf a[1]; for(i=2; i< length(a); i++){printf "_"a[i]}; print "\t"a[length(a)]"\t"$0}'; done |
#  # sed 's/^>//1' |

# NG_ponA 1       >NG_ponA_1      AAAAACAACGGCGGGCGTTGGGCGGTGGTTCAAGAGCCGTTGCCGCAGGGGGCTTTGGTTTCGCTGGATGCAAAA
# NG_ponA 2       >NG_ponA_2      AAAAACAACGGCGGGCGTTGGGCGGGGGTTCAAGAGCCGTTGCTGCAGGGGGCTTTGGTTTCGCTGGATGCAAAA
# [...]

# import sys
# import json
# import requests

# for line in sys.stdin:
#    # NG_gyrA 1       >NG_gyrA_1      CTGTACGCGAT
#    INFO locus (column 1)
#    INFO allele_id (column 2)
#    locus, allele_id, fasta_header, sequence = line.rstrip().split("\t")
#    recordurl = 'https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/' + locus + '/alleles/' + allele_id
#    pubmlst = ""

#    try:
#       pubmlst = requests.get(recordurl).text
#    except requests.exceptions.RequestException as e:
#       raise SystemExit(e)

#    pub_dict = json.loads(pubmlst)
#    comments = ""

#    try:
#       comments = pub_dict["comments"]
#    except KeyError as e:
#       comments = ""

#    output = '\t'.join([locus, allele_id, sequence, str(len(sequence)), comments])

#    # get this into same format as PubMLST online UI output
#    # output = record



# Match NGSTAR with NGMAST porB sequences
def match_porb(ngmast_porb, ngstar_porb):
    '''
    A function that matches short porB sequences from NG-STAR (30nt) to long porB sequences from NG-MAST
    Has two inputs:
    ngmast_porb: Path to ngmast porB fasta file
    ngstar_porb: Path to ngstar porB fasta file
    Returns a dict with ngmast porB IDs as keys and NGSTAR porB IDs as values: ttable
    '''
    # msg(ngmast_porb)
    # msg(ngstar_porb)

    if os.path.isfile(ngmast_porb) and os.path.isfile(ngstar_porb):
        ttable = {}
        ttable["-"] = "-"

        for mast in SeqIO.parse(ngmast_porb, "fasta"):
            mastid = re.sub(r"^porB_", "", str(mast.id))
            ttable[mastid] = "-"

            for star in SeqIO.parse(ngstar_porb, "fasta"):
                starid = re.sub(r"^porB_", "", str(star.id))

                if not (mast.seq.find(star.seq) == -1):
                    ttable[mastid] = starid
        return ttable
        
    else:
        err('ERROR: porB files for NG-STAR or NG-MAST could not be found. Run ngmaster.py --update to update PubMLST databases.')


def read_ngstar(ngsfile):
    '''
    Function that reads the ngstar.txt profile file and returns
    a dict of NGSTAR 7-allele tuples as keys and STs as value
    ST  penA    mtrR    porB    ponA    gyrA    parC    23S
    1   20  10  13  100 100 2   100
    2   22  14  14  100 100 7   100
    3   166 1   13  1   5   1   100
    '''
    ngstable = {}    

    with open(ngsfile,'r') as f:
        for line in f:
            cols = line.rstrip().split("\t")
            tpl_key = tuple(cols[1:])  
            st = cols[0]
            ngstable[tpl_key] = st

    return ngstable


def convert_ngstar(ngstable, mlstngstar):
    
    alleles = tuple(mlstngstar.alleles)
    st = ngstable.get(alleles, "-")

    converted = MlstRecord(mlstngstar.fname, mlstngstar.scheme, st, mlstngstar.alleles)

    return converted


def collate_results(ngmast_res, ngstar_res, ttable, ngstartbl):

    '''
    A function that collates results from running ngmast mlst and ngstar mlst and returns a combined output
    Has three inputs:
    ngmast_res: results from running mlst --legacy --scheme ngmast, a dict with filenames as keys and MlstRecord objects as values
    ngstar_res: results from running mlst --legacy --scheme ngstar, a dict with filenames as keys and MlstRecord objects as values
    ttable: translation table that matches up short porB sequences from NG-STAR (30nt) to long porB sequences from NG-MAST
            created by match_porb(), a dict with ngmast / key -> ngstar / value pairs
    '''

    combined_res = []
    # check that keys for both dicts are the same
    if ngmast_res.keys() == ngstar_res.keys():
        for file in ngmast_res:

            ngstar_res[file].alleles[2] = ttable[ngmast_res[file].porb]
            conv_ngs = convert_ngstar(ngstartbl, ngstar_res[file])

            # check = [a for a in ngmast_res[file].alleles] + [a for a in conv_ngs.alleles]
            # msg(check)
            # msg(type(check))

            ngmastar = MlstRecord(ngmast_res[file].fname,
                "ngmaSTar",
                ngmast_res[file].st + "/" + conv_ngs.st,
                # FIXME a list when otherwise it is a dict ['porB', 'tbpB', 'penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S']
                [a for a in ngmast_res[file].alleles] + [a for a in conv_ngs.alleles]
            )

            combined_res.append(ngmastar)

    else:
        raise KeyError("Not all files have been successfully processed by mlst (ngstar or ngmast scheme). Cannot collate results from two runs. Exiting.")

    return combined_res
