# Modules and Functions
import sys
import os
import os.path
import shutil
import re
from Bio import SeqIO
import subprocess
import requests
import json
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
        self.simil = "" # similar "~"
        self.part = "" # partial "?"

        if self.scheme == 'ngstar':
            self.porb = self.alleles[2]
        else:
            self.porb = self.alleles[0]

        # FIXME how to handle multi allele porB?
        # INFO if there is a multi hit then the approx hits won't be present
        # INFO if there are approximate hits then no multi hits will be present
        # TODO make self.porb a list (so we can do multi match to ngmast via ttable)
        # TODO to get comments for others we might also want to allow multiple alleles
        # INFO for all other alleles this is not needed
        # FIXME only allow comments to be pulled in when allele is an exact match

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

    def get_record(self, sep, comments):
        '''
        Function that returns a comma-separated or tab-separated string of allele IDs
        (and optionally associated PubMLST comments) for each record
        '''

        ngmast = self.alleles[:2]
        ngstar = self.alleles[2:]

        if comments:
            rec_comms = []

            #Retrieve record-specific comments
            for locus, allele in enumerate(ngstar):
                try:
                    rec_comms.append(comments[locus].get(allele,""))
                except IndexError as e:
                    raise SystemExit(e)
                    err(f'error {e}: locus is {locus} when comments is of length {len(comments)} and last index {comments[-1]} and allele {allele} for file {self.fname} and scheme {self.scheme}')

                
            # Interleave comments with allele IDs
            ngstar = list(sum(zip(ngstar,rec_comms), ()))

        record = ngmast + ngstar

        if sep == ',':
            csv_record = []
            for col in record:
                if re.search(r',', col):
                    col = '"' + col + '"'
                csv_record.append(col)
            record = csv_record

        joined_record = sep.join([self.fname, self.scheme, self.st] + record)

        return joined_record


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
    for dir in ["/blast","/pubmlst","/pubmlst/ngmast","/pubmlst/ngstar"]:
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

    with open(db['db'],'w') as f:
        f.write(new_db)
        
    msg(db['db'] + ' ... Done.')


def download_comments(DBpath, db_list):
    '''
    Function to download the comments for individual NG-STAR alleles
    '''

    msg("Downloading NG-STAR information for individual alleles. This may take a while.")
    comms_file = []

    for db in db_list:
        if db["comments"]:
            for seq in SeqIO.parse(db['db'], "fasta"):

                locus, allele = str(seq.id).split("_") 
                recordurl = 'https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/' + db["comments"] + '/alleles/' + allele

                try:
                    comms = requests.get(recordurl).text
                except requests.exceptions.RequestException as e:
                    raise SystemExit(e)

                comm_dict = json.loads(comms)
                allele_comm = comm_dict.get("comments", "")

                comms_file.append("\t".join([locus, allele, allele_comm]))

    with open(DBpath + "/pubmlst/ngstar/allele_comments.tsv",'w') as f:
        for line in comms_file:
            f.write(line + "\n")

def load_ngstar_comments(DBpath):
    '''
    Function that returns a list of comment dicts when --comments is set
    '''

    with open(DBpath + "/pubmlst/ngstar/allele_comments.tsv",'r') as f:
        d_comms = {"penA":{}, "mtrR":{}, "porB":{}, "ponA":{}, "gyrA":{}, "parC":{}, "23S":{}}

        ngstar_comments = []

        for line in f:
            cols = line.rstrip("\n").split("\t")
            d_comms[cols[0]][cols[1]]=cols[2]

        for a in ["penA", "mtrR", "porB", "ponA", "gyrA", "parC", "23S"]:
            ngstar_comments.append(d_comms[a])        
    
    return ngstar_comments


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

        cpmbdb = resource_filename(__name__, 'scripts/mlst-make_blast_db')

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


# Match NGSTAR with NGMAST porB sequences
def match_porb(ngmast_porb, ngstar_porb):
    '''
    A function that matches short porB sequences from NG-STAR (30nt) to long porB sequences from NG-MAST
    Has two inputs:
    ngmast_porb: Path to ngmast porB fasta file
    ngstar_porb: Path to ngstar porB fasta file
    Returns a dict with ngmast porB IDs as keys and NGSTAR porB IDs as values: ttable
    '''

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
        err('ERROR: porB files for NG-STAR or NG-MAST could not be found. Run ngmaster --update to update PubMLST databases.')


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
            cols = line.rstrip('\n').split("\t")
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

            #msg(f"Post-processing file {file} and creating output ...")

            try:
                ngstar_res[file].alleles[2] = "".join([ngmast_res[file].simil, ttable[ngmast_res[file].porb], ngmast_res[file].part])
            except TypeError as e:
                #This gets triggered when porb is a list rather than a single value
                # Go through all elements and create a multi-allele entry using ttable
                if isinstance(ngmast_res[file].porb, list):
                    porblist = []
                    for pb in ngmast_res[file].porb:
                        porblist.append(ttable[pb])
                    ngstar_res[file].alleles[2] = ",".join(porblist)
                else:
                    err(f"Porb allele for record {file} is not list but {type(ngmast_res[file].porb)}. Exiting ...")
                    raise SystemExit

            conv_ngs = convert_ngstar(ngstartbl, ngstar_res[file])

            ngmastar = MlstRecord(ngmast_res[file].fname,
                "ngmaSTar",
                ngmast_res[file].st + "/" + conv_ngs.st,
                ngmast_res[file].alleles + conv_ngs.alleles
            )

            combined_res.append(ngmastar)

    else:
        raise KeyError("Not all files have been successfully processed by mlst (ngstar or ngmast scheme). Cannot collate results from two runs. Exiting.")

    return combined_res
