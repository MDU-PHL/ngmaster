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
import logging
import tempfile
from mlstdb.core.auth import get_client_credentials, retrieve_session_token
from rauth import OAuth1Session, OAuth1Service
from tqdm import tqdm

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
        self.cc = ''  # Add CC field

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

 
    def get_record(self, sep='\t', comments=None, header=None): # FLAG: modified function to deal with wrong number of fields in output
        '''
        Function that returns a comma-separated or tab-separated string of allele IDs
        (and optionally associated PubMLST comments) for each record
        '''
        
        ngmast_loci = self.alleles[:2]
        ngstar_loci = self.alleles[2:]
        # print(f"DEBUG: ngmast_loci: {ngmast_loci}, ngstar_loci: {ngstar_loci}") # FLAG: for debugging
        if comments:
            out = [self.fname, self.scheme, self.st] + self.alleles
            # Interleave alleles and comments
            out_with_comments = []
            # print(f"DEBUG: out: {out}") # FLAG: for debugging
            for i, allele in enumerate(ngstar_loci):  # Assuming 7 loci
                # Get the corresponding header if provided
                column_header = header[3 + 2 * i] if header else f"Allele_{i + 1}"
                comment_header = header[4 + 2 * i] if header else f"Comment_{i + 1}"
                comment = comments[i].get(allele, "")
                # Debugging: Log the header, allele, and comment
                # print(f"DEBUG: {column_header}: {allele}, {comment_header}: {comment}")
                out_with_comments.append(allele)
                out_with_comments.append(comment)
            out = out[:5] + out_with_comments
            out.append(self.cc)  # Add CC to output
        else:
            out = [self.fname, self.scheme, self.st] + self.alleles + [self.cc]  # Add CC to output

        # **Handle special characters for CSV output**
        if sep == ',':
            csv_record = []
            for col in out:
                if isinstance(col, str) and (',' in col or ';' in col):
                    col = f'"{col}"'  # Quote fields containing special characters
                csv_record.append(col)
            out = csv_record

        # Debugging: Log the final row before returning
        # print(f"DEBUG: Final row: {out}")

        return sep.join(map(str, out))

# Log a message to stderr
def msg(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
    msg(*args, **kwargs)
    sys.exit(1)


class PubMLSTAuth: 
    """Class to handle PubMLST authentication state and requests"""
    
    def __init__(self):
        self.auth_tested = False
        self.use_auth = False
    
    def get_response(self, url):
        """
        Get response from URL, using OAuth if possible
        Falls back to unauthenticated request if OAuth fails
        """
        # Test auth only once
        if not self.auth_tested:
            try:
                # Get OAuth credentials from mlstdb
                client_key, client_secret = get_client_credentials('pubmlst')
                session_token, session_secret = retrieve_session_token('pubmlst')
                
                # Initialize OAuth session
                self.session = OAuth1Session(
                    consumer_key=client_key,
                    consumer_secret=client_secret,
                    access_token=session_token,
                    access_token_secret=session_secret,
                )
                self.session.headers.update({"User-Agent": "BIGSdb downloader"})
                self.use_auth = True
                
            except Exception as e:
                msg(f"OAuth setup failed: {e}")
                msg("Please run 'mlstdb fetch' to set up/refresh OAuth credentials")
                msg("Falling back to unauthenticated requests")
                
                self.use_auth = False
                
            self.auth_tested = True
        
        try:
            if self.use_auth:
                response = self.session.get(url)
            else:
                response = requests.get(url)
                
            # response.raise_for_status()
            return response
            
        except Exception as e:
            msg(f"Request failed: {e}")
            if self.use_auth:
                msg("Retrying without authentication")
                response = requests.get(url)
                return response
            raise

# Update DB function FLAG: fix oauth response
def update_db(db_folder, db):
    '''
    A function to update the database with OAuth authentication via mlstdb
    '''
    
    auth_handler = PubMLSTAuth()
    
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
            # print debugging info
            # print("DEBUGGING")
            print(f"Dowloading from URL: {db['url']}")
            
            # Download the main data
            pubmlst = auth_handler.get_response(db['url']).text
            
            # if `schemes` string present in db['url'], then it is a profile URL, get database_version info as well
            if 'schemes' in db['url']:
                print("Downloading database version info")
                scheme_url = db['url'].replace('/profiles_csv', '')
                print(f"Scheme URL: {scheme_url}")
                
                version_response = auth_handler.get_response(scheme_url)
                # version_response.raise_for_status()
                scheme_data = version_response.json()
                
                db_version = scheme_data.get('last_added', 'Not found')
                if db_version is None:
                    db_version = 'No version information available'
                if auth_handler.use_auth == False:
                    db_version = f"{db_version} (Unauthenticated)"
                    
                print(f"Database version: {db_version}")
                
                # Save database version to a file in the appropriate directory
                scheme_type = 'ngstar' if '67' in scheme_url else 'ngmast'
                db_version_path = os.path.join(db_folder, 'pubmlst', scheme_type, 'database_version.txt')
                with open(db_version_path, 'w') as version_file:
                    version_file.write(db_version)
                print(f"Saved version info to: {db_version_path}")  
                
        except requests.exceptions.RequestException as e:
            raise SystemExit(e)
        
        # print(pubmlst)
        
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
    msg("Starting download_comments...This may take a while (few hours) depending on the number of alleles.")
    msg(f"Base directory: {DBpath}")
    
    auth_handler = PubMLSTAuth()

    comments_file_path = os.path.join(DBpath, "pubmlst/ngstar/allele_comments.tsv")
    msg(f"Output file: {comments_file_path}")

    comms_file = []
    
    # Count total sequences first for progress bar
    total_sequences = sum(len(list(SeqIO.parse(db['db'], "fasta"))) 
                         for db in db_list if db["comments"])

    progress_bar = tqdm(total=total_sequences, desc="Downloading comments")

    for db in db_list:
        if db["comments"]:
            msg(f"Processing database: {db['db']} with comments URL base: {db['comments']}")

            for seq in SeqIO.parse(db['db'], "fasta"):
                locus, allele = str(seq.id).split("_")
                recordurl = f'https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/{db["comments"]}/alleles/{allele}'
                msg(f"Fetching URL: {recordurl}")

                try:
                    comms = auth_handler.get_response(recordurl).text
                except requests.exceptions.RequestException as e:
                    err(f"Failed to fetch URL: {recordurl} with error: {e}")

                comm_dict = json.loads(comms)
                allele_comm = comm_dict.get("comments", "")

                comms_file.append("\t".join([locus, allele, allele_comm]))
                progress_bar.update(1)

    progress_bar.close()

    msg(f"Writing comments to file: {comments_file_path}")
    os.makedirs(os.path.dirname(comments_file_path), exist_ok=True)

    with open(comments_file_path, 'w') as f:
        for line in comms_file:
            f.write(line + "\n")

    msg("Finished download_comments.")
    
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


def read_ngstar(ngstar_profile):
    """Read NG-STAR profiles"""
    ngstar = {}
    with open(ngstar_profile) as f:
        header = f.readline().rstrip().split('\t')
        for line in f:
            tokens = line.rstrip().split('\t')
            if len(tokens) >= 8:  # Require at least the base columns (ST + 7 alleles)
                st = tokens[0]
                # Get the alleles
                alleles = tokens[1:8]
                # Get CC if present, otherwise empty string
                cc = tokens[8] if len(tokens) > 8 else ''
                ngstar["/".join(alleles)] = (st, cc)  # Store CC with ST
    return ngstar

def convert_ngstar(ngstartbl, rec):
    """Convert NG-STAR record using profile table"""
    key = "/".join(rec.alleles)
    if key in ngstartbl:
        rec.st = ngstartbl[key][0]  # ST is first element
        rec.cc = ngstartbl[key][1]  # CC is second element
    else:
        rec.st = "-"
        rec.cc = "-"
    return rec

def collate_results(ngmast_res, ngstar_res, ttable, ngstartbl):
    
    '''
    A function that collates results from running ngmast mlst and ngstar mlst and returns a combined output
    Has three inputs:
    ngmast_res: results from running mlst --legacy --scheme ngmast
    ngstar_res: results from running mlst --legacy --scheme ngstar
    ttable: translation table that matches porB sequences
    ngstartbl: NG-STAR profile table for ST and CC lookup
    '''
    
    combined_res = []
    
    if ngmast_res.keys() == ngstar_res.keys():
        for file in ngmast_res:
            
            #msg(f"Post-processing file {file} and creating output ...")
            
            try:
                # Handle porB modification
                if isinstance(ngmast_res[file].porb, list):
                    porblist = []
                    for pb in ngmast_res[file].porb:
                        porblist.append(ttable[pb])
                    ngstar_res[file].alleles[2] = ",".join(porblist)
                else:
                    ngstar_res[file].alleles[2] = "".join([
                        ngmast_res[file].simil,
                        ttable[ngmast_res[file].porb],
                        ngmast_res[file].part
                    ])

                # Convert NG-STAR record and get both ST and CC
                conv_ngs = convert_ngstar(ngstartbl, ngstar_res[file])
                
                # Create combined record with CC included
                ngmastar = MlstRecord(
                    ngmast_res[file].fname,
                    "ngmaSTar",
                    f"{ngmast_res[file].st}/{conv_ngs.st}",
                    ngmast_res[file].alleles + conv_ngs.alleles
                )
                
                # Set the CC value from the converted NG-STAR record
                ngmastar.cc = conv_ngs.cc
                
                combined_res.append(ngmastar)
                
            except TypeError as e:
                if not isinstance(ngmast_res[file].porb, list):
                    print(f"Porb allele for record {file} is not list but {type(ngmast_res[file].porb)}. Exiting ...")
                    raise SystemExit
                raise e
    else:
        raise KeyError("Not all files have been successfully processed by mlst (ngstar or ngmast scheme). Cannot collate results from two runs. Exiting.")
        
    return combined_res

def get_db_version(DBpath):
    """
    Get database versions for NG-MAST and NG-STAR schemes
    Returns a formatted string combining both versions
    """
    try:
        # Read NG-MAST version
        ngmast_version_path = os.path.join(DBpath, 'pubmlst', 'ngmast', 'database_version.txt')
        with open(ngmast_version_path, 'r') as f:
            ngmast_version = f.read().strip()

        # Read NG-STAR version
        ngstar_version_path = os.path.join(DBpath, 'pubmlst', 'ngstar', 'database_version.txt')
        with open(ngstar_version_path, 'r') as f:
            ngstar_version = f.read().strip()

        # Format the combined version string
        return f"ngmast_{ngmast_version}_ngstar_{ngstar_version}"
    except FileNotFoundError:
        return "Database version not available. Run --updatedb first"
    except Exception as e:
        return f"Error reading database versions: {str(e)}"