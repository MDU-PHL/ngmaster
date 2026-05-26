# Modules and Functions
import ngmaster
import sys
import os
import os.path
import shutil
import re
from Bio import SeqIO
import subprocess
import requests
import json
from pathlib import Path
import logging
import tempfile
from mlstdb.core.auth import get_client_credentials, retrieve_api_key, retrieve_session_token
from rauth import OAuth1Session
from tqdm import tqdm
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

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
    """Class to handle PubMLST authentication state and requests.

    Auth priority (determined once at construction, no deferred logic):
      1. Personal API key stored by mlstdb (~/.config/mlstdb/api_keys)
      2. OAuth session tokens stored by mlstdb
      3. Unauthenticated (fallback)
    """

    def __init__(self):
        self.auth_mode = None  # 'api_key' | 'oauth' | None

        # 1. Try API key first (BIGSdb >= v1.53.0, recommended)
        api_key = retrieve_api_key('pubmlst')
        if api_key:
            self.session = requests.Session()
            self.session.headers.update({
                "User-Agent": f"ngmaster/{ngmaster.__version__}",
                "X-API-Key": api_key,
            })
            self.auth_mode = 'api_key'
            return

        # 2. Try OAuth credentials stored by mlstdb
        try:
            client_key, client_secret = get_client_credentials('pubmlst')
            session_token, session_secret = retrieve_session_token('pubmlst')
            if session_token and session_secret:
                self.session = OAuth1Session(
                    consumer_key=client_key,
                    consumer_secret=client_secret,
                    access_token=session_token,
                    access_token_secret=session_secret,
                )
                self.session.headers.update({"User-Agent": f"ngmaster/{ngmaster.__version__}"})
                self.auth_mode = 'oauth'
                return
        except Exception:
            pass

        # 3. Unauthenticated fallback
        msg("No authentication credentials found. Falling back to unauthenticated requests.")
        msg("For authenticated access, run: mlstdb connect --db pubmlst --api-key  (recommended)")
        msg("                           or: mlstdb connect --db pubmlst  (OAuth)")

    def get_response(self, url):
        """
        Get response from URL using the best available authentication.
        Falls back to unauthenticated on auth failure.
        """
        try:
            if self.auth_mode == 'api_key':
                response = self.session.get(url)
            elif self.auth_mode == 'oauth':
                # params={} is required to avoid a rauth 0.7.3 bug where
                # _parse_optional_params raises TypeError when params is None
                response = self.session.get(url, params={})
            else:
                response = requests.get(url)
            return response

        except Exception as e:
            msg(f"Request failed: {e}")
            if self.auth_mode in ('api_key', 'oauth'):
                msg("Retrying without authentication")
                return requests.get(url)
            raise


_thread_local = threading.local()


def _get_thread_auth():
    """Return a PubMLSTAuth instance local to the current thread."""
    if not hasattr(_thread_local, 'auth'):
        _thread_local.auth = PubMLSTAuth()
    return _thread_local.auth

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
                if auth_handler.auth_mode is None:
                    db_version = "2024-12-31_Unauthenticated"
                    
                print(f"Database version: {db_version}")
                
                # Save database version to a file in the appropriate directory
                scheme_type = 'ngstar' if '67' in scheme_url else 'ngmast'
                db_version_path = os.path.join(db_folder, 'pubmlst', scheme_type, 'database_version.txt')
                with open(db_version_path, 'w') as version_file:
                    version_file.write(f"{db_version}\n")
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
    
def download_comments(DBpath, db_list, threads=1):
    '''
    Function to download the comments for individual NG-STAR alleles.
    `threads` controls the number of concurrent HTTP requests to PubMLST
    (default 1; recommended 4 to balance speed and server load).
    '''
    msg(f"Starting download_comments with {threads} thread(s)...This may take a while depending on the number of alleles.")
    msg(f"Base directory: {DBpath}")

    comments_file_path = os.path.join(DBpath, "pubmlst/ngstar/allele_comments.tsv")
    msg(f"Output file: {comments_file_path}")

    # Pre-collect all work items in original order: (index, locus, allele, url)
    work_items = []
    for db in db_list:
        if db["comments"]:
            msg(f"Processing database: {db['db']} with comments URL base: {db['comments']}")
            for seq in SeqIO.parse(db['db'], "fasta"):
                locus, allele = str(seq.id).split("_")
                recordurl = (
                    f'https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef'
                    f'/loci/{db["comments"]}/alleles/{allele}'
                )
                work_items.append((len(work_items), locus, allele, recordurl))

    progress_bar = tqdm(total=len(work_items), desc="Downloading comments")

    def _fetch(item):
        idx, locus, allele, url = item
        auth = _get_thread_auth()
        try:
            comms = auth.get_response(url).text
        except requests.exceptions.RequestException as e:
            msg(f"Failed to fetch URL: {url} with error: {e}")
            comms = '{}'
        comm_dict = json.loads(comms)
        allele_comm = comm_dict.get("comments", "")
        return idx, locus, allele, allele_comm

    results = {}
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_fetch, item): item for item in work_items}
        for future in as_completed(futures):
            idx, locus, allele, allele_comm = future.result()
            results[idx] = (locus, allele, allele_comm)
            progress_bar.update(1)

    progress_bar.close()

    msg(f"Writing comments to file: {comments_file_path}")
    os.makedirs(os.path.dirname(comments_file_path), exist_ok=True)

    with open(comments_file_path, 'w') as f:
        for idx in range(len(work_items)):
            locus, allele, allele_comm = results[idx]
            f.write("\t".join([locus, allele, allele_comm]) + "\n")

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

def process_duplicate_23s_alleles(rlist):
    """
    Process 23S alleles in ngstar scheme to handle duplicates.
    """
    if not rlist or len(rlist) < 2:  # No data or just headers
        return rlist
        
    # Get header and find the 23S column
    header = rlist[0].split("\t")
    
    # If no 23S column, nothing to process
    if "23S" not in header:
        return rlist
    
    # Find the column index for 23S
    s23_index = header.index("23S")
    
    # Process each result row (skip header)
    for i in range(1, len(rlist)):
        fields = rlist[i].split("\t")
        
        # Skip if row doesn't have enough fields
        if len(fields) <= s23_index:
            continue
            
        # Get filename for warning message
        filename = fields[0]
        
        # Check if 23S field contains commas (duplicates)
        s23_value = fields[s23_index]
        if "," in s23_value:
            s23_alleles = s23_value.split(",")
            unique_alleles = list(set(s23_alleles))
            
            # If we found duplicates
            if len(unique_alleles) < len(s23_alleles):
                # If only one unique value, use that
                if len(unique_alleles) == 1:
                    new_value = unique_alleles[0]
                else:
                    # Otherwise, join unique values with commas
                    new_value = ",".join(unique_alleles)
                
                # Replace value in the row
                fields[s23_index] = new_value
                rlist[i] = "\t".join(fields)
                
                # Print warning message
                msg(f"WARNING: Duplicate 23S alleles detected in {filename} ({s23_value}). "
                      f"Squashed to {new_value}.")
    
    return rlist