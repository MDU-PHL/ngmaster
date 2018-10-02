# Modules and Functions
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import os.path
import io
import urllib.request, urllib.error, urllib.parse
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
import subprocess
import shutil
import re
from sys import argv
from subprocess import Popen
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests
from bs4 import BeautifulSoup

# Log a message to stderr
def msg(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
    msg(*args, **kwargs)
    sys.exit(1);

# Import allele database as dictionary
def fasta_to_dict(db):
    dict = {}
    dbFILE = open((db), 'rU')
    for allele in SeqIO.parse(dbFILE, 'fasta'):
        dict[str(allele.seq)] = (allele.id)
        dict[str(allele.seq.reverse_complement())] = (allele.id)
    return dict

# Update DB function
def update_db(db_folder, db_file, db_url, allele_db = False):
    '''
    A function to update the database
    Has four inputs:
     db_folder: a path to save the db
     db_file: a path to save new db files
     db_url: the url to obtain the new db from
     allele_db: True/False. If downloading the sequence db this should be false. If downloading the allele_db then this should be True
    '''
    #check if db folder exists and try to create it if not
    try:
        if not os.path.exists(db_folder):
            os.makedirs(db_folder)
    except:
        err("Could not find/create db folder:'{}'".format( db_folder))
    #download and process db information
    try:
        if os.path.isfile(db_file):
            shutil.copy(db_file, db_file+'.old')
        page = requests.get(db_url)
        if not allele_db:
            # all alleles are in a single tag <textarea>, and lines are
            # separated by '\r\n'
            # the first line will decode the tag area, splitting on the
            # separator, creating a list'
            # the second line joins the elements of the list to create
            # a single string for saving to file.
            soup = BeautifulSoup(page.text, "html.parser")
            new_db = soup.select("textarea")[0].text.split("\r\n")
            new_db = '\n'.join(new_db[0:-1]) 
            # the -1 is to remove the last empty field
        else:
            # page is just a comma separated text file, with <br> tags
            # for newlines
            # need to translate the </br> tags to newline character, and then
            # the text comes out ready to save to file
            new_db = page.text.split("<br>")
            new_db = '\n'.join(new_db)
    except:
        err("Unable to download/process URL:'{}'".format(db_url))
    # save new db to db_file
    try:
        f = open(db_file,'w')
        if not allele_db:
            f.write(new_db)
        else:
            f.write('ST'+','+'POR'+','+'TBPB'+'\n'+new_db)
        f.close()
    except:
        err("Could not save db file: '{}'".format(db_file))
    msg(db_file + ' ... Done.')