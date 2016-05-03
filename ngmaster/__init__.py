#!/usr/bin/env python
# Script by Jason Kwong & Torsten Seemann
# In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)

# Use modern print function from python 3.x
from __future__ import print_function

# Modules and Functions
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import os.path
import StringIO
import urllib2
from urllib2 import urlopen, HTTPError, URLError
import subprocess
import shutil
import re
import tempfile
from sys import argv
from subprocess import Popen
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

# Perform the pure-Python equivalent of in-place `sed` substitution: e.g.,
# `sed -i -e 's/'${pattern}'/'${repl}' "${filename}"`.
def sed_inplace(filename, pattern, repl):
	pattern_compiled = re.compile(pattern)
	with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
		with open(filename) as src_file:
			for line in src_file:
				tmp_file.write(pattern_compiled.sub(repl, line))
	shutil.copystat(filename, tmp_file.name)
	shutil.move(tmp_file.name, filename)

# sed with . matching all characters including newline
def sed_all(filename, pattern, repl):
	pattern_compiled = re.compile(pattern, re.DOTALL)
	with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
		with open(filename, 'r') as src_file:
			d = src_file.read()
			tmp_file.write(pattern_compiled.sub(repl, d))
	shutil.copystat(filename, tmp_file.name)
	shutil.move(tmp_file.name, filename)

# Format database
def format(file):
	sed_all(file, '^.*<textarea name="concatenation".*?>', '')
	sed_all(file, '<\/textarea>.*$', '')
	sed_inplace(file, '&gt;', '>')

# Import allele database as dictionary
def fasta_to_dict(db):
	dict = {}
	dbFILE = open((db), 'rU')
	for allele in SeqIO.parse(dbFILE, 'fasta'):
		dict[str(allele.seq)] = (allele.id)
		dict[str(allele.seq.reverse_complement())] = (allele.id)
	return dict
