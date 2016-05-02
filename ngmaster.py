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
import urllib
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

# Set target lengths for isPcr amplicons and trimmed sequences
porAMPLEN = 737
porTRIMLEN = 490
tbpbAMPLEN = 580
tbpbTRIMLEN = 390


# Usage
parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)\n'
		'\nRef: Martin et al. J Infect Dis, 2004 Apr 15;189(8):1497-1505.\n'
		'See also http://www.ng-mast.net/',
	usage='\n  %(prog)s [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>')
parser.add_argument('fasta', metavar='FASTA', nargs='*', help='input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN')
parser.add_argument('--db', metavar='DB', help='specify custom directory containing allele databases\n'
	'directory must contain database files "POR.tfa", "TBPB.tfa", and "ng_mast.txt"')
parser.add_argument('--csv', action='store_true', default=False, help='output comma-separated format (CSV) rather than tab-separated')
parser.add_argument('--printseq', metavar='FILE', nargs=1, help='specify filename to save allele sequences to (default=off)')
parser.add_argument('--updatedb', action='store_true', default=False, help='update allele database from <www.ng-mast.net>')
parser.add_argument('--test', action='store_true', default=False, help='run test example')
parser.add_argument('--version', action='version', version=
	'=====================================\n'
	'%(prog)s v0.3\n'
	'Updated 1-May-2016 by Jason Kwong\n'
	'Dependencies: isPcr, BioPython\n'
	'=====================================\n'
	'Please cite as:\n'
    'Kwong JC, Goncalves da Silva A, Howden BP and Seemann T.\n'
    'NGMASTER: in silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)\n'
    'GitHub: https://github.com/MDU-PHL/ngmaster\n')
args = parser.parse_args()

# Path to database files
if args.db:
	DBpath = str(args.db).rstrip('/')
else:
	DBpath = os.path.dirname(os.path.realpath(sys.argv[0])) + "/db"
if not os.path.exists(DBpath):
    os.makedirs(DBpath)

porDB = DBpath + '/POR.tfa'
tbpbDB = DBpath + '/TBPB.tfa'
alleleDB = DBpath + '/ng_mast.txt'
tempFILE = DBpath + '/temp'

# Update DB
if args.updatedb:
	msg('WARNING: Updating DB will overwrite existing DB files.')
	yn = raw_input('Continue? [y/n]: ')
	if yn == 'y':
		msg("Updating DB files ... ")
		# Update POR DB
		dbFILE = urllib.URLopener()
		dbFILE.retrieve("http://www.ng-mast.net/sql/fasta.asp?allele=POR", tempFILE)
		format(tempFILE)
		os.rename(tempFILE, porDB)
		msg(porDB + ' ... Done.')
		# Update TBPB DB
		dbFILE = urllib.URLopener()
		dbFILE.retrieve("http://www.ng-mast.net/sql/fasta.asp?allele=TBPB", tempFILE)
		format(tempFILE)
		os.rename(tempFILE, tbpbDB)
		msg(tbpbDB + ' ... Done.')
		# Update allele DB
		dbFILE = urllib.URLopener()
		dbFILE.retrieve("http://www.ng-mast.net/sql/st_comma.asp", tempFILE)
		sed_inplace(tempFILE, ',', '\t')
		sed_inplace(tempFILE, '<br>', '\n')
		with open(alleleDB, 'w') as f:
			f.write('ST' + '\t' + 'POR' + '\t' + 'TBPB' + '\n')
			with open(tempFILE, 'r') as t:
				for line in t:
					f.write(line)
		os.remove(tempFILE)
		msg(alleleDB + ' ... Done.')
	sys.exit(0)

# Check isPcr installed and running correctly
devnull = open(os.devnull, 'w')
checkdep = subprocess.Popen(['which', 'isPcr'], stdout=devnull, stderr=subprocess.PIPE, close_fds=True)
output, errcode = checkdep.communicate()
if checkdep.returncode != 0:
	err('ERROR: Check isPcr is installed correctly and in $PATH.')

# Set separator
if args.csv:
	SEP = ','
else:
	SEP = '\t'

# Run test example
if args.test:
	testSEQ = os.path.dirname(os.path.realpath(sys.argv[0])) + "/test/test.fa"
	msg('Running ngmaster.py on test example (NG-MAST 10699) ...')
	msg('$ ngmaster.py test/test.fa')
	args.fasta = [testSEQ]

# Check if positional arguments
if not args.fasta:
	parser.print_help()
	err("ERROR: too few arguments")

# Import allele profiles as dictionary
NGMAST = {}
with open(alleleDB) as f:
	for line in f:
		if line.strip():
			lines = line.split('\t')
			ST = lines[0]
			porALLELE = lines[1].rstrip('\n')
			tbpbALLELE = lines[2].rstrip('\n')
			alleles = str(porALLELE) + '-' + str(tbpbALLELE)
			NGMAST[alleles] = ST

# Import allele databases
porDICT = fasta_to_dict(porDB)
tbpbDICT = fasta_to_dict(tbpbDB)

# Set up primer database
primerDB = [['por', 'CAAGAAGACCTCGGCAA', 'CCGACAACCACTTGGT'], ['tbpB', 'CGTTGTCGGCAGCGCGAAAAC', 'TTCATCGGTGCGCTCGCCTTG']]
NGprimers = "\n".join(" ".join(map(str,l)) for l in primerDB) + "\n"

# Check queries are in FASTA format
# Run Jim Kent's isPcr to identify amplicon
alleleSEQS = []
print('ID' + SEP + 'NG-MAST' + SEP + 'POR' + SEP + 'TBPB')
for f in args.fasta:
	if os.path.isfile(f) == False:
		msg( 'ERROR: Cannot find "%(f)s". Check file exists.' % globals() )
		continue
	s = open(f, 'r')
	if s.read(1) != '>':
		msg( 'ERROR: "%(f)s" does not appear to be in FASTA format.' % globals() )
		continue
	s.close()

	# Setup lists in case there are multiple hits
	por = None
	porCOUNT = set()
	porKEY = set()
	tbpb = None
	tbpbCOUNT = set()
	tbpbKEY = set()
	# Run isPcr by Jim Kent
	proc = subprocess.Popen(['isPcr', f, 'stdin', 'stdout', '-tileSize=6', '-minPerfect=5', '-stepSize=3', '-maxSize=900'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	PCRout = proc.communicate(input=NGprimers)[0]
	alleleSEQ = StringIO.StringIO()
	alleleSEQ.write(PCRout)
	alleleSEQ.seek(0)
	# Check amplicon length and starting key motif
	for amplicon in SeqIO.parse(alleleSEQ, "fasta"):
		product = amplicon.description.split()
		ampID = product[1]
		ampLEN = product[2]
		if ampID == "por":
			if int(ampLEN[:-2]) > (porAMPLEN-100) and int(ampLEN[:-2]) < (porAMPLEN+100):	# Check por amplicon length
				porSEQ = amplicon.seq.upper()
				start = porSEQ.find('TTGAA')
				if start != -1:											# Check for starting key motif
					newporSEQ = str(porSEQ[start:(start+porTRIMLEN)])	# Trim sequence from starting key motif
					if len(newporSEQ) == porTRIMLEN:
						# Add sequences to print later
						porSEQR = Seq(newporSEQ)
						porRECR = SeqRecord(porSEQR, id=f, description='POR')
						alleleSEQS.append(porRECR)
						# Search trimmed sequence against database dictionary
						try:
							porRESULT = (porDICT[str(porSEQR)])
							por = porRESULT.split('R')[1]
						except KeyError:
							por = 'new'
							continue
						if por not in porCOUNT:
							porCOUNT.add(por)
				else:
					porKEY.add('no_key')
		if ampID == "tbpB":
			if int(ampLEN[:-2]) > (tbpbAMPLEN-100) and int(ampLEN[:-2]) < (tbpbAMPLEN+100):	# Check tbpB amplicon length
				tbpbSEQ = amplicon.seq.upper()
				start = tbpbSEQ.find('CGTCTGAA')
				if start != -1:												# Check for starting key motif
					newtbpbSEQ = str(tbpbSEQ[start:(start+tbpbTRIMLEN)])	# Trim sequence from starting key motif
					if len(newtbpbSEQ) == tbpbTRIMLEN:
							# Add sequences to print later
						tbpbSEQR = Seq(newtbpbSEQ)
						tbpbRECR = SeqRecord(tbpbSEQR, id=f, description='TBPB')
						alleleSEQS.append(tbpbRECR)
						# Search trimmed sequence against database dictionary
						try:
							tbpbRESULT = (tbpbDICT[str(tbpbSEQR)])
							tbpb = tbpbRESULT.split('PB')[1]
						except KeyError:
							tbpb = 'new'
							continue
						if tbpb not in tbpbCOUNT:
							tbpbCOUNT.add(tbpb)
				else:
					tbpbKEY.add('no_key')
	alleleSEQ.close()

	if not por:
		por = '-'
	if not tbpb:
		tbpb = '-'

	# If multiple hits with trimmed sequence (eg. duplicated genes, multiple starting key motifs etc.) print multiple results
	if len(porCOUNT) > 1 or len(tbpbCOUNT) > 1:
		print( f + SEP + 'multiple' + SEP + str(porCOUNT) + SEP + str(tbpbCOUNT) )
	else:
	# Report if starting key motifs present
		if not porCOUNT:
			if porKEY:
				por = 'no_key'
		if not tbpbCOUNT:
			if tbpbKEY:
				tbpb = 'no_key'
		portbpb = str(por) + '-' + str(tbpb)	
	# Print results to screen
		if portbpb in NGMAST:
			type = NGMAST[portbpb]
		else:
			type = "-"
		if not args.test:
			print( f + SEP + type + SEP + por + SEP + tbpb )
		else:
			print( 'test.fa' + SEP + type + SEP + por + SEP + tbpb )
			if type != '10699':
				err('ERROR: Test unsucessful. Check allele database is updated: ngmaster.py --updatedb')
			else:
				msg('... Test successful.')

# Print allele sequences to file
if args.printseq:
	allelesOUT = "".join(args.printseq)
	with open(allelesOUT, "w") as output:
		SeqIO.write(alleleSEQS, output, 'fasta')

