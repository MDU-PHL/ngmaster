#!/usr/bin/env python
# Script by Jason Kwong
# In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)

# Usage
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)\n'
		'\nRef: Martin et al. J Infect Dis, 2004 Apr 15;189(8):1497-1505.\n'
		'See also http://www.ng-mast.net/',
	usage='\n  %(prog)s [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>')
parser.add_argument('fasta', metavar='FASTA', nargs='+', help='input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN')
parser.add_argument('--db', metavar='DB', help='Specify custom directory containing allele databases\n'
	'Directory must contain database files "POR.tfa", "TBPB.tfa", and "ng_mast.txt"')
parser.add_argument('--printseq', metavar='FILE', nargs=1, help='Specify filename to save allele sequences to (default=off)')
parser.add_argument('--version', action='version', version=
	'=====================================\n'
	'%(prog)s v0.1\n'
	'Updated 1-Sept-2015 by Jason Kwong\n'
	'Dependencies: isPcr, BLAST, BioPython\n'
	'=====================================')
args = parser.parse_args()

# Modules and Functions
import sys
import os
import os.path
import StringIO
import subprocess
from subprocess import Popen
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline

# Path to database files
if args.db:
	DBpath = str(args.db).rstrip('/')
else:
	DBpath = os.path.dirname(os.path.realpath(sys.argv[0])) + "/db"

porDB = DBpath + '/POR.tfa'
tbpbDB = DBpath + '/TBPB.tfa'
alleleDB = DBpath + '/ng_mast.txt'

def progexit(n):
	sys.exit(n)

# Need to fix - running isPcr without options produces error code
# Check isPcr installed and running correctly
#def progcheck(isPcr):
#	devnull = open(os.devnull, 'w')
#	checkdep = subprocess.Popen(['isPcr'], stdout=devnull, stderr=subprocess.PIPE, close_fds=True)
#	output, err = checkdep.communicate()
#	if checkdep.returncode != 0:
#		print 'ERROR: Check isPcr is installed correctly.'
#		progexit(1)
#progcheck('isPcr')

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

# Set up primer database
primerDB = [['por', 'CAAGAAGACCTCGGCAA', 'CCGACAACCACTTGGT'], ['tbpB', 'CGTTGTCGGCAGCGCGAAAAC', 'TTCATCGGTGCGCTCGCCTTG']]
NGprimers = "\n".join(" ".join(map(str,l)) for l in primerDB) + "\n"

# Run Jim Kent's isPcr to identify amplicon
# Check queries are in FASTA format
alleleSEQS = []
print 'ID' + '\t' + 'NG-MAST' + '\t' + 'POR' + '\t' + 'TBPB'
for f in args.fasta:
	if os.path.isfile(f) == False:
		print 'ERROR: Cannot find "%(f)s". Check file exists.' % globals()
		continue
	s = open(f, 'r')
	if s.read(1) != '>':
		print 'ERROR: "%(f)s" does not appear to be in FASTA format.' % globals()
		continue
	s.close()

	# Setup lists in case there are multiple hits
	por = None
	porCOUNT = []
	porKEY = []
	tbpb = None
	tbpbCOUNT = []
	tbpbKEY = []
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
			if int(ampLEN[:-2]) > 637 and int(ampLEN[:-2]) < 837:	# Check por amplicon length
				porSEQ = amplicon.seq.upper()
				start = porSEQ.find('TTGAA')
				if start != -1:										# Check for starting key motif
					newporSEQ = str(porSEQ[start:(start+490)])		# Trim sequence from starting key motif
					if len(newporSEQ) == 490:
						if 'por' not in porCOUNT:
							porCOUNT.append('por')
						# Add sequences to print later
						porSEQR = Seq(newporSEQ)
						porRECR = SeqRecord(porSEQR, id=f, description='POR')
						alleleSEQS.append(porRECR)
						# BLAST trimmed sequence against database
						porBLAST = NcbiblastnCommandline(query='-', subject=porDB, evalue='0.001', perc_identity='100', outfmt='"6 qseqid sseqid pident length"')
						stdout, stderr = porBLAST(stdin=porRECR.format('fasta'))
						if stdout:
							BLASTout = stdout.split('\n')
							for line in BLASTout:
								if line.strip():
									lines = line.split('\t')
									if lines[3] == '490':
										porRESULT = lines[1]
										por = porRESULT.split('-')[1]
									elif not por:
										por = 'new'
						elif not por:
							por = 'new'
						porCOUNT.append(por)
					else:
						porKEY.append('no_key')

		if ampID == "tbpB":
			if int(ampLEN[:-2]) > 480 and int(ampLEN[:-2]) < 680:	# Check tbpB amplicon length
				tbpbSEQ = amplicon.seq.upper()
				start = tbpbSEQ.find('CGTCTGAA')
				if start != -1:										# Check for starting key motif
					newtbpbSEQ = str(tbpbSEQ[start:(start+390)])	# Trim sequence from starting key motif
					if len(newtbpbSEQ) == 390:
						if 'tbpB' not in tbpbCOUNT:
							tbpbCOUNT.append('tbpB')
						# Add sequences to print later
						tbpbSEQR = Seq(newtbpbSEQ)
						tbpbRECR = SeqRecord(tbpbSEQR, id=f, description='TBPB')
						alleleSEQS.append(tbpbRECR)
						# BLAST trimmed sequence against database
						tbpbBLAST = NcbiblastnCommandline(query='-', subject=tbpbDB, evalue='0.001', perc_identity='100', outfmt='"6 qseqid sseqid pident length"')
						stdout, stderr = tbpbBLAST(stdin=tbpbRECR.format('fasta'))
						if stdout:
							BLASTout = stdout.split('\n')
							for line in BLASTout:
								if line.strip():
									lines = line.split('\t')
									if lines[3] == '390':
										tbpbRESULT = lines[1]
										tbpb = tbpbRESULT.split('-')[1]
									elif not tbpb:
										tbpb = 'new'
						elif not tbpb:
							tbpb = 'new'
						tbpbCOUNT.append(tbpb)
					else:
						tbpbKEY.append('no_key')
	alleleSEQ.close()

	if not por:
		por = '-'
	if not tbpb:
		tbpb = '-'

	# If multiple hits with trimmed sequence (eg. duplicated genes, multiple starting key motifs etc.) print multiple results
	if len(porCOUNT) > 2 or len(tbpbCOUNT) > 2:
		print f + '\t' + 'multiple' + '\t' + str(porCOUNT) + '\t' + str(tbpbCOUNT)
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
		print f + '\t' + type + '\t' + por + '\t' + tbpb

# Print allele sequences to file
if args.printseq:
	allelesOUT = "".join(args.printseq)
	with open(allelesOUT, "w") as output:
		SeqIO.write(alleleSEQS, output, 'fasta')

progexit(0)
