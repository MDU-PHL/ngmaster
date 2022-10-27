# Script by Jason Kwong & Torsten Seemann
# Re-factored and extended (NG-STAR) by Andreas Stroehlein
# In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)
# and Neisseria gonorrhoeae Sequence Typing for Antimicrobial Resistance (NG-STAR)

# import ngmaster functions
import ngmaster
from ngmaster.utils import *

#imports
from argparse import ArgumentParser, RawTextHelpFormatter
import sys
import os
import os.path
from io import StringIO
from subprocess import run, CalledProcessError
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pkg_resources import resource_filename

# Define REST API URLs from PubMLST
ngm_porb = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NG-MAST_porB/alleles_fasta", "comments":""}
ngm_tbpb = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NG-MAST_tbpB/alleles_fasta", "comments":""}
ngm_profiles = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/71/profiles_csv", "comments":""}
ngs_pena = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NEIS1753/alleles_fasta", "comments":"NEIS1753"} 
ngs_mtrr = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/'mtrR/alleles_fasta", "comments":"\'mtrR"}
ngs_porb = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NG_porB/alleles_fasta", "comments":"NG_porB"}
ngs_pona = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NG_ponA/alleles_fasta", "comments":"NG_ponA"}
ngs_gyra = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NG_gyrA/alleles_fasta", "comments":"NG_gyrA"}
ngs_parc = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NG_parC/alleles_fasta", "comments":"NG_parC"}
ngs_23S = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NG_23S/alleles_fasta", "comments":"NG_23S"}
ngs_profiles = {"url": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/67/profiles_csv", "comments":""}


def main():
    # Usage
    parser = ArgumentParser(
        prog="ngmaster",
        formatter_class=RawTextHelpFormatter,
        description='In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)\n'
            'and Neisseria gonorrhoeae Sequence Typing for Antimicrobial Resistance (NG-STAR)\n'
            '\nPlease cite as:\n'
            '  Kwong JC, Gonçalves da Silva A, Dyet K, Williamson DA, Stinear TP, Howden BP and Seemann T.\n'
            '  NGMASTER: in silico multi-antigen sequence typing for Neisseria gonorrhoeae.\n'
            '  Microbial Genomics 2016; doi: 10.1099/mgen.0.000076\n'
            '  GitHub: https://github.com/MDU-PHL/ngmaster\n',
        usage='\n  %(prog)s [OPTIONS] <fasta1> <fasta2> <fasta3> ... <fastaN>')
    parser.add_argument('fasta', metavar='FASTA', nargs='*', help='input FASTA files eg. fasta1, fasta2, fasta3 ... fastaN')
    parser.add_argument('--db', metavar='DB', help='specify custom directory containing allele databases\n'
        'directory must contain database sequence files (.tfa) and allele profile files (ngmast.txt / ngstar.txt)\n'
        'in mlst format (see <https://github.com/tseemann/mlst#adding-a-new-scheme>)')
    parser.add_argument('--csv', action='store_true', default=False, help='output comma-separated format (CSV) rather than tab-separated')
    parser.add_argument('--printseq', metavar='FILE', nargs=1, help='specify filename to save allele sequences to')
    parser.add_argument('--minid', metavar='MINID', nargs=1, default=95, help='DNA percent identity of full allele to consider \'similar\' [~]')
    parser.add_argument('--mincov', metavar='MINCOV', nargs=1, default=10, help='DNA percent coverage to report partial allele at [?]')
    parser.add_argument('--updatedb', action='store_true', default=False, help='update NG-MAST and NG-STAR allele databases from <https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef>')
    parser.add_argument('--assumeyes', action='store_true', default=False, help='assume you are certain you wish to update db')
    parser.add_argument('--test', action='store_true', default=False, help='run test example')
    parser.add_argument('--comments', action='store_true', default=False, help='Include NG-STAR comments for each allele in output')
    parser.add_argument('--version', action='version', version=f'%(prog)s {ngmaster.__version__}')
    args = parser.parse_args()
    # TODO make sure --db works as expected for running mlst (should be ok), once it has been packaged


    idcov = ['--minid', str(args.minid), '--mincov', str(args.mincov)]

    # Set separator
    if args.csv:
        SEP = ','
    else:
        SEP = '\t'

    # Path to database files
    if args.db:
        DBpath = str(args.db).rstrip('/')
    else:
        DBpath = resource_filename(__name__, 'db')

    ngstar_comments = []
    if args.comments:
        ngstar_comments = load_ngstar_comments(DBpath)
        

    ngm_porb['db'] = DBpath + '/pubmlst/ngmast/porB.tfa'
    ngm_tbpb['db'] = DBpath + '/pubmlst/ngmast/tbpB.tfa'
    ngm_profiles['db'] = DBpath + '/pubmlst/ngmast/ngmast.txt'
    ngs_pena['db'] = DBpath + '/pubmlst/ngstar/penA.tfa'
    ngs_mtrr['db'] = DBpath + '/pubmlst/ngstar/mtrR.tfa'
    ngs_porb['db'] = DBpath + '/pubmlst/ngstar/porB.tfa'
    ngs_pona['db'] = DBpath + '/pubmlst/ngstar/ponA.tfa'
    ngs_gyra['db'] = DBpath + '/pubmlst/ngstar/gyrA.tfa'
    ngs_parc['db'] = DBpath + '/pubmlst/ngstar/parC.tfa'
    ngs_23S['db'] = DBpath + '/pubmlst/ngstar/23S.tfa'
    ngs_profiles['db'] = DBpath + '/pubmlst/ngstar/ngstar.txt'

    db_list = [
    ngm_porb,
    ngm_tbpb,
    ngm_profiles,
    ngs_pena,
    ngs_mtrr,
    ngs_porb,
    ngs_pona,
    ngs_gyra,
    ngs_parc,
    ngs_23S,
    ngs_profiles
    ]

    # Check mlst installed and running correctly
    try:
        checkdep = run(['which', 'mlst'], capture_output=True, text=True, check=True)
        mlstpath = checkdep.stdout.strip()
        mkblastdbpath = "/".join(mlstpath.split('/')[:-2]) + "/scripts/mlst-make_blast_db"  
    except CalledProcessError as e:
        err('ERROR: Could not find mlst executable. Check mlst (https://github.com/tseemann/mlst) is installed correctly and in $PATH.')
        raise SystemExit(e)
        

    # Update DB
    if args.updatedb:
        msg('WARNING: Updating DB will overwrite existing DB files.')
        if not args.assumeyes:
            yn = input('Continue? [y/n]: ')
        else:
            yn = 'y'
        if yn == 'y':
            msg("Updating DB files ... ")
            for db in db_list:
                update_db(DBpath, db)
                if not os.path.isfile(db['db']):
                    err('ERROR: Cannot locate database: "{}"'.format(db['db']))
                    raise SystemExit(e)
            download_comments(DBpath, db_list)
            make_mlst_db(DBpath, mkblastdbpath)
        sys.exit(0)

    # Translation table to match NG-STAR and NG-MAST results
    ttable = match_porb(ngm_porb['db'], ngs_porb['db'])

    # Read in NGSTAR profile table for conversion
    ngstartbl = read_ngstar(ngs_profiles['db'])

    # Run test example
    if args.test:
        testSEQ = resource_filename(__name__, "/test/test.fa")
        msg('\033[94mRunning ngmaster.py on test example (NG-MAST 4186 / NG-STAR 231) ...\033[0m')
        args.fasta = [testSEQ]

    # Check if positional arguments
    if not args.fasta:
        parser.print_help()
        err("ERROR: No FASTA file provided")

################################################
    
    output = {"ngmast" : {}, "ngstar" : {}}
    for scheme in output:

        printseq = []
        if args.printseq:
            printseq = ['--novel', scheme.upper() + "__" + args.printseq[0]]

        try:
            result = subprocess.run([mlstpath, '--legacy', '-q', '--threads', '16', '--datadir', DBpath + '/pubmlst', '--blastdb', DBpath + '/blast/mlst.fa', '--scheme', scheme] + idcov + printseq + args.fasta,  capture_output=True, check=True, text=True)
            rlist = result.stdout.split("\n")[:-1] # drop last empty line

            # INFO Checking number of alleles to catch mlst error
            # Issue #125
            # https://github.com/tseemann/mlst/issues/125
            n_allele = len(rlist[0].split("\t")[3:])

            for rec in rlist[1:]: # drop header
                # allele_dct = {}
                reclist = rec.split("\t")
                fname, st = reclist[0], reclist[2]
                alleles = reclist[3:][:n_allele]

                # Add this record to output dict with filename as key
                output[scheme][fname] = MlstRecord(fname,scheme,st,alleles)

        except CalledProcessError as e:
            raise SystemExit(e)

        
    # Collate results from two runs
    collate_out = collate_results(output['ngmast'], output['ngstar'], ttable, ngstartbl)


################################################

    header = ''
    if args.comments:
        header = ['FILE',
                'SCHEME',
                'NG-MAST/NG-STAR',
                'porB_NG-MAST',
                'tbpB',
                'penA',
                'penA_comments',
                'mtrR',
                'mtrR_comments',
                'porB_NG-STAR',
                'porB_NG-STAR_comments',
                'ponA',
                'ponA_comments',
                'gyrA',
                'gyrA_comments',
                'parC',
                'parC_comments',
                '23S',
                '23S_comments'
                ]

    else:
        header = ['FILE', 'SCHEME', 'NG-MAST/NG-STAR', 'porB_NG-MAST', 'tbpB', 'penA', 'mtrR', 'porB_NG-STAR', 'ponA', 'gyrA', 'parC', '23S']

    print(SEP.join(header))

    for out in collate_out:
        print(out.get_record(sep = SEP, comments = ngstar_comments))  

    if args.test:
        if collate_out[0].st != '4186/231':
            err('ERROR: Test unsucessful. Check allele database is updated: ngmaster.py --updatedb')
        else:
            msg('\033[92m... Test successful.\033[0m')

if __name__ == "__main__":
    main()
