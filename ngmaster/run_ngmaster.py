# Script by Jason Kwong & Torsten Seemann
# Re-factored and extended (NG-STAR) by Andreas Stroehlein
# In silico multi-antigen sequence typing for Neisseria gonorrhoeae (NG-MAST)
# and Neisseria gonorrhoeae Sequence Typing for Antimicrobial Resistance (NG-STAR)

# import ngmaster functions
import ngmaster
from ngmaster.utils import * 
# from utils import * # FLAG: import all functions from the utils.py not the compiled `.pyc` file

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
from pathlib import Path
from mlstdb.core.download import create_blast_db
from concurrent.futures import ThreadPoolExecutor, as_completed


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


def _run_scheme(scheme, mlstpath, DBpath, idcov, printseq_arg, fasta, threads):
    """Run mlst for a single scheme; returns (scheme, dict[fname -> MlstRecord])."""
    printseq = []
    novel_output = None
    if printseq_arg:
        novel_output = scheme.upper() + "__" + printseq_arg[0]
        printseq = ['--novel', novel_output]

    result = subprocess.run(
        [mlstpath, '--legacy', '-q', '--threads', str(threads),
         '--datadir', DBpath + '/pubmlst',
         '--blastdb', DBpath + '/blast/mlst.fa',
         '--scheme', scheme] + idcov + printseq + fasta,
        capture_output=True, check=True, text=True
    )
    if novel_output:
        novel_output_path = Path(novel_output)
        if novel_output_path.exists() and novel_output_path.stat().st_size == 0:
            novel_output_path.unlink()

    rlist = result.stdout.split("\n")[:-1]  # drop last empty line

    if scheme == 'ngstar':
        rlist = process_duplicate_23s_alleles(rlist)

    # INFO Checking number of alleles to catch mlst error
    # Issue #125 https://github.com/tseemann/mlst/issues/125
    n_allele = len(rlist[0].split("\t")[3:])

    scheme_output = {}
    for rec in rlist[1:]:  # drop header
        reclist = rec.split("\t")
        fname, st = reclist[0], reclist[2]
        alleles = reclist[3:][:n_allele]
        scheme_output[fname] = MlstRecord(fname, scheme, st, alleles)

    return scheme, scheme_output


def _validate_fasta_inputs(fasta):
    for fname in fasta:
        path = Path(fname)
        if not path.is_file() or path.stat().st_size == 0:
            err(f"ERROR: Input FASTA is empty or contains no readable sequence: {fname}")

        try:
            has_sequence = any(len(record.seq) > 0 for record in SeqIO.parse(fname, "fasta"))
        except Exception:
            has_sequence = False

        if not has_sequence:
            err(f"ERROR: Input FASTA is empty or contains no readable sequence: {fname}")


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
        'in mlst format (see <https://github.com/tseemann/mlst#adding-a-new-scheme>)\n'
        'overrides $NGMASTER_DB environment variable if both are set\n'
        f'default: $NGMASTER_DB if set, otherwise {Path(__file__).parent / "db"}\n')
    parser.add_argument('--csv', action='store_true', default=False, help='output comma-separated format (CSV) rather than tab-separated')
    parser.add_argument('--json', action='store_true', default=False, help='output JSON format rather than tab-separated')
    parser.add_argument('--printseq', metavar='FILE', nargs=1, help='specify filename to save novel allele sequences to\n'
        '(only alleles marked ~n are written; no file created if all alleles are exact matches)')
    parser.add_argument('--minid', metavar='MINID', type=int, default=95, help='DNA percent identity of full allele to consider \'similar\' [~] (default: 95, range: 0-100)')
    parser.add_argument('--mincov', metavar='MINCOV', type=int, default=50, help='DNA percent coverage to report partial allele at [?] (default: 50, range: 0-100)')
    parser.add_argument('--updatedb', action='store_true', default=False, help='update NG-MAST and NG-STAR allele databases from <https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef>')
    parser.add_argument('--assumeyes', action='store_true', default=False, help='assume you are certain you wish to update db')
    parser.add_argument('--test', action='store_true', default=False, help='run test example')
    parser.add_argument('--comments', action='store_true', default=False, help='Include NG-STAR comments for each allele in output')
    parser.add_argument('--threads', metavar='THREADS', type=int, default=1,
        help='number of threads to use\n'
             '  --updatedb: concurrent HTTP requests to PubMLST (default 1, recommended 4)\n'
             '  main analysis: runs NG-MAST and NG-STAR schemes in parallel\n'
             '                 and passes thread count to the mlst subprocess\n'
             'default: 1\n')
    parser.add_argument("--version", action="store_true", help="Show version information")
    args = parser.parse_args()

    if not (0 <= args.minid <= 100):
        err('ERROR: --minid value must be between 0 and 100')
    if not (0 <= args.mincov <= 100):
        err('ERROR: --mincov value must be between 0 and 100')
    if args.printseq:
        msg('NOTE: --printseq saves only novel (~n) allele sequences. '
            'No output file is created if all alleles are exact matches.')

    idcov = ['--minid', str(args.minid), '--mincov', str(args.mincov)]
    threads = args.threads

    # Set separator
    if args.csv:
        SEP = ','
    else:
        SEP = '\t'

    # Path to database files
    if args.db:
        DBpath = str(args.db).rstrip('/')
    elif os.environ.get('NGMASTER_DB'):
        DBpath = os.environ['NGMASTER_DB'].rstrip('/')
    else:
        DBpath = str(Path(__file__).parent / 'db')

    if args.version:
        print(f"ngmaster {ngmaster.__version__}")
        db_version = get_db_version(DBpath)
        print(f"Database version: {db_version}")
        sys.exit(0)

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
        # mkblastdbpath = "/".join(mlstpath.split('/')[:-2]) + "/scripts/mlst-make_blast_db"  
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
            msg(f"Updating DB files using {threads} thread(s) ... ")
            with ThreadPoolExecutor(max_workers=threads) as executor:
                futures = {executor.submit(update_db, DBpath, db): db for db in db_list}
                for future in as_completed(futures):
                    db = futures[future]
                    future.result()  # re-raise any exception from the worker
                    if not os.path.isfile(db['db']):
                        err('ERROR: Cannot locate database: "{}"'.format(db['db']))
            download_comments(DBpath, db_list, threads=threads)
            msg("\nCreating BLAST database from downloaded pubMLST schemes...")
            create_blast_db(DBpath + '/pubmlst', DBpath + '/blast')
            
        sys.exit(0)

    # Translation table to match NG-STAR and NG-MAST results
    ttable = match_porb(ngm_porb['db'], ngs_porb['db'])

    # Read in NGSTAR profile table for conversion
    ngstartbl = read_ngstar(ngs_profiles['db'])

    # Run test example
    if args.test:
        testSEQ = str(Path(__file__).parent / 'test' / 'test.fa')
        msg('\033[94mRunning ngmaster on test example (NG-MAST 4186 / NG-STAR 231) ...\033[0m')
        args.fasta = [testSEQ]
        print('Test example FASTA file: {}'.format(testSEQ))

    # Check if positional arguments
    if not args.fasta:
        parser.print_help()
        err("ERROR: No FASTA file provided")

    _validate_fasta_inputs(args.fasta)

################################################

    output = {"ngmast": {}, "ngstar": {}}
    # Always parallelise the two schemes (ngmast, ngstar) with max_workers=2
    # These are I/O-bound subprocess calls that benefit from parallel execution regardless of user's --threads setting
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = {
            executor.submit(_run_scheme, scheme, mlstpath, DBpath, idcov, args.printseq, args.fasta, threads): scheme
            for scheme in output
        }
        for future in as_completed(futures):
            try:
                scheme, scheme_output = future.result()
                output[scheme] = scheme_output
            except CalledProcessError as e:
                if e.stderr and e.stderr.strip():
                    msg(e.stderr.strip())
                err(f"ERROR: mlst failed while running {futures[future]}")

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
                '23S_comments',
                'CC'  # Add CC column
                ]

    else:
        header = ['FILE', 'SCHEME', 'NG-MAST/NG-STAR', 'porB_NG-MAST', 'tbpB', 'penA', 'mtrR', 'porB_NG-STAR', 'ponA', 'gyrA', 'parC', '23S', 'CC']  # Add CC column

    if args.json:
        import json
        json_out = []
        for out in collate_out:
            record_list = out.get_list(comments=ngstar_comments, header=header)
            json_out.append(dict(zip(header, record_list)))
        print(json.dumps(json_out, indent=2))
    else:
        print(SEP.join(header))
        for out in collate_out: 
            print(out.get_record(sep = SEP, comments = ngstar_comments))
    
    if args.test:
        if collate_out[0].st != '4186/231':
            err('ERROR: Test unsuccessful. Check allele database is updated: ngmaster --updatedb')
        else:
            msg('\033[92m... Test successful.\033[0m')
            
if __name__ == "__main__":
    main()
