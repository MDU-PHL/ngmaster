#!/bin/bash
set -e

NGMAST="../ngmaster.py"
DEBUG=1
INVARIANT=1
STOP=0

# load test engine
. assert.sh

# no params
assert_raises "$NGMAST" 1

# just help
assert "$NGMAST -h | head -n 1" "usage: "
assert_raises "$NGMAST -h" 0

# test on real fasta
assert "$NGMAST test.fa" "ID\tNG-MAST\tPOR\tTBPB\ntest.fa\t10699\t6277\t4"
assert_raises "$NGMAST test.fa" 0

# test on two files
assert "$NGMAST test.fa test.fa" "ID\tNG-MAST\tPOR\tTBPB\ntest.fa\t10699\t6277\t4\ntest.fa\t10699\t6277\t4"
assert_raises "$NGMAST test.fa test.fa" 0

# test in zero length file
assert "$NGMAST null.fa 2>&1 | head -n 1 | cut -c1-5" "ERROR"

# test in dud fasta file
#assert "$NGMAST noseq.fa"

# test on empty fasta
#assert_raises "$NGMAST /dev/null" 1
#assert "$NGMAST /dev/null | cut -c1-5" "ERROR"

# end of test suite
assert_end ngmaster
