#!/usr/bin/env python

from __future__ import print_function
import os,sys,argparse

a=1
b=3
q=5
r=2

RC=dict(zip('ATGCNatcgn','TACGNtagcn'))

def sam_seq_dict( fai ):
    d = []
    for line in fai:
        chrom, length = line.split()[:2]
        d.append( '\t'.join( [ '@SQ', 'SN:%s'%chrom, 'LN:%s'%length ] ) )
    return d

def rcseq( seq ):
    ''' returns the reverse complement of a DNA sequence '''
    return ''.join([ RC[n] for n in seq[::-1] ])

def get_seq_fasta( fna_file_ptr ):
    ''' returns a dictionary of sequences with sequence ID as key and sequence as value '''
    seqs = {}
    lc = 0
    for line in fna_file_ptr:
        line = line.strip()
        if lc%2==0:
            header = line[1:]
        else:
            seqs[header] = line
        lc+=1
    return seqs

def cigar_list_to_expanded_string( cigar_list ):
    return ''.join([ l*op for l, op in cigar_list ])

def cigar_expanded_string_to_list(cigar_expanded_string):
    cigar = [ [1,cigar_expanded_string[0]] ]
    for idx in range(1,len(cigar_expanded_string)):
        if cigar_expanded_string[idx] != cigar[-1][1]:
            cigar.append( [1,cigar_expanded_string[idx]] )
        else:
            cigar[-1][0] += 1
    return [ (l, o) for l,o in cigar ]

def psl2sam( args ):
    sequences = get_seq_fasta( args.fasta )

    seq_dict = sam_seq_dict( args.fai )
    for s in seq_dict:
        print(s)

    for line in args.psl_file:
        # parse PSL line
        toks = line.strip().split('\t')
        qname = toks[9]
        rname = toks[13]
        qsize = int(toks[10])
        qstart = int(toks[11])
        qend = int(toks[12])
        tstart = int(toks[15])
        strand = toks[8]
        
        if strand=='-':
            tmp = qstart
            qstart = qsize - qend
            qend = qsize - tmp

        # construct SAM record
        sam = [ qname,                                                          # QUERY_NAME
                0 if strand=='+' else 16,                                       # FLAG
                rname,                                                          # TARGET_NAME
                tstart + 1,                                                     # TARGET_POS
                0,                                                              # MAPQ
                '',                                                             # CIGAR, temporarily empty
                '*', 0, 0,                                                      # MATE_TARGET_NAME, MATE_TARGET_POS, INSERT_SIZE
                sequences[ qname ] if strand=='+' else rcseq(sequences[qname]), # SEQ
                len(sequences[ qname ])*'I' ]                                   # QUAL

        # 5'-end clipping
        cigar = [ (qstart,'S') ] if qstart else []
        block_sizes = map(int, [t for t in toks[18].split(',') if t])
        qstarts = map(int, [t for t in toks[19].split(',') if t])
        tstarts = map(int, [t for t in toks[20].split(',') if t])

        qstart0, tstart0 = qstarts[0], tstarts[0]

        block_count = int(toks[17])
        gap_open, gap_ext = 0, 0
        for i in range(1,block_count):
            lqstart = qstarts[i] - qstarts[i-1] - block_sizes[i-1]
            ltstart = tstarts[i] - tstarts[i-1] - block_sizes[i-1]

            #print('Block-%d'%i, qstarts[i], qstarts[i-1], block_sizes[i-1], lqstart, tstarts[i], ltstart)
            if ltstart < lqstart: # ins: the query gap is longer
                gap_open += 1
                gap_ext += lqstart 
                cigar.append( (block_sizes[i-1], 'M') )
                cigar.append( (lqstart,'I') )
                if ltstart != 0:
                    cigar.append( (ltstart,'D') )
                    gap_open += 1
                    gap_ext += ltstart
                qstart0, tstart0 = qstarts[i], tstarts[i]
            if lqstart < ltstart: # del: the reference gap is longer
                gap_open += 1
                gap_ext += ltstart
                cigar.append( (block_sizes[i-1], 'M') )
                if lqstart != 0:
                    cigar.append( (lqstart, 'I') )
                    gap_open += 1
                    gap_ext += lqstart
                cigar.append( (ltstart, 'D') )
                qstart0, tstart0 = qstarts[i], tstarts[i]

        cigar.append( ((qend - qstart0), 'M') )
        # 3'-end clipping
        if qsize != qend:
            cigar.append( ((qsize - qend), 'S') )

        # remove primer alignment
        if args.primer_length:
            offset_head = args.primer_length
            offset_tail = args.primer_length
            # shift alignment start position
            sam[3] = tstart+1+(args.primer_length-qstart) if qstart < args.primer_length else tstart+1
            #print(qstart, qend, qsize, offset_head, offset_tail, sam[3], tstart+1, sep='\t')
            # shorten sequence and quality length
            # EDIT: switch to soft clipping
            if not args.soft_clip_primer:
                sam[9] = sam[9][offset_head:-offset_tail]
                sam[10] = sam[10][offset_head:-offset_tail]

            # remove cigar operations in primer regions that are NOT DELETIONS
            # convert to string form
            cigar_str = cigar_list_to_expanded_string( cigar )
            # remove from 5' end
            cigar_str_new = ''
            removed = 0
            for op in cigar_str:
                if removed < offset_head:
                    # EDIT: switch to soft clipping
                    if args.soft_clip_primer:
                        if op!='D': cigar_str_new += 'S'

                    if op!='D': removed += 1
                    else: sam[3] += 1
                    continue
                cigar_str_new += op
            # remove from 3' end
            cigar_str_new_rev = cigar_str_new[::-1]
            cigar_str_new = ''
            removed = 0
            for op in cigar_str_new_rev:
                if removed < offset_tail:
                    # EDIT: switch to soft clipping
                    if args.soft_clip_primer:
                        if op!='D': cigar_str_new += 'S'
                    if op!='D': removed += 1
                    continue
                cigar_str_new += op
            cigar_str_new = cigar_str_new[::-1]

            if not args.allow_end_deletion:
            # ensure no D at ends
                for idx, c in enumerate(cigar_str_new):
                    if c=='D':
                        sam[3] += 1
                        continue
                    else: break
                cigar_str_new = cigar_str_new[idx:]
                for idx, c in enumerate(cigar_str_new[::-1]):
                    if c=='D': continue
                    else: break
                cigar_str_new = cigar_str_new[:len(cigar_str_new)-idx]

            # convert back to list form
            cigar = cigar_expanded_string_to_list( cigar_str_new )

        sam[5] = ''.join( ['%d%s'%v for v in cigar] )
        matches, mismatches = int(toks[0]), int(toks[1])
        score = a * matches - b * mismatches - q * gap_open - r * gap_ext
        if score < 0: score = 0
        sam.append( "AS:i:%d"%score )

        if mismatches > args.max_variant:
            sam[4] = 0
        elif gap_open > args.max_gap:
            sam[4] = 0
        else:
            sam[4] = int(max((1.0 - (mismatches/float(args.max_variant) + gap_open/float(args.max_gap))/2),0)*40)
            #sam[4] = int(max(matches-mismatches,0)/float(qsize)*40)

        print('\t'.join( map(str,sam) ), file=args.out)


def command_line_interface(*args,**kw):
    ''' command line parsing function '''
    parser = argparse.ArgumentParser(description="A script for converting BLAT's PSL format to SAM \
                                                  - ADAPTED FROM HENG LI's psl2sam.pl in Samtools")

    parser.add_argument('psl_file',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help='input PSL file, or "-" for stdin')

    parser.add_argument('--max_variant',
                        type=int,
                        default=3,
                        help='maximum number of variants expected per amplicon, mapping score = 0 if amplicon contains more than N SNVs')

    parser.add_argument('--max_gap',
                        type=int,
                        default=2,
                        help='maximum number of gaps expected per amplicon, mapping score = 0 if amplicon contains more than N gaps')

    parser.add_argument('--soft_clip_primer',
                        action="store_true",
                        help='soft clip primer sequences rather than hard clip')

    parser.add_argument('--allow_end_deletion',
                        action="store_true",
                        help='allow deletion operations at ends of the alignment')

    parser.add_argument('--out',
                        type=argparse.FileType('w'),
                        help='Output file, or "-" for stdout (default)')

    parser.add_argument('--primer_length',
                        type=int,
                        default=0,
                        help='size of flanking primer to be removed (intentionally kept for better BLAT alignment)')

    parser.add_argument('--fasta',
                        type=argparse.FileType('r'),
                        required=True,
                        default=None,
                        help='FASTA file for sequences used in BLAT alignment')

    parser.add_argument('--fai',
                        type=argparse.FileType('r'),
                        required=True,
                        default=None,
                        help='Reference FASTA file index (fai) for sequence dict of SAM')

    return parser.parse_args(*args,**kw)


if __name__=='__main__':
    args = command_line_interface()
    psl2sam( args )
