#!/usr/bin/env python
'''
Inflates a BAM file from collapsed FASTQ sequences.
- collapsed FASTA sequences: 
    X sequences with identical nucleotide content
    X is encoded in FASTA header in the form, e.g.:
        1_945_TP53_exon11_175_chr17:7573924-7574040_reads_40_MS0312_667FFPE_100%_length_109
    as "reads_40"
- Inflation process simply stutters the SAM record X times
    This preserves amplicon information in TAGS where:
        XC  amplicon chromosome
        XB  amplicon begin
        XE  amplicon end
        XP  percentage
        XL  amplicon lenth
        XR  identical sequences

Usage: [SAM STREAM] | inflate_alignment.py | samtools view -Sbu - | samtools sort - OUTPUT

e.g. bwa mem -Mt THREADS INDEX FASTA | inflate_alignment.py | samtools view -Sbu - | samtools sort - OUTPUT
'''
from __future__ import print_function
import sys, csv, re, argparse

cigar_matcher=re.compile('\*|([0-9]+)([MIDNSHPX=])')

def cigar_str_to_cigar_list( cigar_str ):
    return [ (int(l), op) for l,op in cigar_matcher.findall( cigar_str ) ]

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

def inflate( args ):
    input = args.sam_file
    
    records = 0
    raw = 0
    for row in csv.reader( input, delimiter='\t'):
        # reconstruct header section
        if row[0].startswith('@'):
            print('\t'.join(row),file=args.out)
            continue
    
        # split header and extract fields
        toks = row[0].split('_')
        #toks = re.split('_+', row[0])
        n_length = toks[-1]
        #length = toks[-2]
        #perc = toks[-3]
        n_reads = toks[-4]
        #reads = toks[-5]
        #coord = toks[-6]
        #assert perc.endswith('%')
        #assert reads=='reads'
        #assert length=='length'
        #chrom, grange = coord.split(':')
        #gstart, gend = map(int, grange.split('-'))
        n_reads, n_length = map(int, [n_reads, n_length])
        #perc = float(perc[:-1])/100
    
        # put fields back into SAM format order
        SAM = [ '_'.join(toks) ] + row[1:]
        # add tags for extracted info
        #SAM += ['XC:A:%s'%chrom, 'XB:i:%d'%gstart, 'XE:i:%d'%gend, 'XP:f:%f'%perc, 'XL:i:%d'%n_length, 'XR:i:%d'%n_reads]

        # remove primer alignment
        if args.primer_length:
            offset_head = args.primer_length
            offset_tail = args.primer_length
            # shorten sequence and quality length
            # do hard clipping if not soft clipping
            if not args.soft_clip_primer:
                SAM[9] = SAM[9][offset_head:-offset_tail]
                SAM[10] = SAM[10][offset_head:-offset_tail]

            # remove cigar operations in primer regions that are NOT DELETIONS
            # convert to string form
            cigar = cigar_str_to_cigar_list( SAM[5] )
            cigar_str = cigar_list_to_expanded_string( cigar )
            # remove from 5' end
            cigar_str_new = ''
            removed = 0
            idx = 0
            for op in cigar_str:
                if args.soft_clip_primer:
                    if op=='D' or op=='I' or idx==offset_head:
                        break
                    idx+=1
                #if removed < offset_head:
                    # EDIT: switch to soft clipping
                    #if args.soft_clip_primer:
                    #    if op!='D': cigar_str_new += 'S'

                    #if op!='D': removed += 1
                    #else: SAM[3] += 1
                    #continue
                #cigar_str_new += op
            cigar_str_new = idx*'S' + cigar_str[idx:]
            # shift alignment start position
            SAM[3] = int(SAM[3])+idx #args.primer_length

            # remove from 3' end
            cigar_str_new_rev = cigar_str_new[::-1]
            cigar_str_new = ''
            removed = 0
            idx = 0
            for op in cigar_str_new_rev:
                if args.soft_clip_primer:
                    if op=='D' or op=='I' or idx==offset_tail:
                        break
                    idx+=1
                #if removed < offset_tail:
                #    # EDIT: switch to soft clipping
                #    if args.soft_clip_primer:
                #        if op!='D': cigar_str_new += 'S'
                #    if op!='D': removed += 1
                #    continue
                #cigar_str_new += op
            cigar_str_new = idx*'S' + cigar_str_new_rev[idx:]
            cigar_str_new = cigar_str_new[::-1]

            #if not args.allow_end_deletion:
            # ensure no D at ends
            #    for idx, c in enumerate(cigar_str_new):
            #        if c=='D':
            #            sam[3] += 1
            #            continue
            #        else: break
            #    cigar_str_new = cigar_str_new[idx:]
            #    for idx, c in enumerate(cigar_str_new[::-1]):
            #        if c=='D': continue
            #        else: break
            #    cigar_str_new = cigar_str_new[:len(cigar_str_new)-idx]

            # convert back to list form
            cigar = cigar_expanded_string_to_list( cigar_str_new )

            SAM[3] = str(SAM[3])
            SAM[5] = ''.join( ['%d%s'%v for v in cigar] )
    
        raw += 1
        records += n_reads
        for n in range(n_reads):
            SAM_line = '\t'.join( [SAM[0]+'_%d'%(n+1)] + SAM[1:] )
            print(SAM_line, file=args.out)
    
    print('%d raw records'%raw, file=sys.stderr)
    print('%d SAM records written'%records, file=sys.stderr)


def command_line_interface(*args,**kw):
    ''' command line parsing function '''
    parser = argparse.ArgumentParser(description="A script for inflating grouped reads into BAM file.")

    parser.add_argument('sam_file',
                        type=argparse.FileType('r'),
                        help='input SAM file, or "-" for stdin')

    parser.add_argument('--max_variant',
                        type=int,
                        default=3,
                        help='maximum number of variants expected per amplicon')

    parser.add_argument('--soft_clip_primer',
                        action="store_true",
                        help='soft clip primer sequences if --primer_length is given. Default=False (i.e. hard clip)')

    #parser.add_argument('--allow_end_deletion',
    #                    action="store_true",
    #                    help='allow deletion operations at ends of the alignment')

    parser.add_argument('--out',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Output file, or "-" for stdout (default)')

    parser.add_argument('--primer_length',
                        type=int,
                        default=0,
                        help='size of flanking primer to be removed (intentionally kept for better BLAT alignment)')

    return parser.parse_args(*args,**kw)


if __name__=='__main__':
    args = command_line_interface()
    inflate( args )

