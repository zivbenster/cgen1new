# your ID here
import time
import sys
import argparse
from collections import defaultdict, Counter

import numpy as np

DNA = dict(A='T',T='A',G='C',C='G',N='N')
DEFAULT_COST = {(k1,k2): 1 if k2 != k1 else -1 for k1 in DNA for k2 in DNA}
DEFAULT_COST.update({('N',k): 0 for k in DNA})
DEFAULT_COST.update({(k,'N'): 0 for k in DNA})
GAP = '-'


def fasta_iter(file):
    with open(file) as fasta:
        hdr, seq = '', ''
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                if seq: yield hdr, seq
                hdr = line[1:]
                seq = ''
            else:
                seq += line
        if seq: yield hdr, seq  # yield last seuqence


def revcomp(seq):
    return ''.join(DNA[c] for c in seq[::-1])


def align(a, b, pcost=None, gcost=1, overlap=True):
    """
    :param a: a DNA string (with potentially unknown bases "N"s)
    :param b: a DNA string (with potentially unknown bases "N"s)
    :param pcost: a cost associated with every base pairing, default is 1 per mismatch (N's excluded) and -1 per match
    :param gcost: gap cost
    :param overlap: whether looking for an overlap alignment (i.e. no penalty for flanking gaps)
    :return: the alignment score, and a gapped string pair
    """
    if pcost is None: pcost = DEFAULT_COST
    A, B = len(a) + 1, len(b) + 1
    s, p = np.zeros((A, B)), {}
    s[0, :], s[:, 0] = np.arange(B) * gcost * (1 - overlap), np.arange(A) * gcost * (1 - overlap)
    for ai in range(1,A):
        for bi in range(1,B):
            up = s[ai - 1, bi] + gcost
            left = s[ai, bi - 1] + gcost
            diag = s[ai - 1, bi - 1] + pcost[(a[ai-1],b[bi-1])]
            if diag <= up and diag <= left:
                s[ai, bi], p[ai, bi] = diag, ((ai - 1, bi - 1), 'diag')
            elif left <= diag and left <= up:
                s[ai, bi], p[ai, bi] = left, ((ai, bi - 1), 'left')
            elif up <= diag and up <= left:
                s[ai, bi], p[ai, bi] = up, ((ai - 1, bi), 'up')

    if not overlap:
        b_, a_, (next, step), s_ = [], [], p[ai, bi], s[-1,-1]
    else:
        min_bottom = B - np.argmin(s[-1, -1:0:-1]) - 1
        min_right = A - np.argmin(s[-1:0:-1, -1]) - 1
        s_ = min(s[-1,min_bottom],s[min_right,-1])
        if s_ == s[-1,min_bottom]:
            # minimum was obtained when b is shifted B-min_bottom to the right
            ext = B - min_bottom
            b_, a_ = [b[-ext:][1:]], [GAP * (ext - 1)]
            next, step = p[A - 1, min_bottom]
        else:
            # minimum was obtained when a is shifted A-min_right to the right
            ext = A - min_right
            b_, a_ = [GAP * (ext - 1)], [a[-ext:][1:]]
            next, step = p[min_right, B - 1]

    # retrace solution from optimal point
    while True:
        if step == 'diag': _a, _b,  = a[ai-1], b[bi-1]
        elif step == 'left': _a, _b = GAP, b[bi-1]
        elif step == 'up': _a, _b = a[ai-1], GAP
        a_.append(_a)
        b_.append(_b)
        if next not in p: break
        ai, bi = next
        next, step = p[next]

    a_, b_ = ''.join(a_[::-1]), ''.join(b_[::-1])
    if next[1] == 0:
        b_, a_ = GAP * next[0] + b_, a[:next[0]] + a_
    elif next[0] == 0:
        b_, a_ = b[:next[1]] + b_, GAP * next[1] + a_
    return s_, (a_, b_)


def read_sampler(genome, n_reads=1000, l=30, pct_fwd=0.5, error_rate=0.01):
    """
    A generator for reads. Not efficient, but should work.
    :param genome: a map of names -> sequences
    :param n_reads: number of reads to sample
    :param l: read length
    :param pct_fwd: percent of reads from the forward (reference/+) strand
    :param error_rate: probability of an error
    :return: Yielding their sequence, with the reference name,
             position, strand, and number of errors introduced.
    """
    keys = list(genome.keys())
    if len(keys) == 1: keys = keys * 2
    dna = list(DNA.keys())
    for ni in range(n_reads):
        ref = np.random.choice(keys)
        pos = np.random.randint(0, len(genome[ref])-l)
        seq = genome[ref][pos:pos+l]
        is_fwd = np.random.rand(1) < pct_fwd
        if not is_fwd: seq = revcomp(seq)
        errs = np.nonzero(np.random.rand(l) < error_rate)[0]
        for i in errs:
            seq = seq[:i] + np.random.choice(dna) + seq[i+1:]
        yield seq, (ref, pos, int(is_fwd), len(errs))




class Aligner(object):

    def __init__(self, ref, k=15, seed_checks=3, max_edist=10):
        """
        build an index for k-mers of the reference sequence
        :param ref: an iterator over pairs of (name, iterator over DNA) for a genome. e.g. ("chr1", iter('AGAGGCGC...')), ("chr2",
        :param k: the size of seeds that will be used to search the index
        :param seed_checks: number of seeds tested before a read is discarded
        :param max_edist: maximum allowed edit distance between a read and the reference
        """
        self.k = k
        self.seed_checks = seed_checks
        self.max_edist = max_edist
        self.map = {}
        self.ref = {k:v for (k,v) in ref}
        for r in ref:
            name = r[0]
            seq = r[1]
            for i in range(0, len(seq)-k, 3):
                seed = ''.join(seq[i:i+k])
                if seed not in self.map:
                    self.map[seed] = [(name, i)]
                else:
                    self.map[seed].append((name, i))

    def __getitem__(self, kmer):
        """
        find occurences of the given kmer
        :param kmer: kmer to lookup in underlying reference
        :return: a sorted list of matches (from best to worst), each is a 4-tuple with:
                 - reference name
                 - left-most position of alignment (0-based indexing)
                 - strand: 1 for reference, 0 for reverse
                 - alignment score / edit distance
                 e.g.:
                 [('chr1', 644, 0, -30.0), ('chr2', 818, 0, -30.0)]
        """
        kmer = ''.join(kmer)
        krevc = revcomp(kmer)
        k = self.k
        results = []
        results += self.searchKmer(kmer, k, 1)
        results += self.searchKmer(krevc, k, 0)
        results = sorted(results, key=lambda res: res[3])
        return results

    def searchKmer(self, kmer, k, strand):
        results = []
        for j in range(len(kmer) - k):
            seed = kmer[j:j + k]
            try:
                idxs = self.map[seed]
            except:
                continue
            for idx in idxs:
                name = idx[0]
                i = idx[1]
                pref = max(0, i - 20)
                suff = min(i + len(kmer) + 20, len(self.ref[name]))
                ref_snip = self.ref[name][pref:suff]
                res = align(kmer, ref_snip)
                sub = 0
                while res[1][0][sub] == '-':
                    sub += 1
                index = pref + sub
                results.append((name, index, strand, res[0]))
        return results


def index_and_align(args):
    args.__dict__['output'] = sys.stdout if args.output is None else open(args.output,'w')
    idx = Aligner(fasta_iter(args.reference), k=args.k, max_edist=args.max_edist, seed_checks=args.seed_checks)

    args.output.write('\t'.join(['read_header','read_seq','refseq','position','is_forward', 'edist\n']))
    for h, q in fasta_iter(args.reads):
        alignments = idx[q]
        if alignments:
            for i, a in enumerate(alignments):
                if i < args.report_n:
                    args.output.write('\t'.join([h, q] + [str(x) for x in a]) + '\n')
        else:
            args.output.write('\t'.join([h, q] + ['*'] * 4) + '\n')


def sample_reads(args):
    args.__dict__['output'] = sys.stdout if args.output is None else open(args.output, 'w')
    genome = dict(fasta_iter(args.reference))
    for s, meta in read_sampler(genome, args.n_reads, args.read_len, args.pct_fwd, args.error_rate):
        args.output.write('>%s %i %i %i\n' % meta)
        args.output.write('%s\n' % s)


def test(args):
    if args.what == 'align': test_align()
    if args.what == 'index': test_index()
    if args.what == 'sample': test_sample()
    if args.what == 'all': test_all()


def test_align():
    def align_print_pair(a, b, ol=True):
        s, al = align(a, b, overlap=ol)
        print(s)
        print(al[0])
        print(al[1])

    while True:
        print("-" * 35)
        print("type any non-DNA letter to exit...")
        print("-"*35)
        a = input('string 1: ')
        if not all(i in DNA for i in a): break
        b = input('string 2: ')
        if not all(i in DNA for i in b): break
        print("==== overlap ====")
        align_print_pair(a, b, True)
        print("==== global  ====")
        align_print_pair(a, b, False)


def test_index():
    L = int(1e6)
    a = np.asarray(list(set(DNA.keys()) - set(['N'])))
    ref = [('random', np.random.choice(a, L))]
    mn = int(np.ceil(np.log(L) / 2))
    for k in range(mn+3, mn * 2):
        start = time.time()
        idx = Aligner(ref, k)
        end = time.time()
        hist = Counter()
        for v in idx.map.values(): hist[len(v)] += 1
        sys.stderr.write("K = %i (%f sec)\n" % (k, end - start))
        sys.stderr.write('   collisions\tcount\n')
        for occ, cnt in sorted(tuple(hist.items()), key=lambda k: k[0]):
            sys.stderr.write(''.join(['  ', str(occ - 1), '\t\t', str(cnt), '\n']))


def test_sample():
    a = np.asarray(list(set(DNA.keys()) - set(['N'])))
    genome = {'test1': ''.join(np.random.choice(a, 1000)),
              'test2': ''.join(np.random.choice(a, 1000))}
    sys.stderr.write('=' * 80 + '\n')
    sys.stderr.write('both strands\n')
    sys.stderr.write('=' * 80 + '\n')
    for s, meta in read_sampler(genome, n_reads=20):
        sys.stderr.write('> %s %i %i %i\n' % meta)
        sys.stderr.write(s+'\n')
    sys.stderr.write('=' * 80 + '\n')
    sys.stderr.write('only reference (forward) reads\n')
    sys.stderr.write('=' * 80 + '\n')
    for s, meta in read_sampler(genome, n_reads=20, pct_fwd=1):
        sys.stderr.write('> %s %i %i %i\n' % meta)
        sys.stderr.write(s + '\n')
    sys.stderr.write('=' * 80 + '\n')
    sys.stderr.write('20% errors:\n')
    sys.stderr.write('=' * 80 + '\n')
    for s, meta in read_sampler(genome, n_reads=20, error_rate=.2):
        sys.stderr.write('> %s %i %i %i\n' % meta)
        sys.stderr.write(s + '\n')


def test_all():
    a = np.asarray(list(set(DNA.keys()) - set(['N'])))
    genome = {'chr1': ''.join(np.random.choice(a, 1000)),
              'chr2': ''.join(np.random.choice(a, 1000))}
    idx = Aligner(genome.items())

    def sample_and_print(sampler, header):
        no_aln, wrong_aln = 0, 0
        for i, (read, meta) in enumerate(sampler):
            alns = idx[read]
            if not alns:
                no_aln += 1
                continue
            alns = alns[0]  # just check top alignment
            if alns[:3] != meta[:3]: wrong_aln += 1
        sys.stderr.write('=' * 80 + '\n')
        sys.stderr.write(header + "\n")
        sys.stderr.write('Total: %i, No alignments: %i, Erroneous alignments: %i\n' % (i, no_aln, wrong_aln))

    sample_and_print(read_sampler(genome, error_rate=.001), 'error rate = 0.001, read length = 30')
    sample_and_print(read_sampler(genome, error_rate=.01), 'error rate = 0.01, read length = 30')
    sample_and_print(read_sampler(genome, error_rate=.05), 'error rate = 0.1, read length = 30')
    sample_and_print(read_sampler(genome, error_rate=.05, l=60), 'error rate = 0.1, read length = 60')
    sample_and_print(read_sampler(genome, error_rate=.1), 'error rate = 0.25, read length = 30')
    sample_and_print(read_sampler(genome, error_rate=.1, l=60), 'error rate = 0.25, read length = 60')
    sample_and_print(read_sampler(genome, error_rate=.1, l=90), 'error rate = 0.25, read length = 90')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sps = p.add_subparsers()

    sp = sps.add_parser('test', help='test sub functionality')
    sp.add_argument('what', help='which functionality to test?', choices=['align','index','sample', 'all'])
    sp.set_defaults(func=test)

    sp = sps.add_parser('sample', help='sample reads from reference into fasta format. '
                                      'headers contain ref name, position, strand, and #errors')
    sp.add_argument('reference', help='reference fasta file')
    sp.add_argument('--n_reads', '-n', help='number of reads to sample', default=1000, type=int)
    sp.add_argument('--read_len', '-l', help='read length', default=50, type=int)
    sp.add_argument('--error_rate', '-e', help='base error rate', default=.01, type=float)
    sp.add_argument('--output', '-o', default=None, help='file path for output. default is stdout')
    sp.add_argument('--pct_fwd', '-p', help='fraction of reads that will come from reference strand',
                    default=0.5, type=float)
    sp.set_defaults(func=sample_reads)

    sp = sps.add_parser('align', help='align reads to reference')
    sp.add_argument('reference', help='reference fasta file')
    sp.add_argument('reads', help='a fasta file to align')
    sp.add_argument('--k', help='seed length', default=12)
    sp.add_argument('--max_edist', '-e', help='maximum edit distance allowed for an alignment', default=3)
    sp.add_argument('--error_rate', '-er', help='estimated error rate', default=.01)
    sp.add_argument('--seed_checks', '-s', help='number of seeds to check per read', default=3)
    sp.add_argument('--report_n', '-n', help='how many alignments to report (at most)', default=1)
    sp.add_argument('--output', '-o', default=None, help='file path for output. default is stdout')
    sp.set_defaults(func=index_and_align)

    args = p.parse_args()
    args.func(args)