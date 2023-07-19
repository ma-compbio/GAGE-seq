import os, sys, time, argparse, itertools
from util import print_datetime


# 1101:24434:1344	mm10_chr10	63155823	mm10_chr10	65318943	-	+	UU	G	G	B10,E1	F


def parseArgument():
	parser = argparse.ArgumentParser('Usage')
	args = parser.parse_args()
	sys.stderr.write(str(args))
	sys.stderr.write('\n')
	sys.stderr.flush()
	return args


def process():
	q = list()
	n_tot = 0
	l_q = 0
	# reportSize, reportUnit = 1000, 'k'
	reportSize, reportUnit = 1000000, 'M'

	def write(readid, wellid):
		nonlocal q
		# print('='*20 + ' write ' + '='*20)
		read_pairs = {(s, t) for s, t in q if s[0] != '!' and t[0] != '!' and s != t}
		read_pairs |= {(t, s) for s, t in read_pairs}
		reads = {r for r in itertools.chain.from_iterable(q) if r[0] != '!'}
		q.clear()
		n = len(read_pairs) // 2
		for (i, s), (j, t) in itertools.combinations(enumerate(reads), 2):
			if (s, t) in read_pairs or (t, s) in read_pairs: continue
			sys.stdout.write('\t'.join([
				readid, s[0], str(s[1]), t[0], str(t[1]), s[2], t[2], s[3] + t[3], s[4], t[4],
				wellid,
			]) + '\n')
			n += 1
		assert n == len(reads) * (len(reads) - 1) // 2, (n, len(reads), readid)

	prev_readid, readid, wellid = '', '', ''
	for iline, line in enumerate(sys.stdin):
		if line.startswith('@') or line.startswith('#'):
			sys.stdout.write(line)
			continue
		if n_tot % reportSize == 0:
			sys.stderr.write(
				f'{print_datetime()}Processed {n_tot//reportSize}{reportUnit} reads \t'
				f'avg queue length = {l_q/reportSize:.2f}\n')
			sys.stderr.flush()
			l_q = 0
		line_str = line
		line = line.strip('\n').split('\t')

		readid = line[0]
		if prev_readid != readid: write(prev_readid, wellid)
		sys.stdout.write(line_str)
		prev_readid = readid
		wellid = line[10]
		r1 = (line[1], int(line[2]), line[5], line[7][0], line[8])
		r2 = (line[3], int(line[4]), line[6], line[7][1], line[9])
		q.append((r1, r2))

		l_q += len(q)
		n_tot += 1

	write(readid, wellid)

	# sys.stderr.write(f'{print_datetime()}Finished\n')


if __name__ == '__main__':
	sys.stderr.write(f'{print_datetime()}Begin\n')

	args = parseArgument()

	process()

	sys.stderr.write(f'{print_datetime()}Finished\n')
