import os, sys, time, argparse, resource
from util import print_datetime
from dedup import Read, UnionFind
import numpy as np
from collections import deque, namedtuple


def parseArgument():
	parser = argparse.ArgumentParser('Usage')
	parser.add_argument('--num_shift', type=int, default=5, help='The maximum shift between duplicates')
	parser.add_argument('--num_mismatches', type=int, default=None, help='The maximum number of mismatches between duplicates')
	parser.add_argument('--num_mismatches_UMI', type=int, default=None, help='The maximum number of mismatches between duplicates')
	args = parser.parse_args()
	sys.stderr.write(str(args))
	sys.stderr.write('\n')
	sys.stderr.flush()
	return args


def process(num_shift, num_mismatches, num_mismatches_UMI):
	q = deque()
	q2 = deque()
	n_tot = 0
	l_q = 0
	UF = UnionFind()
	# reportSize, reportUnit = 1000, 'k'
	reportSize, reportUnit = 1000000, 'M'

	def write(r):
		if r.isWritten(): return True
		t = UF.pop(r.qname)
		if t == -1: return False
		if t > 0: r.setDup()
		s = r.write() + '\t'
		s += 'XD:i:' + ('0' if r.is_dup() else str(UF.get_total_size(r.qname, check_is_root=True)))
		s += '\n'
		sys.stdout.write(s)
		# sys.stdout.flush()
		return True

	for iline, line in enumerate(sys.stdin):
		if line.startswith('@'):
			sys.stdout.write(line)
			sys.stdout.flush()
			continue

		# sys.stderr.write(line)
		line = line.strip('\n').split('\t')
		line[0] = f'{line[0]}:{iline}'
		r = Read(line)
		if r.isSecondary() or r.isUmmapped(): assert write(r)
		else:
			for x in q:
				if r.samePosition(x, checkBC=True, num_shift=num_shift): break
				write(x)
			while len(q) and not r.samePosition(q[0], checkBC=True, num_shift=num_shift) and write(q[0]): q.popleft()
			q = deque(_ for _ in q if not _.isWritten())

			flag_wrote = False
			l_q += len(q)
			for x in q:
				if x.isWritten(): continue
				t = r.checkDup(x, num_shift=num_shift, num_mismatches=num_mismatches, num_mismatches_UMI=num_mismatches_UMI)
				if t == -1: continue
				if t == 0:
					flag_wrote = True
					r.setDup()
					assert write(r)
					break
				UF.merge(r.qname, x.qname)
			if not flag_wrote:
				q.append(r)

		n_tot += 1
		if n_tot % reportSize == 0:
			sys.stderr.write(f'{print_datetime()}Processed {n_tot//reportSize}{reportUnit} reads, avg queue length = {l_q/reportSize:.2f}, UF size reduced = {UF.clean()}, UF size = {len(UF.d)}\t{len(UF.obj)}\n')
			sys.stderr.flush()
			l_q = 0

	for x in q: write(x)
	while len(q) and write(q[0]): q.popleft()
	assert len(q) == 0


if __name__ == '__main__':
	sys.stderr.write(f'{print_datetime()}Begin\n')

	args = parseArgument()

	process(num_shift=args.num_shift, num_mismatches=args.num_mismatches, num_mismatches_UMI=args.num_mismatches_UMI)

	sys.stderr.write(f'{print_datetime()}Finished\n')
