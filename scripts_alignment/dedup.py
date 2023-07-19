import os, sys, time, argparse, resource
from util import print_datetime
import numpy as np
from collections import deque, namedtuple, defaultdict


class UnionFind:
	# key must NOT be negative integers
	def __init__(self):
		self.d = defaultdict(lambda: -1)
		self.obj = {}
		self.total_size = defaultdict(lambda: 1)

	def find(self, x):
		if x not in self.d: return x
		t = []
		while x in self.d and (type(self.d[x]) != int or self.d[x] >= 0):
			t.append(x)
			x = self.d[x]
		for _ in t: self.d[_] = x
		return x

	def merge(self, x, y):
		x = self.find(x)
		y = self.find(y)
		if x == y: return -self.d[x]
		if self.d[y] < self.d[x]: x, y = y, x
		self.d[x] += self.d[y]; self.d[y] = x
		self.total_size[x] += self.total_size[y]
		return -self.d[x]

	def pop(self, x, isRoot=None):
		if x not in self.d:
			self.obj.pop(x, None)
			self.total_size.pop(x, None)
			return 0
		r = self.find(x)
		if x == r:
			if self.d[r] < -1: return -1
			self.obj.pop(r, None)
			self.total_size.pop(r, None)
			del self.d[r]
			return 0
		elif isRoot is not None and not isRoot(self.obj[r], self.obj[x]):
			self.d[x] = self.d[r]
			self.d[r] = x
			self.total_size[x] = self.total_size[r]
			return -1
		# del self.d[x]
		self.obj.pop(x, None)
		self.total_size.pop(x, None)
		self.d[r] += 1
		return -self.d[r]

	def attachObj(self, x, obj):
		self.obj[x] = obj

	def get_total_size(self, x, check_is_root=False):
		assert (not check_is_root) or self.find(x) == x
		return self.total_size.get(x, 1)

	def clean(self):
		keys2remove = []
		for x in self.d.keys():
			r = self.find(x)
			if r not in self.d:
				keys2remove.append(x)
		for k in keys2remove:
			del self.d[k]
			assert k not in self.obj
			assert k not in self.total_size
		return len(keys2remove)


def checkDup(s1, s2, n, l=-1, k=20):
	c = 0
	if l == -1: l = min(len(s1), len(s2))
	for i in range(0, l, k):
		j = min(i+k, l)
		c += sum(_ != __ and _ != 'N' and __ != 'N' for _, __ in zip(s1[i:j], s2[i:j]))
		if c > n: return -1
	c = sum(_ != __ for _, __ in zip(s1[:l], s2[:l]))
	return c


class Read:
	def __init__(self, cols):
		# for i in [1, 3, 4]: cols[i] = int(cols[i])
		cols[1] = int(cols[1])
		# python -c "for i, k in enumerate(['qname', 'flag', 'chrom', 'pos', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'seq', 'qual']): print(f'\t\tself.{k} = cols[{i}]')"
		self.keys = ['qname', 'flag', 'chrom', 'pos', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'seq', 'qual', 'attrs']
		# for k, v in zip(self.keys, cols[:11] + [cols[11:]]): setattr(self, k, v)
		self.qname = cols[0]
		self.flag = cols[1]
		self.chrom = cols[2]
		self.pos = int(cols[3])
		self.MAPQ = cols[4]
		self.CIGAR = cols[5]
		self.RNEXT = cols[6]
		self.PNEXT = cols[7]
		self.TLEN = cols[8]
		self.seq = cols[9]
		self.qual = cols[10]
		self.attrs = cols[11:]
		self.barcode = next((a[5:] for a in self.attrs if a.startswith('BC:Z:')), None)
		self.UMI = next((a[5:] for a in self.attrs if a.startswith('RX:Z:')), None)
		self.flag_isWritten = False

	def isUnique(self):
		if self.attrs[0].startswith('NH:i:'): return int(self.attrs[0][5:]) == 1
		if any(attr.startswith('XA:Z:') for attr in self.attrs): return False
		return True

	def isUmmapped(self): return (self.flag & 0x4) == 0x4

	def isSecondary(self): return (self.flag & 0x100) == 0x100

	def getStrand(self): return '-' if self.flag & 0x10 else '+'

	def samePosition(self, r, checkBC=True, num_shift=0):
		if checkBC: assert self.barcode is not None and r.barcode is not None
		return self.chrom == r.chrom and abs(self.pos - r.pos) <= num_shift and (not checkBC or self.barcode == r.barcode)
		# return self.chrom == r.chrom and self.getStrand() == r.getStrand() and abs(self.pos == r.pos) <= num_shift and \
		# 	(not checkBC or self.barcode == r.barcode)

	def checkDup(self, r, checkBC=True, num_shift=0, num_mismatches=None, num_mismatches_UMI=None):
		ret = 0
		if self.getStrand() != r.getStrand(): return -1
		assert self.chrom == r.chrom
		# assert self.getStrand() == r.getStrand()
		assert self.pos >= r.pos
		if checkBC:
			assert self.barcode is not None and r.barcode is not None
			if self.barcode != r.barcode: return -1
		if num_mismatches_UMI is not None:
			assert self.UMI is not None and r.UMI is not None
			c = checkDup(self.UMI, r.UMI, n=num_mismatches_UMI)
			if c == -1: return -1
			ret = max(ret, c)
		if num_shift is not None:
			c = abs(self.pos - r.pos)
			if c > num_shift: return -1
			ret = max(ret, c)
		if num_mismatches is not None:
			c = checkDup(self.seq, r.seq, n=num_mismatches)
			if c == -1: return -1
			ret = max(ret, c)
		return ret

	def setDup(self): self.flag |= 0x400

	def is_dup(self): return (self.flag & 0x400) != 0

	def isWritten(self): return self.flag_isWritten

	def write(self, delimiter='\t', setFlag=True):
		if setFlag: self.flag_isWritten = True
		return delimiter.join(map(str, [getattr(self, k) for k in self.keys[:-1]] + self.attrs))


class ReadPairDNADNA:
	def __init__(self, cols0, cols1, cols2):
		self.keys = ['chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type']
		self.readID = cols0[0]
		self.chrom1 = cols0[1]
		self.pos1 = cols0[2]
		self.chrom2 = cols0[3]
		self.pos2 = cols0[4]
		self.strand1 = cols0[5]
		self.strand2 = cols0[6]
		self.pair_type = cols0[7]
		self.read1 = Read(cols1)
		self.read2 = Read(cols2)

	def samePosition(self, p, checkBC=False):
		return all(getattr(self, k) == getattr(p, k) for k in self.keys) and \
			(not checkBC or (self.read1.barcode == p.read1.barcode and self.read2.barcode == p.read2.barcode))

	def checkDup(self, p, **kwargs):
		c1 = self.read1.checkDup(p.read1, **kwargs)
		if c1 == -1: return -1
		c2 = self.read2.checkDup(p.read2, **kwargs)
		if c2 == -1: return -1
		self.pair_type = 'DD'
		return max(c1, c2)

	def write(self, delimiter='\t', delimiter_sam='\031'):
		return delimiter.join([getattr(self, k) for k in self.keys] + [self.read1.write(delimiter_sam), self.read2.write(delimiter_sam)])


class ReadPairDNADNAShort:
	def __init__(self, cols):
		try:
			self.keys = ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type', 'barcode']
			self.readID = cols[0]
			self.chrom1 = cols[1]
			self.pos1 = int(cols[2])
			self.chrom2 = cols[3]
			self.pos2 = int(cols[4])
			self.strand1 = cols[5]
			self.strand2 = cols[6]
			self.pair_type = cols[7]
			self.seq1 = cols[8]
			self.seq2 = cols[9]
			self.barcode = cols[10]
			self.flipped = cols[11] == 'F'
			self.flag_isWritten = False
			assert self.chrom1 != self.chrom2 or self.pos2 >= self.pos1
		except:
			sys.stderr.write('error\n')
			sys.stderr.write(str(cols) + '\n')

	def isUpstreamOf(self, p):
		if self.barcode is not None: assert self.barcode == p.barcode
		assert self.chrom1 == p.chrom1 and self.chrom2 == p.chrom2
		return self.chrom1 != self.chrom2 or self.pos2-self.pos1 >= p.pos2-p.pos1

	def samePosition(self, p, checkBC=True, checkR1=True, checkR2=False, num_shift_R1=5, num_shift_R2=5):
		if checkBC: assert self.barcode is not None and p.barcode is not None
		return (not checkBC or self.barcode == p.barcode) and \
			(not checkR1 or (self.chrom1 == p.chrom1 and np.abs(self.pos1-p.pos1) <= num_shift_R1 and self.strand1 == p.strand1)) and \
			(not checkR2 or (self.chrom2 == p.chrom2 and np.abs(self.pos2-p.pos2) <= num_shift_R2 and self.strand2 == p.strand2))

	def checkDup(self, p, num_shift_R1=5, num_shift_R2=5, num_shift_R12=-1, num_mismatches_R1=-1, num_mismatches_R2=-1, num_mismatches_R12=-1):
		assert self.chrom1 == p.chrom1 and self.strand1 == p.strand1
		# if self.chrom1 != p.chrom1 or self.strand1 != p.strand1: return -1
		# assert self.chrom1 == p.chrom1 and self.strand1 == p.strand2
		if self.chrom2 != p.chrom2 or self.strand2 != p.strand2: return -1
		n1 = self.pos1-p.pos1
		n2 = np.abs(self.pos2-p.pos2)
		assert n1 >= 0
		if self.flipped == p.flipped:
			if self.flipped:
				if n1 > num_shift_R2 or n2 > num_shift_R1: return -1
			else:
				if n1 > num_shift_R1 or n2 > num_shift_R2: return -1
		else:
			if num_shift_R12 == -1 or n1 > num_shift_R12 or n2 > num_shift_R12: return -1
		ret = max(n1, n2)

		if num_mismatches_R1 != -1:
			c1 = checkDup(self.seq1, p.seq1, n=num_mismatches_R1)
			if c1 == -1: return -1
			ret = max(ret, c1)
		if num_mismatches_R1 != -1:
			c2 = checkDup(self.seq1, p.seq1, n=num_mismatches_R2)
			if c2 == -1: return -1
			ret = max(ret, c2)

		return ret

	def setDup(self): self.pair_type = 'DD'

	def is_dup(self): return self.pair_type == 'DD'

	def isWritten(self): return self.flag_isWritten

	def write(self, delimiter='\t', setFlag=True):
		if setFlag: self.flag_isWritten = True
		return delimiter.join([str(getattr(self, k)) for k in self.keys])

