import gzip, os, sys, time, subprocess, io, gzip, datetime, itertools, argparse, pickle, re
import pandas as pd
import numpy as np
import multiprocessing as mp
import threading
from util import print_datetime

hash_base = 5
# hash_base = 4
hash_power = hash_base**np.arange(10)
def HASH(b): return (b % hash_base) @ hash_power[:len(b)]
# def HASH(b): return (b % 8 // 2) @ hash_power[:len(b)]

def matchBarcode(seq, bc, mode):
	if isinstance(seq, bytes): seq = np.frombuffer(seq, dtype=np.uint8)
	if bc.ndim == 1:
		return bc[HASH(seq)]
	else:
		d = ((seq[None] != bc) & (bc != ord('N'))).sum(1)
		idx = np.argmin(d)
		if (mode == 'unique' or d[idx] <= mode) and (d == d[idx]).sum() == 1:
			return idx
		else:
			return -1

def preprocessBarcode(barcode, mode):
	N, L = barcode.shape
	d = np.full(hash_base**L, -1)
	# TODO: use multi-source BFS to generalize this
	if mode == 0:
		for i, b in enumerate(barcode):
			print(b)
			d[HASH(b)] = i
	elif mode == 1:
		for i, b in enumerate(barcode):
			bb = np.copy(b)
			for idx in range(L):
				# for c in [65, 67, 71, 84]:  # np.frombuffer(bytes('ACGT','utf-8'), dtype=np.uint8)
				for c in [65, 67, 71, 78, 84]:	# np.frombuffer(bytes('ACGNT','utf-8'), dtype=np.uint8)
					bb[idx] = c
					d[HASH(bb)] = i
				bb[idx] = b[idx]
	elif mode == 'unique' or mode == 2:
		# TODO: optimize this by calculating the difference in a cumulative and categorical manner
		b = np.zeros(L, dtype=np.int64)

		def f(i):
			if i == L:
				d[HASH(b)] = matchBarcode(b, barcode, mode)
			else:
				for c in [65, 67, 71, 78, 84]:	# np.frombuffer(bytes('ACGNT','utf-8'), dtype=np.uint8)
					b[i] = c
					f(i + 1)

		f(0)
	else:
		raise NotImplementedError
	return d

def loadBarcode(filename, mode, preprocessed_barcode_file_prefix=None):
	with open(filename, 'r') as f: ret = [line.rstrip('\n').split('\t') for line in f]
	N = len(ret)
	wellID, barcode1, barcode2 = zip(*ret)
	barcode1 = np.array([np.frombuffer(bytes(_, 'utf-8'), dtype=np.uint8) for _ in barcode1])
	barcode2 = np.array([np.frombuffer(bytes(_, 'utf-8'), dtype=np.uint8) for _ in barcode2])
	barcodes = []
	for i, bc in [(1, barcode1), (2, barcode2)]:
		while True:
			if preprocessed_barcode_file_prefix is not None: filename = preprocessed_barcode_file_prefix + f'{i}.pkl'
			else: filename = None
			if filename is not None and os.path.exists(filename):
				sys.stderr.write(f'{print_datetime()}Loading preprocessed barcodes from {filename}\n')
				with open(filename, 'rb') as f: bc = pickle.load(f)
				break
			if N > 2:
				bc = preprocessBarcode(bc, mode)
				if filename is not None:
					with open(filename, 'wb') as f: pickle.dump(bc, f)
					sys.stderr.write(f'{print_datetime()}Saving preprocessed barcodes to {filename}\n')
				break
			else:
				break
		barcodes.append(bc)

		if bc.ndim == 1:
			sys.stderr.write(f'{print_datetime()}R{i}: {(bc != -1).sum()} out of {len(bc)} seq are valid\n')
			sys.stderr.flush()
	barcode1, barcode2 = barcodes
	return N, wellID, barcode1, barcode2


def process(
		wellID_list, barcode1_list, barcode2_list, pos1_list, pos2_list, mode_list,
		in_handles, out_handles, out_handles_u, barcode2keep,
):
	K = len(mode_list)
	num_total = 0
	num_matched = np.zeros([K+1, 3], dtype=np.uint)	# row: R1, R2, both; columns: K individual flags, AND of all flags
	while True:
		try:
			lines_read1, lines_read2 = [next(in_handles[0]) for _ in range(4)], [next(in_handles[1]) for _ in range(4)]
		except StopIteration:
			break
		flag_list = []
		idx_list = []
		for wellID, barcode1, barcode2, pos1, pos2, mode in zip(wellID_list, barcode1_list, barcode2_list, pos1_list, pos2_list, mode_list):
			seq1 = np.frombuffer(lines_read1[1], dtype=np.uint8)[pos1]
			seq2 = np.frombuffer(lines_read2[1], dtype=np.uint8)[pos2]
			idx1 = matchBarcode(seq1, barcode1, mode) if len(seq1) else None
			idx2 = matchBarcode(seq2, barcode2, mode) if len(seq2) else None
			# print(seq1.tostring().decode(), idx1)
			# print(seq2.tostring().decode(), idx2)
			"""
			>= 0: matched
			== -1: unmatched
			np.nan: not required
			"""

			flag_list.append((idx1 != -1, idx2 != -1, idx1 != -1 and idx2 != -1 and (idx1 == idx2 or idx1 is None or idx2 is None)))
			if flag_list[-1][-1]:
				if   idx1 is not None: idx = idx1
				elif idx2 is not None: idx = idx2
				else: raise NotImplementedError
			else: idx = -1
			idx_list.append(idx)

		flag_list = np.array(flag_list, dtype=bool)
		flag_list_all = flag_list.all(0)
		if flag_list_all[-1]:
			wellID = '_'.join(_[__] for _, __ in zip(wellID_list, idx_list))
			# print(wellID in barcode2keep)
			if barcode2keep is not None and wellID not in barcode2keep: flag_list_all[-1] = False
			# print(flag_list_all[-1])
			# assert False
		else: wellID = ''
		if flag_list_all[-1]:
			assert all(_ >= 0 for _ in idx_list)
			t = ('_' + wellID + '\n').encode('utf-8')
			out1 = b''.join([lines_read1[0].rstrip(b'\n').split(b' ')[0] + t] + lines_read1[1:])
			out2 = b''.join([lines_read2[0].rstrip(b'\n').split(b' ')[0] + t] + lines_read2[1:])
			out_handles[0].write(out1)
			out_handles[1].write(out2)
		else:
			out_handles[0].write(b'\n'*4)
			out_handles[1].write(b'\n'*4)
			if out_handles_u is not None:
				out1 = b''.join(lines_read1)
				out2 = b''.join(lines_read2)
				out_handles_u[0].write(out1)
				out_handles_u[1].write(out2)

		num_total += 1
		num_matched[:-1] += flag_list
		num_matched[ -1] += flag_list_all
		if num_total % 10000000 == 0:
		# if num_total % 100000 == 0:
		# if num_total % 1000 == 0:
		# if num_total % 1 == 0:
			sys.stderr.write(f'{print_datetime()}# = {num_total}\n')
			sys.stderr.write(re.sub(r'[\[\] ]', '', np.array2string(num_matched.T/num_total, separator='\t', formatter={'all':'{:.2e}'.format})) + '\n')
			sys.stderr.flush()
		# if num_total == 10000: break
	sys.stderr.write(f'{print_datetime()}# = {num_total}\n')
	sys.stderr.write(re.sub(r'[\[\] ]', '', np.array2string(num_matched.T, separator='\t')) + '\n')
	sys.stderr.flush()


def parseArgument():
	parser = argparse.ArgumentParser('Usage')
	parser.add_argument('--infile1', type=str, help='path/to/read1.fastq to be demultiplexed')
	parser.add_argument('--infile2', type=str, help='path/to/read2.fastq to be demultiplexed')
	parser.add_argument('--outfile1', type=str, help='path/to/read1.fastq, demultiplexed file')
	parser.add_argument('--outfile2', type=str, help='path/to/read2.fastq, demultiplexed file')
	parser.add_argument('--barcode', action='append', type=str, help='path/to/barcode file')
	# parser.add_argument('--preprocessed_barcode_file_prefix', action='append', type=str, help='preprocessed barcode')
	parser.add_argument('--pos1', action='append', type=eval, help='Position of barcode in read 1, e.g., slice(8,11)')
	parser.add_argument('--pos2', action='append', type=eval, help='Position of barcode in read 2')
	parser.add_argument('--mode', type=lambda _: 'unique' if _ == 'unique' else int(_), action='append', help='unique | <INT>')
	parser.add_argument('--unknown1', default=None, type=str, help='path/to/file of mismatched read pairs')
	parser.add_argument('--unknown2', default=None, type=str, help='path/to/file of mismatched read pairs')
	parser.add_argument('--barcode2keep', type=str, default=None, help='path/to/file that contains a list of barcodes to keep')
	args = parser.parse_args()
	sys.stderr.write(str(args))
	sys.stderr.write('\n')
	sys.stderr.flush()
	return args

if __name__ == '__main__':
	sys.stderr.write(f'{print_datetime()}Begin\n')

	args = parseArgument()

	wellID_list = []
	barcode1_list = []
	barcode2_list = []
	for barcode_file, mode in zip(args.barcode, args.mode):
		nbarcode, wellID, barcode1, barcode2 = loadBarcode(barcode_file, mode)
		wellID_list.append(wellID)
		barcode1_list.append(barcode1)
		barcode2_list.append(barcode2)
		sys.stderr.write(f'{print_datetime()}loaded {nbarcode} barcodes\n')
		sys.stderr.flush()

	if args.barcode2keep:
		with open(args.barcode2keep, 'r') as f: barcode2keep = {_.strip().replace('\t', '_') for _ in f}
		sys.stderr.write(f'{print_datetime()}Keep {len(barcode2keep)} barcodes\n')
	else:
		barcode2keep = None
		sys.stderr.write(f'{print_datetime()}Keep all barcodes\n')
	sys.stderr.flush()

	in_handles = np.array([open(f, 'rb') for f in [args.infile1, args.infile2]])
	out_handles = np.array([open(f, 'wb') for f in [args.outfile1, args.outfile2]])
	if args.unknown1 and args.unknown2:
		out_handles_u = np.array([open(f, 'wb') for f in [args.unknown1, args.unknown2]])
	else:
		out_handles_u = np.array([])

	sys.stderr.write(f'{print_datetime()}Launched subprocesses\n')
	sys.stderr.flush()

	process(
		wellID_list=wellID_list, barcode1_list=barcode1_list, barcode2_list=barcode2_list,
		pos1_list=args.pos1, pos2_list=args.pos2, mode_list=args.mode,
		in_handles=in_handles, out_handles=out_handles, out_handles_u=out_handles_u if len(out_handles_u) else None,
		barcode2keep=barcode2keep,
	)

	sys.stderr.write(f'{print_datetime()}Closing file handles\n')
	sys.stderr.flush()

	# for h in itertools.chain(in_handles.flat, out_handles.flat, out_handles_u.flat):
	for h in itertools.chain(out_handles.flat, out_handles_u.flat):
		h.flush()
		h.close()

	sys.stderr.write(f'{print_datetime()}Done\n')
	sys.stderr.flush()
