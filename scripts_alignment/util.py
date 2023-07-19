import os, sys, time, argparse, resource, datetime
import numpy as np
from collections import deque, namedtuple


def print_datetime():
	return datetime.datetime.now().strftime('[%Y/%m/%d %H:%M:%S]\t')


def estimate_lib_size(n, c):
	"""
	from https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/DuplicationMetrics.java#L132
	C / X = 1 - exp( -N / X )
	f(x) = c x - 1 + exp(-n x)
	f'(x) = c - n exp(-n x)
	:param n: # of observed reads
	:param c: # of observed unique reads
	:return: # of unique reads
	"""
	f = lambda n, c, x: c / x - 1 + np.exp(-n / x)
	# decreasing in x
	l = c
	r = c
	while f(n, c, r) > 0:
		r *= 2
	while r - l > 1:
		m = (l + r) // 2
		if f(n, c, m) > 0:
			l = m
		else:
			r = m
	return l


if __name__ == '__main__':
	print(estimate_lib_size(960538255, 522019650))
