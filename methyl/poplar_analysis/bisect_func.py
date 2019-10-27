import bisect

def bisect_ge(a, x):
	# leftmost (first) INDEX greater than or equal to x
	# return None if not found
	i = bisect.bisect_left(a, x)
	if i != len(a):
		return i
	return None

def bisect_lt(a, x):
	# rightmost (last) INDEX less than x
	# return None if not found
	i = bisect.bisect_left(a, x)
	if i:
		return i-1
	return None

def bisect_index(a, x):
	try:
		i = bisect.bisect_left(a, x)
		if i != len(a) and a[i] == x:
			return i
		return None
	except TypeError:
		if x in a:
			return a.index(x)
		else:
			return None
