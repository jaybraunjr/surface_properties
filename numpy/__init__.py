import math
import random as _random

nan = float('nan')

class ndarray:
    def __init__(self, data):
        if isinstance(data, ndarray):
            self.data = data.data
        else:
            self.data = self._convert(data)
    def _convert(self, d):
        if isinstance(d, ndarray):
            return d.data
        if isinstance(d, (list, tuple)):
            return [self._convert(x) for x in d]
        return d
    def __len__(self):
        return len(self.data)
    def __iter__(self):
        for x in self.data:
            yield x
    def _proc_index(self, idx, length):
        if isinstance(idx, slice):
            return list(range(*idx.indices(length)))
        if isinstance(idx, (list, tuple)):
            return list(idx)
        if isinstance(idx, ndarray):
            if len(idx) and isinstance(idx.data[0], bool):
                return [i for i, v in enumerate(idx.data) if v]
            return [int(v) for v in idx.data]
        return [idx]
    def __getitem__(self, key):
        if isinstance(key, tuple):
            rows, cols = key
            row_idx = self._proc_index(rows, len(self.data))
            if isinstance(cols, slice):
                col_range = list(range(*cols.indices(len(self.data[0]))))
                return ndarray([[self.data[r][c] for c in col_range] for r in row_idx])
            else:
                col_idx = self._proc_index(cols, len(self.data[0]))
                if len(col_idx) == 1:
                    return ndarray([self.data[r][col_idx[0]] for r in row_idx])
                return ndarray([[self.data[r][c] for c in col_idx] for r in row_idx])
        result = self.data[key]
        if isinstance(result, list):
            return ndarray(result)
        return result
    def _binary_op(self, other, op):
        if isinstance(other, ndarray):
            other = other.data
        if isinstance(self.data, list):
            return ndarray([op(a, b) for a, b in zip(self.data, other)])
        return op(self.data, other)
    def __sub__(self, other):
        return self._binary_op(other, lambda a,b: a-b)
    def __add__(self, other):
        return self._binary_op(other, lambda a,b: a+b)
    def __mul__(self, other):
        return self._binary_op(other, lambda a,b: a*b)
    def __truediv__(self, other):
        return self._binary_op(other, lambda a,b: a/b)
    __rsub__ = lambda self, other: ndarray(other).__sub__(self)
    __radd__ = lambda self, other: ndarray(other).__add__(self)
    __rmul__ = lambda self, other: ndarray(other).__mul__(self)
    __rtruediv__ = lambda self, other: ndarray(other).__truediv__(self)
    def __gt__(self, other):
        return self._binary_op(other, lambda a,b: a>b)
    def __lt__(self, other):
        return self._binary_op(other, lambda a,b: a<b)
    def __ge__(self, other):
        return self._binary_op(other, lambda a,b: a>=b)
    def __le__(self, other):
        return self._binary_op(other, lambda a,b: a<=b)
    @property
    def shape(self):
        if len(self.data)==0:
            return (0,)
        if isinstance(self.data[0], list):
            return (len(self.data), len(self.data[0]))
        return (len(self.data),)
    def astype(self, typ):
        return ndarray([typ(x) for x in self.data])
    def tolist(self):
        return [x.tolist() if isinstance(x, ndarray) else x for x in self.data]
    def __repr__(self):
        return f"ndarray({self.data})"

def array(obj):
    return ndarray(obj)

def _flatten(obj):
    if isinstance(obj, ndarray):
        obj = obj.data
    if isinstance(obj, (list, tuple)):
        for x in obj:
            yield from _flatten(x)
    else:
        yield float(obj)

def mean(a, axis=None):
    if isinstance(a, ndarray):
        data = a.data
    else:
        data = a
    if axis is None:
        vals = list(_flatten(data))
        return sum(vals)/len(vals) if vals else nan
    if axis == 0:
        n = len(data)
        if n==0:
            return ndarray([])
        m = len(data[0])
        return ndarray([sum(float(row[i]) for row in data)/n for i in range(m)])
    raise NotImplementedError

def std(a, axis=None):
    if isinstance(a, ndarray):
        data = a.data
    else:
        data = a
    if axis is None:
        vals = list(_flatten(data))
        m = mean(vals)
        return math.sqrt(sum((x-m)**2 for x in vals)/len(vals)) if vals else nan
    if axis == 0:
        n = len(data)
        if n==0:
            return ndarray([])
        mvals = [sum(float(row[i]) for row in data)/n for i in range(len(data[0]))]
        return ndarray([math.sqrt(sum((float(row[i])-mvals[i])**2 for row in data)/n) for i in range(len(mvals))])
    raise NotImplementedError

def zeros(n):
    return ndarray([0]*n)

def linspace(start, stop, num):
    if num==1:
        return ndarray([float(start)])
    step = (float(stop)-float(start))/(num-1)
    return ndarray([start + i*step for i in range(num)])

class _Random:
    def rand(self, *shape):
        if len(shape)==0:
            return _random.random()
        if len(shape)==1:
            return ndarray([_random.random() for _ in range(shape[0])])
        if len(shape)==2:
            return ndarray([[_random.random() for _ in range(shape[1])] for __ in range(shape[0])])
        raise NotImplementedError
random = _Random()

class _Linalg:
    def norm(self, x, axis=None):
        if isinstance(x, ndarray):
            data = x.data
        else:
            data = x
        if axis is None:
            return math.sqrt(sum(float(i)**2 for i in data))
        if axis in (1,-1):
            return ndarray([math.sqrt(sum(float(v)**2 for v in row)) for row in data])
        raise NotImplementedError
linalg = _Linalg()

def dot(a,b):
    if isinstance(a, ndarray):
        a = a.data
    if isinstance(b, ndarray):
        b = b.data
    return sum(float(x)*float(y) for x,y in zip(a,b))

def repeat(a, repeats, axis=0):
    if isinstance(a, ndarray):
        data = a.data
    else:
        data = a
    if isinstance(repeats, int):
        repeats = [repeats]*len(data)
    out = []
    for item, r in zip(data, repeats):
        for _ in range(r):
            out.append(item)
    return ndarray(out)

def unique(a, return_counts=False):
    if isinstance(a, ndarray):
        data = a.data
    else:
        data = a
    counts = {}
    for x in data:
        counts[x] = counts.get(x,0)+1
    uniq = list(counts.keys())
    if return_counts:
        return ndarray(uniq), ndarray([counts[x] for x in uniq])
    return ndarray(uniq)

def histogram(pos, weights=None, bins=None):
    n = len(bins)-1
    hist = [0.0]*n
    for i,val in enumerate(pos):
        w = weights[i] if weights else 1
        for j in range(n):
            if (val >= bins[j] and val < bins[j+1]) or (val==bins[-1] and j==n-1):
                hist[j] += w
                break
    return ndarray(hist), ndarray(bins)

class _Char:
    def startswith(self, arr, prefix):
        if isinstance(arr, ndarray):
            arr = arr.data
        return ndarray([str(x).startswith(prefix) for x in arr])
char = _Char()

def isin(elements, test_elements, invert=False):
    if isinstance(elements, ndarray):
        elements = elements.data
    if isinstance(test_elements, ndarray):
        test_elements = test_elements.data
    test_set = set(test_elements)
    mask = [el in test_set for el in elements]
    if invert:
        mask = [not x for x in mask]
    return ndarray(mask)

def transpose(arr):
    if isinstance(arr, ndarray):
        arr = arr.data
    return ndarray([list(x) for x in zip(*arr)])

def savetxt(fname, X):
    if isinstance(X, ndarray):
        X = X.tolist()
    with open(fname, 'w') as f:
        for row in X:
            if isinstance(row, (list, tuple)):
                f.write(" ".join(str(v) for v in row)+"\n")
            else:
                f.write(str(row)+"\n")

def isscalar(obj):
    return not isinstance(obj, (list, tuple, ndarray))
