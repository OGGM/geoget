import six
from distutils.version import LooseVersion
import os
import sys
import unittest
import logging
import matplotlib
import numpy as np
from six.moves.urllib.request import urlopen
from six.moves.urllib.error import URLError


# Defaults
logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

# test dirs
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
TESTDIR_BASE = os.path.join(CURRENT_DIR, 'tmp')

# Some logic to see which environment we are running on

# Some control on which tests to run (useful to avoid too long tests)
# defaults everywhere else than travis
ON_TRAVIS = False
RUN_SLOW_TESTS = False
RUN_DOWNLOAD_TESTS = True

if sys.version_info < (3, 5):
    # Minimal tests
    RUN_SLOW_TESTS = False

# quick n dirty method to see if internet is on
try:
    _ = urlopen('http://www.google.com', timeout=1)
    HAS_INTERNET = True
except URLError:
    HAS_INTERNET = False


def requires_internet(test):
    # Test decorator
    msg = 'requires internet'
    return test if HAS_INTERNET else unittest.skip(msg)(test)


def requires_py3(test):
    # Test decorator
    msg = "requires python3"
    return unittest.skip(msg)(test) if six.PY2 else test


def is_slow(test):
    # Test decorator
    msg = "requires explicit environment for slow tests"
    return test if RUN_SLOW_TESTS else unittest.skip(msg)(test)


def is_download(test):
    # Test decorator
    msg = "requires explicit environment for download tests"
    return test if RUN_DOWNLOAD_TESTS else unittest.skip(msg)(test)


# the code below is copy/pasted from xarray
# TODO: go back to xarray when https://github.com/pydata/xarray/issues/754
def assertEqual(a1, a2):
    assert a1 == a2 or (a1 != a1 and a2 != a2)


def decode_string_data(data):
    if data.dtype.kind == 'S':
        return np.core.defchararray.decode(data, 'utf-8', 'replace')


def data_allclose_or_equiv(arr1, arr2, rtol=1e-05, atol=1e-08):
    from xarray.core import ops

    if any(arr.dtype.kind == 'S' for arr in [arr1, arr2]):
        arr1 = decode_string_data(arr1)
        arr2 = decode_string_data(arr2)
    exact_dtypes = ['M', 'm', 'O', 'U']
    if any(arr.dtype.kind in exact_dtypes for arr in [arr1, arr2]):
        return ops.array_equiv(arr1, arr2)
    else:
        return ops.allclose_or_equiv(arr1, arr2, rtol=rtol, atol=atol)


def assertVariableAllClose(v1, v2, rtol=1e-05, atol=1e-08):
    assertEqual(v1.dims, v2.dims)
    allclose = data_allclose_or_equiv(
        v1.values, v2.values, rtol=rtol, atol=atol)
    assert allclose, (v1.values, v2.values)


def assertDatasetAllClose(d1, d2, rtol=1e-05, atol=1e-08):
    assertEqual(sorted(d1, key=str), sorted(d2, key=str))
    for k in d1:
        v1 = d1.variables[k]
        v2 = d2.variables[k]
        assertVariableAllClose(v1, v2, rtol=rtol, atol=atol)
