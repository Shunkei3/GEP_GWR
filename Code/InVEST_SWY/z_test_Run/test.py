# myutils.py
import numpy as np
def add(a, b):
    return a + b

def row_sums(x):
    """x: 2D array-like"""
    arr = np.asarray(x)
    return arr.sum(axis=1)