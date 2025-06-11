from types import ModuleType
import sys

pyplot = ModuleType('matplotlib.pyplot')

def _func(*args, **kwargs):
    return None

pyplot.__getattr__ = lambda name: _func
sys.modules[__name__ + '.pyplot'] = pyplot
