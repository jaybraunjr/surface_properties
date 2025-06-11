import os
import sys

# Add stubs directory so tests can import placeholder packages
STUBS_DIR = os.path.join(os.path.dirname(__file__), 'stubs')
if os.path.isdir(STUBS_DIR) and STUBS_DIR not in sys.path:
    sys.path.insert(0, STUBS_DIR)

