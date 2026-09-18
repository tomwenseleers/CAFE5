#!/usr/bin/env python3
"""Compatibility entry point for the current nominal refitted-bootstrap workflow.

Output is bootstrap/branch_tests.tsv (signed posterior mean change is primary).
See docs/focal_bootstrap_methods.md for the output contract and assumptions.
"""
from refitted_bootstrap import main

if __name__ == '__main__':
    main()
