#!/usr/bin/env python3
"""Compatibility entry for running this source checkout."""

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

# The checkout command delegates to the same CLI as an installed wheel.
# pylint: disable=wrong-import-position
from omics_workflows.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
