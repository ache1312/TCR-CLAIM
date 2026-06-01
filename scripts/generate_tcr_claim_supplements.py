#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tcrbio.cli import supplements_main


if __name__ == "__main__":
    supplements_main()
