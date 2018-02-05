#!/usr/bin/python3
import os

import pytest

from app import nussinov
from defs import ROOT_DIR


def test_simple_strand():
    test_strand = 'gGGACCUUCCCGGUCUC'.upper()

    base_pairs = nussinov.nussinov_algorithm(test_strand)

    assert (0, 16) == base_pairs[0]
