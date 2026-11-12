from __future__ import annotations

from importlib.metadata import version

import numpy as np
from pyteomics import mzml, mzxml

import eic_smoothing_simulator as sim


def test_installed_versions_match_release_metadata() -> None:
    assert version("eicsmooth") == sim.__version__ == "1.0"
    assert version("numpy") == np.__version__ == "1.26.4"
    assert version("matplotlib") == "3.9.3"
    assert version("pyteomics") == "5.0.1"


def test_pyteomics_xml_readers_import() -> None:
    assert mzml.MzML is not None
    assert mzxml.MzXML is not None
