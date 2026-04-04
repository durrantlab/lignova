r"""Test unichem API."""

import os

from lignova.APIs import UniChemAPI

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_unichem = {
    "write_dir": "./tmp/unichem",
}


def prep_dirs():
    os.makedirs(context_unichem["write_dir"])


if not os.path.exists(context_unichem["write_dir"]):
    prep_dirs()


def test_default_sources():
    api = UniChemAPI()
    assert len(api._SOURCES) == len(UniChemAPI._SOURCES)

    assert api._SOURCES.get("chembl") == 1
    assert api._SOURCES.get("drugbank") == 2
    assert api._SOURCES.get("pubchem") == 22


def test_sources_update():
    api = UniChemAPI()

    original_id = api._SOURCES["chembl"]
    api._SOURCES["chembl"] = -999

    assert api._SOURCES["chembl"] == -999

    api._update_sources()

    assert api._SOURCES["chembl"] == original_id
    assert api._SOURCES["chembl"] != -999
