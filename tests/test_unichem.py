# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Test unichem API."""

import json
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
    assert UniChemAPI._SOURCES.items() <= api._SOURCES.items()
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


def test_mapping_url():
    url_a = UniChemAPI._mapping_url(1, 2)
    url_b = UniChemAPI._mapping_url(2, 1)
    assert url_a == url_b
    assert "src1src2" in url_a
    url = UniChemAPI._mapping_url(1, 22)
    assert url.startswith(UniChemAPI._FTP_BASE)
    assert url.endswith("src1src22.txt.gz")
    assert "/src_id1/" in url
    url = UniChemAPI._mapping_url(1, 1)
    assert "src1src1" in url


def test_get_remote_etag():
    api = UniChemAPI(task="mapping")
    url = UniChemAPI._mapping_url(1, 2)
    etag = api._get_remote_etag(url)
    assert isinstance(etag, str)
    etag = api._get_remote_etag("https://ftp.ebi.ac.uk/this/does/not/exist.gz")
    assert etag == ""
