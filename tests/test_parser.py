import os

import h5py
import numpy as np
from loguru import logger

from lignova.hdf5.parser import HDF5Parser

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_hdf5_parser = {
    "write_dir": "./tmp/pubchem",
    "test_file": "test.hdf5",
}
file_path = os.path.join(
    context_hdf5_parser["write_dir"], context_hdf5_parser["test_file"]
)


def prep_dirs():
    os.makedirs(context_hdf5_parser["write_dir"])


if not os.path.exists(context_hdf5_parser["write_dir"]):
    prep_dirs()


# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))


def test_create():
    parser = HDF5Parser(file_path)
    parser.create()
    assert os.path.exists(file_path)


def test_write():
    parser = HDF5Parser(file_path)
    data = np.array([1, 2, 3])
    parser.write("dataset", data)
    assert os.path.exists(file_path)
    assert np.array_equal(parser.read("dataset"), data)


def test_read():
    parser = HDF5Parser(file_path)
    data = np.array([1, 2, 3])
    result = parser.read("dataset")
    assert np.array_equal(result, data)


def test_read_attributes():
    parser = HDF5Parser(file_path)
    attributes = {"attr1": "value1", "attr2": "value2"}
    parser.write_attributes("dataset_1", attributes)
    result = parser.read_attributes("dataset_1")
    assert result == attributes


def test_write_attributes():
    parser = HDF5Parser(file_path)
    attributes = {"attr7": "value1", "attr8": "value2"}
    parser.write_attributes("dataset_2", attributes)
    assert os.path.exists(file_path)
    assert parser.read_attributes("dataset_2") == attributes


def test_find_file_stats():
    parser = HDF5Parser(file_path)
    parser.find_file_stats()
    statfile = file_path.replace(".hdf5", "_stats.txt")
    assert os.path.exists(statfile)
    with open(statfile, "r") as f:
        stats = f.read()
    assert "File Path: " in stats
