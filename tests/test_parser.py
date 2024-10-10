import os

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from loguru import logger

from lignova.hdf5.parquet import ParquetParser
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

parquet_file_path = os.path.join(context_hdf5_parser["write_dir"], "test.parquet")


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
    parser.write(data, "dataset")
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


def test_parquet_create():
    parser = ParquetParser(parquet_file_path, pa.schema([]))
    parser.create()
    assert os.path.exists(parquet_file_path)


def test_parquet_write():
    # Define the schema for the nested structure
    schema = pa.schema(
        [
            ("id", pa.int64()),
            ("name", pa.string()),
            (
                "attributes",
                pa.struct(
                    [
                        ("age", pa.int64()),
                        (
                            "address",
                            pa.struct(
                                [
                                    ("street", pa.string()),
                                    ("city", pa.string()),
                                    ("zip", pa.int64()),
                                ]
                            ),
                        ),
                    ]
                ),
            ),
        ]
    )
    logger.info(f"Schema: {schema.names}")
    parser = ParquetParser(parquet_file_path, schema)
    data = [
        {
            "id": 2,
            "name": "Jane Smith",
            "attributes": {
                "age": 25,
                "address": {
                    "street": "456 Maple Ave",
                    "city": "Othertown",
                    "zip": 67890,
                },
            },
        },
    ]
    parser.write(data, parser.schema)
    assert os.path.exists(parquet_file_path)
    result = parser.read()
    assert result.schema.names == ["id", "name", "attributes"]
    assert parser.convert_to_pandas().equals(pd.DataFrame(data))


def test_parquet_read():
    schema = pa.schema(
        [
            ("id", pa.int64()),
            ("name", pa.string()),
            (
                "attributes",
                pa.struct(
                    [
                        ("age", pa.int64()),
                        (
                            "address",
                            pa.struct(
                                [
                                    ("street", pa.string()),
                                    ("city", pa.string()),
                                    ("zip", pa.int64()),
                                ]
                            ),
                        ),
                    ]
                ),
            ),
        ]
    )
    parser = ParquetParser(parquet_file_path, schema)
    data = [
        {
            "id": 2,
            "name": "Jane Smith",
            "attributes": {
                "age": 25,
                "address": {
                    "street": "456 Maple Ave",
                    "city": "Othertown",
                    "zip": 67890,
                },
            },
        },
    ]
    result = parser.read()
    assert result.schema.names == ["id", "name", "attributes"]
    assert parser.convert_to_pandas().equals(pd.DataFrame(data))


def test_filter_data():
    schema = pa.schema(
        [
            ("id", pa.int64()),
            ("name", pa.string()),
            (
                "attributes",
                pa.struct(
                    [
                        ("age", pa.int64()),
                        (
                            "address",
                            pa.struct(
                                [
                                    ("street", pa.string()),
                                    ("city", pa.string()),
                                    ("zip", pa.int64()),
                                ]
                            ),
                        ),
                    ]
                ),
            ),
        ]
    )
    parser = ParquetParser(parquet_file_path, schema)
    data = [
        {
            "id": 2,
            "name": "Jane Smith",
            "attributes": {
                "age": 25,
                "address": {
                    "street": "456 Maple Ave",
                    "city": "Othertown",
                    "zip": 67890,
                },
            },
        },
        {
            "id": 3,
            "name": "John Doe",
            "attributes": {
                "age": 30,
                "address": {
                    "street": "123 Elm St",
                    "city": "Anytown",
                    "zip": 12345,
                },
            },
        },
    ]
    parser.write(data, parser.schema)
    result = parser.filter_data(lambda x: (x == "John Doe"), "name")
    assert result.equals(pd.DataFrame([data[1]]))
