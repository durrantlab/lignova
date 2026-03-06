import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
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
    """Prepare directories for writing files."""
    os.makedirs(context_hdf5_parser["write_dir"])


if not os.path.exists(context_hdf5_parser["write_dir"]):
    prep_dirs()


# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))


def _parquet_schema() -> pa.Schema:
    """Return the shared nested schema used by all Parquet tests."""
    return pa.schema(
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


def _sample_data() -> list[dict]:
    """Return a single-row sample dataset used by most Parquet tests."""
    return [
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


def test_create():
    r"""Test the creation of HDF5 files."""
    parser = HDF5Parser(file_path)
    parser.create()
    assert os.path.exists(file_path)


def test_write():
    r"""Test writing data to HDF5 files."""
    parser = HDF5Parser(file_path)
    data = np.array([1, 2, 3])
    parser.write(data, "dataset")
    assert os.path.exists(file_path)
    assert np.array_equal(parser.read("dataset"), data)


def test_read():
    r"""Test reading data from HDF5 files."""
    parser = HDF5Parser(file_path)
    data = np.array([1, 2, 3])
    result = parser.read("dataset")
    assert np.array_equal(result, data)


def test_write_attributes():
    r"""Test writing attributes to HDF5 files."""
    parser = HDF5Parser(file_path)
    attributes = {"attr7": "value1", "attr8": "value2"}
    parser.write_attributes("dataset_2", attributes)
    assert os.path.exists(file_path)
    assert parser.read_attributes("dataset_2") == attributes


def test_read_attributes():
    r"""Test reading attributes from HDF5 files."""
    parser = HDF5Parser(file_path)
    attributes = {"attr1": "value1", "attr2": "value2"}
    parser.write_attributes("dataset_1", attributes)
    result = parser.read_attributes("dataset_1")
    assert result == attributes


def test_find_file_stats():
    r"""Test finding file stats."""
    parser = HDF5Parser(file_path)
    parser.find_file_stats()
    statfile = file_path.replace(".hdf5", "_stats.txt")
    assert os.path.exists(statfile)
    with open(statfile, "r", encoding="utf-8") as f:
        stats = f.read()
    assert "File Path: " in stats
    assert stats.split("\n")[0] == f"File Path: {file_path}"


def test_parquet_create():
    r"""Test the creation of Parquet files."""
    parser = ParquetParser(parquet_file_path, pa.schema([]))
    parser.create()
    assert os.path.exists(parquet_file_path)


def test_parquet_write():
    r"""Test writing data to Parquet files."""
    schema = _parquet_schema()
    data = _sample_data()
    parser = ParquetParser(parquet_file_path, schema)
    parser.write(data)
    assert os.path.exists(parquet_file_path)
    result = parser.read()
    assert result.schema.names == ["id", "name", "attributes"]
    assert parser.convert_to_pandas().equals(pd.DataFrame(data))


def test_parquet_read_schema():
    r"""Test reading data from Parquet files."""
    schema = _parquet_schema()
    parser = ParquetParser(parquet_file_path, schema)
    result = parser.find_schema()
    assert result == schema


def test_parquet_read():
    r"""Test reading data from Parquet files."""
    schema = _parquet_schema()
    data = _sample_data()
    parser = ParquetParser(parquet_file_path, schema)
    result = parser.read()
    assert result.schema.names == ["id", "name", "attributes"]
    assert parser.convert_to_pandas().equals(pd.DataFrame(data))


def test_filter_data():
    r"""Test filtering data from Parquet files."""
    schema = _parquet_schema()
    data = _sample_data() + [
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
    parser = ParquetParser(parquet_file_path, schema)
    parser.write(data)
    result = parser.filter_data(lambda x: (x == "John Doe"), "name")
    assert result.equals(pd.DataFrame([data[1]]))


def test_parquet_read_metadata():
    r"""Test reading metadata from Parquet files."""
    schema = _parquet_schema()
    data = _sample_data()
    parser = ParquetParser(parquet_file_path, schema)
    parser.write(data)
    metadata = parser.read_metadata()
    assert metadata.num_rows == 1
    assert metadata.num_columns == 6
    assert metadata.num_row_groups == 1
    data = _sample_data() + [
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
    parser.write(data,1)
    metadata = parser.read_metadata()
    assert metadata.num_row_groups == 2
    parser.write(data, group_size=2)
    assert parser.read_metadata().num_row_groups == 1


def test_parquet_open_writer():
    r"""Test streaming writes using open_writer."""
    schema = _parquet_schema()
    writer_path = os.path.join(context_hdf5_parser["write_dir"], "test_writer.parquet")
    parser = ParquetParser(writer_path, schema)
    batch_1 = [
        {
            "id": 1,
            "name": "Alice",
            "attributes": {
                "age": 28,
                "address": {
                    "street": "10 Main St",
                    "city": "Townsville",
                    "zip": 11111,
                },
            },
        },
    ]
    batch_2 = [
        {
            "id": 2,
            "name": "Bob",
            "attributes": {
                "age": 35,
                "address": {
                    "street": "20 Oak Rd",
                    "city": "Villageton",
                    "zip": 22222,
                },
            },
        },
    ]
    with parser.open_writer() as writer:
        writer.write_table(
            pa.Table.from_pandas(pd.DataFrame(batch_1), schema=schema)
        )
        writer.write_table(
            pa.Table.from_pandas(pd.DataFrame(batch_2), schema=schema)
        )
    assert os.path.exists(writer_path)
    result = parser.convert_to_pandas()
    assert len(result) == 2
    assert list(result["name"]) == ["Alice", "Bob"]


def test_parquet_to_csv():
    r"""Test exporting Parquet data to a CSV file."""
    schema = _parquet_schema()
    data = _sample_data()
    parser = ParquetParser(parquet_file_path, schema)
    parser.write(data)
    csv_output_path = os.path.join(context_hdf5_parser["write_dir"], "test_output.csv")
    parser.to_csv(csv_output_path)
    assert os.path.exists(csv_output_path)


def test_parquet_read_lazy():
    r"""Test lazy reading from Parquet files returns a Scanner."""
    schema = _parquet_schema()
    data = _sample_data()
    parser = ParquetParser(parquet_file_path, schema)
    parser.write(data)
    scanner = parser.read(columns=["id", "name"], lazy=True)
    assert isinstance(scanner, ds.Scanner)
    table = scanner.to_table()
    assert table.num_rows == 1
    assert table.schema.names == ["id", "name"]


def test_parquet_convert_to_table():
    r"""Test converting Parquet data to a PyArrow Table."""
    schema = _parquet_schema()
    data = _sample_data()
    parser = ParquetParser(parquet_file_path, schema)
    parser.write(data)
    table = parser.convert_to_table(column_names=["id", "name"])
    assert isinstance(table, pa.Table)
    assert table.schema.names == ["id", "name"]
    assert table.num_rows == 1


def test_parquet_open_dataset():
    r"""Test opening a Parquet file as a dataset."""
    schema = _parquet_schema()
    data = _sample_data()
    parser = ParquetParser(parquet_file_path, schema)
    parser.write(data)
    dataset = ParquetParser.open_dataset(parquet_file_path, schema=schema)
    assert isinstance(dataset, ds.Dataset)
    table = dataset.to_table()
    assert table.num_rows == 1
    assert table.schema.names == ["id", "name", "attributes"]