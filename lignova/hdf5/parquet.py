r"""Implementation of the parquet class ."""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from loguru import logger

from .files import FormatManager


class ParquetParser(FormatManager):
    r"""Class for parsing Parquet files.

    Parameters:
    ----------
        file_path (str): Path to the Parquet file.
    """

    def create(self) -> None:
        r"""Create a Parquet file using PyArrow."""
        # Create an empty parquet nested structure
        schema = pa.schema([])
        # Create an empty table
        table = pa.Table.from_pandas(pd.DataFrame(), schema=schema)
        # Write the table to a Parquet file
        pq.write_table(table, self.file_path)
        logger.info(f"Parquet file created at {self.file_path}")

    def read(self, schema: pa.schema, columns: str | list | None = None) -> pa.Table:
        r"""Read a Parquet file using PyArrow.
        parameters:
        ----------
            schema: pa.schema
                Schema of the data to read.
            columns: str | list | None
                Columns to read from the Parquet file.
        Returns:
        ----------
            table: pyarrow.Table
        """
        # Read the Parquet file
        table = pq.read_table(self.file_path, schema=schema, columns=columns)
        logger.info(f"Data read from Parquet file at {self.file_path}")
        return table

    def write(self, data: pd.DataFrame, data_scheme: pa.Schema) -> None:
        r"""Write data to a Parquet file using PyArrow

        Parameters:
        ----------
            data : pd.DataFrame
                Data to write to the Parquet file.
            schema : pa.Schema
                Schema of the data to write.
        """
        # Create a table from the pandas DataFrame
        table = pa.Table.from_pandas(data, schema=data_scheme)
        # Write the table to a Parquet file
        pq.write_table(table, self.file_path)
        logger.info(f"Data written to Parquet file at {self.file_path}")
