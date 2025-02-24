r"""Implementation of the parquet class ."""

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from loguru import logger

from .files import FormatManager


class ParquetParser(FormatManager):
    r"""Class for parsing Parquet files.

    Args:
        file_path (str): Path to the Parquet file.
    """

    def __init__(self, file_path: str, schema: pa.schema) -> None:
        super().__init__(file_path)
        self.schema = schema

    def create(self) -> None:
        r"""Create a Parquet file using PyArrow."""
        # Create an empty parquet nested structure
        # Create an empty table
        table = pa.Table.from_pandas(pd.DataFrame())
        # Write the table to a Parquet file
        pq.write_table(table, self.file_path)
        logger.info(f"Parquet file created at {self.file_path}")

    def read(self, column: str | list | None = None) -> ds.Dataset:
        r"""Read a Parquet file using PyArrow.

        Args:

            columns: str | list | None
                Columns to read from the Parquet file.

        Returns:
            dataset: pyarrow.Dataset
        """
        # Read the Parquet file
        dataset = ds.dataset(self.file_path, format="parquet", schema=self.schema)
        if column is not None:
            dataset = dataset.scanner(columns=column)
        logger.info(f"Data read from Parquet file at {self.file_path}")
        return dataset

    def write(self, data: pd.DataFrame | list, data_scheme: pa.Schema) -> None:
        r"""Write data to a Parquet file using PyArrow.

        Args:
            data : pd.DataFrame
                Data to write to the Parquet file.
            schema : pa.Schema
                Schema of the data to write.
        """
        if isinstance(data, list):
            logger.debug("Converting list of dictionaries to pandas DataFrame")
            data = pd.DataFrame(data, columns=data_scheme.names)
        elif not isinstance(data, pd.DataFrame):
            raise ValueError(
                "Data must be a pandas DataFrame or a list of dictionaries"
            )
        # Remove duplicate rows
        # data = data.drop_duplicates(inplace=True,keep='first',ignore_index=True)
        # Create a table from the pandas DataFrame
        table = pa.Table.from_pandas(data, schema=data_scheme)
        # Write the table to a Parquet file
        pq.write_table(table, self.file_path)
        logger.info(f"Data written to Parquet file at {self.file_path}")

    def convert_to_table(self, column_names: str | list | None = None) -> pa.Table:
        r"""Convert data to a PyArrow Table.

        Args:
            column_names: str | list | None
                Columns to convert to a PyArrow Table.

        Returns:
            table: pa.Table
        """
        if column_names is None:
            data = self.read()
            data = data.scanner()
        else:
            data = self.read(column=column_names)
        return data.to_table()

    def convert_to_pandas(self, column_names: str | list | None = None) -> pd.DataFrame:
        r"""Convert data to a pandas DataFrame.

        Args:
            column_names: str | list | None
                Columns to convert to a pandas DataFrame.

        Returns:
            df: pd.DataFrame
        """
        if column_names is None:
            data = self.read()
            data = data.scanner().to_table()
        else:
            data = self.read(column=column_names).to_table()
        return data.to_pandas()

    def filter_data(self, condition: callable, column: str) -> pd.DataFrame:
        r"""Filter data based on a column value.

        Args:
            condition: callable
                Condition to filter the data.
            column: str
                Column to filter the data.

        Returns:
            df: pd.DataFrame
        """
        data = self.read()
        field = ds.field(column)
        filter_expr = condition(field)
        logger.debug(filter_expr)
        data = data.scanner(filter=filter_expr).to_table()
        return data.to_pandas()
