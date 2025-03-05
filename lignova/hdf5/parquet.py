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

    def __init__(self, *args, **kwargs):
        r"""Initialize the ParquetParser class."""
        super().__init__(*args, **kwargs)

    def create(self) -> None:
        r"""Create a Parquet file using PyArrow."""
        # Create an empty table
        table = pa.Table.from_pandas(pd.DataFrame())
        # Write the table to a Parquet file
        pq.write_table(table, self.file_path)
        logger.info(f"Parquet file created at {self.file_path}")

    def read(self, columns: str | list | None = None) -> ds.Dataset:
        r"""Read a Parquet file using PyArrow.

        Args:
            columns: Columns to read from the Parquet file.

        Returns:
            dataset: pyarrow.Dataset
        """
        # ensure that the schema is not None
        if self.schema is None:
            raise ValueError("Schema must be provided to read the Parquet file")
        # Read the Parquet file
        dataset = ds.dataset(self.file_path, format="parquet", schema=self.schema)
        if columns is not None:
            dataset = dataset.scanner(columns=columns)
        logger.info(f"Data read from Parquet file at {self.file_path}")
        return dataset

    def write(self, data: pd.DataFrame | list) -> None:
        r"""Write data to a Parquet file using PyArrow.

        Args:
            data : pandas dataframe with the data to write to parquet file.
        """
        if self.schema is None:
            raise ValueError("Schema must be provided to write the Parquet file")
        if isinstance(data, list):
            logger.debug("Converting list of dictionaries to pandas DataFrame")
            data = pd.DataFrame(data, columns=self.schema.names)
        elif not isinstance(data, pd.DataFrame):
            raise ValueError(
                "Data must be a pandas DataFrame or a list of dictionaries"
            )
        # Remove duplicate rows
        # data = data.drop_duplicates(inplace=True,keep='first',ignore_index=True)
        # Create a table from the pandas DataFrame
        table = pa.Table.from_pandas(data, schema=self.schema)
        # Write the table to a Parquet file
        pq.write_table(table, self.file_path)
        logger.info(f"Data written to Parquet file at {self.file_path}")

    def convert_to_table(self, column_names: str | list | None = None) -> pa.Table:
        r"""Convert data to a PyArrow Table.

        Args:
            column_names: Columns to convert to a PyArrow Table.

        Returns:
            A PyArrow Table containing the data from the specified columns in the Parquet file
        """
        if column_names is None:
            data = self.read()
            data = data.scanner()
        else:
            data = self.read(columns=column_names)
        return data.to_table()

    def convert_to_pandas(self, column_names: str | list | None = None) -> pd.DataFrame:
        r"""Convert data to a pandas DataFrame.
        Args:
            column_names: Columns to convert to a pandas DataFrame.
        Returns:
            A pandas DataFrame containing the data from the specified columns in the Parquet file
        """
        if column_names is None:
            data = self.read()
            data = data.scanner().to_table()
        else:
            data = self.read(columns=column_names).to_table()
        return data.to_pandas()

    def filter_data(self, condition: callable, column: str) -> pd.DataFrame:
        r"""Filter data based on a column value.
        Args:
            condition: the Condition to filter the data.
            column: Column to filter the data.
        Returns:
            A pandas DataFrame containing the filtered data.
        """
        data = self.read()
        field = ds.field(column)
        filter_expr = condition(field)
        logger.debug(filter_expr)
        data = data.scanner(filter=filter_expr).to_table()
        return data.to_pandas()
