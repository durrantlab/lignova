r"""Implementation of the parquet class."""

import os

import pandas as pd
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from loguru import logger

from .files import FormatManager


class ParquetParser(FormatManager):
    r"""Class for parsing Parquet files.

    Args:
        file_path (str): Path to the Parquet file.
    """

    VALID_COMPRESSION = ("snappy", "zstd", "gzip", "lz4", "brotli", "none")

    def __init__(self, *args, compression: str = "snappy", **kwargs) -> None:
        r"""Initialize the ParquetParser class."""
        super().__init__(*args, **kwargs)
        if compression not in self.VALID_COMPRESSION:
            raise ValueError(
                f"Invalid compression option: {compression}. "
                f"Valid options are: {self.VALID_COMPRESSION}"
            )
        self.compression = compression
        if self.schema is None:
            self.schema = self.find_schema()

    def find_schema(self) -> pa.Schema | None:
        r"""Find the schema of a Parquet file using PyArrow.

        Returns:
            schema: pyarrow.Schema
        """
        # check if file exists
        if not os.path.exists(self.file_path):
            logger.warning(f"File {self.file_path} does not exist.")
            return None
        return pq.read_schema(self.file_path)

    def create(self) -> None:
        r"""Create a Parquet file using PyArrow."""
        # Create an empty table
        table = pa.Table.from_pandas(pd.DataFrame())
        # Write the table to a Parquet file
        pq.write_table(table, self.file_path)
        logger.info(f"Parquet file created at {self.file_path}")

    def read(
        self,
        columns: str | list | None = None,
        filters: list | None | ds.Expression = None,
        lazy: bool = False,
    ) -> pa.Table | ds.Scanner:
        r"""Read a Parquet file using PyArrow and offers two apporaches
            - memory intensive approach when lazy is false, where the entire file is loaded
                into memory as a PyArrow Table.
            - When lazy is true, it returns a PyArrow Dataset object that
                allows for more efficient querying and filtering
        Args:
            columns: Columns to read from the Parquet file.
            filters: Filters to apply when reading the Parquet file. Accepts either
                a list of tuples (e.g., [('column_name', '==', value)])
                or a PyArrow expression. e.g., ds.field('column_name') == value.
            lazy: if true retuns a Scanner object instead of a table. Default is False.
        Returns:
            pa.Table (lazy=False) or ds.Scanner (lazy=True)
        """
        # ensure that the schema is not None
        if self.schema is None:
            raise ValueError("Schema must be provided to read the Parquet file")
        if lazy:
            dataset = ds.dataset(
                self.file_path,
                format="parquet",
                schema=self.schema,
            )
            scanner_kwargs: dict = {}
            if columns is not None:
                scanner_kwargs["columns"] = columns
            if filters is not None:
                scanner_kwargs["filter"] = filters
            return dataset.scanner(**scanner_kwargs)

        table = pq.read_table(
            self.file_path,
            columns=columns,
            filters=filters or None,
        )
        logger.info(f"Data read from Parquet file at {self.file_path}")
        return table

    def write(
        self,
        data: pa.Table | pd.DataFrame | list,
        group_size: int | None = None,
        compression: str | None = None,
    ) -> None:
        r"""Write data to a Parquet file.
        Args:
            data: Data to write.it accepts a PyArrow Table, a pandas DataFrame, or a list of dictionaries
            group_size: number of rows per group in the Parquet file (default is None, which means no grouping).
            compression: Compression codec (default ``"snappy"``).
        """
        compression = compression or self.compression

        if isinstance(data, pa.Table):
            table = data
            self.schema = table.schema
        else:
            if self.schema is None:
                raise ValueError(
                    "Schema must be provided to write DataFrame / list data"
                )
            if isinstance(data, list):
                logger.debug("Converting list of dictionaries to pandas DataFrame")
                data = pd.DataFrame(data, columns=self.schema.names)
            elif not isinstance(data, pd.DataFrame):
                raise ValueError(
                    "Data must be a pa.Table, pd.DataFrame, or list of dicts"
                )
            table = pa.Table.from_pandas(data, schema=self.schema)

        pq.write_table(
            table,
            self.file_path,
            compression=compression,
            row_group_size=group_size,
        )
        logger.info(f"Data written to Parquet file at {self.file_path}")

    def read_metadata(self) -> pq.FileMetaData:
        r"""Read only Parquet file metadata (schema + row counts).
            No row data is loaded, making this quick stats.
        Returns:
            Parquet FileMetaData object.
        """
        return pq.read_metadata(self.file_path)

    def open_writer(
        self,
        schema: pa.Schema | None = None,
        compression: str | None = None,
    ) -> pq.ParquetWriter:
        r"""Open a ParquetWriter for incremental / streaming writes.

        The caller is responsible for closing the writer (use as a
        context manager or call .close()).

        Args:
            schema: Arrow schema for the file.  Falls back to
                ``self.schema`` if not provided.
            compression: Compression codec.  Falls back to the
                instance default set at init.

        Returns:
            An open pq.ParquetWriter.
        """
        schema = schema or self.schema
        compression = compression or self.compression
        if schema is None:
            raise ValueError(
                "A schema is required to open a ParquetWriter. "
                "Provide one explicitly or ensure the file already exists."
            )
        parent = os.path.dirname(self.file_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        return pq.ParquetWriter(
            self.file_path,
            schema=schema,
            compression=compression,
        )

    def convert_to_table(self, column_names: str | list | None = None) -> pa.Table:
        r"""Convert data to a PyArrow Table.

        Args:
            column_names: Columns to convert to a PyArrow Table.

        Returns:
            A PyArrow Table containing the data from the specified columns in the Parquet file
        """
        return self.read(columns=column_names, lazy=False)

    def convert_to_pandas(self, column_names: str | list | None = None) -> pd.DataFrame:
        r"""Convert data to a pandas DataFrame.
        Args:
            column_names: Columns to convert to a pandas DataFrame.
        Returns:
            A pandas DataFrame containing the data from the specified columns in the Parquet file
        """
        return self.read(columns=column_names).to_pandas()

    def to_csv(self, output_path: str) -> None:
        r"""Convert data to a CSV file.
        Args:
            output_path: Path to the output CSV file.
        """
        table = self.convert_to_table()
        while any(pa.types.is_struct(field.type) for field in table.schema):
            table = table.flatten()
        pcsv.write_csv(table, output_path)
        logger.info(f"Exported Parquet data to CSV file in {output_path}")

    def filter_data(
        self,
        condition: callable,
        column: str,
        lazy: bool = False,
    ) -> pd.DataFrame | ds.Scanner:
        r"""Filter data based on a column value.
        Args:
            condition: the Condition to filter the data.
            column: Column to filter the data.
            lazy: if true retuns a Scanner object instead of a table. Default is False.
        Returns:
            pandas DataFrame (Lazy=False) or ds.Scanner (Lazy=True) containing the filtered data.
        """
        filter_expr = condition(ds.field(column))
        scanner = self.read(filters=filter_expr, lazy=True)
        if lazy:
            return scanner
        return scanner.to_table().to_pandas()

    @staticmethod
    def open_dataset(
        path: str,
        schema: pa.Schema | None = None,
    ) -> ds.Dataset:
        r"""Open a directory of Parquet files as a unified lazy dataset.
        Args:
            path: Directory (or single file) containing .parquet files.
            schema: Optional schema to enforce.

        Returns:
            A PyArrow Dataset object representing the Parquet files at the specified path.
        """
        kwargs = {"format": "parquet"}
        if schema is not None:
            kwargs["schema"] = schema
        return ds.dataset(path, **kwargs)
