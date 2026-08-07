"""Bulk download + Parquet ingest for PubChem BioAssay activity data."""

import gzip
import hashlib
import os

import httpx
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.dataset as ds
from loguru import logger
from wqm.api import BaseAPI
from wqm.api.errors import PermanentAPIError, TransientAPIError

from lignova.hdf5 import ParquetParser

from .model import AssayInfo, _AssayRecord

PUBCHEM_BIOASSAY_FTP = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay"
"""HTTPS mirror of the PubChem BioAssay FTP tree."""

BIOACTIVITIES_PATH = "Extras/bioactivities.tsv.gz"
"""Path (under the bioassay root) to the whole-database bioactivity dump."""

_MIB = 1024 * 1024  # bytes in one mebibyte

# Maps bioactivities columns to the AssayRecord fields.
_TSV_TO_RECORD = {
    "CID": "CID",
    "Activity Outcome": "Activity Outcome",
    "Activity Name": "Activity Name",
    "Activity Value": "Activity Value [uM]",  #
    "Activity Qualifier": "Activity Qualifier",
    "Activity Unit": "Activity Unit",
    "Protein Accession": "Target Accession",
    "Gene ID": "Target GeneID",
    "PMID": "PMID",
}

# Columns read from the Bioactivities dump
_TSV_USE_COLUMNS = [
    "AID",
    "CID",
    "Activity Outcome",
    "Activity Name",
    "Activity Value",
    "Activity Qualifier",
    "Activity Unit",
    "Protein Accession",
    "Gene ID",
    "PMID",
]


class PubChemBulk(BaseAPI):
    """Client for the PubChem BioAssay bulk activity dump.
    Download once the assay info for all AIDs into a single Parquet file, then load individual AIDs on demand.
    """

    def __init__(self, base_url: str = PUBCHEM_BIOASSAY_FTP, **kwargs) -> None:
        """Initialize the bulk client.

        Args:
            base_url: the HTTPS for the PubChem BioAssay FTP root.
        """
        super().__init__(base_url, **kwargs)

    async def remote_signature(self, path: str = BIOACTIVITIES_PATH) -> str:
        """Return a change signature for the remote file via a HEAD request.
        so that we can skip downloading if the file is unchanged.

        Args:
            path: File path under the bioassay root.

        Returns:
            A signature string that changes if the remote file changes
        """
        try:
            resp = await self._request("HEAD", path)
        except (PermanentAPIError, TransientAPIError):
            logger.warning("HEAD failed for {path}; no remote signature.", path=path)
            return ""
        etag = resp.headers.get("ETag", "")
        if etag:
            return etag.strip('"')
        last_mod = resp.headers.get("Last-Modified", "")
        length = resp.headers.get("Content-Length", "")
        return f"{last_mod}|{length}" if (last_mod or length) else ""

    @staticmethod
    def _etag(file_path: str) -> str:
        """Path of the signature etag written alongside the downloaded file.

        Args:
            file_path: Path to the downloaded file.

        Returns:
                The path to the etag file."""
        return str(file_path) + ".etag"

    def _is_current(self, file_path: str, signature: str) -> bool:
        """True if etag file exists and its stored signature matches the remote one

        Args:
            file_path: Path to the downloaded file.
            signature: Remote signature string to compare against.

        Returns:
                True if the file is current, False otherwise."""
        etag = self._etag(file_path)
        if not os.path.exists(file_path) or not os.path.exists(etag) or not signature:
            return False
        with open(etag, "r", encoding="utf-8") as file:
            return file.read().strip() == signature

    async def download(
        self,
        file_path: str,
        path: str = BIOACTIVITIES_PATH,
        force: bool = False,
        verify_md5: bool = False,
    ) -> str:
        """Stream the bulk file to disk, skipping the download if already current.

        Args:
            file_path: Local path for the .tsv.gz file.
            path: Remote path under the bioassay root.
            force: Re-download even if the local file appears current.
            verify_md5: If True, fetch checksum.md5 and verify the file's MD5
                after download (slower, but catches truncation).

        Returns:
            The path to the downloaded file.

        Raises:
            TransientAPIError: If the download stream fails.
        """
        file_path = str(file_path)
        signature = await self.remote_signature(path)

        if not force and self._is_current(file_path, signature):
            logger.info(
                "{file_path} matches the remote signature. Skipping download.",
                file_path=file_path,
            )
            return file_path

        parent = os.path.dirname(file_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        tmp = file_path + ".part"
        url = f"{self.base_url}{path}"

        logger.info("Downloading file to {file_path}", file_path=file_path)
        bytes_written = 0
        try:
            async with self._client.stream("GET", path) as resp:
                resp.raise_for_status()
                with open(tmp, "wb") as file:
                    async for chunk in resp.aiter_bytes(_MIB):
                        file.write(chunk)
                        bytes_written += len(chunk)
        except httpx.HTTPError as exc:
            if os.path.exists(tmp):
                os.remove(tmp)
            logger.error("Download failed because of {exc}", exc=exc)
            raise TransientAPIError(url=url) from exc

        logger.info("Downloaded {n:.2f} GB.", n=bytes_written / (1024**3))
        if verify_md5:
            try:
                await self._verify_md5(tmp, path)
            except Exception:
                if os.path.exists(tmp):
                    os.remove(tmp)
                raise

        os.replace(tmp, file_path)
        if signature:
            with open(self._etag(file_path), "w", encoding="utf-8") as file:
                file.write(signature)
        return file_path

    async def _verify_md5(self, file: str, path: str) -> None:
        """Verify a downloaded file against the directory's checksum.md5.

        Args:
            file: Local file to hash.
            path: Remote path of the file which is used to find its checksum line.

        Raises:
            TransientAPIError or PermanentAPIError: If the checksum can't be fetched.
            ValueError: If no checksum line matches the file, or the MD5 differs.
        """
        checksum_path = f"{os.path.dirname(path)}/checksum.md5"
        resp = await self._get(checksum_path)

        filename = os.path.basename(path)
        expected = None
        for line in resp.text.splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1].lstrip("*").endswith(filename):
                expected = parts[0]
                break
        if expected is None:
            raise ValueError(
                f"No checksum line for {filename} in {checksum_path}; cannot verify."
            )

        h = hashlib.md5()
        with open(file, "rb") as fh:
            for block in iter(lambda: fh.read(_MIB), b""):
                h.update(block)
        actual = h.hexdigest()
        if actual != expected:
            raise ValueError(
                f"MD5 mismatch for {filename}: expected {expected}, got {actual}"
            )
        logger.info("MD5 verified for {f}.", f=filename)

    def to_parquet(
        self,
        tsv_gz_path: str,
        parquet_path: str,
    ) -> str:
        """Convert the downloaded bioactivities tsv file to a single Parquet file.

        Args:
            tsv_gz_path: Path to the downloaded .tsv.gz.
            parquet_path: Path to write the Parquet file to. Will be overwritten if it exists.

        Returns:
            The parquet_path written.
        """
        parquet_path = str(parquet_path)
        read_opts = pcsv.ReadOptions(block_size=_MIB * 64)
        parse_opts = pcsv.ParseOptions(delimiter="\t")
        convert_opts = pcsv.ConvertOptions(
            include_columns=_TSV_USE_COLUMNS,
            column_types={c: pa.string() for c in _TSV_USE_COLUMNS},
            null_values=["", "NULL", "null"],
            strings_can_be_null=True,
        )

        parser = ParquetParser(parquet_path)
        writer = None
        try:
            with gzip.open(tsv_gz_path, "rb") as file:
                reader = pcsv.open_csv(
                    file,
                    read_options=read_opts,
                    parse_options=parse_opts,
                    convert_options=convert_opts,
                )
                for batch in reader:
                    table = pa.Table.from_batches([batch])
                    if writer is None:
                        writer = parser.open_writer(schema=table.schema)
                    writer.write_table(table)
        finally:
            if writer is not None:
                writer.close()
        logger.info("Wrote {src} to {dst}", src=tsv_gz_path, dst=parquet_path)
        return parquet_path

    def load_assay(self, aid: int, parquet_path: str) -> AssayInfo:
        """Build an AssayInfo for one AID from the Parquet file.
        Note that the AssayInfo returned does not have the CIDs properties like SMILES etc.

        Args:
            aid: PubChem Assay ID.
            parquet_path: Path to the Parquet created from the bioactivities tsv.gz file.

        Returns:
            An AssayInfo. records is empty if the AID is absent; per-row PMIDs
            ride along on each record, and AssayInfo.pubmed_id is derived.
        """
        parquet_path = str(parquet_path)
        if not os.path.exists(parquet_path):
            raise FileNotFoundError(
                f"{parquet_path} does not exist. Run download() + to_parquet() first."
            )
        parser = ParquetParser(parquet_path)
        scanner = parser.read(filters=ds.field("AID") == str(aid), lazy=True)
        table = scanner.to_table()
        if table.num_rows == 0:
            return AssayInfo(aid=aid)

        records: list[_AssayRecord] = [
            _AssayRecord.model_validate(
                {alias: row.get(col) for col, alias in _TSV_TO_RECORD.items()}
            )
            for row in table.to_pylist()
        ]
        return AssayInfo(aid=aid, records=records)
