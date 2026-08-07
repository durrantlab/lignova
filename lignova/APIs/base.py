r"""Implementation for a base class for file parsers."""

import random
import time
from abc import ABC
from typing import Any

import requests
from loguru import logger


class BaseAPI(ABC):
    r"""Abstract class for REST/FTP API calls."""

    _BASE_URL: str

    _MAX_RETRIES: int = 3
    _ATTEMPTS_WAIT: float = 5.0
    _RETRYABLE_STATUS_CODES: set[int] = {503, 429}
    _DEFAULT_TIMEOUT: int = 30

    def __init__(
        self,
        base_url: str | None = None,
        response_format: str = "json",
        timeout: int | None = None,
    ):
        """Initialise the API client.

        Args:
            base_url: Override for the base URL. Falls back to _BASE_URL
            response_format: Expected response format json or xml.
            timeout: Request timeout in seconds.
        """
        self.base_url = base_url or self._BASE_URL
        self.response_format = response_format
        self.timeout = timeout or self._DEFAULT_TIMEOUT
        self._session = requests.Session()
        self._session.headers.update(self._default_headers())

    def _default_headers(self) -> dict[str, str]:
        """Return default headers for the session."""
        return {
            "accept": "application/json"
            if self.response_format == "json"
            else "application/xml"
        }

    def _request(
        self,
        method: str,
        url: str,
        params: dict[str, Any] | None = None,
        json: dict[str, Any] | None = None,
        stream: bool = False,
    ) -> requests.Response | None:
        """Make an HTTP request with automatic retries on transient errors.

        Args:
            method: HTTP verb (GET, POST)
            url: Full request URL.
            params: Query-string parameters.
            json: JSON body (for POST/PUT).
            stream: Whether to stream the response.

        Returns:
            requests.Response on success, None if all retries fail
        """
        for attempt in range(1, self._MAX_RETRIES + 1):
            try:
                resp = self._session.request(
                    method,
                    url,
                    params=params,
                    json=json,
                    stream=stream,
                    timeout=self.timeout,
                )

                if resp.status_code in self._RETRYABLE_STATUS_CODES:
                    wait = self._ATTEMPTS_WAIT * attempt + random.uniform(0, 5)
                    logger.warning(
                        f"HTTP {resp.status_code} from {url} (attempt {attempt}/{self._MAX_RETRIES}). Retrying in {wait} seconds."
                    )
                    time.sleep(wait)
                    continue

                if not resp.ok:
                    logger.error(
                        f"HTTP {resp.status_code} from {url} (attempt {attempt}/{self._MAX_RETRIES}). No more retries."
                    )
                    return None

                return resp

            except requests.RequestException as exc:
                logger.error(
                    f"Request error for {url} (attempt {attempt}/{self._MAX_RETRIES}): {exc}"
                )
                if attempt < self._MAX_RETRIES:
                    time.sleep(self._ATTEMPTS_WAIT * attempt)

        logger.error(f"All {self._MAX_RETRIES} attempts failed for {url}.")
        return None

    def _get(
        self, url: str, params: dict[str, Any] | None = None, stream: bool = False
    ) -> requests.Response | None:
        """Shorthand for a GET request."""
        return self._request("GET", url, params=params, stream=stream)

    def _post(
        self,
        url: str,
        json: dict[str, Any] | None = None,
        params: dict[str, Any] | None = None,
    ) -> requests.Response | None:
        """Shorthand for a POST request."""
        return self._request("POST", url, json=json, params=params)

    def _get_json(
        self, url: str, params: dict[str, Any] | None = None
    ) -> dict[str, Any] | None:
        """GET and parse JSON, returning None on failure."""
        resp = self._get(url, params=params)
        if resp is None:
            return None
        try:
            data = resp.json()
            return data
        except ValueError:
            logger.error(f"Failed to parse JSON from {url}")
            return None

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(base_url={self.base_url!r})"
