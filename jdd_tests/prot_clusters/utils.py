import hashlib
import json
import os
import urllib.parse
from datetime import date

import requests


def request_with_cache(url, use_date=False):
    """
    Download the contents of a url, caching the results in a file.

    Args:
        url (str): The url to download
        use_date (bool): If True, append the current date to the filename. If
            you expect data to change frequently, set this to True to pull the
            latest data.
    Returns:
        str: The contents of the url.
    """

    if not os.path.exists("cache"):
        os.mkdir("cache")

    # Hash the url
    m = hashlib.md5()
    url_hash = m.update(url.encode("utf-8"))
    url_hash = m.hexdigest()

    filename = "./cache/" + url_hash

    if use_date:
        today = date.today()
        date_string = today.strftime("%m-%d-%Y")
        filename += "." + date_string

    filename += ".txt"

    if os.path.exists(filename):
        # print("Using cached file: " + filename)
        with open(filename, "r") as f:
            return f.read()

    # print("Downloading file: " + url[:100])

    r = requests.get(url, allow_redirects=True)
    content = r.content.decode("utf-8")

    with open(filename, "w") as f:
        f.write(content)

    return content


def make_url(url, params):
    """
    Make a url with params.

    Args:
        url (str): The url to download
        params (any): The params to add to the url.
    Returns:
        str: The url with params.
    """

    # If params is not a string, json encode it.
    if type(params) != str:
        params = json.dumps(params)
    else:
        # It is a string
        params = params.replace("\n", " ")
        while "  " in params:
            params = params.replace("  ", " ")

    # Regardless, url encode it.
    params = urllib.parse.quote(params)

    return url + params
