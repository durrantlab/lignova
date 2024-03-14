import re

from .utils import request_with_cache


def download_clustered_data():
    """
    Download the latest clustered data from the rcsb:
    https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-90.txt

    Returns:
        list: A list of lists, where each inner list is a cluster of pdbids.
    """
    content = request_with_cache(
        "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-90.txt",
        True,
    )

    content = re.sub(r"AF_AF[A-Z0-9].+?_[0-9]", "", content)
    content = re.sub(r"MA_MA[A-Z0-9].+?_[0-9]", "", content)

    while " \n" in content:
        content = content.replace(" \n", "\n")

    content = content.strip()

    while "\n\n" in content:
        content = content.replace("\n\n", "\n")

    while "  " in content:
        content = content.replace("  ", " ")

    data = [l.split(" ") for l in content.split("\n")]

    return data
