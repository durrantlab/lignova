import json

from .utils import make_url, request_with_cache


def get_pdbids_with_ligands_and_not_mutations():
    """
    Get all PDB ids that have ligands and no mutations.

    Returns:
        set: A set of pdbids.
    """

    total_count = 9999999
    current_idx = 0

    wt_pdbs_with_ligs = set([])

    while current_idx < total_count:
        search_url = "https://search.rcsb.org/rcsbsearch/v2/query?json="
        search_params = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "entity_poly.rcsb_mutation_count",
                            "operator": "equals",
                            "negation": False,
                            "value": 0,
                        },
                    },
                    {
                        "type": "terminal",
                        "service": "text_chem",
                        "parameters": {
                            "attribute": "chem_comp.type",
                            "operator": "exact_match",
                            "negation": False,
                            "value": "non-polymer",
                        },
                    },
                ],
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": current_idx, "rows": 10000},
                "results_content_type": ["experimental"],
                "sort": [{"sort_by": "score", "direction": "desc"}],
                "scoring_strategy": "combined",
            },
        }

        print(current_idx, "/", total_count)

        current_idx += 10000

        # Convert that into json and url encode it.
        url = make_url(search_url, search_params)

        # Make the request
        results = json.loads(request_with_cache(url, True))

        if total_count == 9999999:
            total_count = results["total_count"]

        wt_pdbs_with_ligs.update([x["identifier"] for x in results["result_set"]])

    return wt_pdbs_with_ligs
