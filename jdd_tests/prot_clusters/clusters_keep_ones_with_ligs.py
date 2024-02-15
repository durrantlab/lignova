from tqdm import tqdm

def keep_only_cluster_members_with_ligs(clusters, wt_pdbs_with_ligs):
    """
    Keep only cluster members that are wt proteins and have ligands.

    Args:
        clusters (list): A list of lists, where each inner list is a cluster of pdbids.
        wt_pdbs_with_ligs (set): A set of pdbids, the wt ones with ligands.
    """

    new_data = []
    uniq_pdbids = set([])

    for cluster in tqdm(clusters):
        new_cluster = []
        for pdbid_and_chain in cluster:
            pdbid = pdbid_and_chain.split("_")[0]
            if pdbid in wt_pdbs_with_ligs:
                new_cluster.append(pdbid_and_chain)
                uniq_pdbids.add(pdbid)
        if len(new_cluster) > 0:
            new_data.append(new_cluster)

    # print(len(uniq_pdbids))
    uniq_pdbids = list(uniq_pdbids)
    uniq_pdbids.sort()
    return new_data, uniq_pdbids
