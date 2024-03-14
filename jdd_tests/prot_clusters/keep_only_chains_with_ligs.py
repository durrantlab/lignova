def keep_only_chains_with_ligs(filtered_clusters, entity_id_to_ligs):
    """
    Keep only chains with ligands.

    Args:
        filtered_clusters (list):   A list of lists, where each inner list is a
                                    cluster of pdbids.
        pdbid_intid_to_ligs (dict): A dictionary with entity ids as keys and
                                    ligands as values.

    Returns:
        list: A list of lists, where each inner list is a cluster of pdbids.
    """
    new_filtered_clusters = []
    for i, cluster in enumerate(filtered_clusters):
        new_cluster = []
        for entity_id in cluster:
            if entity_id in entity_id_to_ligs:
                new_cluster.append(entity_id)
        new_filtered_clusters.append(new_cluster)

    # Remove empty clusters
    new_filtered_clusters = [c for c in new_filtered_clusters if len(c) > 0]

    return new_filtered_clusters
