def keep_only_chains_with_ligs(filtered_clusters, pdbid_intid_to_ligs):
    new_filtered_clusters = []
    for i, cluster in enumerate(filtered_clusters):
        new_cluster = []
        for pdbid_intid in cluster:
            if pdbid_intid in pdbid_intid_to_ligs:
                new_cluster.append(pdbid_intid)
        new_filtered_clusters.append(new_cluster)

    # Remove empty clusters
    new_filtered_clusters = [c for c in new_filtered_clusters if len(c) > 0]

    return new_filtered_clusters