# First, download the latest clustered data from the rcsb:
# https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-90.txt

import json
from prot_clusters.keep_only_chains_with_ligs import keep_only_chains_with_ligs
from prot_clusters.clusters_keep_ones_with_ligs import keep_only_cluster_members_with_ligs
from prot_clusters.collect_lig_data import collect_all_lig_data
from prot_clusters.download_pdb_clusters import download_clustered_data
from prot_clusters.utils import make_url, request_with_cache
from prot_clusters.wt_prots_with_ligs import get_pdbids_with_ligands_and_not_mutations

def main():
    """
    Main function.
    """

    # Download the latest clustered data from the rcsb:
    clusters = download_clustered_data()

    # Get all PDB ids that have ligands and no mutations.
    wt_pdbs_with_ligs = get_pdbids_with_ligands_and_not_mutations()

    # Keep only cluster members that are wt proteins and have ligands.
    filtered_clusters, uniq_pdbids = keep_only_cluster_members_with_ligs(clusters, wt_pdbs_with_ligs)
    
    # Collect all ligand data for a list of pdbids.
    entity_id_to_pdbid_chain, pdbid_chain_to_ligs, pdbid_intid_to_ligs = collect_all_lig_data(uniq_pdbids)

    filtered_clusters_with_ligs = keep_only_chains_with_ligs(filtered_clusters, pdbid_intid_to_ligs)

    # Now switch from the entity_id to the pdbid_chain (auth). Multiple chains
    # map to the same entity_id. Let's just pick the first one.
    enhanced_clusters = []
    for cluster in filtered_clusters_with_ligs:
        enhanced_cluster = []
        for entity_id in cluster:
            for pdb_chain in entity_id_to_pdbid_chain[entity_id]:
                for lig in pdbid_chain_to_ligs[pdb_chain]:
                    enhanced_cluster.append((pdb_chain, lig))
                    # print(entity_id, pdb_chain, lig)
        enhanced_clusters.append(enhanced_cluster)

    # Save to json
    with open('enhanced_clusters.json', 'w') as f:
        json.dump(enhanced_clusters, f)

if __name__ == "__main__":
    main()
