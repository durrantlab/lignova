r"Implementation for the wrapper for the docking program Schrodinger GLIDE."
import os

import MDAnalysis as mda
import pandas as pd
import pyarrow as pa
from loguru import logger

from lignova.analysis.rmsd import RMSD
from lignova.analysis.utils import clean_and_standardize_file
from lignova.docking.combind import Combind
from lignova.docking.contexts import GlideContext
from lignova.docking.contexts.combind import CombindContext
from lignova.docking.glide import Glide
from lignova.docking.utils import manipulate_complexes
from lignova.hdf5.parquet import ParquetParser
from lignova.io import write_text
from lignova.structure.editing import convert_cif2pdb, write_mda_universe
from lignova.structure.ligand import DockedLigand, Ligand, PreparedLigand
from lignova.structure.protein import PreparedProtein, Protein
from lignova.structure.utils import (
    chery_pick_ligand,
    get_ligand_names,
    separate_protein_ligand,
    validate_ligands,
    validate_pdb,
)

# OUTLINE
# 1. Importing necessary modules
# 2. Reading the input files to get the ligand and receptor files
# 3. Running the preparation of the ligand and receptor files
# 4. Running the docking program
# 5. Running the post-processing of the docking results to get the top poses using combind score
# 6. Calculating the RMSD of the top poses with the reference ligand pose
# 7. Writing the output files with which proteins passed the docking
# and the RMSD values <= 2.5 Angstroms


def get_pdb_ids_from_parquet(file_path: str, schema: pa.Schema | None = None) -> list:
    r"""
    Get the pdb ids from the parquet file

    Args:
        file_path : The path to the parquet file

    Returns:
        The list of pdb ids.
    """
    if not os.path.exists(file_path):
        logger.error(f"The file {file_path} does not exist")
        raise FileNotFoundError(f"The file {file_path} does not exist")
    if schema is None:
        schema = pa.schema(
            [
                ("Protein Cluster number", pa.int64()),
                ("PDB/Gene ID", pa.string()),
                ("Compound ID", pa.string()),
                ("Smiles", pa.string()),
                ("Ligand Cluster number", pa.int64()),
            ]
        )
    new_parquet = ParquetParser(file_path, schema)
    raw_prot_ids = new_parquet.convert_to_pandas()["PDB/Gene ID"].unique()
    pdb_ids = [
        raw_prot_ids[i]
        for i in range(len(raw_prot_ids))
        if any(char.isalpha() for char in raw_prot_ids[i])
    ]
    return pdb_ids


def extract_parquet_clusters(
    file_path: str, pdb_id: str, same_ligand_cluster: bool = True
) -> pd.DataFrame:
    r"""
    Read the parquet file containing the clustered protein and ligand information

    Args:
        file_path : The path to the parquet file
        pdb_id : The pdb id of the protein of interest
        same_ligand_cluster : If true, only return the members in the same ligand cluster
            as the protein of interest not only the same protein cluster

    Returns:

        The dataframe containing the protein and ligand information
    """
    if not os.path.exists(file_path):
        logger.error(f"The file {file_path} does not exist")
        raise FileNotFoundError(f"The file {file_path} does not exist")
    data = pd.read_parquet(file_path)
    prot_data = data[data["PDB/Gene ID"] == pdb_id]
    # logger.info(f"The protein data is {prot_data}")
    # check if the protein is in the dataframe
    if prot_data.empty:
        logger.error(f"The protein {pdb_id} is not in the dataframe")
        raise ValueError(f"The protein {pdb_id} is not in the dataframe")
    prot_cluster_number = prot_data["Protein Cluster number"].values[0]
    lig_cluster_number = prot_data["Ligand Cluster number"].values[0]
    if same_ligand_cluster:
        members = data[
            (data["Ligand Cluster number"] == lig_cluster_number)
            & (data["Protein Cluster number"] == prot_cluster_number)
        ]
    else:
        members = data[data["Protein Cluster number"] == prot_cluster_number]
    return members


def parse_ligand_members(
    cluster_members: pd.DataFrame,
    pdb_id: str,
    find_pdb_ligand: bool = False,
    input_dir: str | None = None,
    water: bool = True,
    datatype: str = "pdb",
    type_water: str | None = None,
) -> pd.DataFrame | mda.Universe:
    r"""Parse the cluster members to separate real PDB ligands and PubChem ligands.

    Args:
        cluster_members : DataFrame containing ligand cluster information.
        pdb_id : PDB ID of the protein.
        find_pdb_ligand : Whether to find PDB ligands only or all ligands.
        input_dir : Directory containing PDB files.
        water : Whether to remove water molecules when extracting ligands.
        datatype : Output format when extracting ligands ("pdb" or "pandas").
        type_water : Water type to retain if water=False.

    Returns:
        Ligands as a pandas DataFrame or MDAnalysis Universe object.
    """
    if input_dir is None:
        logger.error("Input directory must be provided.")
        raise ValueError("input_dir must be provided to locate the PDB file.")

    pdb_file = os.path.join(input_dir, f"{pdb_id.lower()}.pdb")
    if not os.path.exists(pdb_file):
        logger.error(f"PDB file {pdb_file} does not exist.")
        raise FileNotFoundError(f"PDB file {pdb_file} not found.")

    # Load the PDB file and detect ligand residue names
    ligand_resnames = get_ligand_names(pdb_file)

    logger.info(f"Detected ligand residue names in PDB: {ligand_resnames}")

    if not ligand_resnames:
        logger.error("No ligands detected in the PDB file.")
        raise ValueError("No ligands detected in the PDB file.")

    # Split the cluster_members into pdb ligands and pubchem ligands
    all_pdb_ligands = cluster_members[
        cluster_members["Compound ID"].isin(ligand_resnames)
    ].drop_duplicates()

    pubchem_ligands = cluster_members[
        ~cluster_members["Compound ID"].isin(ligand_resnames)
    ].drop_duplicates()

    if all_pdb_ligands.empty:
        logger.error("No matching cluster members found for PDB ligands.")
        raise ValueError("No matching cluster members found for PDB ligands.")

    # Handle datatype
    if find_pdb_ligand:
        if datatype == "pdb":
            ligand_universe_list = []
            for lig_id in all_pdb_ligands["Compound ID"]:
                _, ligand = chery_pick_ligand(
                    pdb_file,
                    lig_id,
                    remove_water=water,
                    water_selection=type_water,
                )
                ligand_universe_list.append(ligand)

            ligand = (
                ligand_universe_list[0]
                if len(ligand_universe_list) == 1
                else ligand_universe_list
            )

        elif datatype == "pandas":
            ligand = (
                all_pdb_ligands[["Smiles", "Compound ID"]]
                .rename(columns={"Smiles": "SMILES", "Compound ID": "s_m_title"})
                .reset_index(drop=True)
            )

        else:
            logger.error(f"Invalid datatype: {datatype}")
            raise ValueError(f"Invalid datatype: {datatype}")

    else:
        # combine all ligands
        all_ligands = pd.concat([all_pdb_ligands, pubchem_ligands], ignore_index=True)
        all_ligands = all_ligands.drop_duplicates()
        all_ligands = (
            all_ligands[["Smiles", "Compound ID"]]
            .rename(columns={"Smiles": "SMILES", "Compound ID": "s_m_title"})
            .reset_index(drop=True)
        )
        ligand = all_ligands

    return ligand


def get_pdb_coordinates(pdb_id: str, work_dir: str) -> None:
    """
    This function takes a list of PDB IDs and downloads the PDB files to a specified directory
    if they pass the validation test. (no mutation, has ligand, no covalent bond, x-ray structures,)

    Args:
        pdb_id : The PDB ids to be downloaded
        work_dir : The working directory where the PDB file will be downloaded.
    """
    current_dir = os.getcwd()
    protein = Protein()
    # check if the output directory exists
    if not os.path.exists(work_dir):
        logger.info("Output_Dir not found,Creating it in working directory")
        os.mkdir(os.path.join(current_dir, work_dir))
    if (
        not os.path.exists(os.path.join(work_dir, pdb_id.lower() + ".pdb"))
        and validate_pdb(pdb_id)
        and validate_ligands(pdb_id)
    ):
        logger.info(f"Downloading PDB file for {pdb_id}")
        file_ext = (
            "pdb" if protein.get_pdb_from_rcsb(pdb_id).startswith("HEADER") else "cif"
        )
        protein.load(
            pdb_id=pdb_id,
            write=True,
            write_path=os.path.join(work_dir, pdb_id.lower() + "." + file_ext),
        )
        if file_ext == "cif":
            logger.info(f"Converting {pdb_id} to pdb format")
            convert_cif2pdb(
                os.path.join(work_dir, pdb_id.lower() + ".cif"),
                os.path.join(work_dir, pdb_id.lower() + ".pdb"),
            )
    elif os.path.exists(os.path.join(work_dir, pdb_id.lower() + ".pdb")):
        logger.info(f"{pdb_id} already exists in the directory")
    else:
        logger.warning(f"{pdb_id} failed validation test")


def prep_ligands(
    ligand_object: pd.DataFrame | mda.core.groups.AtomGroup,
    context: GlideContext | None,
    file_name: str,
) -> PreparedLigand:
    r"""
    Prepare the ligands for docking

    Args:
        ligand_file : The ligand data to be prepared
        context : The context object with information about the preparation
        file_name : The name of the file to be written

    Returns:
        a PreparedLigand object with the prepared ligand
    """
    glide = Glide()
    if not os.path.exists(context.write_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    if isinstance(ligand_object, pd.DataFrame):
        file_path = write_text(ligand_object, file_ext=".csv")
        file_name = file_name + "_lig_pubchem.csv"
    else:
        file_path = write_text(ligand_object, file_ext=".pdb")
        file_name = file_name + "_lig.pdb"
    temp_path = os.path.dirname(file_path)
    logger.debug(f"Temp location is {temp_path}")
    # Rename the temporary file
    new_path = os.path.join(temp_path, f"{file_name}")
    os.rename(file_path, new_path)
    logger.debug(f"Renamed temporary file to {new_path}")
    ligand = Ligand(new_path)
    try:
        if not os.path.exists(
            os.path.join(context.write_dir, ligand.file_id + "_prepared.mae")
        ):
            glide.PrepLigand(ligand, context)
        prepped_lig = PreparedLigand(
            os.path.join(context.write_dir, ligand.file_id + "_prepared.mae")
        )
    except Exception as exc:
        logger.error(f"Error in preparing ligand {ligand.file_id}")
        raise exc
    os.remove(new_path)
    if os.path.exists(os.path.join(context.write_dir, ligand.file_id + ".mae")):
        os.remove(os.path.join(context.write_dir, ligand.file_id + ".mae"))
    return prepped_lig


def prep_proteins(
    pdb_file: str,
    context: GlideContext,
    lig_asl: str | None,
    type_water: str | None = None,
) -> PreparedProtein:
    """
    Prepare the protein for docking

    Args:
        pdb_file : The path to the pdb file
        context : The context object with information about the preparation
        lig_asl : The ligand asl file to be used in the preparation i.e the ligand residue nam
        type_water : The type of water to be kept when samplewater is True.
            it can be "surface"| "interfacial"| "all" or None. Default is None
    Returns:
        a PreparedProtein object with the prepared protein
    """
    temp_prot = Protein(file_path=pdb_file)
    if not os.path.exists(context.write_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    if not os.path.exists(os.path.join(context.write_dir, temp_prot.file_id + ".pdb")):
        protein, _ligand = separate_protein_ligand(
            pdb_file,
            remove_water=not (context.samplewater),
            keep_het_chain="A",
            water_selection=type_water,
        )
        write_mda_universe(
            protein, os.path.join(context.write_dir, f"{temp_prot.file_id}.pdb")
        )
    input_prot = Protein(
        file_path=os.path.join(context.write_dir, f"{temp_prot.file_id}.pdb")
    )
    glide = Glide()
    try:
        if not os.path.exists(
            os.path.join(context.write_dir, input_prot.file_id + "_grid.zip")
        ):
            glide.PrepProtein(input_prot, context, lig_asl)
        prepped_prot = PreparedProtein(
            os.path.join(context.write_dir, input_prot.file_id + "_grid.zip")
        )
    except Exception as exc:
        logger.error(f"Error in preparing protein {input_prot.file_id}")
        raise exc
    os.remove(os.path.join(context.write_dir, f"{temp_prot.file_id}.pdb"))
    return prepped_prot


def dock_ligands(
    prepped_protein: PreparedProtein,
    prepped_ligand: PreparedLigand,
    context: GlideContext,
) -> DockedLigand:
    r"""
    Dock the ligands to the protein

    Args:
        prepped_protein : The prepared protein object to dock the ligand to
        prepped_ligand :  The prepared ligand object to be docked
        context : The glidecontext object with information about the docking
    Returns:
        a DockedLigand object with the docked ligand
    """
    glide = Glide()
    if not os.path.exists(context.write_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    logger.info(f"Docking {prepped_ligand.file_id} to {prepped_protein.file_id}")
    try:
        if not os.path.exists(
            os.path.join(
                context.write_dir,
                prepped_ligand.file_id.replace("_prepared", "_docking_pv.maegz"),
            )
        ):
            glide.run(prepped_protein, prepped_ligand, context)
        docked_ligand = DockedLigand(
            os.path.join(
                context.write_dir,
                prepped_ligand.file_id.replace("_prepared", "_docking_pv.maegz"),
            )
        )
    except Exception as exc:
        logger.error(f"Error in docking ligand {prepped_ligand.file_id}")
        raise exc
    return docked_ligand


def run_combind(
    docked_ligand: DockedLigand, context: CombindContext, screen: bool = False
) -> str:
    """
    Run the combind program to get the top poses

    Args:
        docked_ligand : The DockedLigand object containing the docking results
        context : CombindContext with information about the combind docking
        screen : If true, we use combindVS to screen the ligands and generate features
            Default is False

    Returns:
        The path to the combind docking results CSV file
    """
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    if not os.path.exists(context.work_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    logger.info(f"Generating features for {docked_ligand.file_id}")
    try:
        if not os.path.exists(
            os.path.join(context.work_dir, docked_ligand.file_id + "_features")
        ):
            combind.featurize(
                docking_filepaths=docked_ligand.file_path,
                file_name=docked_ligand.file_id,
                screen=screen,
            )
        logger.debug(
            f"Fetures generated for {docked_ligand.file_id}.Selecting top poses"
        )
        if not (
            os.path.exists(
                os.path.join(context.work_dir, docked_ligand.file_id + "_poses.csv")
                or os.path.join(context.work_dir, docked_ligand.file_id + "_screen.npy")
            )
        ):
            if screen:
                combind.compute_combind_score(
                    features_dir=os.path.join(
                        context.work_dir, docked_ligand.file_id + "_features"
                    ),
                    filename=docked_ligand.file_id,
                )
                combind.apply_combind_score(
                    docking_filepath=docked_ligand.file_path,
                    combind_score_file=docked_ligand.file_path.replace(
                        ".maegz", "_screen.npy"
                    ),
                    output_filename=docked_ligand.file_id,
                )
                combind.extract_data_csv(
                    docking_file=os.path.join(
                        context.work_dir,
                        docked_ligand.file_id + "_combind_sorted.maegz",
                    ),
                    filename=docked_ligand.file_id + "_poses",
                    top_poses=True,
                )
            else:
                combind.select_pose(
                    docked_ligand.file_id + "_poses",
                    os.path.join(context.work_dir, docked_ligand.file_id + "_features"),
                )
        if os.path.exists(
            os.path.join(context.work_dir, docked_ligand.file_id + "_combind.maegz")
        ):
            os.remove(
                os.path.join(context.work_dir, docked_ligand.file_id + "_combind.maegz")
            )
        logger.debug(f"Top poses selected csv for {docked_ligand.file_id} is written")
    except Exception as exc:
        logger.error(f"Error in scoring ligand {docked_ligand.file_id}")
        raise exc
    return os.path.join(context.work_dir, docked_ligand.file_id + "_poses.csv")


def get_top_combind_pose(
    combind_csv: str, glide_docking_file: DockedLigand, context: CombindContext
) -> DockedLigand:
    r"""Parse the combind docking results to get the top pose

    Args:
        combind_csv : The path to the combind docking results CSV file
        glide_docking_file : DockedLigand object with the glide docking results
        context : CombindContext with information about the combind docking
    Returns:
        The DockedLigand object containing the top pose
    """
    if not os.path.exists(context.work_dir):
        logger.warning(f"The directory {context.work_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    if not os.path.exists(combind_csv):
        logger.error(f"The file {combind_csv} does not exist")
        raise FileNotFoundError(f"The file {combind_csv} does not exist")
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    try:
        if not os.path.exists(
            os.path.join(
                context.work_dir, glide_docking_file.file_id + "_top_poses.maegz"
            )
        ):
            combind.get_3d_top_pose(
                glide_docking_file.file_path,
                combind_csv,
                glide_docking_file.file_id + "_top_poses",
            )
        logger.debug(f"Top poses selected for {glide_docking_file.file_id}")
        top_poses = DockedLigand(
            os.path.join(
                context.work_dir, glide_docking_file.file_id + "_top_poses_pv.maegz"
            )
        )
    except Exception as exc:
        logger.error(f"Error in selecting top poses for {glide_docking_file.file_id}")
        raise exc
    return top_poses


def extract_pdb_top_poses(
    combind_result: DockedLigand, context: CombindContext, pdb_lig: str
) -> DockedLigand:
    r"""Parse the combind docking results top poses to get each complex in a separate file

    Args:
        combind_result : The DockedLigand object containing the combind docking results
        context : CombindContext with information about the combind docking
        pdb_lig : The id of the pdb ligand to be extracted from the complex.
    Returns:
        The DockedLigand object containing the top pose
    """
    pdb_file = None
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    glide_context = GlideContext.get_current()
    glide_context.write_dir = context.work_dir
    glide_context.set_current(glide_context)
    if not os.path.exists(context.work_dir):
        logger.warning(f"The directory {context.work_dir} does not exist.Creating it")
        os.makedirs(context.work_dir)
    if not os.path.exists(combind_result.file_path):
        logger.error(f"The file {combind_result.file_path} does not exist")
        raise FileNotFoundError(f"The file {combind_result.file_path} does not exist")
    try:
        if not os.path.exists(
            os.path.join(context.work_dir, combind_result.file_id + "_split_lig.maegz")
        ):
            manipulate_complexes(
                input_file=combind_result.file_path,
                context=glide_context,
                outfile_name=combind_result.file_id + "_merge.maegz",
                mode="merge",
            )
            logger.debug(f"Complexes file for {combind_result.file_id} is created")
            manipulate_complexes(
                input_file=os.path.join(
                    context.work_dir, combind_result.file_id + "_merge.maegz"
                ),
                context=glide_context,
                outfile_name=combind_result.file_id + "_split_lig.maegz",
                mode="split_ligand",
            )
            if os.path.exists(
                os.path.join(context.work_dir, combind_result.file_id + "_merge.maegz")
            ):
                os.remove(
                    os.path.join(
                        context.work_dir, combind_result.file_id + "_merge.maegz"
                    )
                )
        logger.debug(f"Top poses pdb complex for {combind_result.file_id} is extracted")
    except Exception as exc:
        logger.error(f"Error in selecting top poses for {combind_result.file_id}")
        raise exc
    property_csv = combind.extract_data_csv(
        os.path.join(context.work_dir, combind_result.file_id + "_split_lig.maegz"),
        combind_result.file_id + "_split_lig",
        filter_data=False,
    )
    logger.debug(f"The property csv file is {property_csv}")
    # read the property_csv and find the value of s_lp_Variant that contains the pdb_lig
    combind_data = pd.read_csv(property_csv)
    logger.debug(f"The combind data is {combind_data}")
    matches = combind_data[
        combind_data["s_m_title"].str.contains(pdb_lig, case=False, na=False)
    ]
    logger.debug(f"The matches are {matches} for {pdb_lig}")
    # find the index of the pdb ligand from the s_m_title column and extract the s_lp_Variant
    slpvariant = matches["s_lp_Variant"].values[0]
    parsed_file = Glide().parse_mae_files(
        os.path.join(context.work_dir, combind_result.file_id + "_split_lig.maegz"),
        glide_context,
        slpvariant,
    )
    # run mae to mol2 function then clean and standardize the file
    mol2_file = Glide().mae_to_mol(parsed_file.file_path, glide_context)
    pdb_file = clean_and_standardize_file(mol2_file.file_path)
    if os.path.exists(mol2_file.file_path):
        os.remove(mol2_file.file_path)
    if os.path.exists(parsed_file.file_path):
        os.remove(parsed_file.file_path)
    return DockedLigand(pdb_file)


def calc_rmsd_spyrmsd(
    reference_file: DockedLigand,
    target_file: DockedLigand,
    context: GlideContext,
) -> float:
    """
    Calculate the symmetry corrected RMSD between the reference and target ligands

    Args:
        reference_file : DockedLigand object containing the reference complex
        target_file : DockedLigand object containing the target complex
        context : GlideContext containing the information about the docking
    Returns:
        A float value of the symmetry corrected RMSD between the reference and target ligands
    """
    if not os.path.exists(reference_file.file_path):
        logger.error(f"The file {reference_file.file_path} does not exist")
        raise FileNotFoundError(f"The file {reference_file.file_path} does not exist")
    if not os.path.exists(target_file.file_path):
        logger.error(f"The file {target_file.file_path} does not exist")
        raise FileNotFoundError(f"The file {target_file.file_path} does not exist")
    if not os.path.exists(context.write_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    if reference_file.file_ext == "maegz" or reference_file.file_ext == "mae":
        logger.error(f"The file {reference_file.file_path} is not a pdb file")
        raise ValueError(f"The file {reference_file.file_path} is not a pdb file")
    ref_pdb = clean_and_standardize_file(reference_file.file_path)
    if target_file.file_ext == "maegz" or target_file.file_ext == "mae":
        logger.error(f"The file {target_file.file_path} is not a pdb file")
        raise ValueError(f"The file {target_file.file_path} is not a pdb file")
    target_pdb = clean_and_standardize_file(target_file.file_path)
    ref_lig_object = DockedLigand(ref_pdb)
    target_lig_object = DockedLigand(target_pdb)

    rmsd = RMSD(target_lig_object, ref_lig_object, context)
    res = rmsd.symmetry_rmsd()
    logger.info(f"The RMSD between the reference and target ligands is {res}")
    if os.path.exists(target_lig_object.file_path):
        os.remove(target_lig_object.file_path)
    if os.path.exists(ref_lig_object.file_path):
        os.remove(ref_lig_object.file_path)
    return res


def main(
    csv_file: str,
    output_files_dir: dict,
    save_file: str,
    parquet_file: str,
    pdbfile_dir: str,
    glide_only: bool = False,
):
    r"""
    Main function to run the docking wrapper.

    Args:
        csv_file : Path to the CSV file containing the PDB IDs to dock.
        output_files_dir : Dictionary where key is the input directory,
            value is the RMSD column name.
        save_file : Path to the CSV file where output will be saved.
        parquet_file : Path to the parquet file containing ligand data.
        pdbfile_dir : The directory where the pdb files are located.
        glide_only : If true, we don't run combind and only run glide docking
    """
    # Load or initialize the output DataFrame
    if os.path.exists(save_file):
        save_file_df = pd.read_csv(save_file)
        if save_file_df.empty:
            logger.warning(f"{save_file} is empty. Initializing new DataFrame.")
            save_file_df = pd.DataFrame(
                columns=["PDB_ID"] + list(output_files_dir.values())
            )
    else:
        logger.info(f"{save_file} does not exist. Creating new output file.")
        save_file_df = pd.DataFrame(
            columns=["PDB_ID"] + list(output_files_dir.values())
        )

    # Read PDB IDs to process
    pdb_df = pd.read_csv(csv_file)
    # limit the pdb_df to the one with 5ix0 pdb id
    if pdb_df.empty or "PDB_ID" not in pdb_df.columns:
        logger.error("CSV file contains no data or missing 'PDB_ID' column.")
        return

    # sort the rows based on the Number of PubChem Ligands column from smallest to largest
    pdb_df = pdb_df.sort_values(by=["Number of PubChem Ligands"])
    # exclude rows with 0 in the Total Water column
    pdb_df = pdb_df[pdb_df["Total Water"] != 0]
    pdbs = pdb_df["PDB_ID"].dropna().unique()
    logger.info(f"Number of PDB IDs to process: {len(pdbs)}")
    for pdbid in pdbs:
        logger.info(f"--- Processing {pdbid} ---")
        # chck if the pdb file exists in the pdbfile_dir
        if not os.path.exists(os.path.join(pdbfile_dir, f"{pdbid.lower()}.pdb")):
            logger.info(f"Downloading PDB file for {pdbid}")
            get_pdb_coordinates(pdbid, pdbfile_dir)
        pdbfile = os.path.join(pdbfile_dir, f"{pdbid.lower()}.pdb")
        ref_lig_n_list = (
            []
        )  # Initialize with a default value to avoid unassigned variable usage
        for directory, column_name in output_files_dir.items():
            # check the directory to determine water
            if "./water" in directory:
                water = "all"
            elif "./nt_water" in directory:
                water = None
            else:
                water = "interfacial"
            # check if the pdb_df["Bridge Water"] of the current row is 0 and if so skip
            if (
                pdb_df[pdb_df["PDB_ID"] == pdbid]["Bridge Water"].values[0] == 0
                and water == "interfacial"
            ):
                logger.info(f"Skipping {pdbid} due to no bridge water.")
                continue
            if (
                pdbid in save_file_df["PDB_ID"].values
                and not save_file_df.loc[save_file_df["PDB_ID"] == pdbid, column_name]
                .isnull()
                .all()
            ):
                logger.info(f"{pdbid} already has value in {column_name}. Skipping.")
                continue
            logger.info(f"Directory: {directory} -> Column: {column_name}")
            combind_context = CombindContext.get_current()
            combind_context.work_dir = directory
            combind_context.set_current(combind_context)
            glide = Glide()
            context = GlideContext.get_current()
            context.write_dir = directory
            context.samplewater = False if water is None else True
            context.set_current(context)
            ligand_members = extract_parquet_clusters(parquet_file, pdbid.upper())
            ligand_info = parse_ligand_members(
                ligand_members,
                pdbid.upper(),
                input_dir=pdbfile_dir,
                find_pdb_ligand=False,
                water=True,
                datatype="pandas",
            )
            pdb_ligands = parse_ligand_members(
                ligand_members,
                pdbid.upper(),
                input_dir=pdbfile_dir,
                find_pdb_ligand=True,
                water=True,
                datatype="pandas",
            )
            print(f"The ligand info is {ligand_info}")
            if pdb_ligands.empty:
                logger.warning(f"No valid ligand found for {pdbid}. Skipping.")
                continue
            pdb_lig = pdb_ligands["s_m_title"].values[0]
            logger.info(f"PDB ligand: {pdb_lig}")
            # check if the prepared protein file exists
            if not os.path.exists(
                os.path.join(directory, f"{pdbid.lower()}_protein_prepared.mae")
            ):
                # prepare the protein
                grid_file = prep_proteins(pdbfile, context, pdb_lig, type_water=water)
            else:
                logger.info(
                    f"Protein prepared file already exists: {pdbid.lower()}_protein_prepared.mae"
                )
                grid_file = PreparedProtein(
                    os.path.join(directory, f"{pdbid.lower()}_grid.zip")
                )
            prepared_path = PreparedLigand(
                os.path.join(directory, f"{pdbid.lower()}_protein_prepared.mae")
            )
            if os.path.exists(prepared_path.file_path) and not os.path.exists(
                prepared_path.file_path.replace(".mae", "_merged.mae")
            ):
                new_reflig_obj = glide.split_prepared_prot(
                    prepared_path.file_path, context, pdb_lig
                )

                ref_lig_n_list = [
                    glide.mae_to_mol(obj.file_path, context, n_structure=1)
                    for obj in new_reflig_obj
                ]
            # check if the prepared ligand file exists
            if not os.path.exists(
                os.path.join(directory, f"{pdbid.lower()}_lig_pubchem_prepared.mae")
            ):
                prepared_lig = prep_ligands(ligand_info, context, pdbid)
            else:
                logger.info(
                    f"Ligand prepared file already exists: {pdbid.lower()}_lig_pubchem_prepared.mae"
                )
                prepared_lig = PreparedLigand(
                    os.path.join(directory, f"{pdbid.lower()}_lig_pubchem_prepared.mae")
                )
            # check if the docking file exists
            if not os.path.exists(
                os.path.join(
                    directory,
                    f"{pdbid.lower()}_lig_pubchem_docking_pv.maegz",
                )
            ):
                try:
                    # dock the ligand
                    docked_ligands = dock_ligands(
                        grid_file,
                        prepared_lig,
                        context,
                    )
                except Exception as e:
                    logger.error(f"Docking failed for {pdbid} in {directory}: {e}")
                    continue
            else:
                logger.info(f"Docking result file already exists for {pdbid.lower()}")
                docked_ligands = DockedLigand(
                    os.path.join(
                        directory,
                        f"{pdbid.lower()}_lig_pubchem_docking_pv.maegz",
                    )
                )
            # check if the combind file exists
            if not glide_only:
                if not os.path.exists(
                    os.path.join(
                        directory,
                        f"{pdbid.upper()}_lig_pubchem_docking_pv_top_poses_pv.maegz",
                    )
                ):
                    if water is "interfacial":
                        continue
                    combind_context = CombindContext.get_current()
                    combind_context.work_dir = directory
                    combind_context.set_current(combind_context)
                    raw_combind_csv = run_combind(
                        docked_ligands, combind_context, screen=False
                    )
                    docking_result_path = get_top_combind_pose(
                        raw_combind_csv,
                        docked_ligands,
                        combind_context,
                    )

                else:
                    logger.info(
                        f"Combind result file already exists for {pdbid.upper()}"
                    )
                    docking_result_path = DockedLigand(
                        os.path.join(
                            directory,
                            f"{pdbid.upper()}_lig_pubchem_docking_pv_top_poses_pv.maegz",
                        )
                    )
                srodinger = extract_pdb_top_poses(
                    docking_result_path, combind_context, pdb_lig
                )
            else:
                docking_result_path = docked_ligands
                glide.sort_docking_results(
                    docking_result_path.file_path, context, best=1
                )

                sorted_path = os.path.join(
                    directory,
                    f"{pdbid.upper()}_lig_pubchem_docking_pv_sorted.maegz",
                )
                file_sorted = DockedLigand(sorted_path)

                srodinger = extract_pdb_top_poses(file_sorted, combind_context, pdb_lig)
            rmsd_list = []
            for ref_lig_n in ref_lig_n_list:
                try:
                    ref_object = clean_and_standardize_file(ref_lig_n.file_path)
                    rmsd = calc_rmsd_spyrmsd(
                        DockedLigand(file_path=ref_object),
                        srodinger,
                        context,
                    )
                    rmsd_list.append(rmsd)
                    logger.info(f"RMSD for {pdbid} in {directory} is {rmsd}")
                except Exception as e:
                    logger.error(
                        f"Failed to calculate RMSD for {pdbid} in {directory}: {e}"
                    )
                    continue

            if len(rmsd_list) == 1:
                rmsd_value = rmsd_list[0]
            else:
                rmsd_value = ",".join(map(str, rmsd_list))

            if pdbid not in save_file_df["PDB_ID"].values:
                new_row = {"PDB_ID": pdbid, column_name: rmsd_value}
                save_file_df = pd.concat(
                    [save_file_df, pd.DataFrame([new_row])], ignore_index=True
                )
            else:
                row_index = save_file_df.index[save_file_df["PDB_ID"] == pdbid][0]
                save_file_df.at[row_index, column_name] = rmsd_value

            save_file_df.to_csv(save_file, index=False)


if __name__ == "__main__":
    PARQUET_FILENAME = "final_ligand_cluster_0.7_Tc.parquet"
    CSV_FILE = "data_summary_v2.csv"
    SAVE_FILE = "merged_df_v2c.csv"
    FILE_PDB = "./raw"
    OUTPUT_FILES_DIR = {
        "./nt_water": "Combind RMSD Water Deletion",
        "./bridg_water": "Combind RMSD Selective Retention",
        "./water": "Combind RMSD Water Retention",
    }
    main(
        CSV_FILE,
        OUTPUT_FILES_DIR,
        SAVE_FILE,
        PARQUET_FILENAME,
        FILE_PDB,
        glide_only=False,
    )
