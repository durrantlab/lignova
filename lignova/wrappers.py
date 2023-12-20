r" Implemtnation for a wrapper to use lignova for validation"
from typing import TextIO, Union

import glob
import os

from line_profiler import profile
from loguru import logger

from lignova.docking import Glide
from lignova.docking.combind import Combind
from lignova.docking.contexts import GlideContext
from lignova.docking.contexts.combind import CombindContext
from lignova.structure.editing import write_mda_universe
from lignova.structure.ligand import DockedLigand, Ligand
from lignova.structure.protein import Protein
from lignova.structure.utils import is_xray_structure, separate_protein_ligand


@profile
def clean_cluster_files(file_path: str, delim: list = ["[", "]"]):
    """
    This function takes a file containing a list of PDB IDs and returns a list of PDB IDs.
    Parameters
    ----------
    file_path : str
        The path to the file containing the PDB IDs.
    delim : str, optional
        The delimiter used in the file. The default is ":".
    Returns
    -------
    pdb_ids : list
        A list of PDB IDs.
    """
    with open(file_path, "r", encoding="utf-8") as file:
        pdb_ids = []
        for line in file:
            if line.startswith("Cluster: []"):
                continue
            if line.startswith("Cluster:"):
                cluster_id = (
                    line.split(":")[1]
                    .strip()
                    .replace(delim[0], "")
                    .replace(delim[1], "")
                    .replace("\n", "")
                )
                pdb_ids.extend(cluster_id.split(", "))
        pdb_ids = [x for x in pdb_ids if x not in (" ", "''")]
        return pdb_ids


@profile
def get_coordinates(
    pdb_ids: Union[str, TextIO, list], work_dir: str, limit: int = 1000
):
    """
    This function takes a list of PDB IDs and downloads the PDB files to a specified directory.
    Parameters
    ----------
    pdb_ids : Union[str,TextIO,list]
        A list of PDB IDs to be downloaded.
    work_dir : str
        The working directory where the PDB files will be downloaded.
    limit : int, optional
        The number of PDB IDs to be downloaded. The default is 50.
    Returns
    -------
    None.
    """
    current_dir = os.getcwd()
    protein = Protein()
    # check if the output directory exists
    if not os.path.exists(work_dir):
        logger.info("Output_Dir not found,Creating it in working directory")
        os.mkdir(os.path.join(current_dir, work_dir))
    # check if the input is a list or a file
    if isinstance(pdb_ids, list):
        logger.info("Input is a list")
        for pdb_id in pdb_ids:
            logger.info(f"Downloading PDB file for {pdb_id}")
            file_ext = (
                "pdb"
                if protein.get_pdb_from_rcsb(pdb_id).startswith("HEADER")
                else "cif"
            )
            protein.load(
                pdb_id,
                write=True,
                write_path=os.path.join(work_dir, pdb_id + "." + file_ext),
            )
    elif isinstance(pdb_ids, str):
        if os.path.exists(pdb_ids):
            logger.info("Input is a file")
            pdb_ids = clean_cluster_files(
                pdb_ids, delim=["[", "]"]
            )  # NOTE:change this delim when dealing with filtered to ( and )
            logger.info(f"Found a total of {len(pdb_ids)} PDB IDs in the file")
            limit = len(pdb_ids) if len(pdb_ids) < limit else limit
            for i in range(limit):
                if "|" in pdb_ids[i]:
                    pdb_id = pdb_ids[i].split("_")[0]
                    tmp = pdb_id.strip("'")
                else:
                    pdb_id = pdb_ids[i]
                    tmp = pdb_id.lstrip("'").rstrip("'")
                if os.path.exists(os.path.join(work_dir, tmp.lower() + ".pdb")):
                    logger.info(f"PDB file for {tmp} already exists")
                    continue

                logger.info(f"Downloading PDB file for {tmp}")
                file_ext = (
                    "pdb"
                    if protein.get_pdb_from_rcsb(tmp).startswith("HEADER")
                    else "cif"
                )
                protein.load(
                    pdb_id=tmp,
                    write=True,
                    write_path=os.path.join(work_dir, tmp.lower() + "." + file_ext),
                )
        else:
            # check if the input is a single pdb id
            logger.info("Input is a single PDB ID")
            pdb_id = pdb_ids
            logger.info(f"Downloading PDB file for {pdb_id}")
            file_ext = (
                "pdb"
                if protein.get_pdb_from_rcsb(pdb_id).startswith("HEADER")
                else "cif"
            )
            protein.load(
                pdb_id=pdb_id,
                write=True,
                write_path=os.path.join(work_dir, pdb_id.lower() + "." + file_ext),
            )


@profile
def prep_structure(
    input_dir: str, output_dir: str, pdb_id: Union[str, list, None], limit: int = 500
):
    """
    This function takes a PDB ID and a ligand ID and returns a Structure object.
    Parameters
    ----------
    directory : str
        The directory where the PDB file is located.
    pdb_id : Optional[str]
        The PDB ID of the protein. The default is None.
    limit : int
        The number of PDB IDs to prep. The default is 50.
    """
    glide = Glide()
    failed = []
    if os.path.exists(os.path.join(output_dir, "failed.txt")):
        with open(
            os.path.join(output_dir, "failed.txt"), "r", encoding="utf-8"
        ) as file:
            file_line = file.read().splitlines()
        failed = list(set(failed).union(set(file_line)))
    context = GlideContext.get_current()
    context.write_dir = output_dir
    context.set_current(context)
    if not os.path.exists(input_dir):
        raise FileNotFoundError("Input directory not found")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if pdb_id is not None:
        if isinstance(pdb_id, str):
            pdb_id = [pdb_id]
    else:
        # NOTE:REMEMBER TO DEAL WITH THE CIF FILES
        pdb_id = glob.glob(os.path.join(input_dir, "*.pdb"))
        logger.info(f"Found a total of {len(pdb_id)} PDB files in the directory")
        if len(pdb_id) == 0:
            raise FileNotFoundError("No PDB files found in the directory")
    limit = min(len(pdb_id), limit)
    for pdb_file in pdb_id[:limit]:
        ids = os.path.basename(pdb_file).split(".")[0]
        if not os.path.exists(
            os.path.join(input_dir, ids.lower() + ".pdb")
        ) and not os.path.exists(os.path.join(input_dir, ids.lower() + ".cif")):
            logger.warning(f"{ids} not found in the directory. Downloading it now")
            get_coordinates(ids, input_dir)
        # check if the pdb_file has an extention and if so find the file in the input dir
        if "." not in pdb_file:
            pdb_file = glob.glob(os.path.join(input_dir, ids.lower() + ".*"))[0]
        prot = Protein(os.path.join(input_dir, pdb_file))
        logger.debug(prot.file_path)
        if "_lig" not in pdb_file and not os.path.exists(
            os.path.join(input_dir, pdb_file.replace(".pdb", "_lig.pdb"))
        ):
            if is_xray_structure(prot.file_path):
                logger.info(f"Separating {pdb_file}")
                protein, ligand = separate_protein_ligand(prot.file_path)
                if len(ligand.atoms) == 0:
                    logger.warning(f"No ligand found in {prot.file_name}")
                    os.remove(pdb_file)
                    # add the pdb file to the failed list
                    failed.append(prot.file_name)
                    continue
                write_mda_universe(protein, os.path.join(input_dir, pdb_file))
                # NOTE:I can have this named with the lig name should i do it?
                write_mda_universe(
                    ligand,
                    os.path.join(input_dir, pdb_file.replace(".pdb", "_lig.pdb")),
                )
            else:
                logger.warning(f"{pdb_file} is not an X-ray structure")
                os.remove(pdb_file)
                continue
        if (
            len(
                glob.glob(os.path.join(output_dir, prot.file_name.replace(".pdb", "*")))
            )
            == 2
            or "_lig" in pdb_file
            or prot.file_name in failed
        ):
            logger.info(f"{prot.file_name} already prepped")
            continue
        raw_prot = Protein(os.path.join(input_dir, pdb_file))
        lig_file = os.path.join(input_dir, raw_prot.file_id + "_lig.pdb")
        raw_lig = Ligand(file_path=lig_file)
        glide.convert_to_mae(raw_lig, context)
        lig_mae = Ligand(os.path.join(context.write_dir, raw_lig.file_id + ".mae"))
        glide.convert_to_mae(raw_prot, context)
        prot_mae = Protein(
            file_path=os.path.join(context.write_dir, raw_prot.file_id + ".mae")
        )
        try:
            glide.PrepLigand(lig_mae, context)
            glide.PrepProtein(prot_mae, context)
            os.remove(os.path.join(context.write_dir, raw_prot.file_id + "_grid.log"))
            os.remove(os.path.join(context.write_dir, prot_mae.file_name))
            os.remove(os.path.join(context.write_dir, raw_prot.file_id + "_lig.mae"))
            os.remove(
                os.path.join(
                    context.write_dir,
                    raw_prot.file_id + "_protein_prepared.mae",
                )
            )
        except Exception:
            logger.warning(f"Could not prepare {prot.file_name}")
            # save the pdb file in a folder called failed
            failed.append(prot.file_name)
            # remove any files in the output directory With the same name
            for file in glob.glob(
                os.path.join(output_dir, prot.file_name.replace(".pdb", "*"))
            ):
                os.remove(file)
            continue
    if len(failed) > 0:
        with open(
            os.path.join(output_dir, "failed.txt"), "w", encoding="utf-8"
        ) as file:
            file.write("\n".join(failed))


@profile
def dock_ligand(
    input_dir: str, output_dir: str, pdb: [str, list, None], limit: int = 10
):
    """
    Dock the ligand to the protein
    Parameters
    ----------
    input_dir : str
        The directory containing the prepped protein and ligand files
    output_dir : str
        The directory to save the docked files
    pdb : Optional[str,list]
        The pdb file to dock. If None, all the pdb files in the input directory will be docked
    limit : int, optional
        The number of files to dock, by default 10
    """
    failed = []
    if os.path.exists(os.path.join(output_dir, "failed.txt")):
        with open(
            os.path.join(output_dir, "failed.txt"), "r", encoding="utf-8"
        ) as file:
            file_line = file.read().splitlines()
        failed = list(set(failed).union(set(file_line)))
    # check if the output directory exists and create it if it does not
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # check if the input directory exists and if not raise an error
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"{input_dir} does not exist")
    # check if the pdb is not none and if it is not a list, raise an error
    if pdb is not None:
        if isinstance(pdb, str):
            pdb = [pdb]
    else:
        pdb = glob.glob(os.path.join(input_dir, "*_grid.zip"))
        logger.info(f"Found {len(pdb)} docking files")
        if len(pdb) == 0:
            raise FileNotFoundError(f"No docking files found in {input_dir}")
    glide = Glide()
    context = GlideContext.get_current()
    context.write_dir = output_dir
    context.set_current(context)
    limit = min(limit, len(pdb))
    pdb = pdb[:limit]
    for pdb_file in pdb:
        logger.info(f"Docking {pdb_file}")
        prep_prot = Protein(pdb_file)
        prep_lig = Ligand(
            os.path.join(
                input_dir, prep_prot.file_id.replace("grid", "lig_prepared.mae")
            )
        )
        logger.info(f"Docking {prep_lig.file_name} to {prep_prot.file_name}")
        if os.path.exists(
            os.path.join(
                context.write_dir,
                prep_prot.file_id.replace("grid", "lig_docking_pv.maegz"),
            )
        ) and os.path.exists(
            os.path.join(
                context.write_dir, prep_prot.file_id.replace("grid", "lig_docking.csv")
            )
            or prep_lig.file_name in failed
        ):
            logger.info(f"{prep_lig.file_name} already docked")
            continue
        try:
            glide.run(prep_prot, prep_lig, context)
            os.remove(
                os.path.join(
                    context.write_dir,
                    prep_prot.file_id.replace("grid", "lig_docking.in"),
                )
            )
            logger.info("Docking is completed")
            os.remove(
                os.path.join(
                    context.write_dir,
                    prep_prot.file_id.replace("grid", "lig_docking.log"),
                )
            )
            os.remove(
                os.path.join(
                    context.write_dir,
                    prep_prot.file_id.replace("grid", "lig_docking_skip.csv"),
                )
            )
        except Exception:
            logger.warning(
                f"Could not dock {prep_lig.file_name} to {prep_prot.file_name}"
            )
            # save the pdb file in a folder called failed
            failed.append(prep_lig.file_name)
            # remove any files in the output directory With the same name
            continue
    if len(failed) > 0:
        with open(
            os.path.join(output_dir, "failed.txt"), "w", encoding="utf-8"
        ) as file:
            file.write("\n".join(failed))


def combind_pose_selction(
    input_dir: str, output_dir: str, pdb: Union[str, list, None], limit: int = 50
):
    """
    Combine the poses from the docking
    Parameters
    ----------
    input_dir : str
        The directory containing the docked files
    output_dir : str
        The directory to save the combined files
    pdb : Optional[str,list]
        The pdb file to combine. If None, all the pdb files in the input directory will be combined
    limit : int, optional
        The number of files to combine, by default 50
    """
    context = CombindContext.get_current()
    context.work_dir = output_dir
    context.set_current(context)
    # check if the output directory exists and create it if it does not
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # check if the input directory exists and if not raise an error
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"{input_dir} does not exist")
    # check if the pdb is not none and if it is not a list, raise an error
    if pdb is not None:
        if isinstance(pdb, str):
            pdb = [pdb]
    else:
        pdb = glob.glob(os.path.join(input_dir, "*_docking_pv.maegz"))
        logger.info(f"Found {len(pdb)} docking files")
        if len(pdb) == 0:
            raise FileNotFoundError(f"No docking files found in {input_dir}")
    limit = min(limit, len(pdb))
    pdb = pdb[:limit]
    for pdb_file in pdb:
        docking_file = DockedLigand(pdb_file)
        logger.info(f"Working on {docking_file.file_name}")
        combind = Combind(
            command=context.command,
            work_dir=context.work_dir,
            schrodinger=context.schrodinger,
            schrodinger_env=context.schrodinger_env,
        )
        if not os.path.exists(
            os.path.join(
                context.work_dir, docking_file.file_name.split("_pv")[0] + "_features"
            )
        ):
            combind.featurize(
                docking_file.file_path, docking_file.file_name.split("_pv")[0]
            )
            files = glob.glob(
                os.path.join(input_dir, docking_file.file_name.split("_pv")[0] + "*")
            )
            if len(files) > 2:
                for i in files:
                    if i != docking_file.file_path and i != os.path.join(
                        input_dir, docking_file.file_name.split("_pv")[0] + ".csv"
                    ):
                        os.remove(i)
            os.remove(
                os.path.join(
                    context.work_dir,
                    docking_file.file_name.split("_pv")[0] + "_features.log",
                )
            )
        if not os.path.exists(
            os.path.join(
                context.work_dir, docking_file.file_name.split("_pv")[0] + "_screen.npy"
            )
        ):
            combind.compute_combind_score(
                features_dir=os.path.join(
                    context.work_dir,
                    docking_file.file_name.split("_pv")[0] + "_features",
                ),
                filename=docking_file.file_name.split("_pv")[0],
            )
        if not os.path.exists(
            os.path.join(
                context.work_dir,
                docking_file.file_name.split("_pv")[0] + "_combind_sorted.maegz",
            )
        ):
            combind.apply_combind_score(
                docking_filepath=docking_file.file_path,
                combind_score_file=os.path.join(
                    context.work_dir,
                    docking_file.file_name.split("_pv")[0] + "_screen.npy",
                ),
                output_filename=docking_file.file_name.split("_pv")[0],
            )
            os.remove(
                os.path.join(
                    context.work_dir,
                    docking_file.file_name.split("_pv")[0] + "_combind_sort.log",
                )
            )
            os.remove(
                os.path.join(
                    context.work_dir,
                    docking_file.file_name.split("_pv")[0] + "_combind.maegz",
                )
            )
        if not os.path.exists(
            os.path.join(
                context.work_dir,
                docking_file.file_name.split("_pv")[0] + "_combind.csv",
            )
        ):
            combind.extract_data_csv(
                docking_file=os.path.join(
                    context.work_dir,
                    docking_file.file_name.split("_pv")[0] + "_combind_sorted.maegz",
                ),
                filename=docking_file.file_name.split("_pv")[0] + "_combind",
            )


if __name__ == "__main__":
    RAW_FILE = "/home/mma121/PubChem_small/try_schrodinger/clusters.csv"
    FILTERED_FILE = (
        "/home/mma121/PubChem_small/try_schrodinger/new_clusters_pdb_postfilter.csv"
    )
    RAW_INPUT_DIR = "/home/mma121/PubChem_small/representatives"
    PREPPED_DIR = "./prepped"
    DOCKED_DIR = "./docked"
    COMBIND_DIR = "./combind"
    # get_coordinates(RAW_FILE, RAW_INPUT_DIR, limit=10000)
    pdb_files = glob.glob(os.path.join(RAW_INPUT_DIR, "*.pdb"))
    logger.info(f"{len(pdb_files)} were found in the input directory")
    cif_files = glob.glob(os.path.join(RAW_INPUT_DIR, "*.cif"))
    prep_structure(
        RAW_INPUT_DIR,
        PREPPED_DIR,
        pdb_id=pdb_files,
        limit=500,
    )
    dock_ligand(PREPPED_DIR, DOCKED_DIR, pdb=None, limit=500)
    combind_pose_selction(DOCKED_DIR, COMBIND_DIR, pdb=None, limit=500)
