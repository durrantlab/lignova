r" Implemtnation for a wrapper to use lignova for validation"
from typing import TextIO

import ast
import glob
import os
import shutil

import pandas as pd
from line_profiler import profile
from loguru import logger

from lignova.analysis.rmsd import RMSD
from lignova.analysis.utils import obabel_result_parser
from lignova.docking import Glide
from lignova.docking.combind import Combind
from lignova.docking.contexts import GlideContext
from lignova.docking.contexts.combind import CombindContext
from lignova.docking.utils import convert_to_pdb, manipulate_complexes
from lignova.structure.editing import write_mda_universe
from lignova.structure.ligand import DockedLigand, Ligand
from lignova.structure.protein import Protein
from lignova.structure.utils import (
    is_xray_structure,
    separate_protein_ligand,
    validate_ligands,
    validate_pdb,
)


@profile
def find_cluster_reps(file_path: str, delim: str = ":"):
    """
    This function takes the mmseqs cluster file and returns the representative PDB IDs.
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
            if line.startswith("Cluster "):
                logger.debug(line.split(delim)[1].rstrip())
                if "|" in line.split(delim)[1].rstrip():
                    # split the line by _ and get the first element
                    cluster_id = line.split(delim)[1].rstrip().split("_")[0]
                    cluster_id = cluster_id.split("|")[0]
                    cluster_id = cluster_id.strip()
                    pdb_ids.append(cluster_id)
                else:
                    cluster_id = line.split(delim)[1].rstrip()
                    cluster_id = ast.literal_eval(cluster_id)
                    pdb_ids.extend(cluster_id)
            else:
                continue
        pdb_ids = [x for x in pdb_ids if x not in (" ", "''")]
        return pdb_ids


@profile
def get_coordinates(pdb_ids: str | list, work_dir: str, limit: None | int = 500):
    """
    This function takes a list of PDB IDs and downloads the PDB files to a specified directory
    if they pass the validation test. (no mutation, has ligand, no covalent bond, x-ray structures,)
    Parameters
    ----------
    pdb_ids : str|list
        A list of PDB IDs or a single PDB ids to be downloaded
    work_dir : str
        The working directory where the PDB files will be downloaded.
    limit : int, optional
        The number of PDB IDs to be downloaded. The default is 500.
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
    if isinstance(pdb_ids, str):
        pdb_ids = [pdb_ids]
    limit = min(len(pdb_ids), limit)
    for i in range(limit):
        if "|" in pdb_ids[i]:
            pdb_id = pdb_ids[i].split("_")[0]
            tmp = pdb_id.strip("'")
        else:
            pdb_id = pdb_ids[i]
            tmp = pdb_id.lstrip("'").rstrip("'")
        logger.info(f"Working on {tmp}")
        if os.path.exists(os.path.join(work_dir, tmp.lower() + ".pdb")):
            logger.info(f"PDB file for {tmp} already exists")
            continue
        if validate_pdb(tmp) and validate_ligands(tmp):
            logger.info(f"Downloading PDB file for {tmp}")
            file_ext = (
                "pdb" if protein.get_pdb_from_rcsb(tmp).startswith("HEADER") else "cif"
            )
            protein.load(
                pdb_id=tmp,
                write=True,
                write_path=os.path.join(work_dir, tmp.lower() + "." + file_ext),
            )
        else:
            logger.warning(f"{tmp} failed validation test")


@profile
def prep_structure(
    input_dir: str, output_dir: str, pdb_id: str | list | None, limit: int = 50
):
    """
    This function takes a PDB ID and a ligand ID and returns a Structure object.
    Parameters
    ----------
    input_dir : str
        The directory containing the PDB files.
    output_dir : str
        The directory to save the prepped files.
    pdb_id : str|list
        The PDB ID to be prepped. If None, all the PDB files in the input directory will be prepped.
    limit : int, optional
        The number of PDB files to be prepped. The default is 50.
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
    context.samplewater = False
    context.prot_watdist = "0"
    context.set_current(context)
    if not os.path.exists(input_dir):
        raise FileNotFoundError("Input directory not found")
    if not os.path.exists(output_dir):
        logger.info("Output_Dir not found,Creating it in working directory")
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
        if not validate_ligands(ids) and not validate_pdb(ids):
            logger.warning(f"{ids} failed validation test")
            failed.append(ids)
            continue
        if not os.path.exists(
            os.path.join(input_dir, ids.lower() + ".pdb")
        ) and not os.path.exists(os.path.join(input_dir, ids.lower() + ".cif")):
            logger.warning(f"{ids} not found in the directory. Downloading it now")
            get_coordinates(ids, input_dir)
        pdb_file = glob.glob(os.path.join(input_dir, ids.lower() + ".pdb"))[0]
        prot = Protein(pdb_file)
        logger.debug(f"Protein file path:{prot.file_path}")
        logger.debug(prot.file_name.replace(".pdb", "_lig.pdb"))
        if not os.path.exists(
            os.path.join(input_dir, prot.file_name.replace(".pdb", "_lig.pdb"))
        ):
            logger.info(f"Separating {prot.file_path} into protein and ligand")
            protein, ligand = separate_protein_ligand(
                prot.file_path, remove_water=True, keep_het_chain="A"
            )
            write_mda_universe(protein, pdb_file)
            write_mda_universe(
                ligand,
                prot.file_path.replace(".pdb", "_lig.pdb"),
            )
        if (
            (
                os.path.exists(
                    os.path.join(
                        output_dir, prot.file_name.replace(".pdb", "_lig_prepared.mae")
                    )
                )
                and os.path.exists(
                    os.path.join(
                        output_dir,
                        prot.file_name.replace(".pdb", "_protein_prepared.mae"),
                    )
                )
            )
            or "_lig" in pdb_file
            or prot.file_name in failed
        ):
            logger.info(f"{prot.file_name} already prepped")
            continue
        logger.info(
            f"Prepping {prot.file_name} and {prot.file_name.replace('.pdb','_lig.pdb')}"
        )
        lig_file = os.path.join(input_dir, prot.file_id + "_lig.pdb")
        raw_lig = Ligand(file_path=lig_file)
        if not os.path.exists(raw_lig.file_path.replace(".pdb", ".mae")):
            glide.convert_to_mae(raw_lig, context)
        lig_mae = Ligand(os.path.join(context.write_dir, raw_lig.file_id + ".mae"))
        if not os.path.exists(prot.file_path.replace(".pdb", ".mae")):
            glide.convert_to_mae(prot, context)
        prot_mae = Protein(
            file_path=os.path.join(context.write_dir, prot.file_id + ".mae")
        )
        try:
            if not os.path.exists(lig_mae.file_path.replace(".mae", "_prepared.mae")):
                glide.PrepLigand(lig_mae, context)
            if not os.path.exists(
                prot_mae.file_path.replace(".mae", "_protein_prepared.mae")
            ):
                glide.PrepProtein(prot_mae, context)
            os.remove(os.path.join(context.write_dir, prot.file_id + "_grid.log"))
        except Exception:
            logger.warning(f"Could not prepare {prot.file_name}")
            # save the pdb file in a folder called failed
            # logger.debug(f"Removing {prot.file_name} from the directory")
            failed.append(prot.file_name)
            """
            # remove any files in the output directory With the same name
            for file in glob.glob(
                os.path.join(output_dir, prot.file_name.replace(".pdb", "*"))
            ):
                logger.debug(f"Removing {file}")
                os.remove(file)
            continue
        """
    if len(failed) > 0:
        with open(
            os.path.join(output_dir, "failed.txt"), "w", encoding="utf-8"
        ) as file:
            file.write("\n".join(failed))


@profile
def dock_ligand(
    input_dir: str, output_dir: str, pdb: str | list | None = None, limit: int = 10
):
    """
    Dock the ligand to the protein
    Parameters
    ----------
    input_dir : str
        The directory containing the prepped protein and ligand files
    output_dir : str
        The directory to save the docked files
    pdb : Optional[str,list,None]
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
            raise FileNotFoundError(
                f"No prepared structures files found in {input_dir}"
            )
    glide = Glide()
    context = GlideContext.get_current()
    context.write_dir = output_dir
    context.set_current(context)
    limit = min(limit, len(pdb))
    pdb = pdb[:limit]
    for pdb_file in pdb:
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
            if not os.path.exists(
                os.path.join(
                    context.write_dir,
                    prep_prot.file_id.replace("grid", "lig_docking_pv_sorted.maegz"),
                )
            ):
                logger.debug(
                    os.path.join(
                        context.write_dir,
                        prep_prot.file_id.replace("grid", "lig_docking_pv.maegz"),
                    )
                )
                results = os.path.join(
                    context.write_dir,
                    prep_prot.file_id.replace("grid", "lig_docking_pv.maegz"),
                )
                glide.sort_docking_results(results, context)
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
    input_dir: str, output_dir: str, pdb: str | list | None, limit: int = 50
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
    # check if the failed file exists and if so, read it into failed_list
    failed_list = []
    if os.path.exists(os.path.join(output_dir, "failed.txt")):
        with open(
            os.path.join(output_dir, "failed.txt"), "r", encoding="utf-8"
        ) as file:
            failed_list = file.read().splitlines()
    for pdb_file in pdb:
        if "." not in pdb_file:
            pdb_file = glob.glob(
                os.path.join(input_dir, pdb_file.lower() + "_lig_docking_pv.maegz")
            )[0]
            # if not found, and not in the failed.txt file in the input dir,
            # then run the dock_ligand function after logging a warning
            while not os.path.exists(pdb_file) and pdb_file not in os.path.join(
                input_dir, "failed.txt"
            ):
                logger.warning(
                    f"{pdb_file} not found in the directory. Prepping it now"
                )
                dock_ligand(
                    input_dir, output_dir, pdb_id=pdb_file.split("_docking_pv.maegz")[0]
                )
        docking_file = DockedLigand(pdb_file)
        logger.info(f"Working on {docking_file.file_name}")
        combind = Combind(
            command=context.command,
            work_dir=context.work_dir,
            schrodinger=context.schrodinger,
            schrodinger_env=context.schrodinger_env,
        )
        try:
            if (
                not os.path.exists(
                    os.path.join(
                        context.work_dir,
                        docking_file.file_name.split("_pv")[0] + "_features",
                    )
                )
                and docking_file.file_name not in failed_list
            ):
                combind.featurize(
                    docking_file.file_path, docking_file.file_name.split("_pv")[0]
                )
                files = glob.glob(
                    os.path.join(
                        input_dir, docking_file.file_name.split("_pv")[0] + "*"
                    )
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
            if (
                not os.path.exists(
                    os.path.join(
                        context.work_dir,
                        docking_file.file_name.split("_pv")[0] + "_screen.npy",
                    )
                )
                and docking_file.file_name not in failed_list
            ):
                combind.compute_combind_score(
                    features_dir=os.path.join(
                        context.work_dir,
                        docking_file.file_name.split("_pv")[0] + "_features",
                    ),
                    filename=docking_file.file_name.split("_pv")[0],
                )
            if (
                not os.path.exists(
                    os.path.join(
                        context.work_dir,
                        docking_file.file_name.split("_pv")[0]
                        + "_combind_sorted.maegz",
                    )
                )
                and docking_file.file_name not in failed_list
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
            if (
                not os.path.exists(
                    os.path.join(
                        context.work_dir,
                        docking_file.file_name.split("_pv")[0] + "_combind.csv",
                    )
                )
                and docking_file.file_name not in failed_list
            ):
                combind.extract_data_csv(
                    docking_file=os.path.join(
                        context.work_dir,
                        docking_file.file_name.split("_pv")[0]
                        + "_combind_sorted.maegz",
                    ),
                    filename=docking_file.file_name.split("_pv")[0] + "_combind",
                )
        except Exception as e:
            failed_list.append(docking_file.file_name)
            logger.error(f"Failed to process {docking_file.file_name}: {str(e)}")
            # Delete files related to the failed docking file
            files = glob.glob(
                os.path.join(
                    context.work_dir, docking_file.file_name.split("_pv")[0] + "*"
                )
            )
            for file in files:
                if not file.endswith(".log"):
                    # if it is a directory, delete it and its contents
                    if os.path.isdir(file):
                        for subfile in glob.glob(os.path.join(file, "*")):
                            os.remove(subfile)
                        os.rmdir(file)
                    else:
                        os.remove(file)
            files = glob.glob(
                os.path.join(input_dir, docking_file.file_name.split("_pv")[0] + "*")
            )
            if len(files) > 2:
                for i in files:
                    if i != docking_file.file_path and i != os.path.join(
                        input_dir, docking_file.file_name.split("_pv")[0] + ".csv"
                    ):
                        os.remove(i)
    failed_list = list(set(failed_list))
    # Write failed list to a file
    with open(os.path.join(output_dir, "failed.txt"), "w", encoding="utf-8") as file:
        file.write("\n".join(failed_list))


def calculate_rmsd_mda(
    input_dir: str,
    output_dir: str,
    pdb: str | list | None,
    reference_dir: str,
    limit: int,
):
    """
    Parse the docking files to calculate the rmsd using mdanalysis
    Parameters
    ----------
    input_dir : str
        The directory containing the docked files
    output_dir : str
        The file to save the parsed files
    pdb : Optional[str,list]
        The pdb file to parse. If None, all the pdb files in the input directory will be parsed
    reference_dir : str
        The directory containing the reference files
    limit : int, optional
        The number of files to parse, by default 50
    """
    context = GlideContext.get_current()
    context.write_dir = output_dir
    context.set_current(context)
    # check if the input directory exists and if not raise an error
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"{input_dir} does not exist")
    # check if the output directory exists and create it if it does not
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # check if pdb is not none and if it is not a list, raise an error
    if pdb is not None:
        if isinstance(pdb, str):
            pdb = [pdb]
    # find all the docking files in the input directory
    else:
        pdb = glob.glob(os.path.join(input_dir, "*.maegz"))
        logger.info(f"Found {len(pdb)} docking files")
        if len(pdb) == 0:
            raise FileNotFoundError(f"No docking files found in {input_dir}")
    limit = min(limit, len(pdb))
    pdb = pdb[:limit]
    for docked_pdb_file in pdb:
        logger.info(f"Working on {docked_pdb_file}")
        # get the pdb_id from the file name
        pdb_id = os.path.basename(docked_pdb_file).split("_")[0]
        # get the reference file from the reference directory
        reference_file = glob.glob(os.path.join(reference_dir, pdb_id + ".pdb"))[0]
        # get the docking file
        docking_file = DockedLigand(docked_pdb_file)
        ref = Protein(reference_file)
        try:
            rmsd = RMSD(docking_file, ref, context)
            val = rmsd.rmsd_mda(
                ref_lig=os.path.join(reference_dir, pdb_id + "_lig.pdb")
            )
            logger.info(f"RMSD for {docking_file.file_name} is {val}")
            # save the rmsd value with the pdb name to a file in the output directory as a csv file
            with open(os.path.join(output_dir, pdb_id + ".csv"), "a") as file:
                file.write("PDB_ID,RMSD\n")
                file.write(f"{pdb_id},{val}")
            # remove the .maegz and .pdb files from the output directory using glob
            for file in glob.glob(os.path.join(output_dir, "*.maegz")) and glob.glob(
                os.path.join(output_dir, pdb_id + "*.pdb")
            ):
                os.remove(file)
        except Exception as e:
            logger.error(f"Could not calculate RMSD for {docking_file.file_name}")
            logger.error(str(e))
            continue


def calculate_rmsd_schrodinger(
    dock_ligand_dir: str,
    reference_dir: str,
    csv_file_name: str,
    number: str | int = 3,
):
    """
    Calculate the rmsd between the docked ligand and the reference ligand using rmsd_schrodinger
    Parameters
    ----------
    dock_ligand_dir : str
        The directory containing the docked ligand files
    reference_dir : str
        The directory containing the reference ligand files
    number : str|int
        The number of ligand to calculate the rmsd for, by default 3
    csv_file_name : str
        The name of the csv file to save the rmsd values
    """
    done = []
    # check if the dock_ligand_dir exists
    if not os.path.exists(dock_ligand_dir):
        raise FileNotFoundError(f"{dock_ligand_dir} does not exist")
    # check if the reference_dir exists
    if not os.path.exists(reference_dir):
        raise FileNotFoundError(f"{reference_dir} does not exist")
    # check if the csv_file_name exists in the dock_ligand_dir
    if os.path.exists(os.path.join(dock_ligand_dir, csv_file_name)):
        logger.info(f"{csv_file_name} already exists in {dock_ligand_dir}")
        # read the csv file using pandas
        rmsd_df = pd.read_csv(os.path.join(dock_ligand_dir, csv_file_name))
        # get the 2nd column of the csv file
        column = rmsd_df.iloc[:, 1]
        done.extend(list(set(column)))
    logger.info(done)
    logger.info(f"Found {len(done)} PDB IDs already done")
    # make a pandas dataframe to save the rmsd values
    rmsd_df = pd.DataFrame(
        columns=[
            "Index",
            "Title",
            "Mode",
            "RMSD",
            "Max dist.",
            "Max dist atom index pair",
            "ASL",
        ]
    )
    # get all the docked ligand files in the dock_ligand_dir
    docked_ligands = glob.glob(os.path.join(dock_ligand_dir, "*.maegz"))
    # loop over the docked ligand files
    for docked_ligand in docked_ligands:
        logger.info(f"Working on {docked_ligand}")
        # get the pdb_id from the file name
        pdb_id = os.path.basename(docked_ligand).split("_")[0]
        if pdb_id.upper() in done:
            logger.info(f"{pdb_id} already done")
            continue
        # get the reference ligand file from the reference_dir
        # check if the reference ligand file exists and if not continue
        if not os.path.exists(
            os.path.join(reference_dir, pdb_id + "_protein_prepared.mae")
        ):
            logger.warning(f"{pdb_id} not found in the reference directory")
            continue
        reference_ligand = glob.glob(
            os.path.join(reference_dir, pdb_id + "*_protein_prepared.mae")
        )[0]
        docked_ligand = DockedLigand(docked_ligand)
        reference_ligand = Ligand(reference_ligand)
        context = GlideContext.get_current()
        context.write_dir = dock_ligand_dir
        context.set_current(context)
        rmsd = RMSD(docked_ligand, reference_ligand, context=context)
        logger.info(f"Calculating RMSD for {docked_ligand.file_name}")
        rmsd.rmsd_schrodinger("rmsd.csv")
        # load the rmsd.csv file using pandas and get the first 3 rows excluding the header
        new_rmsd_df = pd.read_csv("rmsd.csv")
        new_rmsd_df = new_rmsd_df.iloc[1 : number + 1]
        rmsd_df = pd.concat([rmsd_df, new_rmsd_df], ignore_index=True)
    # save the rmsd_df to a csv file in the dock_ligand_dir with
    # the name csv_file_name and the columns names as the header
    rmsd_df.to_csv(os.path.join(dock_ligand_dir, csv_file_name), index=False)


def calc_rmsd_obabel(
    dock_ligand_dir: str,
    reference_dir: str,
    csv_file_name: str,
    number: str | int = 3,
):
    """
    Calculate the rmsd between the docked ligand and the reference ligand using rmsd_obabel
    Parameters
    ----------
    dock_ligand_dir : str
        The directory containing the docked ligand files
    reference_dir : str
        The directory containing the reference ligand files
    number : str or int
        The number of ligand to calculate the rmsd for, by default 3
    csv_file_name : str
        The name of the csv file to save the rmsd values
    """
    done = []
    # check if the dock_ligand_dir exists
    if not os.path.exists(dock_ligand_dir):
        raise FileNotFoundError(f"{dock_ligand_dir} does not exist")
    # check if the reference_dir exists
    if not os.path.exists(reference_dir):
        raise FileNotFoundError(f"{reference_dir} does not exist")
    # check if the csv_file_name exists in the dock_ligand_dir
    if os.path.exists(os.path.join(dock_ligand_dir, csv_file_name)):
        logger.info(f"{csv_file_name} already exists in {dock_ligand_dir}")
        # read the csv file using pandas
        rmsd_df = pd.read_csv(os.path.join(dock_ligand_dir, csv_file_name))
        # get the 2nd column of the csv file
        column = rmsd_df.iloc[:, 1]
        done.extend(list(set(column)))
    else:
        logger.info(f"{csv_file_name} does not exist in {dock_ligand_dir}")
        # make the csv file
        with open(
            os.path.join(dock_ligand_dir, csv_file_name), "w", encoding="utf-8"
        ) as file:
            file.write("file,RMSD\n")
    logger.info(done)
    logger.info(f"Found {len(done)} PDB IDs already done")
    # make a pandas dataframe to save the rmsd values
    rmsd_df = pd.DataFrame(
        columns=[
            "file",
            "RMSD",
        ]
    )
    # get all the docked ligand files in the dock_ligand_dir
    docked_ligands = glob.glob(os.path.join(dock_ligand_dir, "*.maegz"))
    # loop over the docked ligand files
    for docked_ligand in docked_ligands:
        logger.info(f"Working on {docked_ligand}")
        # get the pdb_id from the file name
        pdb_id = os.path.basename(docked_ligand).split("_")[0]
        if pdb_id.upper() in done:
            logger.info(f"{pdb_id} already done")
            continue
        # get the reference ligand file from the reference_dir
        # check if the reference ligand file exists and if not continue
        if not os.path.exists(
            os.path.join(reference_dir, pdb_id + "_protein_prepared.mae")
        ):
            logger.warning(f"{pdb_id} not found in the reference directory")
            continue
        full_reference_ligand = glob.glob(
            os.path.join(reference_dir, pdb_id + "*_protein_prepared.mae")
        )[0]
        # run manipulate to get the ligand file
        docked_ligand = DockedLigand(docked_ligand)
        context = GlideContext.get_current()
        context.write_dir = dock_ligand_dir
        context.set_current(context)
        if not os.path.exists(
            os.path.join(context.write_dir, docked_ligand.file_id + "_split_lig.pdb")
        ):
            logger.info(f"Splitting {docked_ligand.file_name}")
            manipulate_complexes(
                docked_ligand.file_path,
                context=context,
                mode="merge",
                outfile_name=docked_ligand.file_id + "_merge.maegz",
            )
            manipulate_complexes(
                os.path.join(context.write_dir, docked_ligand.file_id + "_merge.maegz"),
                context=context,
                mode="split_ligand",
                outfile_name=docked_ligand.file_id + "_split_lig.maegz",
            )
            convert_to_pdb(
                os.path.join(
                    context.write_dir, docked_ligand.file_id + "_split_lig.maegz"
                ),
                context=context,
            )
        logger.debug(f"Reference ligand: {full_reference_ligand}")
        reference_ligand = Ligand(full_reference_ligand)
        logger.debug(reference_ligand.file_id)
        if not os.path.exists(
            os.path.join(context.write_dir, reference_ligand.file_id + "_lig.pdb")
        ):
            logger.info(f"Splitting {reference_ligand.file_name}")
            # convert mae to pdb using convert_to_pdb
            convert_to_pdb(full_reference_ligand, context=context)
            protein, ligand = separate_protein_ligand(
                os.path.join(context.write_dir, reference_ligand.file_id + ".pdb"),
                remove_water=False,
            )
            print(ligand)
            write_mda_universe(
                ligand,
                os.path.join(context.write_dir, reference_ligand.file_id + "_lig.pdb"),
            )
            """
            manipulate_complexes(
                full_reference_ligand,
                context=context.get_current(),
                mode="merge",
                outfile_name=reference_ligand.file_id+"_merge.maegz",
            )
            manipulate_complexes(
                full_reference_ligand,
                context=context.get_current(),
                mode="split_ligand",
                outfile_name=reference_ligand.file_id+"_split_lig.maegz",
            )
            convert_to_pdb(os.path.join(context.write_dir,reference_ligand.file_id+ "_split_lig.maegz")
            ,context=context)
            """
        dock_lig_pdb = DockedLigand(
            os.path.join(context.write_dir, docked_ligand.file_id + "_split_lig.pdb")
        )
        ref_lig_pdb = DockedLigand(
            os.path.join(context.write_dir, reference_ligand.file_id + ".pdb")
        )
        rmsd = RMSD(dock_lig_pdb, ref_lig_pdb, context=context)
        logger.info(f"Calculating RMSD for {docked_ligand.file_name}")
        res = rmsd.rmsd_obabel(save=False)
        values = obabel_result_parser(res)
        # load the values dictionary to a dataframe
        current_rmsd_df = pd.DataFrame(values.items(), columns=["Index", "Value"])
        print(current_rmsd_df)
        # Set the index explicitly
        current_rmsd_df.set_index("file", inplace=True)
        # rename the index to the file name
        current_rmsd_df.index = [f"{docked_ligand.file_id}_{i}" for i in range(1, 4)]
        # load the rmsd.csv file using pandas and append the current_rmsd_df to it excluding the header
        rmsd_df = pd.concat([rmsd_df, current_rmsd_df])
        # clean up the directory
        for file in (
            glob.glob(os.path.join(context.write_dir, "*_split_lig.pdb"))
            and glob.glob(os.path.join(context.write_dir, "*_split_lig.maegz"))
            and glob.glob(os.path.join(context.write_dir, "*_merge.maegz"))
        ):
            os.remove(file)
        # save the rmsd_df to a csv file in the dock_ligand_dir with the name csv_file_name
        # and the columns names as the header
        rmsd_df.to_csv(os.path.join(dock_ligand_dir, csv_file_name), index=False)


def calc_rmsd_spyrmsd(
    dock_ligand_dir: str,
    reference_dir: str,
    csv_file_name: str,
    number: str | int = 3,
):
    """
    Calculate the rmsd between the docked ligand and the reference ligand using rmsd_obabel
    Parameters
    ----------
    dock_ligand_dir : str
        The directory containing the docked ligand files
    reference_dir : str
        The directory containing the reference ligand files
    number : str |int
        The number of ligand to calculate the rmsd for, by default 3
    csv_file_name : str
        The Path of the csv file to save the rmsd values,
        it will be created in the dock_ligand_dir if the path not exists
    """
    done = []
    failed = []
    count = 0
    context = GlideContext.get_current()
    context.write_dir = dock_ligand_dir
    context.set_current(context)
    # check if the dock_ligand_dir exists
    if not os.path.exists(dock_ligand_dir):
        raise FileNotFoundError(f"{dock_ligand_dir} does not exist")
    # check if the reference_dir exists
    if not os.path.exists(reference_dir):
        raise FileNotFoundError(f"{reference_dir} does not exist")
    # check if the csv_file_name is a path
    if "/" not in csv_file_name:
        csv_file_name = os.path.join(dock_ligand_dir, csv_file_name)
    # check if failed.txt exists in the dock_ligand_dir
    if os.path.exists(os.path.join(dock_ligand_dir, "failed.txt")):
        with open(os.path.join(dock_ligand_dir, "failed.txt"), "r") as file:
            failed = file.read().splitlines()
        done.extend([x.split("_")[0] for x in failed])
    # check if the csv_file_name exists in the dock_ligand_dir
    if os.path.exists(csv_file_name):
        logger.info(
            f"{csv_file_name} already exists in {os.path.dirname(csv_file_name)}"
        )
        rmsd_df = pd.read_csv(csv_file_name)
        column = rmsd_df.iloc[:, 0]
        done.extend(list(set([i.split("_")[0] for i in list(set(column))])))
    else:
        logger.info(
            f"{csv_file_name} does not exist in {os.path.dirname(csv_file_name)}"
        )
        with open(os.path.join(csv_file_name), "w") as file:
            file.write("file,RMSD\n")
        rmsd_df = pd.read_csv(os.path.join(csv_file_name))
    logger.info(f"Found {len(done)} PDB IDs already done")
    # get all the docked ligand files in the dock_ligand_dir
    docked_ligands = glob.glob(os.path.join(context.write_dir, "*_sorted.maegz"))
    # loop over the docked ligand files
    for docked_ligand in docked_ligands:
        flag = 0
        # get the pdb_id from the file name
        pdb_id = os.path.basename(docked_ligand).split("_")[0]
        if pdb_id in done:
            logger.info(f"{pdb_id} already done")
            continue
        if not os.path.exists(
            os.path.join(reference_dir, pdb_id + "_protein_prepared.mae")
        ):
            logger.error(f"Prepared {pdb_id} file not found in the reference directory")
            raise FileNotFoundError(
                f"Prepared {pdb_id} file not found in the reference directory"
            )
        full_reference_ligand = glob.glob(
            os.path.join(reference_dir, pdb_id + "_protein_prepared.mae")
        )[0]
        logger.debug(f"Reference ligand: {full_reference_ligand}")
        reference_ligand = Ligand(full_reference_ligand)
        logger.debug(f"The reference file id: {reference_ligand.file_id}")
        docked_ligand = DockedLigand(docked_ligand)
        # check if the docked_ligand file_id had _merge or split_lig in it and if so skip it
        if "_merge" in docked_ligand.file_id or "_split_lig" in docked_ligand.file_id:
            continue
        logger.info(f"Working on {docked_ligand.file_path}")
        # run manipulate to get the ligand file
        if not os.path.exists(
            os.path.join(context.write_dir, docked_ligand.file_id + "_merge.maegz")
        ):
            logger.debug(
                f"Complex generation for {os.path.join(context.write_dir, docked_ligand.file_id)}"
            )
            manipulate_complexes(
                docked_ligand.file_path,
                context=context,
                mode="merge",
                outfile_name=docked_ligand.file_id + "_merge.maegz",
            )
        if not os.path.exists(
            os.path.join(context.write_dir, docked_ligand.file_id + "_split_lig.pdb")
        ):
            logger.info(f"Splitting {docked_ligand.file_name} to get the ligand")

            if not os.path.exists(
                os.path.join(
                    context.write_dir, docked_ligand.file_id + "_split_lig.maegz"
                )
            ):
                manipulate_complexes(
                    os.path.join(
                        context.write_dir, docked_ligand.file_id + "_merge.maegz"
                    ),
                    context=context,
                    mode="split_ligand",
                    outfile_name=docked_ligand.file_id + "_split_lig.maegz",
                )
            if convert_to_pdb(
                os.path.join(
                    context.write_dir, docked_ligand.file_id + "_split_lig.maegz"
                ),
                context=context,
            ):
                flag = 1
        if not os.path.exists(
            os.path.join(context.write_dir, reference_ligand.file_id + "_lig.pdb")
        ):
            logger.info(f"Splitting {reference_ligand.file_name}")
            # convert mae to pdb using convert_to_pdb
            if not os.path.exists(
                os.path.join(context.write_dir, reference_ligand.file_id + ".pdb")
            ):
                convert_to_pdb(full_reference_ligand, context=context)
            if flag == 1:
                protein, ligand = separate_protein_ligand(
                    os.path.join(context.write_dir, reference_ligand.file_id + ".pdb"),
                    remove_water=True,
                    reference=os.path.join(
                        context.write_dir, docked_ligand.file_id + "_split_lig-1.pdb"
                    ),
                )
            else:
                logger.debug(
                    f"PDB for reference ligand: {os.path.join(context.write_dir, reference_ligand.file_id + '.pdb')}"
                )
                protein, ligand = separate_protein_ligand(
                    os.path.join(context.write_dir, reference_ligand.file_id + ".pdb"),
                    remove_water=True,
                    reference=os.path.join(
                        context.write_dir, docked_ligand.file_id + "_split_lig.pdb"
                    ),
                )
            write_mda_universe(
                ligand,
                os.path.join(context.write_dir, reference_ligand.file_id + "_lig.pdb"),
            )
        ref_lig_pdb = DockedLigand(
            os.path.join(context.write_dir, reference_ligand.file_id + "_lig.pdb")
        )
        try:
            if flag == 1:
                values = {}
                for i in range(1, number + 1):
                    dock_lig_pdb = DockedLigand(
                        os.path.join(
                            context.write_dir,
                            docked_ligand.file_id + f"_split_lig-{str(i)}.pdb",
                        )
                    )
                    rmsd = RMSD(dock_lig_pdb, ref_lig_pdb, context=context)
                    logger.info(f"Calculating RMSD for {docked_ligand.file_name}")
                    res = rmsd.symmetry_rmsd()
                    values[i] = str(list(obabel_result_parser(res).values()))
            else:
                dock_lig_pdb = DockedLigand(
                    os.path.join(
                        context.write_dir, docked_ligand.file_id + "_split_lig.pdb"
                    )
                )
                rmsd = RMSD(dock_lig_pdb, ref_lig_pdb, context=context)
                logger.info(f"Calculating RMSD for {docked_ligand.file_name}")
                res = rmsd.symmetry_rmsd()
                values = obabel_result_parser(res)
        except Exception as e:
            logger.error(
                f"Could not calculate RMSD for {docked_ligand.file_name} because {str(e)}"
            )
            failed.append(docked_ligand.file_id)
            continue

        print(values)
        # load the values dictionary to a dataframe
        current_rmsd_df = pd.DataFrame(values.items(), columns=["file", "RMSD"])
        # Set the index explicitly
        current_rmsd_df.set_index("file", inplace=True)
        # rename the index to the file name
        current_rmsd_df.index = [
            f"{docked_ligand.file_id}_{i}" for i in range(len(values))
        ]
        # make the index the first column of the dataframe and name it file
        current_rmsd_df.reset_index(inplace=True)
        current_rmsd_df.rename(columns={"index": "file"}, inplace=True)
        # load the rmsd.csv file using pandas and append the current_rmsd_df to it excluding the header
        rmsd_df = pd.concat([rmsd_df, current_rmsd_df])
        # clean up the directory
        # load the values dictionary to a dataframe
        current_rmsd_df = pd.DataFrame(values.items(), columns=["Index", "Value"])
        print(current_rmsd_df)
        # clean up the directory
        for file in (
            glob.glob(os.path.join(context.write_dir, "*.pdb"))
            + glob.glob(os.path.join(context.write_dir, "*_split_lig.maegz"))
            + glob.glob(os.path.join(context.write_dir, "*_merge.maegz"))
        ):
            os.remove(file)
        # save the rmsd_df to a csv file in the dock_ligand_dir with the name csv_file_name and the same columns
        rmsd_df.to_csv(csv_file_name, index=False)
        # save the failed list to a file
        with open(os.path.join(dock_ligand_dir, "failed.txt"), "w") as file:
            file.write("\n".join(failed))


if __name__ == "__main__":
    CLUSTER_FILE = "../new_clusters_cluster_parsed.csv"
    reps = find_cluster_reps(CLUSTER_FILE)
    logger.info(f"Found {len(reps)} representatives")
    RAW_INPUT_DIR = "../../representatives"
    PREPPED_DIR = "./prepped"
    DOCKED_DIR = "./docked"
    COMBIND_DIR = "./combind"
    get_coordinates(reps, RAW_INPUT_DIR, limit=6000)
    # prep_structure("./trial", "./prepped", ["6n2w.pdb", "1xoi.pdb"])
    # dock_ligand("./prepped", "./trial")
    # calc_rmsd_spyrmsd("./trial", "./prepped", "rmsd.csv")
