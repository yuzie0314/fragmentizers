#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import FragmentCatalog
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from typing import Any, Dict, List
import sqlite3
import pickle
import os
import glob
import utilities
import constants

"""
In this script, we will first generate fingerprints catalogue if you don't have your custom library.
"""


class FpCatalogue(object):
    def __init__(self, logger, sdf_folder_to_glob: str, outdir: str, threads: int):
        self.logger = logger
        self.sdf_folder_to_glob = sdf_folder_to_glob
        self.outdir = outdir
        self.threads = threads
        self.catalog: FragmentCatalog = None

    def build_sdf_catalogue(self):
        self.logger.info(f"building sdf catalogue inside {self.outdir}")
        sdf_list = glob.glob(
            os.path.join(self.sdf_folder_to_glob, "*.sdf")
        )  # TODO checking path!!
        suppl = Chem.SDMolSupplier()
        for sdf_file in sdf_list:
            suppl.SetData(sdf_file)

        params = FragmentCatalog.FragCatParams()
        params.maxMolSize = constants.maxMolSize
        params.preprocessMols = constants.preprocessMols
        params.useHs = constants.useHs
        params.keepFragIdsPickle = constants.keepFragIdsPickle
        params.threads = self.threads

        fcgen = FragmentCatalog.FragCatGenerator()
        self.catalog = fcgen.Build(suppl, params=params)
        FragmentCatalog.SaveCatalog(
            self.catalog, f"{self.outdir}/{self.sdf_folder_to_glob}_catalog.pkl"
        )
        self.reverse_fp_catalogue(suppl)

    def generate_molecule_fragments(self):
        """
        Generate fragments for each molecule in the given catalog and save them as separate SDF files in the specified directory.
        """
        out_folder: str = ""
        # Loop over each molecule in the catalog
        for i, mol in enumerate(self.catalog.GetMols()):
            # Get the molecule ID
            mol_id = mol.GetProp("_Name")
            # Generate fragments for the molecule
            fragments = FragmentCatalog.FragmentsFromMol(mol, self.catalog)
            out_folder = f"{self.outdir}/fragments/{mol_id}"
            os.makedirs(out_folder, exist_ok=True)
            # Save each fragment as a separate SDF file
            for j, frag in enumerate(fragments):
                frag_id = f"{mol_id}_frag{j}"
                writer = Chem.SDWriter(f"{out_folder}/{frag_id}.sdf")
                writer.write(frag)
                writer.close()

    def reverse_fp_catalogue(self, suppl):
        """
        Create reverse dictionary mapping molecule IDs to fragment IDs and SMILES strings
        """
        reverse_frag_cat: Dict[int, Dict[Any]] = dict()
        for frag_smiles, mol_ids in self.catalog.items():
            for mol in suppl:
                mol_id = mol.GetProp("_Name")
                mol_smiles = Chem.MolToSmiles(mol)
                if mol_id in mol_ids:
                    if mol_id not in reverse_frag_cat:
                        reverse_frag_cat[mol_id] = {
                            "smiles": mol_smiles,
                            "fragments": list(),
                        }
                    reverse_frag_cat[mol_id]["fragments"].append(frag_smiles)
        self.store_to_sqlite(reverse_frag_cat)

    def store_to_sqlite(self, reverse_frag_cat):
        # Connect to database and create new table
        conn = sqlite3.connect(f"{self.outdir}/mol_to_fragments.db")
        c = conn.cursor()
        # Create table to store molecule information
        c.execute(
            "CREATE TABLE IF NOT EXISTS molecules (id TEXT PRIMARY KEY, smiles TEXT)"
        )

        # Create table to store fragment information
        c.execute(
            "CREATE TABLE IF NOT EXISTS fragments (id INTEGER PRIMARY KEY AUTOINCREMENT, smiles TEXT, mol_id TEXT, FOREIGN KEY (mol_id) REFERENCES molecules(id))"
        )
        # Insert molecule information into database
        for mol_id, mol_info in reverse_frag_cat.items():
            c.execute(
                "INSERT INTO molecules (id, smiles) VALUES (?, ?)",
                (mol_id, mol_info["smiles"]),
            )

        # Get ID of last inserted row to use as foreign key for fragments table
        mol_rowid = c.lastrowid

        # Insert fragment information into database
        for frag_smiles in mol_info["fragments"]:
            c.execute(
                "INSERT INTO fragments (smiles, mol_id) VALUES (?, ?)",
                (frag_smiles, mol_id),
            )

        # Commit changes and close connection
        conn.commit()
        conn.close()


def main():
    sdf, sdf_folder, build_fingerprint_db, outdir, threads = utilities.arguments()
    logger = utilities.setup_logger("Logger")
    # Prepare outputs folder
    os.makedirs(f"{outdir}/catalogue", exist_ok=True)

    if build_fingerprint_db:
        db_building = FpCatalogue(logger, sdf_folder, f"{outdir}/catalogue", threads)
    else:
        pass


if __name__ == "__main__":
    main()
