#!/usr/bin/env python3
## Made by Joonhyeok Choi
## For converting ZINC SMILES to MOL, SDF, PDB, or PDBQT format
## Using ETKDGv3 version


import os
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse
from multiprocessing import Pool, cpu_count
from oddt.toolkits.extras.rdkit import MolToPDBQTBlock
from random import *

def smiles_to_format(args):
    smiles, output_dir, zinc_id, output_format = args
    print(f"Processing {zinc_id}")
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useRandomCoords = True
    params.useMacrocyclesTorsions = True
    params.useMachineLearning = True
    params.numThreads = 4
    params.randomSeed = randrange(10000)

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, 10, params)
    #AllChem.EmbedMolecule(mol, params)
    #AllChem.UFFOptimizeMolecule(mol)

    # Save the molecule in format
    if output_format == 'mol':
        output_file = os.path.join(output_dir, f"{zinc_id}.mol")
        Chem.MolToMolFile(mol, output_file)
    elif output_format == 'sdf':
        output_file = os.path.join(output_dir, f"{zinc_id}.sdf")
        writer = Chem.SDWriter(output_file)
        writer.write(mol)
        writer.close()
    elif output_format == 'pdb':
        output_file = os.path.join(output_dir, f"{zinc_id}.pdb")
        Chem.MolToPDBFile(mol, output_file)
    elif output_format == 'pdbqt':
        output_file = os.path.join(output_dir, f"{zinc_id}.pdbqt")
        with open(output_file, 'w') as f:
            f.write(f"REMARK  Converted from SMILES to PDBQT format\n")
            f.write(MolToPDBQTBlock(mol, computeCharges=True))
    else:
        print("Invalid output format. Please choose 'mol', 'sdf', 'pdb', or 'pdbqt'.")
        
def process_data_file(data_file, output_dir, output_format):
    zinc_dict = {}

    with open(data_file, "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader, None)
        for row in reader:
            zinc_id, smiles = row
            zinc_dict[zinc_id] = smiles

    with Pool(cpu_count()) as pool:
        args_list = [(smiles, output_dir, zinc_id, output_format) for zinc_id, smiles in zinc_dict.items()]

        pool.map(smiles_to_format, args_list)
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert ZINC SMILES to MOL, SDF, PDB, or PDBQT format.")
    parser.add_argument("data_file", help="CSV file containing ZINC IDs and SMILES.")
    parser.add_argument("-o", "--output_dir", default="output", help="Output directory for files.")
    parser.add_argument("-f", "--format", choices=['mol', 'sdf', 'pdb', 'pdbqt'], default='mol', help="Output format (mol, sdf, pdb, pdbqt).")
    
    args = parser.parse_args()

    data_file = args.data_file
    output_dir = args.output_dir
    output_format = args.format

    # Output directory check and creation
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else: #os.path.exists(output_dir)
        print(f"Output directory '{output_dir}' already exists. Please choose a different directory.")
        exit(0)
    
    
    process_data_file(data_file, output_dir, output_format)
