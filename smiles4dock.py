#!/usr/bin/env python3
## Made by Joonhyeok Choi
## For converting ZINC SMILES to MOL, SDF, PDB, or PDBQT format

import os
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse
from multiprocessing import Pool, cpu_count
from oddt.toolkits.extras.rdkit import MolToPDBQTBlock

def smiles_to_mol(args):
    smiles, output_dir, zinc_id = args
    print(f"Processing {zinc_id}")
    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
    #AllChem.MMFFOptimizeMolecule(mol_h)

    # Save the molecule in mol format
    output_file = os.path.join(output_dir, f"{zinc_id}.mol")
    Chem.MolToMolFile(mol_h, output_file)

def smiles_to_sdf(args):
    smiles, output_dir, zinc_id = args
    print(f"Processing {zinc_id}")
    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
    #AllChem.MMFFOptimizeMolecule(mol_h)

    # Save the molecule in sdf format
    output_file = os.path.join(output_dir, f"{zinc_id}.sdf")
    writer = Chem.SDWriter(output_file)
    writer.write(mol_h)
    writer.close()

def smiles_to_pdb(args):
    smiles, output_dir, zinc_id = args
    print(f"Processing {zinc_id}")
    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
    #AllChem.MMFFOptimizeMolecule(mol_h)

    # Save the molecule in pdb format
    output_file = os.path.join(output_dir, f"{zinc_id}.pdb")
    Chem.MolToPDBFile(mol_h, output_file)

def smiles_to_pdbqt(args):
    smiles, output_dir, zinc_id = args
    print(f"Processing {zinc_id}")
    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
    #AllChem.MMFFOptimizeMolecule(mol_h)

    # Save the molecule in pdbqt format
    output_file = os.path.join(output_dir, f"{zinc_id}.pdbqt")
    with open(output_file, 'w') as f:
        f.write(f"REMARK  Converted from SMILES to PDBQT format\n")
        f.write(MolToPDBQTBlock(mol_h, computeCharges=True))

def process_data_file(data_file, output_dir, output_format):
    zinc_dict = {}

    with open(data_file, "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader, None)
        for row in reader:
            zinc_id, smiles = row
            zinc_dict[zinc_id] = smiles

    with Pool(cpu_count()) as pool:
        args_list = [(smiles, output_dir, zinc_id) for zinc_id, smiles in zinc_dict.items()]

        if output_format == 'mol':
            pool.map(smiles_to_mol, args_list)
        elif output_format == 'sdf':
            pool.map(smiles_to_sdf, args_list)
        elif output_format == 'pdb':
            pool.map(smiles_to_pdb, args_list)
        elif output_format == 'pdbqt':
            pool.map(smiles_to_pdbqt, args_list)
        else:
            print("Invalid output format. Please choose 'mol', 'sdf', 'pdb', or 'pdbqt'.")

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
    
    
    process_data_file(data_file, output_dir, output_format)
