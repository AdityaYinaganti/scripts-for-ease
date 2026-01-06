import pandas as pd
import argparse
import random
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_smiles(max_atoms):
    start_atoms = ['C', 'N', 'O']
    mol = Chem.MolFromSmiles(random.choice(start_atoms))
    

    mol.UpdatePropertyCache()
    emol = Chem.EditableMol(mol)
    
    elements = [6, 6, 6, 7, 8, 9, 17, 35] # C, C, C, N, O, F, Cl, Br
    
    num_atoms_to_add = random.randint(3, max_atoms)
    
    for _ in range(num_atoms_to_add):
        current_mol = emol.GetMol()
        try:
            current_mol.UpdatePropertyCache()
        except:
            break
            
        valid_indices = [a.GetIdx() for a in current_mol.GetAtoms() if a.GetNumImplicitHs() > 0]
        
        if not valid_indices:
            break
            
        parent_idx = random.choice(valid_indices)
        new_element = random.choice(elements)
        
        new_idx = emol.AddAtom(Chem.Atom(new_element))
        emol.AddBond(parent_idx, new_idx, Chem.rdchem.BondType.SINGLE)
        
    final_mol = emol.GetMol()
    try:
        Chem.SanitizeMol(final_mol)
        return Chem.MolToSmiles(final_mol)
    except:
        return generate_smiles(max_atoms)

def main():
    parser = argparse.ArgumentParser(description="Generate a CSV of unique, scientifically accurate random SMILES.")
    parser.add_argument('--rows', '-r', type=int, required=True, help="Number of rows (entities) to generate.")
    parser.add_argument('--output', '-o', type=str, required=True, help="Name of the output CSV file.")
    parser.add_argument('--size', '-s', type=int, default=12, help="Max number of atoms per molecule.")

    args = parser.parse_args()

    unique_smiles = set()
    data = []
    
    counter = 1
    while len(data) < args.rows:
        smi = generate_smiles(args.size)
        if smi not in unique_smiles:
            unique_smiles.add(smi)
            data.append({'ID': counter, 'SMILES': smi})
            counter += 1

    df = pd.DataFrame(data)
    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()