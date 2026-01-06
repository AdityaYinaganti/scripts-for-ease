import pandas as pd
import argparse
import random
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_smiles(max_atoms):
    noble_gases = {2, 10, 18, 36, 54, 86, 118}
    all_elements = [i for i in range(1, 119) if i not in noble_gases]
    
    start_atoms = ['C', 'N', 'O', 'P', 'S']
    mol = Chem.MolFromSmiles(random.choice(start_atoms))
    mol.UpdatePropertyCache()
    emol = Chem.EditableMol(mol)
    
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
        new_atomic_num = random.choice(all_elements)
        
        try:
            new_idx = emol.AddAtom(Chem.Atom(new_atomic_num))
            emol.AddBond(parent_idx, new_idx, Chem.rdchem.BondType.SINGLE)
        except:
            continue
        
    final_mol = emol.GetMol()
    try:
        Chem.SanitizeMol(final_mol)
        return Chem.MolToSmiles(final_mol)
    except:
        return generate_smiles(max_atoms)

def main():
    parser = argparse.ArgumentParser(description="Generate random SMILES using the full Periodic Table.")
    parser.add_argument('--rows', '-r', type=int, required=True, help="Number of rows to generate.")
    parser.add_argument('--output', '-o', type=str, required=True, help="Output CSV filename.")
    parser.add_argument('--size', '-s', type=int, default=10, help="Max atoms per molecule.")

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
