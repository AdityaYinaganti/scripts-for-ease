import argparse
import pandas as pd
import random
from rdkit import Chem
from rdkit.Chem import Descriptors

def get_diverse_fragment():
    fragments = [
        'F', 'Cl', 'Br', 'I', 'C', 'CC', 'C(C)C', 'C#N',
        'OC', 'OCC', 'OC(C)C', 'S(=O)(=O)C', 'S(=O)(=O)N', 
        'N', 'N(C)C', 'NC(=O)C', 'C(=O)NC', 'C(=O)N(C)C',
        'C(=O)O', 'C(=O)C', 'CC(=O)C',
        'c1ccccc1', 'c1ccncc1', 'c1cc(F)ccc1', 'C1CC1', 'C1CCCC1', 'C1CCOCC1'
    ]
    return random.choice(fragments)

def is_drug_like(mol):
    if mol is None: return False
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        return (150 < mw < 500 and -1 < logp < 4.5 and hbd <= 5 and hba <= 10)
    except:
        return False

def generate_split_isomers(n_rows, e_file, z_file):
    e_data = []
    z_data = []
    seen_smiles = set()
    
    while len(e_data) < n_rows:
        r = [get_diverse_fragment() for _ in range(4)]
        if r[0] == r[1] or r[2] == r[3]:
            continue
            
        # Geometric construction
        s_z = f"{r[0]}/C({r[1]})=C(\\{r[2]}){r[3]}"
        s_e = f"{r[0]}/C({r[1]})=C(/{r[2]}){r[3]}"
        
        m_z = Chem.MolFromSmiles(s_z)
        m_e = Chem.MolFromSmiles(s_e)
        
        if m_z and m_e and is_drug_like(m_z):
            # Normalize SMILES
            smiles_z = Chem.MolToSmiles(m_z, isomericSmiles=True)
            smiles_e = Chem.MolToSmiles(m_e, isomericSmiles=True)
            
            if smiles_z not in seen_smiles:
                z_data.append({"SMILES": smiles_z})
                e_data.append({"SMILES": smiles_e})
                seen_smiles.add(smiles_z)

    # Save to separate files
    pd.DataFrame(e_data).to_csv(e_file, index=False)
    pd.DataFrame(z_data).to_csv(z_file, index=False)
    print(f"Success! Generated {n_rows} rows in {e_file} and {z_file}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate separate E and Z isomer CSVs.")
    parser.add_argument("-r", "--rows", type=int, required=True, help="Number of rows")
    parser.add_argument("-e", "--e_out", type=str, default="e_isomers.csv", help="E isomers filename")
    parser.add_argument("-z", "--z_out", type=str, default="z_isomers.csv", help="Z isomers filename")
    
    args = parser.parse_args()
    generate_split_isomers(args.rows, args.e_out, args.z_out)
