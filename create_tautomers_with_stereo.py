import random
import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

def get_chiral_tautomer_fragments():
    """
    This function will return the fragments which has both tautomerism
    and stereochemistry. 
    """
    return [
        'CC(F)C(=O)C',    # Chiral alpha-carbon
        'CC(Cl)C=O',      # Chiral alpha-carbon
        'NC(C)C(=O)O',    # Amino acid style
        'C1[C@H](C)C(=O)CC1', # Cyclic chiral ketone
        'C[C@@H](O)C(=N)N',   # Chiral imine/amine
        'CC(N)C(=S)N',    # Thioamide
        'OC[C@H](N)C=O'   # Hydroxy-aldehyde
    ]

def generate_stereo_tautomers(n_parents, output_file):
    enumerator = rdMolStandardize.TautomerEnumerator()
    # Adding this function to not create too many tautomers (Optional)
    enumerator.SetMaxTautomers(50) 
    
    data = []
    frag_pool = get_chiral_tautomer_fragments()
    found_parents = 0
    
    while found_parents < n_parents:
        # Joining multiple fragments to create a more complex molecules. 
        f1 = random.choice(frag_pool)
        f2 = random.choice(frag_pool)
        combined_smiles = f"{f1}.{f2}"
        mol = Chem.MolFromSmiles(combined_smiles)
        
        if not mol: continue

        all_tautomers = enumerator.Enumerate(mol)
        
        # This will keep tautomers which has atleast one stereoisomers. 
        stereo_tautomers = [
            t for t in all_tautomers 
            if len(Chem.FindMolChiralCenters(t, includeUnassigned=True)) > 0
        ]

        if not stereo_tautomers:
            continue

        for idx, t_mol in enumerate(stereo_tautomers):
            # This will make sure that stereochemistry is intact. 
            Chem.AssignStereochemistry(t_mol, force=True, cleanIt=True)
            
            data.append({
                "Parent_ID": f"MOL-{found_parents + 1}",
                "Tautomer_ID": f"{found_parents + 1}_{idx + 1}",
                "SMILES": Chem.MolToSmiles(t_mol, isomericSmiles=True),
                "Chiral_Centers": len(Chem.FindMolChiralCenters(t_mol, includeUnassigned=True)),
                "MW": round(Descriptors.MolWt(t_mol), 2),
                "LogP": round(Descriptors.MolLogP(t_mol), 2)
            })
        
        found_parents += 1

    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"Saved {len(df)} stereo-retaining tautomers to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--number", type=int, default=10)
    parser.add_argument("-o", "--output", type=str, default="stereo_tautomers.csv")
    args = parser.parse_args()
    
    generate_stereo_tautomers(args.number, args.output)