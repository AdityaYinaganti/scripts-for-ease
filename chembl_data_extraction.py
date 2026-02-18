import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from tqdm import tqdm

def get_all_helm_biologics():
    output_file = "chembl_full_helm_dataset.csv"
    print("Connecting to ChEMBL...")
    molecule = new_client.molecule
    biotherapeutic = new_client.biotherapeutic

    # 1. Fetch all biotherapeutics that have a HELM notation
    print("Searching for biotherapeutics with HELM notation...")
    query = biotherapeutic.filter(helm_notation__isnull=False).only(
        ['molecule_chembl_id', 'helm_notation']
    )

    data = []
    try:
        for entry in tqdm(query):
            chembl_id = entry.get('molecule_chembl_id')
            helm = entry.get('helm_notation')
            
            # 2. Getting the SMILES for this specific ID
            mol_data = molecule.get(chembl_id)
            structures = mol_data.get('molecule_structures')
            smiles = structures.get('canonical_smiles') if structures else None
            
            # 3. Calculating the NUM_RINGS
            num_rings = 0
            if smiles:
                try:
                    rdkit_mol = Chem.MolFromSmiles(smiles)
                    if rdkit_mol:
                        num_rings = rdkit_mol.GetRingInfo().NumRings()
                except:
                    num_rings = "Check SMILES"
            
            data.append({
                "SMILES": smiles,
                "HELM": helm,
                "CHEMBL_ID": chembl_id,
                "NUM_RINGS": num_rings
            })

            # The API connection can get lost sometimes, therefore autosaving the file after every 20 seconds. 
            if len(data) % 20 == 0:
                pd.DataFrame(data).to_csv(output_file, index=False)

    except Exception as e:
        print(f"\nError: {e}")

    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"\nComplete! Found {len(df)} records. Saved to {output_file}")

if __name__ == "__main__":
    get_all_helm_biologics()
