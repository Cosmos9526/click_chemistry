from rdkit import Chem
from rdkit.Chem import AllChem
import random
import json
from typing import List, Dict, Tuple
from rdkit.Chem import rdMolDescriptors

# --- Main parameters (DBCO-focused, cleavable) ---
unique_smiles: Dict[str, str] = {
    # Enzyme-cleavable peptides (NH2 left, COOH right)
    'Val-Cit': 'CC(C)[C@@H](C(=O)N[C@@H](CCCNC(=O)N)C(=O)O)N',
    'Phe-Lys': 'c1ccc(cc1)C[C@@H](C(=O)N[C@@H](CCCCN)C(=O)O)N',
    'Ala-Ala': 'C[C@@H](C(=O)N[C@@H](C)C(=O)O)N',
    'Gly-Gly-Phe-Gly': 'NCC(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)NCC(=O)O',

    # Self-immolative spacers (NH2 for attachment)
    'PABC': 'NC1=CC=C(COC(=O)N)C=C1',
    'PAB-OH': 'NC1=CC=C(CO)C=C1',
    'PABC-PNP': 'NC1=CC=C(COC(=O)OC2=CC=C([N+](=O)[O-])C=C2)C=C1',

    # Disulfide-cleavable (bifunctional)
    'Disulfide-SS': 'NCCSSCC(=O)O',

    # pH/Hydrazone-cleavable
    'Hydrazone-linker': 'NNC(=O)CCO',
    'AcBut-hydrazone': 'CC(=O)OC1=CC=C(CCCC(=O)NN)C=C1',

    # Glycosidase-cleavable (with COOH)
    'β-glucuronide': 'O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1C(=O)O',
    'Me-triacetyl-β-glucuronate': 'CC(=O)O[C@H]1[C@@H](OC(=O)C)[C@H](OC(=O)C)[C@@H](OC(=O)C)O[C@@H]1C(=O)O',

    # DBCO-acid for amide (from standard derivative)
    'DBCO-acid': 'O=C(CCC(O)=O)N1CC2=C(C=CC=C2)C#CC3=C1C=CC=C3',

    # PEG spacers (NH2-PEG-COOH)
    'PEG3': 'NCCOCCC(=O)O',  # Matches example
    'PEG4': 'NCCOCCOCCOCC(=O)O'
}

# Component categories (DBCO-only for click)
component_categories = {
    'click': ['DBCO-acid'],
    'spacer': ['PEG3', 'PEG4'],
    'cleavable_enzyme': ['Val-Cit', 'Phe-Lys', 'Ala-Ala', 'Gly-Gly-Phe-Gly'],
    'cleavable_disulfide': ['Disulfide-SS'],
    'cleavable_ph': ['Hydrazone-linker', 'AcBut-hydrazone'],
    'cleavable_glycosidase': ['β-glucuronide', 'Me-triacetyl-β-glucuronate'],
    'self_immolative': ['PABC', 'PAB-OH', 'PABC-PNP']
}

# Cleavage type mapping
cleavage_mapping = {
    'Val-Cit': 'enzyme',
    'Phe-Lys': 'enzyme',
    'Ala-Ala': 'enzyme',
    'Gly-Gly-Phe-Gly': 'enzyme',
    'Disulfide-SS': 'disulfide',
    'Hydrazone-linker': 'pH',
    'AcBut-hydrazone': 'pH',
    'β-glucuronide': 'glycosidase',
    'Me-triacetyl-β-glucuronate': 'glycosidase'
}

# SMIRKS reactions (specific for building)
reaction_smarts: List[str] = [
    '[O:1][C:2](=[O:3])[Cl].[NH2:4]>>[O:1][C:2](=[O:3])[NH:4].[Cl]',  # Carbamate (kept for possible use)
    '[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]',  # Amide (COOH + NH2) - now primary for DBCO
    '[C:1](=[O:2])[OH].[OH:3]>>[C:1](=[O:2])[O:3].O',  # Ester (COOH + OH, specific [OH])
    '[S:1].[S:2]>>[S:1][S:2]',  # Disulfide
    '[C:1]=O.[N:3][N:4]>>[C:1]=[N:3][N:4]',  # Hydrazone
    '[O:1]C(=O)N[C:2]>>[O:1][C:2]=O',  # Carbamate self-immolative
    '[C:1](O)O>>[C:1]=O.O'  # Glucuronide
]

reactions = [AllChem.ReactionFromSmarts(s) for s in reaction_smarts]

def generate_random_pattern() -> Tuple[List[str], str]:
    """Generate random pattern: DBCO + spacer + cleavable + opt spacer + self-immolative"""
    cleavage_types = ['enzyme', 'disulfide', 'pH', 'glycosidase']
    cleavage_type = random.choice(cleavage_types)
    
    pattern = []
    
    # DBCO-acid + spacer (amide)
    pattern.append('DBCO-acid')
    pattern.append(random.choice(component_categories['spacer']))
    
    # Cleavable
    if cleavage_type == 'enzyme':
        pattern.append(random.choice(component_categories['cleavable_enzyme']))
    elif cleavage_type == 'disulfide':
        pattern.append(random.choice(component_categories['cleavable_disulfide']))
    elif cleavage_type == 'pH':
        pattern.append(random.choice(component_categories['cleavable_ph']))
    elif cleavage_type == 'glycosidase':
        pattern.append(random.choice(component_categories['cleavable_glycosidase']))
    
    # Optional spacer
    if random.random() > 0.5:
        pattern.append(random.choice(component_categories['spacer']))
    
    # Self-immolative
    pattern.append(random.choice(component_categories['self_immolative']))
    
    # Cleavage type
    cleavable_comp = next((comp for comp in pattern if comp in cleavage_mapping), None)
    final_cleavage = cleavage_mapping.get(cleavable_comp, 'enzyme')
    
    return pattern, final_cleavage

def generate_linker(components: List[str]) -> str:
    """Connect components using reactions (amide first for DBCO)"""
    try:
        mols = [Chem.MolFromSmiles(unique_smiles[comp]) for comp in components]
        if any(m is None for m in mols):
            return "Error: Invalid SMILES"

        if len(mols) < 2:
            return "Error: At least 2 components"

        current = mols[0]
        for next_mol in mols[1:]:
            connected = False
            for rxn in reactions:
                for order in [(current, next_mol), (next_mol, current)]:
                    products = rxn.RunReactants(order)
                    if products:
                        prod = products[0][0]
                        if Chem.SanitizeMol(prod) == 0:
                            temp_smiles = Chem.MolToSmiles(prod, isomericSmiles=True)
                            if '?' not in temp_smiles and Chem.MolFromSmiles(temp_smiles):
                                current = prod
                                connected = True
                                break
                if connected:
                    break
            if not connected:
                return "Error: Connection failed"
        
        smiles = Chem.MolToSmiles(current, isomericSmiles=True)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Error: Invalid SMILES"
        if rdMolDescriptors.CalcExactMolWt(mol) > 900:
            return "Error: MW too high"
        # Check for DBCO: exactly 1 triple bond
        triple_count = sum(1 for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.TRIPLE)
        if triple_count != 1:
            return "Error: No DBCO triple bond"
        return smiles
    except Exception as e:
        return f"Error: {str(e)}"

# Generate 1000 unique linkers
random.seed(42)
successful = []
unique_smiles_set = set()
target = 1000
attempts = 0
max_attempts = 50000

print(f"Starting generation of {target} unique DBCO linkers...")

while len(successful) < target and attempts < max_attempts:
    pattern, cleavage_type = generate_random_pattern()
    smiles = generate_linker(pattern)
    if smiles.startswith("Error"):
        attempts += 1
        continue
    
    if smiles not in unique_smiles_set:
        unique_smiles_set.add(smiles)
        name = '-'.join(pattern)
        successful.append({
            "name": name,
            "smiles": smiles,
            "cleavage_type": cleavage_type,
            "components": pattern
        })
        if len(successful) % 100 == 0:
            print(f"{len(successful)} linkers generated...")

    attempts += 1

print(f"✅ {len(successful)} unique linkers generated (attempts: {attempts}).")

# Save
with open('dbco_cleavable_linkers.json', 'w', encoding='utf-8') as f:
    json.dump(successful, f, ensure_ascii=False, indent=2)

print("Saved to dbco_cleavable_linkers.json!")

# Test with example pattern
test_pattern = ['DBCO-acid', 'PEG3', 'Val-Cit', 'PAB-OH']
test_smiles = generate_linker(test_pattern)
print("\nExample SMILES (DBCO-PEG3-ValCit-PABOH):", test_smiles)
if not test_smiles.startswith("Error"):
    test_mol = Chem.MolFromSmiles(test_smiles)
    print("Valid:", test_mol is not None)
    print("MW:", rdMolDescriptors.CalcExactMolWt(test_mol))
    print("Triple bonds:", sum(1 for b in test_mol.GetBonds() if b.GetBondType() == Chem.rdchem.BondType.TRIPLE))
