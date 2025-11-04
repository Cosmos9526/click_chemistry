from typing import Dict, Any
from rdkit import Chem
from rdkit.Chem import AllChem
from config import AZIDE, SMART_REACTION, LINKER_CORE
from utils import validate_smiles

def build_linker_azide(core_smiles: str, azide_smiles: str) -> str:
    """
    Combine linker core (nanobody_tag with *) with azide for the click reaction.
    
    :param core_smiles: SMILES of the starred linker core (nanobody_tag)
    :param azide_smiles: SMILES of the azide
    :return: Combined SMILES for linker-azide
    :raises ValueError: If invalid SMILES
    """
    combined = core_smiles + '.' + azide_smiles
    combined = combined.replace('(*)','9')  # Replace placeholder to connect via 'ring' trick
    mol = Chem.MolFromSmiles(combined)
    if mol is None:
        raise ValueError("❌ Invalid linker azide SMILES: " + combined)
    final_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=False)
    print("🔗 Built Linker-Azide SMILES (short):", final_smiles[:50] + "...")
    return final_smiles

def run_click_reaction(drug_smiles: str, linker_core_smiles: str, azide_smiles: str, reaction_smarts: str) -> str:
    """
    Run the click reaction between the drug (best SMILES with alkyne from linker) and linker-azide (nanobody + azide).
    
    :param drug_smiles: SMILES of the drug (best from evaluation, with alkyne in starred_linker)
    :param linker_core_smiles: SMILES of the starred nanobody_tag (attachment for azide)
    :param azide_smiles: SMILES of the azide
    :param reaction_smarts: SMARTS for the reaction
    :return: Final product SMILES
    :raises ValueError: If reaction fails or invalid input
    """
    if drug_smiles is None:
        raise ValueError("❌ No drug SMILES provided.")
    
    linker_azide = build_linker_azide(linker_core_smiles, azide_smiles)
    drug_mol = Chem.MolFromSmiles(drug_smiles)
    if drug_mol is None:
        raise ValueError("❌ Invalid drug SMILES: " + drug_smiles)
    
    linker_mol = Chem.MolFromSmiles(linker_azide)
    if linker_mol is None:
        raise ValueError("❌ Invalid linker-azide SMILES: " + linker_azide)
    
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    products = rxn.RunReactants([drug_mol, linker_mol])
    if not products or not products[0]:
        raise ValueError("❌ Reaction failed: No product generated.")
    
    product_mol = products[0][0]
    try:
        Chem.SanitizeMol(product_mol)
    except Exception as e:
        print(f"⚠️ Warning: Product sanitization failed: {str(e)}")
    
    product_smiles = Chem.MolToSmiles(product_mol, isomericSmiles=True, canonical=False)
    print("✅ Click reaction successful. Final SMILES length:", len(product_smiles))
    return product_smiles

if __name__ == "__main__":
    # Test with proper alkyne sample (phenylacetylene for click reaction)
    sample_drug = "C#CC1=CC=CC=C1"  # Simple alkyne with phenyl
    final = run_click_reaction(sample_drug, LINKER_CORE, AZIDE, SMART_REACTION)
    print("Test Final SMILES (short):", final[:50] + "...")
    print("OK: reaction imported successfully.")