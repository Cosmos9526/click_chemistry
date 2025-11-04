from typing import List, Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import QED, Draw
import matplotlib.pyplot as plt
from utils import generate_variants
from config import DEFAULT_MANUAL_SMILES

def evaluate_clickable_variants(drug_variants: List[Dict[str, float]], scaffold_smiles: str) -> List[Dict[str, Any]]:
    """
    Evaluate combined molecules by QED after combining variants with scaffold.
    
    :param drug_variants: List of variants (e.g., treatment)
    :param scaffold_smiles: SMILES of the scaffold/linker (starred_smiles from JSON)
    :return: Sorted list of combined molecules by QED (descending)
    """
    combined_molecules = []
    for variant in drug_variants:
        combined_smiles = scaffold_smiles + '.' + variant['SMILES']
        combined_smiles = combined_smiles.replace('(*)','9')  # Replace placeholder for combination (RDKit handles as bond)
        mol = Chem.MolFromSmiles(combined_smiles)
        if mol:
            try:
                qed_value = QED.qed(mol)
                combined_molecules.append({
                    "SMILES": Chem.MolToSmiles(mol, isomericSmiles=True, canonical=False),
                    "QED": qed_value,
                    "PartialCharge": variant['PartialCharge'],
                    "MOL": mol
                })
                print(f"✅ Combined SMILES (short): {combined_smiles[:50]}..., QED: {qed_value:.3f}")
            except Exception as e:
                print(f"⚠️ Error evaluating variant: {str(e)}")
                continue
    # Sort by QED descending
    sorted_combined = sorted(combined_molecules, key=lambda x: x["QED"], reverse=True)
    return sorted_combined

def show_best_molecule(molecule_data: List[Dict[str, Any]]) -> Optional[str]:
    """
    Display the best molecule based on QED and return its SMILES.
    
    :param molecule_data: List of evaluated molecules
    :return: SMILES of the best molecule
    """
    if not molecule_data:
        print("❌ No valid molecule data found.")
        return None
    
    best = molecule_data[0]
    print(f"✅ Best QED: {best['QED']:.3f}, PartialCharge: {best['PartialCharge']:.3f}")
    
    # Display image
    plt.figure(figsize=(6, 6))
    img = Draw.MolToImage(best["MOL"], size=(600, 600))
    plt.imshow(img)
    plt.axis('off')
    plt.title("Best Molecule Structure")
    plt.show()
    
    print(f"✅ Best SMILES (short): {best['SMILES'][:100]}...")
    return best['SMILES']

def select_best_linker(linkers: List[Dict[str, Any]], treatment_variants: List[Dict[str, float]]) -> Dict[str, Any]:
    """
    Loop over linkers, evaluate each, and select the best based on max QED.
    (In full pipeline, use OverallScore from properties.py for better selection.)
    
    :param linkers: List from config.load_linkers()
    :param treatment_variants: Variants from generate_variants
    :return: Best linker dict with added 'best_smiles' and 'qed'
    """
    best_linker = None
    best_qed = -1.0
    for linker in linkers:
        scaffold = linker["starred_smiles"]
        print(f"\n🔍 Evaluating linker: {linker['name']}")
        sorted_combined = evaluate_clickable_variants(treatment_variants, scaffold)
        if sorted_combined:
            this_qed = sorted_combined[0]["QED"]
            if this_qed > best_qed:
                best_qed = this_qed
                best_linker = linker.copy()
                best_linker["best_smiles"] = sorted_combined[0]["SMILES"]
                best_linker["qed"] = this_qed
                print(f"📈 New best: {linker['name']} (QED: {this_qed:.3f})")
    
    if best_linker:
        print(f"\n🎯 Overall Best Linker: {best_linker['name']} with QED {best_qed:.3f}")
        return best_linker
    else:
        raise ValueError("❌ No valid linker found.")

if __name__ == "__main__":
    # Test with default
    print("Testing evaluate.py...")
    from config import load_linkers
    linkers = load_linkers()
    variants = generate_variants(DEFAULT_MANUAL_SMILES, "Test Treatment")
    # Test single evaluate (first linker)
    first_scaffold = linkers[0]["starred_smiles"]
    sorted_combined = evaluate_clickable_variants(variants, first_scaffold)
    best_smiles = show_best_molecule(sorted_combined)
    print("OK: evaluate imported successfully.")