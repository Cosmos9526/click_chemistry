import json
import os
from typing import Dict, List, Any
from rdkit import Chem

# Default inputs (manual_smiles changeable via main args)
DEFAULT_MANUAL_SMILES = "[H][C@]12C[C@@H](O)C=C[C@]11CCN(C)CC3=C1C(O2)=C(OC)C=C3"

# Fixed SMILES (nanobody tag with star at end)
NANOBODY_TAG_SMILES = "N[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CCCCN(*))C(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)O"

# Linker core = nanobody tag
LINKER_CORE = NANOBODY_TAG_SMILES

# Azide and reaction (fixed)
AZIDE = "N=N-N(*)"
SMART_REACTION = "[#6:7][C:6]#[C:5].[A:4][#7:3]-[N:2]=[#7:1]>>[A:4]-[#7:3]-1-[#6:5]=[#6:6](-[#6:7])-[#7:1]=[#7:2]-1"

# Functional groups (for attach_star)
FUNCTIONAL_GROUPS = {
    "amine": "[N;H2,H3][C]",       # Primary/secondary amines
    "hydroxyl": "[O;H1][C]",      # Hydroxyl group
    "carboxyl": "[C](=[O])[OH,O-]",  # Carboxyl group (neutral or ionized)
    "thiol": "[S;H1,H0][C]"       # Thiol group
}

# Load linkers from JSON
def load_linkers(json_path: str = "merged_starred_linkers.json") -> List[Dict[str, Any]]:
    """
    Load starred linkers from JSON file.
    Returns list of dicts with 'name', 'starred_smiles', etc. (full info for report).
    Validates only 'starred_smiles'. If file not found, raises error (no samples).
    """
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"❌ JSON file '{json_path}' not found. Please provide the full file for complete linkers.")
    
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Validate each starred_smiles (keep full dict for report if valid)
    valid_linkers = []
    for linker in data:
        starred_smiles = linker.get("starred_smiles", "")
        mol = Chem.MolFromSmiles(starred_smiles)
        if mol is not None:
            valid_linkers.append(linker)  # Keep full info (name, cleavage_type, etc.) for report
        else:
            print(f"⚠️ Invalid starred_smiles for {linker.get('name', 'Unknown')}: Skipping.")
    
    if not valid_linkers:
        raise ValueError("❌ No valid starred_smiles found in JSON. Check the file.")
    
    print(f"✅ Loaded {len(valid_linkers)} valid linkers from JSON. (Only starred_smiles used for processing; full info for report)")
    return valid_linkers

# Validate core SMILES
def validate_core_smiles():
    """Validate nanobody tag and linker core."""
    mol = Chem.MolFromSmiles(NANOBODY_TAG_SMILES)
    if mol is None:
        raise ValueError("❌ Invalid nanobody_tag_smiles.")
    print("✅ Valid Nanobody Tag SMILES (Linker Core):", NANOBODY_TAG_SMILES[:50] + "...")

if __name__ == "__main__":
    validate_core_smiles()
    linkers = load_linkers()
    print("Sample linker names:", [l['name'] for l in linkers[:2]])
    print("Focus: Only 'starred_smiles' processed; others (name, cleavage_type) for report.")