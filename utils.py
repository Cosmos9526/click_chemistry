from typing import List, Dict, Any
from rdkit import Chem
from rdkit.Chem import AllChem
from config import FUNCTIONAL_GROUPS  # Import from config

def validate_smiles(smiles: str, name: str) -> Chem.Mol:
    """Validate SMILES and return Mol object."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"❌ Invalid SMILES string for {name}: {smiles}")
    print(f"✅ Valid {name} SMILES: {smiles[:50]}..." if len(smiles) > 50 else f"✅ Valid {name} SMILES: {smiles}")
    return mol

def attach_star_to_treatment(smiles: str) -> str:
    """
    Attach a star (*) to the appropriate functional group in the treatment molecule.
    If it's a peptide, attach to N-terminal; otherwise, to the priority functional group.
    
    :param smiles: Input SMILES of the molecule
    :return: Starred SMILES
    :raises ValueError: If SMILES is invalid or no functional group is found
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"❌ Invalid SMILES: {smiles}")
    
    rw_mol = Chem.RWMol(mol)
    
    # Detect peptide by counting amide bonds (>2 for short peptides)
    amide_pattern = Chem.MolFromSmarts('[C:1](=[O:2])[N:3][C:4]')
    is_pep = len(mol.GetSubstructMatches(amide_pattern)) > 2
    
    if is_pep:
        print("✅ Detected peptide. Attempting to attach star to N-terminal NH2...")
        pattern = Chem.MolFromSmarts('[N;H2][C]')  # Primary amine connected to carbon
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            idx = matches[0][0]  # First match
            star_idx = rw_mol.AddAtom(Chem.Atom("*"))
            rw_mol.AddBond(idx, star_idx, Chem.BondType.SINGLE)
            print("✅ Star attached to N-terminal NH2.")
        else:
            raise ValueError("❌ N-terminal NH2 not found.")
    else:
        print("✅ Detected non-peptide drug. Searching for functional group...")
        found = False
        for name, smarts in FUNCTIONAL_GROUPS.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                # For carboxyl, select the OH/O- oxygen atom
                if name == "carboxyl":
                    idx = None
                    for match in matches:
                        for atom_idx in match:
                            atom = mol.GetAtomWithIdx(atom_idx)
                            if atom.GetSymbol() == 'O' and (atom.GetNumExplicitHs() == 1 or atom.GetFormalCharge() < 0):
                                idx = atom_idx
                                break
                        if idx is not None:
                            break
                    if idx is None:
                        continue
                else:
                    idx = matches[0][0]
                
                star_idx = rw_mol.AddAtom(Chem.Atom("*"))
                rw_mol.AddBond(idx, star_idx, Chem.BondType.SINGLE)
                found = True
                print(f"✅ Star attached to {name} group.")
                break
        
        if not found:
            raise ValueError("❌ No suitable functional group (N, O, C, S) found.")
    
    # Sanitize the molecule to ensure valid chemistry
    try:
        Chem.SanitizeMol(rw_mol)
    except Exception as e:
        raise ValueError(f"❌ Error in molecule sanitization: {str(e)}")
    
    starred_smiles = Chem.MolToSmiles(rw_mol, isomericSmiles=True, canonical=False).replace("*", "(*)").replace("((*))", "(*)")
    print("⭐️ Initial Starred SMILES for Treatment:", starred_smiles)
    return starred_smiles

def attach_star_to_carrier(smiles: str) -> str:
    """
    Attach a star (*) to the C-terminal of the carrier peptide by removing the oxygen and attaching to carbon.
    (Adapted for nanobody_tag_smiles, which already has * at lysine side chain, but this adds/ensures C-terminal attach if needed.
    Since linker_core has *, this may be used for variants.)
    
    :param smiles: Input SMILES of the carrier peptide
    :return: Starred SMILES
    :raises ValueError: If SMILES is invalid or no C-terminal found
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"❌ Invalid SMILES for carrier: {smiles}")
    
    rw_mol = Chem.RWMol(mol)
    
    # Find C-terminal carboxyl group
    pattern = Chem.MolFromSmarts('[C](=[O])[O;H1,H0,-1]')  # COOH or COO-
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        raise ValueError("❌ C-terminal carboxyl group not found in carrier peptide.")
    
    target_c_idx = None
    oxygen_idx = None
    for match in matches:
        c_idx = match[0]  # Carbon of COOH/COO-
        for atom_idx in match[1:]:  # Check oxygen atoms
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O' and (atom.GetNumExplicitHs() >= 0 or atom.GetFormalCharge() <= 0):
                oxygen_idx = atom_idx
                target_c_idx = c_idx
                break
        if target_c_idx is not None and oxygen_idx is not None:
            break
    
    if target_c_idx is None or oxygen_idx is None:
        raise ValueError("❌ No suitable C-terminal oxygen found in carrier peptide.")
    
    # Remove the oxygen (OH or O-)
    rw_mol.RemoveAtom(oxygen_idx)
    
    # Attach star to the carbon
    star_idx = rw_mol.AddAtom(Chem.Atom("*"))
    rw_mol.AddBond(target_c_idx, star_idx, Chem.BondType.SINGLE)
    
    # Sanitize the molecule
    try:
        Chem.SanitizeMol(rw_mol)
    except Exception as e:
        raise ValueError(f"❌ Error in carrier sanitization: {str(e)}")
    
    starred_smiles = Chem.MolToSmiles(rw_mol, isomericSmiles=True, canonical=False).replace("*", "(*)").replace("((*))", "(*)")
    print("⭐️ Initial Starred SMILES for Carrier:", starred_smiles)
    return starred_smiles

def generate_variants(smiles: str, entity_type: str = "Treatment") -> List[Dict[str, float]]:
    """
    Generate SMILES variants with star attached to terminal O/N atoms with negative charge.
    Merged function for both treatment and carrier.
    
    :param smiles: Input SMILES (without star)
    :param entity_type: "Treatment" or "Carrier" for logging
    :return: List of variant dictionaries with SMILES and PartialCharge
    :raises ValueError: If SMILES is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"❌ Invalid SMILES for {entity_type}: {smiles}")
    
    AllChem.ComputeGasteigerCharges(mol)
    terminal_atoms = []
    for atom in mol.GetAtoms():
        if atom.HasProp('_GasteigerCharge'):
            charge = float(atom.GetProp('_GasteigerCharge'))
            if charge < 0:  # Only negative charges
                symbol = atom.GetSymbol()
                degree = atom.GetDegree()
                if (symbol == 'O' and degree == 1) or (symbol == 'N' and degree <= 2):
                    terminal_atoms.append({"idx": atom.GetIdx(), "charge": charge, "symbol": symbol})
    
    variants = []
    for atom in terminal_atoms:
        rw_mol = Chem.RWMol(mol)
        star_idx = rw_mol.AddAtom(Chem.Atom('*'))
        rw_mol.AddBond(atom["idx"], star_idx, Chem.BondType.SINGLE)
        
        # Sanitize the variant
        try:
            Chem.SanitizeMol(rw_mol)
        except Exception as e:
            print(f"⚠️ Warning: Sanitize of variant failed for {entity_type} atom {atom['symbol']} (idx: {atom['idx']}): {str(e)}")
            continue
        
        new_smiles = Chem.MolToSmiles(rw_mol, isomericSmiles=True, canonical=False).replace("*", "(*)").replace("((*))", "(*)")
        variants.append({"SMILES": new_smiles, "PartialCharge": atom["charge"]})
    
    print(f"⭐️ {entity_type} Variants ({len(variants)} found):")
    for v in variants:
        print(f"  - SMILES: {v['SMILES'][:50]}..., PartialCharge: {v['PartialCharge']:.3f}")
    
    return variants

if __name__ == "__main__":
    # Test imports and basic functions
    from config import DEFAULT_MANUAL_SMILES
    print("Testing utils...")
    validate_smiles(DEFAULT_MANUAL_SMILES, "Test Manual")
    print("OK: utils imported successfully.")