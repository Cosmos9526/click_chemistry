from typing import Dict
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED, Draw, rdMolDescriptors
import matplotlib.pyplot as plt
import pandas as pd

def calculate_properties(smiles: str) -> Dict[str, float]:
    """
    Calculate molecular properties for a given SMILES string.
    Returns a dictionary with properties and an OverallScore.
    
    :param smiles: Input SMILES
    :return: Dict of properties
    :raises ValueError: If SMILES invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("❌ Invalid SMILES string.")
    AllChem.ComputeGasteigerCharges(mol)
    
    properties = {}
    properties['MolecularWeight'] = Descriptors.MolWt(mol)
    properties['LogP'] = Descriptors.MolLogP(mol)
    properties['QED'] = QED.qed(mol)
    properties['TPSA'] = Descriptors.TPSA(mol)
    properties['HeavyAtomCount'] = Descriptors.HeavyAtomCount(mol)
    properties['SAScore'] = Descriptors.HeavyAtomCount(mol) / 10.0 + Descriptors.RingCount(mol)
    partial_charges = [float(atom.GetProp('_GasteigerCharge')) for atom in mol.GetAtoms() if atom.HasProp('_GasteigerCharge')]
    properties['AvgPartialCharge'] = sum(partial_charges) / len(partial_charges) if partial_charges else 0.0
    properties['NumRotatableBonds'] = Descriptors.NumRotatableBonds(mol)
    properties['NumHDonors'] = Descriptors.NumHDonors(mol)
    properties['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
    properties['RingCount'] = Descriptors.RingCount(mol)
    properties['FractionCSP3'] = rdMolDescriptors.CalcFractionCSP3(mol)
    properties['ExactMolWt'] = Descriptors.ExactMolWt(mol)
    properties['MolMR'] = Descriptors.MolMR(mol)
    properties['NumValenceElectrons'] = Descriptors.NumValenceElectrons(mol)
    properties['NumAromaticRings'] = Descriptors.NumAromaticRings(mol)
    properties['BioavailabilityScore'] = 1.0 if (properties['NumRotatableBonds'] <= 10 and properties['TPSA'] <= 140) else 0.0
    properties['StabilityScore'] = 1.0 / (1 + properties['NumRotatableBonds']) if properties['NumRotatableBonds'] > 0 else 1.0
    properties['NormLogP'] = max(0.0, 1 - abs(properties['LogP'] - 1.5) / 3.5)
    properties['NormTPSA'] = max(0.0, min(1.0, (200 - properties['TPSA']) / 60)) if properties['TPSA'] > 140 else 1.0
    properties['NormSAScore'] = max(0.0, 1 - properties['SAScore'] / 10)
    properties['NormFractionCSP3'] = properties['FractionCSP3']
    properties['NormMolecularWeight'] = max(0.0, min(1.0, (1000 - properties['MolecularWeight']) / 500)) if properties['MolecularWeight'] > 500 else 1.0
    properties['OverallScore'] = (
        properties['QED'] +
        properties['BioavailabilityScore'] +
        properties['StabilityScore'] +
        properties['NormLogP'] +
        properties['NormTPSA'] +
        properties['NormSAScore'] +
        properties['NormFractionCSP3'] +
        properties['NormMolecularWeight']
    ) / 8.0  # Normalized average for better scaling
    
    print("✅ Calculated properties (key samples): QED={:.3f}, OverallScore={:.3f}".format(properties['QED'], properties['OverallScore']))
    return properties

def display_properties_table(properties: Dict[str, float]) -> None:
    """Display properties in a pandas table."""
    props_df = pd.DataFrame(list(properties.items()), columns=['Property', 'Value'])
    props_df['Value'] = props_df['Value'].apply(lambda x: f"{x:.3f}" if isinstance(x, float) else x)
    print("\n✅ Chemical Properties Table:")
    print(props_df.to_string(index=False))

def draw_final_molecule(smiles: str) -> None:
    """Draw the final molecule image."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        plt.figure(figsize=(8, 8))
        img = Draw.MolToImage(mol, size=(800, 800))
        plt.imshow(img)
        plt.axis('off')
        plt.title("Final Molecule Structure")
        plt.show()
        print("🖼️ Final molecule image displayed.")
    else:
        print("❌ Could not draw the molecule: Invalid final SMILES.")

if __name__ == "__main__":
    # Test with sample (from reaction test)
    sample_final = "c1ccc(cc1)C#Cn2nnc3ccccc3n2NCCCC[C@H](NC(=O)[C@H](Cc4cnc[nH]4)NC(=O)[C@H](Cc5cnc[nH]5)NC(=O)[C@H](Cc6cnc[nH]6)NC(=O)[C@H](CCCCN)NC(=O)[C@H](Cc7cnc[nH]7)NC(=O)[C@H](Cc8cnc[nH]8)NC(=O)[C@H](Cc9cnc[nH]9)NC(=O)[C@H](Cc%10cnc[nH]%10)NC(=O)[C@H](Cc%11cnc[nH]%11)NC(=O)O)C(=O)O"  # Shortened sample for test
    props = calculate_properties(sample_final)
    display_properties_table(props)
    draw_final_molecule(sample_final)
    print("OK: properties imported successfully.")