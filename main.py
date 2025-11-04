#!/usr/bin/env python3
"""
Main pipeline: Drug-linker-nanobody conjugation via click chemistry.
- Input: manual_smiles (treatment drug)
- Process: Generate variants → Evaluate with each linker (from JSON) → Select best by OverallScore
- Output: Click reaction → Final SMILES + properties table + image + report + JSON export
"""

import argparse
from typing import Dict, Any
import json
import os
import pandas as pd  # For table in report
from rdkit import Chem  # Added for type hints
from config import DEFAULT_MANUAL_SMILES, load_linkers, LINKER_CORE, AZIDE, SMART_REACTION, validate_core_smiles
from utils import generate_variants, attach_star_to_treatment
from evaluate import evaluate_clickable_variants
from reaction import run_click_reaction
from properties import calculate_properties, display_properties_table, draw_final_molecule

def main_pipeline(manual_smiles: str) -> None:
    """
    Run full pipeline.
    
    :param manual_smiles: Treatment drug SMILES
    """
    print("🚀 Starting Click Chemistry Pipeline...")
    
    # Validate and prepare
    validate_core_smiles()
    treatment_mol = validate_smiles(manual_smiles, "Manual (Treatment)")
    
    # Attach star to treatment (for logging, but variants use original)
    starred_treatment = attach_star_to_treatment(manual_smiles)
    
    # Generate treatment variants (once, shared for all linkers)
    treatment_variants = generate_variants(manual_smiles, "Treatment")
    if not treatment_variants:
        raise ValueError("❌ No treatment variants generated.")
    
    # Load linkers and select best by OverallScore
    linkers = load_linkers()
    best_linker = select_best_linker_with_overall(linkers, treatment_variants)
    
    # Run click reaction for best
    best_combined_smiles = best_linker["best_smiles"]
    final_smiles = run_click_reaction(
        drug_smiles=best_combined_smiles,
        linker_core_smiles=LINKER_CORE,
        azide_smiles=AZIDE,
        reaction_smarts=SMART_REACTION
    )
    
    # Calculate final properties and display
    final_props = calculate_properties(final_smiles)
    display_properties_table(final_props)
    draw_final_molecule(final_smiles)
    
    # Report: Best linker details (from JSON)
    print("\n📊 FINAL REPORT:")
    print(f"🎯 Best Linker: {best_linker['name']}")
    print(f"   Cleavage Type: {best_linker.get('cleavage_type', 'N/A')}")
    print(f"   Components: {', '.join(best_linker.get('components', []))}")
    print(f"   OverallScore (pre-reaction): {best_linker['overall_score']:.3f}")
    print(f"   QED (pre-reaction): {best_linker['qed']:.3f}")
    print(f"✅ Final SMILES saved to 'outputs/final_best_smiles.txt' (length: {len(final_smiles)})")
    
    # Save final SMILES to outputs dir
    os.makedirs("outputs", exist_ok=True)
    with open('outputs/final_best_smiles.txt', 'w') as f:
        f.write(final_smiles)
    
    # Export final output to JSON in outputs
    final_report = {
        "best_linker": {
            "name": best_linker['name'],
            "cleavage_type": best_linker.get('cleavage_type', 'N/A'),
            "components": best_linker.get('components', []),
            "overall_score": round(best_linker['overall_score'], 3),
            "qed": round(best_linker['qed'], 3),
            "best_combined_smiles": best_combined_smiles[:200] + "..." if len(best_combined_smiles) > 200 else best_combined_smiles  # Short for JSON
        },
        "final_smiles": final_smiles,
        "final_properties": {k: round(v, 3) if isinstance(v, float) else v for k, v in final_props.items()},
        "treatment_input": manual_smiles,
        "nanobody_tag": LINKER_CORE[:200] + "..." if len(LINKER_CORE) > 200 else LINKER_CORE
    }
    with open('outputs/final_report.json', 'w') as f:
        json.dump(final_report, f, indent=4)
    print("💾 Final report saved to 'outputs/final_report.json' (includes best linker, final SMILES, properties).")
    print("💾 Project complete! Check mounted 'outputs/' for details.")

def select_best_linker_with_overall(linkers: list[Dict[str, Any]], treatment_variants: list[Dict[str, float]]) -> Dict[str, Any]:
    """
    Modified: Loop over linkers, evaluate, calculate OverallScore on best_combined, select max OverallScore.
    Skip linkers without 'starred_smiles'.
    """
    best_linker = None
    best_overall = -1.0
    evaluated_count = 0
    for linker in linkers:
        if "starred_smiles" not in linker:
            print(f"⚠️ Skipping linker '{linker.get('name', 'Unknown')}' - missing 'starred_smiles'.")
            continue
        scaffold = linker["starred_smiles"]
        print(f"\n🔍 Evaluating linker: {linker['name']} ({len(linkers)} total)")
        sorted_combined = evaluate_clickable_variants(treatment_variants, scaffold)
        if sorted_combined:
            evaluated_count += 1
            best_combined_smiles = sorted_combined[0]["SMILES"]
            this_qed = sorted_combined[0]["QED"]
            # Calculate OverallScore on combined
            this_props = calculate_properties(best_combined_smiles)
            this_overall = this_props['OverallScore']
            if this_overall > best_overall:
                best_overall = this_overall
                best_linker = linker.copy()
                best_linker["best_smiles"] = best_combined_smiles
                best_linker["qed"] = this_qed
                best_linker["overall_score"] = this_overall
                print(f"📈 New best: {linker['name']} (OverallScore: {this_overall:.3f})")
    
    print(f"\n📊 Evaluated {evaluated_count}/{len(linkers)} linkers successfully.")
    if best_linker:
        print(f"🎯 Overall Best Linker Selected: {best_linker['name']} with OverallScore {best_overall:.3f}")
        return best_linker
    else:
        raise ValueError("❌ No valid linker found after evaluation.")

def validate_smiles(smiles: str, name: str) -> Chem.Mol:  # Now Chem imported at top
    """Validate SMILES and return Mol object (inline from utils)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"❌ Invalid SMILES string for {name}: {smiles}")
    print(f"✅ Valid {name} SMILES: {smiles[:50]}..." if len(smiles) > 50 else f"✅ Valid {name} SMILES: {smiles}")
    return mol

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Click Chemistry Pipeline")
    parser.add_argument("--manual_smiles", type=str, default=DEFAULT_MANUAL_SMILES,
                        help="Treatment drug SMILES (default: " + DEFAULT_MANUAL_SMILES[:30] + "...)")
    args = parser.parse_args()
    main_pipeline(args.manual_smiles)