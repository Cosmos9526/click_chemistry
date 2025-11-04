Click Chemistry Pipeline
A simple Python pipeline for drug-linker-nanobody conjugation using RDKit for SMILES processing, click reaction simulation, and property evaluation. Selects the best linker from a JSON library and generates final conjugated SMILES.
Requirements

Docker (for containerized run)
Python 3.12+ with RDKit, pandas, matplotlib (handled in Docker)

Quick Start

Build Docker Image:
textdocker build -t click-chemistry .

Run Pipeline (outputs saved to ./outputs/):
textmkdir -p outputs
docker run -v $(pwd)/outputs:/app/outputs -it click-chemistry python main.py --manual_smiles "YOUR_SMILES_HERE"

Example SMILES: [H][C@]12C[C@@H](O)C=C[C@]11CCN(C)CC3=C1C(O2)=C(OC)C=C3
Takes ~5-10 min for 106 linkers.



Outputs

outputs/final_best_smiles.txt: Final conjugated SMILES
outputs/final_report.json: Best linker details, properties, and scores
Console: Properties table and logs

Files

main.py: Entry point (arg: --manual_smiles)
config.py, utils.py, etc.: Core logic
merged_starred_linkers.json: Linker library (106 entries)

Notes

Edit DEFAULT_MANUAL_SMILES in config.py for defaults.
For local run (no Docker): pip install rdkit pandas matplotlib then python main.py.

MIT License. Built with RDKit for chemistry sim.
