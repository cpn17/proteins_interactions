# Presentation

The program allows you to:

- download a structure using its PDB ID;

- read the ATOM records from a PDB file;

- represent the structure using Protein, Chain, Residue, and Atom objects;

- select two protein chains;

- detect interface residue pairs based on a configurable threshold, set by default to 6 Å;

- detect four categories of interactions: hydrophobic contacts, hydrogen bonds, salt bridges, aromatic π–π interactions;

- generate a CSV file containing the residues in contact, their minimum distance, and the types of interactions detected;

- visualize the interface using a PyMOL script.


# Installation

## With uv :

1. Install uv by following astral documentation : 

https://docs.astral.sh/uv/getting-started/installation/

2.  Synchronize the dependancies

uv sync

3. Clone this repository to your local space 

git clone https://github.com/cpn17/proteins_interactions.git

4. Launch the program :

- PDB ID 2XA0, chaîns A and C :

uv run proteins-contacts 2xa0 --chains A,C --threshold 6.0

- PDB ID 2DQJ, chains H and Y :

uv run proteins-contacts 2dqj --chains H,Y --threshold 6.0

5. Personnalize the command with your PDB structures choices with the syntax : 

proteins-contacts PDB_ID --chains CHAIN1,CHAIN2 [options]

Where the options are : 

--chains CHAIN1,CHAIN2 : chaîns to analyze
--threshold 6.0 : threshol to détect the interface in Angstrom (Å)









