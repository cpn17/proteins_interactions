# Presentation

The program allows you to:

- download a structure using its PDB ID;

- read the ATOM records from a PDB file;

- represent the structure using Protein, Chain, Residue, and Atom objects;

- select two protein chains;

- detect interface residue pairs based on a configurable threshold, set by default to 6 Å;

- detect four categories of interactions: hydrophobic contacts, hydrogen bonds, salt bridges, aromatic π–π interactions;

- generate a CSV file containing the residues in contact, their minimum distance, and the types of interactions detected.


# Installation


**1. Install uv by following astral documentation** 

https://docs.astral.sh/uv/getting-started/installation/

Please **restart your terminal** after uv successfulled installation.


**2. Clone this repository to your local space** 

```bash
git clone https://github.com/cpn17/proteins_interactions.git
```

**3.  Synchronize the environement**

- Go to the local repository : 

```bash
cd proteins_interactions
```

- Synchronize :

```bash
uv sync
```


**4. Launch the program with examples**

- PDB ID 2XA0, chaîns A and C :

```bash
uv run proteins-contacts 2xa0 --chains A,C --threshold 6.0
```

- PDB ID 2DQJ, chains H and Y :

```bash
uv run proteins-contacts 2dqj --chains H,Y --threshold 6.0
```


**5. Check the results**

For each PDB code, the program will 

- create directory named results/ 

- create a sub-directory named : Year-Month-Day-Hour-Minute-Second 

- In this sub directory, a .csv file contains all the contacts detected with classification into 4 types.

Examples : 

- results/2026-09-10_09-05-23/2xa0_A_C_contacts.csv

- results/2026-09-10_09-10-06/2dqj_H_Y_contacts.csv


**6. Personnalize the command with your PDB structures choices with the syntax**

proteins-contacts PDB_ID --chains CHAIN1,CHAIN2 [options]

Where the options are : 

--chains CHAIN1,CHAIN2 : chaîns to analyze

--threshold 6.0 : threshold to détect the interface, in Angstrom (Å)


**7. Visualization the contacts with Pymol (in developpement)**

These function is not completed and must be improved

```bash
pymol visualization/2xa0.pml
```

```bash
pymol visualization/2dqj.pml
```









