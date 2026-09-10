# 2DQJ : chain H / chain Y
# Geometric visualization of the protein-protein interface

# Load structure
fetch 2DQJ, async=0
remove solvent
hide everything
bg_color white

# Define chains
select chain_H, chain H
select chain_Y, chain Y

# General representation
show cartoon, chain_H
show cartoon, chain_Y
color marine, chain_H
color orange, chain_Y

# Interface residues: heavy atoms within 6 A
select interface_H, byres ((chain H and not elem H) within 6.0 of (chain Y and not elem H))
select interface_Y, byres ((chain Y and not elem H) within 6.0 of (chain H and not elem H))

# Show interface residues
show sticks, interface_H
show sticks, interface_Y
color marine, interface_H
color orange, interface_Y

# Standard atom colors
color red, (interface_H or interface_Y) and elem O
color blue, (interface_H or interface_Y) and elem N
color yellow, (interface_H or interface_Y) and elem S

# Residue labels
label interface_H and name CA, "%s%s" % (resn, resi)
label interface_Y and name CA, "%s%s" % (resn, resi)

set label_size, 16
set label_font_id, 7
set label_outline_color, white

# General contacts <= 6 A

distance contacts_6A, interface_H and not elem H, interface_Y and not elem H, 6.0
set dash_color, grey70, contacts_6A
hide labels, contacts_6A

# Hydrophobic approximation: C...C <= 4 A
select carbon_H, interface_H and elem C
select carbon_Y, interface_Y and elem C

distance hydrophobic, carbon_H, carbon_Y, 4.0
set dash_color, grey50, hydrophobic
hide labels, hydrophobic


# Hydrogen-bond geometric approximation: N/O <= 4.1 A
select polar_H, interface_H and elem N+O
select polar_Y, interface_Y and elem N+O

distance hbonds, polar_H, polar_Y, 4.1
set dash_color, blue, hbonds
hide labels, hbonds


# Ionic interaction approximation
# ARG/LYS positive, ASP/GLU negative
select positive_H, interface_H and resn ARG+LYS and elem N
select negative_H, interface_H and resn ASP+GLU and elem O
select positive_Y, interface_Y and resn ARG+LYS and elem N
select negative_Y, interface_Y and resn ASP+GLU and elem O

distance salt_Hpos_Yneg, positive_H, negative_Y, 5.5
distance salt_Hneg_Ypos, negative_H, positive_Y, 5.5

set dash_color, yellow, salt_Hpos_Yneg
set dash_color, yellow, salt_Hneg_Ypos
hide labels, salt_Hpos_Yneg
hide labels, salt_Hneg_Ypos


# Appearance

set stick_radius, 0.18
set dash_radius, 0.05
set dash_gap, 0.25
set dash_length, 0.20
set orthoscopic, on
set depth_cue, 0
set antialias, 2

# Hide general 6 A contacts by default because they clutter view
disable contacts_6A

# Center on interface
orient interface_H or interface_Y
zoom interface_H or interface_Y, 6

# ============================================================
# Print residue pairs in terminal
# ============================================================

python
from pymol import cmd
import math

def collect_pairs(selection1, selection2, cutoff):
    """Collect unique residue pairs within a distance cutoff."""
    atoms1 = cmd.get_model(selection1).atom
    atoms2 = cmd.get_model(selection2).atom
    pairs = {}
    for atom1 in atoms1:
        for atom2 in atoms2:
            dx = atom1.coord[0] - atom2.coord[0]
            dy = atom1.coord[1] - atom2.coord[1]
            dz = atom1.coord[2] - atom2.coord[2]
            distance = math.sqrt(dx**2 + dy**2 + dz**2)
            if distance <= cutoff:
                key = (atom1.chain, atom1.resn, atom1.resi,
                       atom2.chain, atom2.resn, atom2.resi)
                if key not in pairs or distance < pairs[key]:
                    pairs[key] = distance
    return pairs

def print_pairs(title, pairs):
    """Print residue pairs and their minimum distance."""
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)
    print("Number of residue pairs:", len(pairs))
    for pair, distance in sorted(pairs.items()):
        chain1, resn1, resi1, chain2, resn2, resi2 = pair
        print(f"{chain1} {resn1}{resi1} - {chain2} {resn2}{resi2} : {distance:.2f} A")

# General interface contacts
general_pairs = collect_pairs(
    "interface_H and not elem H",
    "interface_Y and not elem H",
    6.0
)
print_pairs("GENERAL CONTACTS <= 6.0 A", general_pairs)

# Hydrophobic approximation
hydrophobic_pairs = collect_pairs(
    "carbon_H",
    "carbon_Y",
    4.0
)
print_pairs("HYDROPHOBIC CONTACTS: C...C <= 4.0 A", hydrophobic_pairs)

# Hydrogen-bond geometric approximation
hbond_pairs = collect_pairs(
    "polar_H",
    "polar_Y",
    4.1
)
print_pairs("HYDROGEN-BOND GEOMETRIC CONTACTS: N/O...N/O <= 4.1 A", hbond_pairs)

# Salt-bridge approximation
salt_pairs = collect_pairs(
    "positive_H",
    "negative_Y",
    5.5
)
salt_pairs.update(
    collect_pairs(
        "negative_H",
        "positive_Y",
        5.5
    )
)
print_pairs("IONIC CONTACTS <= 5.5 A", salt_pairs)
python end