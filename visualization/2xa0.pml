
# 2XA0 : BCL-2 chain A / BAX chain C
# Visualization of the protein-protein interface

# Load structure
fetch 2XA0, async=0
remove solvent
hide everything
bg_color white

# Define protein chains
select BCL2, chain A
select BAX, chain C

# General representation
show cartoon, BAX
color marine, BCL2
color orange, BAX

# Interface residues: heavy atoms within 6 A
select interface_A, chain A and not elem H within 6.0 of (chain C and not elem H)
select interface_C, chain C and not elem H within 6.0 of (chain A and not elem H)

# Show interface residues
show sticks, byres interface_A
show sticks, byres interface_C
color marine, byres interface_A
color orange, byres interface_C

# Standard atom colors
color red, (interface_A or interface_C) and elem O
color blue, (interface_A or interface_C) and elem N
color yellow, (interface_A or interface_C) and elem S


# Residue labels
# One label per interface residue, placed on the alpha carbon
label (byres interface_A) and name CA, "%s%s" % (resn, resi)
label (byres interface_C) and name CA, "%s%s" % (resn, resi)

# Label appearance
set label_size, 16
set label_font_id, 7
set label_outline_color, white


# General contacts <= 6 A
distance contacts_6A, interface_A, interface_C, 6.0
set dash_color, grey70, contacts_6A
hide labels, contacts_6A


# Hydrophobic approximation: C...C <= 4 A
select carbon_A, interface_A and elem C
select carbon_C, interface_C and elem C

distance hydrophobic, carbon_A, carbon_C, 4.0
set dash_color, grey50, hydrophobic
hide labels, hydrophobic

# Hydrogen-bond geometric approximation: N/O <= 4.1 A
select polar_A, interface_A and elem N+O
select polar_C, interface_C and elem N+O

distance hbonds, polar_A, polar_C, 4.1
set dash_color, blue, hbonds
hide labels, hbonds

# Ionic interaction approximation
# ARG/LYS positive, ASP/GLU negative

select positive_A, chain A and resn ARG+LYS and elem N
select negative_A, chain A and resn ASP+GLU and elem O
select positive_C, chain C and resn ARG+LYS and elem N
select negative_C, chain C and resn ASP+GLU and elem O

distance salt_Apos_Cneg, positive_A, negative_C, 5.5
distance salt_Aneg_Cpos, negative_A, positive_C, 5.5

set dash_color, yellow, salt_Apos_Cneg
set dash_color, yellow, salt_Aneg_Cpos
hide labels, salt_Apos_Cneg
hide labels, salt_Aneg_Cpos

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
orient interface_A or interface_C
zoom interface_A or interface_C, 6