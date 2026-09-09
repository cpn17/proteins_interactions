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