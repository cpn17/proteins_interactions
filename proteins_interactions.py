from optparse import OptionParser

from Bio.PDB import PDBList
import math
from openbabel import openbabel

AA_STANDARD = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
HYDROPHOBIC_DISTANCE = 4.0
HBOND_DISTANCE = 4.1
HBOND_ANGLE = 100.0
MIN_DISTANCE = 0.5
SALT_BRIDGE_DISTANCE = 5.5
PISTACK_DISTANCE = 5.5
PISTACK_ANGLE_DEVIATION = 30.0
PISTACK_OFFSET = 2.0

class Atom:
    """Represent an atom extracted from a PDB file."""
    def __init__(self, atom_serial_number, atom_name, 
                 alternate_location_indicator, x, y, z,
                 occupancy, temperature_factor, element_symbol, charge_on_the_atom):
        self.atom_serial_number = atom_serial_number
        self.atom_name = atom_name
        self.alternate_location_indicator = alternate_location_indicator
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.temperature_factor = temperature_factor
        self.element_symbol = element_symbol
        self.charge_on_the_atom = charge_on_the_atom

    def is_hydrogen(self):
        """Return True if the atom is a hydrogen atom, otherwise False."""
        return self.element_symbol == "H"

    def distance_to(self, other_atom):
        """Calculate euclidean distance to other atom."""
        dx = self.x - other_atom.x
        dy = self.y - other_atom.y
        dz = self.z - other_atom.z
        distance = math.sqrt(dx**2 + dy**2 + dz**2)
        return distance

    def __str__(self):
        """Return a string describing the atom and its coordinates."""
        return f"Atom name : {self.atom_name}, coords : {self.x:.3f} {self.y:.3f} {self.z:.3f}"

class Residue:
    """Represent a residue containing many atoms."""
    def __init__(self, residue_name, residue_sequence_number, code_for_insertion_of_residues):
        self.residue_name = residue_name
        self.residue_sequence_number = residue_sequence_number
        self.code_for_insertion_of_residues = code_for_insertion_of_residues
        self.atoms = []
        self.hydrophobic_atoms = []
        self.hbond_donors = []
        self.hbond_acceptors = []
        self.charged_groups = []
        self.aromatic_rings = []

    def find_atom(self, atom_name):
        """Find an atom in the residue from its name."""
        for atom in self.atoms:
            if atom.atom_name == atom_name:
                return atom
        return None

    def find_atom_by_serial(self, atom_serial_number):
        """Find an atom in the residue from its serial number."""
        for atom in self.atoms:
            if atom.atom_serial_number == atom_serial_number:
                return atom
        return None

    def find_atom_by_coordinates(self, x, y, z, tolerance=0.01):
        """Find an atom in the residue from its coordinates."""
        for atom in self.atoms:
            if abs(atom.x - x) < tolerance and abs(atom.y - y) < tolerance and abs(atom.z - z) < tolerance:
                return atom
        return None

    def add_atom(self, atom):
        """Add an atom to the residue."""
        if atom.is_hydrogen():
            self.atoms.append(atom)
            return
        for index, existing_atom in enumerate(self.atoms):
            if existing_atom.atom_name == atom.atom_name:
                if atom.occupancy > existing_atom.occupancy:
                    self.atoms[index] = atom
                return
        self.atoms.append(atom)

    def contains_hydrogen(self):
        """Return True if the residue contains at least one hydrogen atom."""
        for atom in self.atoms:
            if atom.is_hydrogen():
                return True
        return False

    def __str__(self):
        """Return a string describing the residue."""
        return f"Residue : {self.residue_name} Residue sequence number : {self.residue_sequence_number} Residue code for insertion : {self.code_for_insertion_of_residues}"

class Chain:
    """Represent a protein chain containing a collection of residues."""
    def __init__(self, chain_identifier):
        """Initialize a protein chain from its PDB chain identifier."""
        self.chain_identifier = chain_identifier
        self.residues = []

    def find_residue(self, residue_sequence_number, code_for_insertion_of_residues):
        """Find a residue in the chain.

        Parameters
        ----------
        residue_sequence_number : int
            Sequence number of the residue.
        code_for_insertion_of_residues : str
            PDB insertion code of the residue.

        Returns
        -------
        Residue or None
            Matching residue if found, otherwise None.
        """
        for residue in self.residues:
            if residue.residue_sequence_number == residue_sequence_number and residue.code_for_insertion_of_residues == code_for_insertion_of_residues:
                return residue
        return None

    def add_atom_to_residue(self, residue_name, residue_sequence_number, code_for_insertion_of_residues, atom):
        """Add an atom to the corresponding residue of the chain.

        A new residue is created if it is not already present.

        Parameters
        ----------
        residue_name : str
            Three-letter residue name.
        residue_sequence_number : int
            Sequence number of the residue.
        code_for_insertion_of_residues : str
            PDB insertion code of the residue.
        atom : Atom
            Atom to add.

        Returns
        -------
        None
        """
        residue = self.find_residue(residue_sequence_number, code_for_insertion_of_residues)
        if residue is None:
            residue = Residue(residue_name, residue_sequence_number, code_for_insertion_of_residues)
            self.residues.append(residue)
        residue.add_atom(atom)

    def number_of_atoms(self):
        """Return the total number of atoms in the chain."""
        number = 0
        for residue in self.residues:
            number += len(residue.atoms)
        return number

class Protein: 
    """Represent a protein structure organized into chains, residues and atoms."""
    def __init__(self, protein_name):
        self.protein_name = protein_name
        self.chains = []
        self.missing_residues = []
        self.missing_atoms = []
        self.excluded_waters = set()
        self.excluded_hetero_molecules = set()
        self.alternative_locations = set()

    def find_chain(self, chain_identifier):
        """Find a chain from its identifier.

        Parameters
        ----------
        chain_identifier : str
            Identifier of the chain.

        Returns
        -------
        Chain or None
            Matching chain if found, otherwise None.
        """
        for chain in self.chains:
            if chain.chain_identifier == chain_identifier:
                return chain
        return None
    
    def add_atom_to_chain(self, chain_identifier, residue_name, residue_sequence_number, code_for_insertion_of_residues, atom):
        """Add an atom to the corresponding chain and residue.

        A new chain is created if it is not already present.

        Parameters
        ----------
        chain_identifier : str
            Identifier of the chain.
        residue_name : str
            Three-letter residue name.
        residue_sequence_number : int
            Sequence number of the residue.
        code_for_insertion_of_residues : str
            PDB insertion code of the residue.
        atom : Atom
            Atom to add.

        Returns
        -------
        None
        """
        chain = self.find_chain(chain_identifier)
        if chain is None : 
            chain = Chain(chain_identifier)
            self.chains.append(chain)
        chain.add_atom_to_residue(residue_name, residue_sequence_number, code_for_insertion_of_residues, atom)
        if atom.alternate_location_indicator != "":
            self.alternative_locations.add((chain_identifier, residue_name, residue_sequence_number, code_for_insertion_of_residues, atom.atom_name, atom.alternate_location_indicator))

    def number_of_residues(self):
        """Return the total number of residues in the protein."""
        number = 0
        for chain in self.chains:
            number += len(chain.residues)
        return number

    def number_of_atoms(self):
        """Return the total number of atoms in the protein."""
        number = 0
        for chain in self.chains:
            number += chain.number_of_atoms()
        return number

    def contains_hydrogen(self):
        """Return True if the protein contains at least one hydrogen atom."""
        for chain in self.chains: 
            for residue in chain.residues: 
                if residue.contains_hydrogen():
                    return True
        return False

    def __str__(self):
        """Return a string summarizing the protein structure."""
        return f"Protein : {self.protein_name}, chains : {len(self.chains)}, {self.number_of_residues()} residues, {self.number_of_atoms()} atoms"

def download_pdb(pdb_id):
    """Download a PDB structure from the Protein Data Bank.

    Parameters
    ----------
    pdb_id : str
        PDB identifier, ex : 2xa0.

    Returns
    -------
    str
        Path of the downloaded PDB file.
    """
    pdb_list = PDBList()
    file_name = pdb_list.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")
    return file_name

def read_missing_residue(line, protein):
    """Store a missing residue reported in a REMARK 465 line."""
    fields = line.split()
    if len(fields) >= 5 and fields[2] in AA_STANDARD:
        residue_name = fields[2]
        chain_identifier = fields[3]
        residue_sequence_number = fields[4]
        protein.missing_residues.append((chain_identifier, residue_name, residue_sequence_number))

def read_missing_atoms(line, protein):
    """Store missing atoms reported in a REMARK 470 record."""
    fields = line.split()
    if len(fields) >= 6 and fields[2] in AA_STANDARD: 
        residue_name = fields[2]
        chain_identifier = fields[3]
        residue_sequence_number = fields[4]
        atom_names = fields[5:] # list of missing atoms of a residue
        protein.missing_atoms.append((chain_identifier, residue_name, residue_sequence_number, atom_names))

def read_pdb(file_name):
    """Read protein atoms from a PDB file.

    Lines startswith "ATOM" are parsed and organized into Atom, Residue, Chain and Protein objects.

    Parameters
    ----------
    file_name : str
        Name or path of the PDB file.

    Returns
    -------
    Protein
        Protein object containing the atoms, residues and chains read from the PDB file.
    """
    protein = Protein(file_name)
    with open(file_name, "r") as file: 
        for line in file:
            if line.startswith("REMARK 465"):
                read_missing_residue(line, protein)
                continue
            if line.startswith("REMARK 470"):
                read_missing_atoms(line, protein)
                continue
            line_keyword = line[0:6].strip()
            if line_keyword == "HETATM":
                residue_name = line[17:20].strip()
                chain_identifier = line[21:22].strip()
                residue_sequence_number = int(line[22:26])
                code_for_insertion_of_residues = line[26:27].strip()
                if residue_name == "HOH":
                    protein.excluded_waters.add((chain_identifier, residue_sequence_number, code_for_insertion_of_residues))
                else:
                    protein.excluded_hetero_molecules.add((chain_identifier, residue_name, residue_sequence_number, code_for_insertion_of_residues))
                continue
            if line_keyword != "ATOM":
                continue
            atom_serial_number = int(line[6:11])
            atom_name = line[12:16].strip()
            alternate_location_indicator = line[16:17].strip()
            residue_name = line[17:20].strip()
            chain_identifier = line[21:22].strip()
            residue_sequence_number = int(line[22:26])
            code_for_insertion_of_residues = line[26:27].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occupancy = float(line[54:60])
            temperature_factor = float(line[60:66])
            element_symbol = line[76:78].strip().upper()
            charge_on_the_atom = line[78:80].strip()
            atom = Atom(atom_serial_number, atom_name, alternate_location_indicator, x, y, z, occupancy, temperature_factor, element_symbol, charge_on_the_atom)
            protein.add_atom_to_chain(chain_identifier, residue_name, residue_sequence_number, code_for_insertion_of_residues, atom)
    return protein

def find_interface_pairs(chain1, chain2, threshold):
    """Find residue pairs in contact and their minimum heavy-atom distance."""
    interface_pairs = []
    for residue1 in chain1.residues:
        for residue2 in chain2.residues:
            minimum_distance = None
            for atom1 in residue1.atoms:
                if atom1.is_hydrogen():
                    continue
                for atom2 in residue2.atoms:
                    if atom2.is_hydrogen():
                        continue
                    distance = atom1.distance_to(atom2)
                    if minimum_distance is None or distance < minimum_distance:
                        minimum_distance = distance
            if minimum_distance is not None and minimum_distance <= threshold:
                interface_pairs.append((residue1, residue2, minimum_distance))
    return interface_pairs

def add_hydrogens(file_name):
    """Add hydrogen atoms to a PDB file using Open Babel."""
    output_file = file_name.rsplit(".", 1)[0] + "_H.pdb"
    conversion = openbabel.OBConversion()
    conversion.SetInAndOutFormats("pdb", "pdb")
    molecule = openbabel.OBMol()
    conversion.ReadFile(molecule, file_name)
    molecule.AddHydrogens()
    conversion.WriteFile(molecule, output_file)
    return output_file

# Detection of hydrophobic contacts
def read_openbabel_molecule(file_name):
    """Read a PDB file with Open Babel.

    Parameters
    ----------
    file_name : str
        Path of the PDB file.

    Returns
    -------
    openbabel.OBMol
        Molecule containing atoms and perceived bonds.
    """
    conversion = openbabel.OBConversion()
    conversion.SetInFormat("pdb")
    molecule = openbabel.OBMol()
    conversion.ReadFile(molecule, file_name)
    return molecule
                
def prepare_hydrophobic_atoms(protein, molecule):
    """Identify hydrophobic carbon atoms using Open Babel."""
    for ob_atom in openbabel.OBMolAtomIter(molecule):
        if ob_atom.GetAtomicNum() != 6:
            continue
        hydrophobic = True
        for neighbor in openbabel.OBAtomAtomIter(ob_atom):
            if neighbor.GetAtomicNum() not in (1, 6):
                hydrophobic = False
                break
        if not hydrophobic:
            continue
        ob_residue = ob_atom.GetResidue()
        if ob_residue is None:
            continue
        chain = protein.find_chain(ob_residue.GetChain())
        if chain is None:
            continue
        residue = chain.find_residue(ob_residue.GetNum(), "")
        if residue is None:
            continue
        atom_name = ob_residue.GetAtomID(ob_atom).strip()
        atom = residue.find_atom(atom_name)
        if atom is not None:
            residue.hydrophobic_atoms.append(atom)

def map_openbabel_atom(protein, ob_atom):
    """Map an Open Babel atom to an Atom object."""
    ob_residue = ob_atom.GetResidue()
    if ob_residue is None:
        return None
    chain = protein.find_chain(ob_residue.GetChain())
    if chain is None:
        return None
    residue = chain.find_residue(ob_residue.GetNum(), "")
    if residue is None:
        return None
    return residue.find_atom_by_coordinates(ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ())

def prepare_hbond_acceptors(protein, molecule):
    """Identify hydrogen-bond acceptor atoms using Open Babel."""
    for ob_atom in openbabel.OBMolAtomIter(molecule):
        if not ob_atom.IsHbondAcceptor():
            continue
        atom = map_openbabel_atom(protein, ob_atom)
        if atom is None:
            continue
        ob_residue = ob_atom.GetResidue()
        chain = protein.find_chain(ob_residue.GetChain())
        residue = chain.find_residue(ob_residue.GetNum(), "")
        residue.hbond_acceptors.append(atom)

def prepare_hbond_donors(protein, molecule):
    """Identify hydrogen-bond donor and hydrogen atom pairs using Open Babel."""
    for ob_donor in openbabel.OBMolAtomIter(molecule):
        if not ob_donor.IsHbondDonor():
            continue
        donor = map_openbabel_atom(protein, ob_donor)
        if donor is None:
            continue
        ob_residue = ob_donor.GetResidue()
        chain = protein.find_chain(ob_residue.GetChain())
        residue = chain.find_residue(ob_residue.GetNum(), "")
        for ob_hydrogen in openbabel.OBAtomAtomIter(ob_donor):
            if not ob_hydrogen.IsHbondDonorH():
                continue
            hydrogen = map_openbabel_atom(protein, ob_hydrogen)
            if hydrogen is not None:
                residue.hbond_donors.append((donor, hydrogen))

def calculate_angle(atom1, vertex, atom2):
    """Calculate the angle between three atoms in degrees."""
    vector1 = (atom1.x - vertex.x, atom1.y - vertex.y, atom1.z - vertex.z)
    vector2 = (atom2.x - vertex.x, atom2.y - vertex.y, atom2.z - vertex.z)
    dot_product = sum(a * b for a, b in zip(vector1, vector2))
    norm1 = math.sqrt(sum(a**2 for a in vector1))
    norm2 = math.sqrt(sum(a**2 for a in vector2))
    cosine = dot_product / (norm1 * norm2)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))

def prepare_protein(protein, file_name):
    """Prepare chemical features used for interaction detection."""
    molecule = read_openbabel_molecule(file_name)
    prepare_hydrophobic_atoms(protein, molecule)
    prepare_hbond_acceptors(protein, molecule)
    prepare_hbond_donors(protein, molecule)
    prepare_charged_groups(protein)
    prepare_aromatic_rings(protein)

def detect_hydrophobic_contact(residue1, residue2):
    """Detect a hydrophobic contact between two residues.

    Parameters
    ----------
    residue1 : Residue
        First residue.
    residue2 : Residue
        Second residue.

    Returns
    -------
    float or None
        Minimum distance between hydrophobic atoms if a contact is found, otherwise None.
    """
    minimum_distance = None
    for atom1 in residue1.hydrophobic_atoms:
        for atom2 in residue2.hydrophobic_atoms:
            distance = atom1.distance_to(atom2)
            if distance <= HYDROPHOBIC_DISTANCE:
                if minimum_distance is None or distance < minimum_distance:
                    minimum_distance = distance
    return minimum_distance

def find_hydrogen_bonds(residue1, residue2):
    """Find hydrogen bonds between two residues."""
    hydrogen_bonds = []
    for donor, hydrogen in residue1.hbond_donors:
        for acceptor in residue2.hbond_acceptors:
            distance = donor.distance_to(acceptor)
            if not MIN_DISTANCE < distance < HBOND_DISTANCE:
                continue
            angle = calculate_angle(donor, hydrogen, acceptor)
            if angle > HBOND_ANGLE:
                hydrogen_bonds.append((donor, hydrogen, acceptor, distance, angle))
    for donor, hydrogen in residue2.hbond_donors:
        for acceptor in residue1.hbond_acceptors:
            distance = donor.distance_to(acceptor)
            if not MIN_DISTANCE < distance < HBOND_DISTANCE:
                continue
            angle = calculate_angle(donor, hydrogen, acceptor)
            if angle > HBOND_ANGLE:
                hydrogen_bonds.append((donor, hydrogen, acceptor, distance, angle))
    return hydrogen_bonds

class ChargedGroup:
    """Represent a charged functional group."""
    def __init__(self, charge_type, atoms):
        """Initialize a charged group."""
        self.charge_type = charge_type
        self.atoms = atoms

    def center(self):
        """Calculate the geometric center of the charged group."""
        x = sum(atom.x for atom in self.atoms) / len(self.atoms)
        y = sum(atom.y for atom in self.atoms) / len(self.atoms)
        z = sum(atom.z for atom in self.atoms) / len(self.atoms)
        return x, y, z

def prepare_charged_groups(protein):
    """Identify charged functional groups in protein residues."""
    group_atoms = {"ARG": ("positive", ["NE", "CZ", "NH1", "NH2"]),
                   "LYS": ("positive", ["NZ"]),
                   "ASP": ("negative", ["CG", "OD1", "OD2"]),
                   "GLU": ("negative", ["CD", "OE1", "OE2"])}
    for chain in protein.chains:
        for residue in chain.residues:
            if residue.residue_name not in group_atoms:
                continue
            charge_type, atom_names = group_atoms[residue.residue_name]
            atoms = []
            for atom_name in atom_names:
                atom = residue.find_atom(atom_name)
                if atom is not None:
                    atoms.append(atom)
            if len(atoms) == len(atom_names):
                residue.charged_groups.append(ChargedGroup(charge_type, atoms))

def distance_between_points(point1, point2):
    """Calculate the Euclidean distance between two points."""
    dx = point1[0] - point2[0]
    dy = point1[1] - point2[1]
    dz = point1[2] - point2[2]
    return math.sqrt(dx**2 + dy**2 + dz**2)

def find_salt_bridges(residue1, residue2):
    """Find salt bridges between two residues."""
    salt_bridges = []
    for group1 in residue1.charged_groups:
        for group2 in residue2.charged_groups:
            if group1.charge_type == group2.charge_type:
                continue
            distance = distance_between_points(group1.center(), group2.center())
            if MIN_DISTANCE < distance < SALT_BRIDGE_DISTANCE:
                salt_bridges.append((group1, group2, distance))
    return salt_bridges

def refine_hydrogen_bonds(hydrogen_bonds, salt_bridges):
    """Refine hydrogen bonds using detected salt bridges."""
    filtered_bonds = []
    for residue1, residue2, bond in hydrogen_bonds:
        donor, hydrogen, acceptor, distance, angle = bond
        salt_bridge_hbond = False
        for salt_residue1, salt_residue2, bridge in salt_bridges:
            group1, group2, salt_distance = bridge
            if donor in group1.atoms and acceptor in group2.atoms:
                salt_bridge_hbond = True
            if donor in group2.atoms and acceptor in group1.atoms:
                salt_bridge_hbond = True
        if not salt_bridge_hbond:
            filtered_bonds.append((residue1, residue2, bond))
    best_bonds = {}
    for residue1, residue2, bond in filtered_bonds:
        donor, hydrogen, acceptor, distance, angle = bond
        if donor not in best_bonds or angle > best_bonds[donor][2][4]:
            best_bonds[donor] = (residue1, residue2, bond)
    return list(best_bonds.values())

class AromaticRing:
    """Represent an aromatic ring."""
    def __init__(self, atoms):
        """Initialize an aromatic ring."""
        self.atoms = atoms

    def center(self):
        """Calculate the geometric center of the ring."""
        x = sum(atom.x for atom in self.atoms) / len(self.atoms)
        y = sum(atom.y for atom in self.atoms) / len(self.atoms)
        z = sum(atom.z for atom in self.atoms) / len(self.atoms)
        return x, y, z

    def normal(self):
        """Calculate a normal vector to the ring plane."""
        atom1 = self.atoms[0]
        atom2 = self.atoms[2]
        atom3 = self.atoms[4]
        vector1 = (atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z)
        vector2 = (atom3.x - atom1.x, atom3.y - atom1.y, atom3.z - atom1.z)
        normal = (vector1[1] * vector2[2] - vector1[2] * vector2[1],
                  vector1[2] * vector2[0] - vector1[0] * vector2[2],
                  vector1[0] * vector2[1] - vector1[1] * vector2[0])
        norm = math.sqrt(sum(value**2 for value in normal))
        return tuple(value / norm for value in normal)

def prepare_aromatic_rings(protein):
    """Identify aromatic rings in protein residues."""
    ring_atoms = {"PHE": [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
                      "TYR": [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
                      "HIS": [["CG", "ND1", "CE1", "NE2", "CD2"]],
                      "TRP": [["CG", "CD1", "NE1", "CE2", "CD2"],
                              ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]]}
    for chain in protein.chains:
        for residue in chain.residues:
            if residue.residue_name not in ring_atoms:
                continue
            for atom_names in ring_atoms[residue.residue_name]:
                atoms = []
                for atom_name in atom_names:
                    atom = residue.find_atom(atom_name)
                    if atom is not None:
                        atoms.append(atom)
                if len(atoms) == len(atom_names):
                    residue.aromatic_rings.append(AromaticRing(atoms))

def angle_between_vectors(vector1, vector2):
    """Calculate the angle between two vectors in degrees."""
    dot_product = sum(a * b for a, b in zip(vector1, vector2))
    norm1 = math.sqrt(sum(a**2 for a in vector1))
    norm2 = math.sqrt(sum(a**2 for a in vector2))
    cosine = dot_product / (norm1 * norm2)
    cosine = max(-1.0, min(1.0, cosine))
    angle = math.degrees(math.acos(cosine))
    return min(angle, 180.0 - angle)

def calculate_ring_offset(center1, center2, normal):
    """Calculate the lateral offset between two ring centers."""
    vector = (center2[0] - center1[0],
              center2[1] - center1[1],
              center2[2] - center1[2])
    projection = sum(vector[i] * normal[i] for i in range(3))
    perpendicular = tuple(vector[i] - projection * normal[i] for i in range(3))
    return math.sqrt(sum(value**2 for value in perpendicular))

def find_aromatic_interactions(residue1, residue2):
    """Find pi-stacking interactions between two residues."""
    interactions = []
    for ring1 in residue1.aromatic_rings:
        for ring2 in residue2.aromatic_rings:
            center1 = ring1.center()
            center2 = ring2.center()
            distance = distance_between_points(center1, center2)
            if not MIN_DISTANCE < distance < PISTACK_DISTANCE:
                continue
            normal1 = ring1.normal()
            normal2 = ring2.normal()
            angle = angle_between_vectors(normal1, normal2)
            offset1 = calculate_ring_offset(center1, center2, normal1)
            offset2 = calculate_ring_offset(center2, center1, normal2)
            offset = min(offset1, offset2)
            interaction_type = None
            if angle < PISTACK_ANGLE_DEVIATION and offset < PISTACK_OFFSET:
                interaction_type = "parallel"
            if 90 - PISTACK_ANGLE_DEVIATION < angle < 90 + PISTACK_ANGLE_DEVIATION and offset < PISTACK_OFFSET:
                interaction_type = "T-shaped"
            if interaction_type is not None:
                interactions.append((ring1, ring2, distance, angle, offset, interaction_type))
    return interactions

def main():
    parser = OptionParser()
    parser.add_option("--chains", dest="chains", help="Two chains to compare, separated by a comma, ex: A,C")
    parser.add_option("--threshold", dest="threshold", type="float", default=6.0, help="Distance threshold in Angstroms, default: 6.0")
    options, args = parser.parse_args()
    chain_identifiers = options.chains.split(",")
    threshold = options.threshold
    if len(args) != 1:
        parser.error("Please give one PDB identifant")
    if options.chains is None:
        parser.error("Please specify 2 chain with --chains, ex: --chains A,C")
    if len(chain_identifiers) != 2:
            parser.error("--chains requires two chain identifiers separated by a comma, ex: A,C")
    pdb_id = args[0]
    file_name = download_pdb(pdb_id)
    protein = read_pdb(file_name)
    print(protein)
    if not protein.contains_hydrogen():
        print("No hydrogen atoms found. Adding hydrogens with Open Babel")
        file_name = add_hydrogens(file_name)
        protein = read_pdb(file_name)
    
    # Prepare to classification
    prepare_protein(protein, file_name)
    lys57 = protein.find_chain("C").find_residue(57, "")
    asp111 = protein.find_chain("A").find_residue(111, "")

    print("LYS57 donors:")
    for donor, hydrogen in lys57.hbond_donors:
        print("donor:", donor.atom_name, "H:", hydrogen.atom_serial_number)

    print("ASP111 acceptors:")
    for acceptor in asp111.hbond_acceptors:
        print("acceptor:", acceptor.atom_name)

    for donor, hydrogen in lys57.hbond_donors:
        for acceptor in asp111.hbond_acceptors:
            distance = donor.distance_to(acceptor)
            angle = calculate_angle(donor, hydrogen, acceptor)
            print("donor:", donor.atom_name,
                "H:", hydrogen.atom_serial_number,
                "acceptor:", acceptor.atom_name,
                "distance:", round(distance, 2),
                "angle:", round(angle, 2))

    chain1 = protein.find_chain(chain_identifiers[0])
    chain2 = protein.find_chain(chain_identifiers[1])
    for chain in protein.chains:
        print("Chain :", chain.chain_identifier, "Residues :", len(chain.residues), "Atoms :", chain.number_of_atoms())
    print("Excluded waters :", len(protein.excluded_waters))
    print("Excluded hetero molecules :", len(protein.excluded_hetero_molecules))
    print("Missing residues :", len(protein.missing_residues))
    print("Missing atoms :", len(protein.missing_atoms))
    print("Alternative locations :", len(protein.alternative_locations))
    for molecule in protein.excluded_hetero_molecules:
        print(molecule)
    if protein.contains_hydrogen():
        print("Hydrogen atoms are present in the structure.")
    else : 
        print("No hydrogen atoms were found in the structure.")
    print("Selected chains :", chain1.chain_identifier, "and", chain2.chain_identifier)
    # Identify residues in interface
    interface_pairs = find_interface_pairs(chain1, chain2, threshold)
    print('Number of interface residue pairs:', len(interface_pairs))
    for residue1, residue2, minimum_distance in interface_pairs:
        print(chain1.chain_identifier, residue1.residue_name, residue1.residue_sequence_number, "-", chain2.chain_identifier, residue2.residue_name, residue2.residue_sequence_number, "minimum distance :", round(minimum_distance,4))
    # Classification into hydrophobic contacts :
    hydrophobic_contacts = []
    for residue1, residue2, interface_distance in interface_pairs:
        hydrophobic_distance = detect_hydrophobic_contact(residue1, residue2)
        if hydrophobic_distance is not None:
            hydrophobic_contacts.append((residue1, residue2, hydrophobic_distance))
    print("Hydrophobic contacts :", len(hydrophobic_contacts))
    for residue1, residue2, distance in hydrophobic_contacts:
        print(chain1.chain_identifier, residue1.residue_name, residue1.residue_sequence_number,
            "-", chain2.chain_identifier, residue2.residue_name, residue2.residue_sequence_number,
            "hydrophobic distance :", round(distance, 2))

    # Ponts salins
    salt_bridges = []
    for residue1, residue2, minimum_distance in interface_pairs:
        for bridge in find_salt_bridges(residue1, residue2):
            salt_bridges.append((residue1, residue2, bridge))

    print("Salt bridges :", len(salt_bridges))
    for residue1, residue2, bridge in salt_bridges:
        group1, group2, distance = bridge
        print(chain1.chain_identifier, residue1.residue_name, residue1.residue_sequence_number,
            "-", chain2.chain_identifier, residue2.residue_name, residue2.residue_sequence_number,
            "salt bridge distance :", round(distance, 2))
    # Liaisons H
    hydrogen_bonds = []
    for residue1, residue2, minimum_distance in interface_pairs:
        for bond in find_hydrogen_bonds(residue1, residue2):
            hydrogen_bonds.append((residue1, residue2, bond))

    # Refine
    hydrogen_bonds = refine_hydrogen_bonds(hydrogen_bonds, salt_bridges)

    print("Hydrogen bonds :", len(hydrogen_bonds))
    for residue1, residue2, bond in hydrogen_bonds:
        donor, hydrogen, acceptor, distance, angle = bond
        print(chain1.chain_identifier, residue1.residue_name, residue1.residue_sequence_number,
              "-", chain2.chain_identifier, residue2.residue_name, residue2.residue_sequence_number,
              "donor:", donor.atom_name,
              "H:", hydrogen.atom_serial_number,
              "acceptor:", acceptor.atom_name,
              "distance:", round(distance, 2),
              "angle:", round(angle, 2))

    # Aromatic
    aromatic_interactions = []
    for residue1, residue2, minimum_distance in interface_pairs:
        for interaction in find_aromatic_interactions(residue1, residue2):
            aromatic_interactions.append((residue1, residue2, interaction))

    print("Aromatic interactions :", len(aromatic_interactions))
    for residue1, residue2, interaction in aromatic_interactions:
        ring1, ring2, distance, angle, offset, interaction_type = interaction
        print(chain1.chain_identifier, residue1.residue_name, residue1.residue_sequence_number,
            "-", chain2.chain_identifier, residue2.residue_name, residue2.residue_sequence_number,
            "type:", interaction_type,
            "distance:", round(distance, 2),
            "angle:", round(angle, 2),
            "offset:", round(offset, 2))
    
if __name__ == "__main__":
    main()