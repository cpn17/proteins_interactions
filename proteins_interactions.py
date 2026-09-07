from optparse import OptionParser
from Bio.PDB import PDBList
import math
AA_STANDARD = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
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

    def add_atom(self, atom):
        """Add an atom to the residue.

        If alternative locations exist for the same atom, keep the atom with the highest occupancy.

        Parameters
        ----------
        atom : Atom
            Atom to add to the residue.

        Returns
        -------
        None
        """
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
    """Find interface residues between 2 proteins chains."""
    interface_pairs = set()
    for residue1 in chain1.residues:
        for atom1 in residue1.atoms:
            for residue2 in chain2.residues:
                for atom2 in residue2.atoms:
                    distance = atom1.distance_to(atom2)
                    if distance <= threshold:
                        interface_pairs.add((residue1, residue2))
    return interface_pairs

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
    for residue1, residue2 in interface_pairs:
        print(chain1.chain_identifier, residue1.residue_name, residue1.residue_sequence_number, "-", chain2.chain_identifier, residue2.residue_name, residue2.residue_sequence_number)

if __name__ == "__main__":
    main()