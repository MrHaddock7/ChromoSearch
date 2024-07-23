# from Bio import PDB

# # Load the PDB files
# parser = PDB.PDBParser(QUIET=True)
# structure_blue = parser.get_structure("blue_protein", "data/Alphafold predictions/Blue/7sws.pdb")
# structure_green = parser.get_structure("green_protein", "data/Alphafold predictions/Green/Green_Protein_9941d/Green_Protein_9941d_unrelaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb")

# # Print basic information about the structures
# blue_info = structure_blue[0]["A"].child_dict.keys()
# green_info = structure_green[0]["A"].child_dict.keys()

# print(f'{blue_info}\n\n') 
# print(green_info)

from Bio.PDB import PDBParser, PDBIO, Superimposer, Selection
from Bio.PDB.PDBIO import Select

class NonHetSelect(Select):
    def accept_residue(self, residue):
        return not residue.id[0].strip()

def add_chromophore_to_mutated(mutated_structure, chromophore_structure):
    """
    Add the chromophore to the mutated structure and save the combined structure.
    """
    # Superimpose the mutated protein onto the wild-type protein to align them
    sup = Superimposer()
    sup.set_atoms(
        Selection.unfold_entities(mutated_structure[0], 'C'),  # C for Calpha atoms
        Selection.unfold_entities(chromophore_structure[0], 'C')
    )
    sup.apply(mutated_structure.get_atoms())

    # Add the chromophore to the mutated structure
    for chain in chromophore_structure[0]:
        if chain.id == 'A':  # Assuming chromophore is in chain A
            mutated_structure[0].add(chain)

    # Save the new structure with the chromophore added
    io = PDBIO()
    io.set_structure(mutated_structure)
    io.save('/mnt/data/mutated_with_chromophore.pdb', NonHetSelect())

# Load the structures
parser = PDBParser(QUIET=True)
structure_blue = parser.get_structure("blue_protein", "data/Alphafold predictions/Blue/7sws.pdb")
structure_green = parser.get_structure("green_protein", "data/Alphafold predictions/Green/Green_Protein_9941d/Green_Protein_9941d_unrelaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb")

# Extract the chromophore
chromophore = structure_blue[0]['A']['CRQ']  # Assuming the chromophore's residue name is CRQ
chromophore_structure = parser.get_structure("chromophore_structure", "data/Alphafold predictions/Blue/7sws.pdb")
chromophore_chain = [residue for residue in chromophore_structure.get_residues() if residue.get_resname() == 'CRQ']

# Add the chromophore to the mutated structure
add_chromophore_to_mutated(structure_green, chromophore_structure)

# Save and re-load the new mutated structure with the chromophore for visualization
structure_mutated_with_chromophore = parser.get_structure("mutated_with_chromophore", "data/Alphafold predictions/Green/Green_Protein_9941d/Green_Protein_9941d_unrelaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb")

# Re-analyze the interactions
green_binding_residues = get_chromophore_binding_residues(structure_mutated_with_chromophore[0], "CRQ")
green_residue_ids = [(res.get_resname(), res.get_id()[1]) for res in green_binding_residues]

print(f'Green protein key binding residues with chromophore: {green_residue_ids}')