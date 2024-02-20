import pandas as pd
import pyrosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint, BoundFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import IncludeCurrent, RestrictToRepacking, NoRepackDisulfides
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover


def modify_residues_in_pdb(input_pdb, output_pdb, sequence, start_res, chain_id):
    aa_codes = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
                'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
    
    with open(input_pdb, 'r') as file:
        lines = file.readlines()

    with open(output_pdb, 'w') as file:
        for line in lines:
            if line.startswith("ATOM") and line[21] == chain_id:
                res_num = int(line[22:26].strip())
                if start_res <= res_num < start_res + len(sequence):
                    aa_three_letter = aa_codes[sequence[res_num - start_res]]
                    line = line[:17] + aa_three_letter + line[20:]
            file.write(line)

# Assuming pose is your loaded Pose object
def add_distance_constraints(pose, substrate_res, serine_res, tolerance=1.5):
    # Use the standard full-atom score function with constraints
    scorefxn = pyrosetta.get_fa_scorefxn()
    
    # Define the atoms for the constraint
    og_atom_id = AtomID(pose.residue(serine_res).atom_index("OG"), serine_res)
    c_atom_id = AtomID(pose.residue(substrate_res).atom_index("C"), substrate_res)
    
    # Calculate the current distance to set as the target for the constraint
    distance = pose.residue(serine_res).xyz("OG").distance(pose.residue(substrate_res).xyz("C"))
    func = BoundFunc(distance - tolerance, distance + tolerance, 1.0, "circularharmonic")
    
    # Create and add the constraint
    constraint = AtomPairConstraint(og_atom_id, c_atom_id, func)
    pose.add_constraint(constraint)
    
    return pose, scorefxn



pyrosetta.init()
#extra_options="-mute all -ex1 -ex2aro"

df = pd.read_csv('sequences.csv')

output_folder = r'/mnt/c/General/Code/Docking/output/'
original_pdb_path = r'/mnt/c/General/Code/Docking/input/inputchainfix.pdb'
chain = 'B'
start_residue = 169
end_residue = 179

for index, row in df.iterrows():
    modified_pdb_path = f'/mnt/c/General/Code/Docking/temp/modified_structure_{index}.pdb'
    
    modify_residues_in_pdb(original_pdb_path, modified_pdb_path, row['Sequence'], start_residue, chain)

    # Load the modified PDB structure
    pose = Pose()
    pyrosetta.rosetta.core.import_pose.pose_from_file(pose, modified_pdb_path)
    pose, scorefxn = add_distance_constraints(pose, substrate_res=2, serine_res=286, tolerance=1.5)

 # Pre-docking setup: TaskFactory for repacking before docking
    tf = TaskFactory()
    tf.push_back(IncludeCurrent())  # Include the current rotamers of the side chains
    tf.push_back(RestrictToRepacking())  # Restrict to repacking to avoid introducing steric clashes
    tf.push_back(NoRepackDisulfides())  # Avoid repacking disulfides

    # Pre-docking repacking: Use a PackRotamersMover for repacking the interface
    pack_mover = PackRotamersMover()
    pack_mover.task_factory(tf)
    pack_mover.score_function(pyrosetta.get_fa_scorefxn())
    pack_mover.apply(pose)

    # Initialize and configure FlexPepDockingProtocol
    flexpep_dock = FlexPepDockingProtocol()

    # Apply FlexPepDocking to the pose
    flexpep_dock.apply(pose)

    # Optional: Post-docking relaxation for further optimization
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.apply(pose)

    # Save the docked and optionally relaxed structure
    output_filename = f'{output_folder}/redocked_structure_{index}.pdb'
    pose.dump_pdb(output_filename)