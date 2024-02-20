import pandas as pd
import pyrosetta
from pyrosetta import Pose, get_fa_scorefxn, init
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint, BoundFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
from pyrosetta.rosetta.protocols.simple_moves import SmallMover, ShearMover


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
def add_distance_constraints(pose, res1, res2, tolerance=1.5):
    scorefxn = pyrosetta.get_fa_scorefxn()  # Make sure the scoring function includes constraints
    ca1 = AtomID(pose.residue(res1).atom_index("CA"), res1)
    ca2 = AtomID(pose.residue(res2).atom_index("N"), res2)
    distance = pose.residue(res1).xyz("CA").distance(pose.residue(res2).xyz("N"))
    
    # Adjust the function call to match the supported signature
    # Using the first constructor signature as an example
    func = BoundFunc(distance - tolerance, distance + tolerance, 1.0, "circularharmonic")

    constraint = AtomPairConstraint(ca1, ca2, func)
    pose.add_constraint(constraint)
    return scorefxn

pyrosetta.init(extra_options="-mute all -ex1 -ex2aro")

df = pd.read_csv('sequences.csv')

output_folder = r'/mnt/c/General/Code/Docking/output/'
original_pdb_path = r'/mnt/c/General/Code/Docking/input/template.pdb'
chain = 'B'
start_residue = 169
end_residue = 179

for index, row in df.iterrows():
    modified_pdb_path = f'/mnt/c/General/Code/Docking/temp/modified_structure_{index}.pdb'
    
    modify_residues_in_pdb(original_pdb_path, modified_pdb_path, row['Sequence'], start_residue, chain)

    # Load the modified PDB structure
    pose = Pose()
    pyrosetta.rosetta.core.import_pose.pose_from_file(pose, modified_pdb_path)
    print("Total number of residues in pose:", pose.total_residue())

    #scorefxn = add_distance_constraints(pose, 171, 467)
    # Initialize the FlexPepDocking protocol
   # flexpep_dock = FlexPepDockingProtocol()

    # Apply FlexPepDocking to the pose
    #flexpep_dock.apply(pose)

    # Save the docked structure
    #output_filename = f'{output_folder}/redocked_structure_{index}.pdb'
    #pose.dump_pdb(output_filename)