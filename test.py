import pandas as pd
import pyrosetta
from pyrosetta import Pose, rosetta
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.loops import Loop, Loops, add_single_cutpoint_variant, remove_cutpoint_variants
from pyrosetta.rosetta.protocols.loops.loop_mover.refine import LoopMover_Refine_CCD
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree, DockMCMProtocol
from pyrosetta.toolbox import mutate_residue

# Initialize PyRosetta
pyrosetta.init()

# Load the peptide sequences from the CSV file
df = pd.read_csv('sequences.csv')

# Load the original PDB structure
original_pose = Pose()
rosetta.core.import_pose.pose_from_file(original_pose, r'/mnt/c/General/Code/Docking/input/inputchainfix.pdb')

# Loop start and end positions
start_residue = 169
end_residue = 179

for index, row in df.iterrows():
    # Create a copy of the original pose
    pose = Pose()
    pose.assign(original_pose)

    # Replace amino acids 169-179 with the current peptide
    for i, aa in enumerate(row['Sequence'], start=start_residue):
        if i > end_residue:
            break
        mutate_residue(pose, i, aa)

    # Define the loop to be remodeled with a cut point
    cut_point = (start_residue + end_residue) // 2  # Midpoint of the loop
    loop = Loop(start_residue, end_residue, cut_point)
    loops = Loops()
    loops.add_loop(loop)

    # Add cutpoint variants to the pose
    add_single_cutpoint_variant(pose, loop)

    # Set up and apply the LoopMover_Refine_CCD
    ccd_mover = LoopMover_Refine_CCD(loops)
    ccd_mover.apply(pose)

    # Remove cutpoint variants after remodeling
    remove_cutpoint_variants(pose)

    # Repack the region to optimize side-chain conformations
    scorefxn = pyrosetta.get_fa_scorefxn()
    task = pyrosetta.standard_packer_task(pose)
    task.restrict_to_repacking()
    task.temporarily_fix_everything()
    for i in range(start_residue, end_residue + 1):
        task.temporarily_set_pack_residue(i, True)
    pack_mover = PackRotamersMover(scorefxn, task)
    pack_mover.apply(pose)

    # Create a MoveMap and freeze chain A
    movemap = MoveMap()
    movemap.set_bb(False)  # Disable backbone movements for all residues
    movemap.set_chi(False)  # Disable sidechain movements for all residues

    # Enable movements for chain B residues
    for res_id in range(1, pose.total_residue() + 1):
        if pose.chain(res_id) == 'B':
            movemap.set_bb(res_id, True)
            movemap.set_chi(res_id, True)

    # Setup for docking with the MoveMap
    jump_nums = rosetta.utility.vector1_int()
    jump_nums.append(1)
    setup_foldtree(pose, "A_B", jump_nums)

    # Create a DockMCMProtocol object
    dock_mover = DockMCMProtocol()
    dock_mover.set_movemap(movemap)
    dock_mover.set_scorefxn(pyrosetta.get_fa_scorefxn())

    # Perform the docking
    dock_mover.apply(pose)

    # Specify the output folder path
    output_folder = r'/mnt/c/General/Code/Docking/output/'

    # Save the output PDB file
    output_filename = f'{output_folder}redocked_structure_{index}.pdb'
    pose.dump_pdb(output_filename)
