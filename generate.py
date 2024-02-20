from pyrosetta.rosetta.protocols.docking import DockMCMProtocol
import multiprocessing
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.docking import DockMCMProtocol
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

def docking_simulation(i, pdb_path):
    # Initialize PyRosetta
    init()

    # Load the PDB file
    pose = pose_from_pdb(pdb_path)

    # Define the MoveMap
    movemap = MoveMap()
    for resi in range(165, 180):
        movemap.set_bb(resi, True)
        movemap.set_chi(resi, True)

    # Setup Minimization Mover with the MoveMap
    min_mover = MinMover()
    min_mover.movemap(movemap)
    # Further setup for the MinMover might be needed

    # Setup the docking protocol
    dock_mcm = DockMCMProtocol()
    dock_mcm.set_partners("A_B")

    # Apply the MinMover and then the Docking Protocol
    min_mover.apply(pose)
    dock_mcm.apply(pose)

    # Save the docked structure
    pose.dump_pdb(f"/mnt/c/General/Code/generatePDBS/output/docked_structure_{i}.pdb")

# The rest of your multiprocessing code remains the same

if __name__ == '__main__':
    pdb_path = "/mnt/c/General/Code/generatePDBS/sourcePDB/P24_renumber.pdb"
    num_simulations = 100

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.starmap(docking_simulation, [(i, pdb_path) for i in range(num_simulations)])
