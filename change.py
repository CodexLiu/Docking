def change_chain_in_pdb(pdb_path, new_pdb_path, start_residue, end_residue, new_chain):
    with open(pdb_path, 'r') as file:
        lines = file.readlines()

    with open(new_pdb_path, 'w') as file:
        for line in lines:
            # Check if line starts with ATOM or HETATM (which contains residue information)
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract residue number and convert it to integer
                res_num = int(line[22:26])

                # If residue number is in the specified range, change the chain
                if start_residue <= res_num <= end_residue:
                    line = line[:21] + new_chain + line[22:]

            file.write(line)

# Usage
pdb_input_path = '/mnt/c/General/Code/Docking/input/P24_renumber.pdb'  # Path to your original PDB file
pdb_output_path = 'inputchainfix.pdb'   # Path for the modified PDB file
change_chain_in_pdb(pdb_input_path, pdb_output_path, 169, 179, 'B')
