#!/bin/bash

# Update and Upgrade the System
sudo apt-get update && sudo apt-get upgrade -y

# Install Python3 and pip if they are not installed
sudo apt-get install python3 python3-pip -y

# It's a good practice to use virtual environments for Python projects
# Install virtualenv if you don't have it
pip3 install virtualenv

# Create a virtual environment for your project (optional)
virtualenv pyrosetta_env
# Activate the virtual environment
source pyrosetta_env/bin/activate

# Install pandas using pip
pip install pandas

# Placeholder for PyRosetta installation
# You need to download PyRosetta from the website after accepting the license agreement
# Assuming you've downloaded PyRosetta to your AWS instance, you can install it with:
pip install pyrosetta.whl

# Deactivate the virtual environment when you're done
deactivate
