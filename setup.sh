#!/bin/bash

# Update and upgrade the system
sudo apt-get update && sudo apt-get upgrade -y

# Install Python3 and pip if they are not installed
sudo apt-get install python3 python3-pip -y

# Install dependencies for PyRosetta
sudo apt-get install -y curl

# Optional: Setup virtual environment for PyRosetta
# Install virtualenv if not already installed
pip3 install virtualenv

# Create a virtual environment named 'pyrosetta_env'
virtualenv pyrosetta_env

# Activate the virtual environment
source pyrosetta_env/bin/activate

# Define PyRosetta credentials and version
USERNAME="levinthal"
PASSWORD="paradox"
PYROSETTA_VERSION='release-295'  # Example version, adjust as needed
PYTHON_VERSION='38'  # Adjust based on your Python version (e.g., 37, 38, 39, 310)
OS='linux'  # Adjust if using a different OS (e.g., mac, ubuntu)
ARCHIVE_FORMAT='tar.bz2'

# Determine the correct PyRosetta build directory name
BUILD_DIR="PyRosetta4.Release.python${PYTHON_VERSION}.${OS}.${PYROSETTA_VERSION}"

# Construct the download URL
DOWNLOAD_URL="https://graylab.jhu.edu/download/PyRosetta4/archive/release/${BUILD_DIR}.${ARCHIVE_FORMAT}"

# Create a temporary directory for the PyRosetta download
TEMP_DIR=$(mktemp -d)
echo "Temporary directory created at ${TEMP_DIR}"

# Navigate to the temporary directory
cd ${TEMP_DIR}

# Download PyRosetta
echo "Downloading PyRosetta from ${DOWNLOAD_URL}"
curl -u ${USERNAME}:${PASSWORD} ${DOWNLOAD_URL} -o pyrosetta.${ARCHIVE_FORMAT}

# Unpack the archive
echo "Unpacking PyRosetta"
tar -xf pyrosetta.${ARCHIVE_FORMAT}

# Install PyRosetta
echo "Installing PyRosetta"
cd ${BUILD_DIR}
pip install .

# Cleanup
echo "Cleaning up temporary files"
cd ..
rm -rf ${TEMP_DIR}

echo "PyRosetta installation completed successfully."