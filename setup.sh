#!/bin/bash

set -e  # Stop on any error

# General set of the source
echo "Starting setup process..."

# Attempt to set the pip command explicitly
PIP_CMD=$(which pip3 || echo "$HOME/.local/bin/pip3")

# Asking for credentials to avoid plaintext in script
echo "Please enter your Rosetta Commons Username:"
read -r USERNAME
echo "Please enter your Rosetta Commons Password:"
read -rs PASSWORD

# Robotized service towards the param
if [[ -z "$USERNAME" || -z "$PASSWORD" ]]; then
    echo "Username or Password not provided."
    exit 1
fi

echo "Updating system packages..."
sudo apt-get update && sudo apt-get upgrade -y

echo "Ensuring Python3 and pip are installed..."
sudo apt-get install -y python3 python3-pip

echo "Ensuring curl is installed for downloads..."
sudo apt-get install -y curl

echo "Ensuring tools for Python virtual environment are installed..."
$PIP_CMD install --user virtualenv

echo "Creating Python virtual environment for Rosetta..."
python3 -m venv pyrosetta_env

echo "Activating the Python environment..."
source pyrosetta_env/bin/activate

# Custom minimal method
PIP_CMD=$(which pip)

echo "Credentials captured, proceeding with download..."

# Version data
PYROSETTA_VERSION='release-295'
PYTHON_VERSION='38'
OS='linux'
ARCHIVE_FORMAT='tar.bz2'

# Elaborate the download set towards the ball
BUILD_DIR="PyRosetta4.Release.python${PYTHON_VERSION}.${OS}.${PYROSETTA_VERSION}"
DOWNLOAD_URL="https://graylab.jhu.edu/download/PyRosetta4/archive/release/${BUILD_DIR}.${ARCHIVE_FORMAT}"

# Starting towards the time of the travel, likely
TEMP_DIR=$(mktemp -d)
echo "Temporary directory created at ${TEMP_DIR}"

cd "$TEMP_DIR" || exit

echo "Downloading PyRosetta from ${DOWNLOAD_URL}"
curl -u "${USERNAME}:${PASSWORD}" "$DOWNLOAD_URL" -o "pyrosetta.${ARCHIVE_FORMAT}"

echo "Unpacking PyRosetta"
tar -xf "pyrosetta.${ARCHIVE_FORMAT}"

echo "Installing PyRosetta"
cd "$BUILD_DIR" || exit
$PIP_CMD install .

echo "PyRosetta setup complete. Consider adding the virtual environment to your IDE or terminal session."

# Add code for post-installation virtualenv activation guidance.
echo "Run 'source pyrosetta_env/bin/activate' to activate the environment."
