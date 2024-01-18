#!/bin/bash

# Define the name of your virtual environment
venv_name="./.venv"

# Create a virtual environment
python -m venv $venv_name

# Activate the virtual environment
source $venv_name/bin/activate

# Install dependencies from requirements.txt
pip install -r requirements.txt

# Display a message indicating the activation of the virtual environment
echo "Virtual environment $venv_name is activated. To deactivate, run 'deactivate'."
