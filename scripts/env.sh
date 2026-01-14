#!/bin/bash

# Script: Setup Genome Repair Tools (GRT) Environment
# Function: Create and configure GRT environment using mamba, install specified packages

# Set environment name
ENV_NAME="GRT"

# 1. Create new environment (if it doesn't exist)
if mamba env list | grep -q "^$ENV_NAME"; then
    echo "Environment $ENV_NAME already exists."
else
    echo "Creating new environment: $ENV_NAME"
    mamba create -n $ENV_NAME python=3.11.14 -y
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to create environment!"
        exit 1
    fi
fi

# 2. Install packages
echo "Installing packages into $ENV_NAME..."

mamba install -n $ENV_NAME \
    python=3.11.14 \
    tgsgapcloser=1.2.1 \
    flye=2.9.6 \
    nextdenovo=2.5.2 \
    verkko=2.2.1 \
    mummer4 \
    hifiasm=0.25.0 \
    bioconda::shasta=0.14.0 \
    -y

if [ $? -ne 0 ]; then
    echo "ERROR: Package installation failed!"
    exit 1
fi

# 3. Activate and show installed software
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

echo ""
echo "=========================================="
echo "Successfully installed:"
echo "=========================================="
echo "1. python=3.11.14"
echo "2. tgsgapcloser=1.2.1"
echo "3. flye=2.9.6"
echo "4. nextdenovo=2.5.2"
echo "5. verkko=2.2.1"
echo "6. mummer4"
echo "7. hifiasm=0.25.0"
echo "8. shasta=0.14.0"
echo "=========================================="
echo ""
echo "Environment '$ENV_NAME' is ready."
echo "To activate: conda activate $ENV_NAME"
