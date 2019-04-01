#!/usr/bin/env bash
# author: github.com/ruxi
# reproducibly create conda env

read -p "Create new conda env (y/n)?" CONT

if [ "$CONT" == "n" ]; then
    echo "exit";
else
    # user chooses to create conda env
    # prompt user for conda env name
    echo "Creating new conda environment, choose name"
    read input_variable
    echo "Name $input_variable was chosen";

    # Create environment.yml or not
    if [ ! -f environment.yml ]; then
        echo "File 'environment.yml' not available. Installing base packages."
        conda create --name $input_variable python=3 bedtools samtools picard
        conda install -n $input_variable -c r r
        conda install -n $input_variable -c bioconda bioconductor-rtracklayer
    else
        echo "Using 'environment.yml' to create environment"
        conda env create

    fi
    echo "to exit: source deactivate"
fi
