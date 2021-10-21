# Substituent Replacement

This codebase is an adaptation of the Pat Walter's [blog](https://patwalters.github.io/practicalcheminformatics/jupyter/chembl/2021/07/05/replace-rgroups.html) 
on using the "R-Group Replacement Database" published by [Takeuchi et al](https://www.future-science.com/doi/10.2144/fsoa-2021-0062). 
By supplying a molecule and it's core scaffold, the code will produce a series
of derivatives in which R-groups have been substituted based on transformations
in the "R-Group Replacement Database"

### Requirements
* python3.8
* RDKit 2021-03-5


### Installation Instructions

#### 1. Create a virtual environment
```bash 
python -m venv env_name
```

##### 2. Install RDKit
RDKit can be easily installed via Cyclica's [rdkit-installer](https://github.com/cyclica/rdkit-installer).
Install RDKit as follows:
```bash
git clone https://github.com/cyclica/rdkit-installer.git

./rdkit-installer/install-rdkit python=3.8 --release=Release_2021_03_5 /path/to/install/dir
```

#### 3. Install Substituent Replacement
```bash
git clone https://github.com/davidkuter/substituent_replacement.git

cd substituent_replacement

pip install -r requirements.txt
pip install -e .
```

### Usage

An example script can be found in *contrib/example/run_script*. Required inputs are:
>* **target_smiles:** SMILES of the molecule you want to perform substituent replacement on
>* **core_smarts:** SMARTS pattern of the core structure.
>* **out_folder:** Path to where results will be stored
>* **level:** The number of layers to be used from the R-Group Replacement Database

Currently, the output is a tsv file with information about the derivatives and
a folder containing SVGs of the parent and derivatives.

