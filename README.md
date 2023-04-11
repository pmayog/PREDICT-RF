
# PREDICT-RF: A ML-based algorithm for ligand binding site prediction of proteins

PREDICT-RF is a ML-based algorithm written in Python 3 for ligand binding sites (LBSs) prediction.

Determining the pocket of proteins and the residues involved in ligand interactions is a major basis in many biostructural areas such as drug design. Moreover, resolving ligand-binding residues can be advantageous in computational prognoses of different physiological and pathological traits.  


## Installation

PREDICT-RF requires some programs, modules, and packages in order to start running:

- **Python 3**: This program has been written using python 3.10. Although no other Python versions have been tested, other versions may also work.
- **Psi Blast**: Psi Blast program must be installed on your computer, taking into account that it must run in any folder and have read, write and execute permissions.
- **DSSP**: DSSP program must be installed on your computer, taking into account that it must run in any folder and have read, write and execute permissions.
- **Protein database**: In order to run psi blast, you need to have access to a protein database. However, the model was trained using **SwissProt**, so this is the recommended one.
- **tables.py**: File is located inside NAME.gz.tar folder.
- **SASA.py**: File is located inside NAME.gz.tar folder.
- **Model**: As NAME is an ML-based algorithm, in order to make predictions, it requires a model, that is passed to NAME as ‘joblib_RandForest_undersampling_model.pkl’.
- **Chimera**: Chimera is an optional program, it is not required in order to run NAME. However, the program creates as output a .cmd file, which is passed to Chimera in order to visualize the LBS-forming residues. Be aware that Chimera can read .cmd files, but its newest version, ChimeraX cannot recognize .cmd extension.
- **AAIndex**
- **Biopython**
- **Numpy**
- **Pandas**
- **Scikit-Learn** (version 1.1.3)
- **Joblib**

AAIndex, Biopython, Numpy, Pandas and Scikit-learn modules must be installed on your computer. Remember that these must be in the correct path so Python can find them. These can be easily installed using:

```bash
pip install <module>
```
For more information on how to install these programs check:
- https://github.com/amckenna41/aaindex 
- https://biopython.org/
- https://numpy.org/install/
- https://pandas.pydata.org/docs/getting_started/install.html
- https://scikit-learn.org/stable/install.html
- https://pypi.org/project/joblib/#downloads


In order to install PREDICT-RF, two downloadable .gz.tar files are available: 

- *predict-rf_self_installation.tar.gz*:

A .tar.gz file containing predict-rf.py, predict-rf.yml and other necessary files. First, untar and unzip the folder, create the conda environment with predict-rf.yml file and activate the predict-rf environment. Then, call python and run predict-rf.py file.

It can be downloaded at https://github.com/pmayog/PREDICT-RF

```bash
tar -xvf predict-rf.tar.gz
```
```bash
conda env create -f predict-rf.yml
```
```bash
conda activate project
```
```bash
python predict-rf.py -i INPUT -db DATABASE
```

*Note*: When creating project environment, check if *“~/bin/miniconda3/envs/predict-rf”* path exists in your computer. Otherwise, modify predict-rf.yml file and change the path.

- *predict-rf.tar.gz*:

A tar.gz file containing a predict-rf executable file and all its dependencies. The user has to untar and unzip the file and write the path to the executable in order to call the program (untared predict-rf.tar.gz is named “dist”).
It is recommended to make an alias of the executable file. 

It can be downloaded at https://drive.google.com/drive/folders/1D3g-f9G_Ji_eyrpA4cUR8V9SGsIVTRdN?usp=share_link. 
 
```bash
tar -xvf predict-rf.tar.gz
```
```bash
alias predict-rf="path_to_/dist/predict-rf/predict-rf"
```
```bash
predict-rf -i INPUT -db DATABASE
```

*Note*: Remember that Psiblast and DSSP programs must run as commands (as dssp and psiblast commands).
## Usage

PREDICT-RF is an algorithm developed for Linux-based operating systems that run in the terminal with the following command and options:

```bash
predict-rf.py -i INPUT -db DATABASE [-e] [-k] [-v]
```



 #### *-i, --input: INPUT* 
 INPUT must be a PDB file or a directory including at least one PDB file. If a non-PDB file is given as input, it will be skipped.

 #### *-db, --database: DATABASE* 
 DATABASE must include the path to the database used when running PsiBlast. Database selection can affect Psi Blast output. If there is no PSSM output, try to change the database.

#### *[-e], [--eval]* 
OPTIONAL. Specify the e-value used to run PsiBlast, 0.001 by default. This selection can affect Psi Blast output. If there is no PSSM output, try to change the threshold value.

#### *[-s], [--show]* 
OPTIONAL. The program will not delete output files created by PsiBlast, DSSP, and other temporary files.

#### *[-k], [--keep]*
OPTIONAL. Will keep files created by DSSP and psi blast, as well as other temporary files that will be commented later on (see Ouput section)

#### *[-v], [--verbose]*
OPTIONAL. Display in the terminal all information when running the program.


## Ouput

After running PREDICT-RF, there will appear two files: 

- **.txt file** named [protein_id]protein_binding_residues.txt

This file contains a list of all binding-related residues that have been predicted. It presents three columns per line representing the residue, the position and the chain:

![](https://live.staticflickr.com/65535/52804913049_3287c6273b_o.png)

*Note:* It must be taken into consideration that positions are an author’s reference and may differ from PDB’s standardized positions. It can be easily seen in PDB’s sequence section:

![](https://live.staticflickr.com/65535/52805096170_65ffd5fdb7_z.jpg)

- **.cmd file** named [protein_id]protein_binding_residues.cmd

This file represents a set of commands executable by Chimera. It can be opened from the terminal with the following command:

```bash
chimera [protein_id]protein_binding_residues.cmd
```
The commands included in the file can be resumed as follows:

![](https://live.staticflickr.com/65535/52804960589_e0aee8aec9_o.png)

As a result we can obtain something similar to the following image:

![](https://live.staticflickr.com/65535/52804157867_7f11daa54f_w.jpg)

*Note:* .cmd files are only supported by classic Chimera software but not by its newest version ChimeraX.

Moreover, if --keep option is used, other files will be created: 

- **.fa file** named [protein_id][chain].fa

A fasta file for each chain given a protein

![](https://live.staticflickr.com/65535/52805011559_0c71e337fa_o.png)

- **.dssp file** named [protein_id]_protein.dssp

A dssp file for each PDB file, containing secondary structures prediction per residue.

![](https://live.staticflickr.com/65535/52805161000_41fb16d873_z.jpg)

- **.pssm file** named [protein_id][chain].pssm

A pssm file for each chain given a protein containing two matrices.

![](https://live.staticflickr.com/65535/52804759266_c76379eb0f_c.jpg)

## Examples

#### Single PDB file without optional arguments

In this example, we will use ‘6uh0.pdb’ as the input file, which contains the carbonic anhydrase 2 protein in complex with SB4-202.

We will run the following command:

```bash
python predict-rf.py -i 6uh0.pdb -db ~/Documents/swissprot
```

In this case, as we are not using other options, we will not have output files from Psi Blast, DSSP, and other temporary files, and it will not show any message in the terminal except warnings.

After running this, we will get two files:

- 6uh0_binding_residues.txt
- 6uh0_binding_residues.cmd

#### Directory with PDB and other files without optional arguments

In this case, we will use as input a directory called ‘random’, in which there are 3 PDB files, 5 Python scripts, and 1 text file.

We will run the following command:

```bash
python predict-rf.py -i random/ -db ~/Documents/swissprot
```

Although -v option is not introduced, the terminal will display a message to warn the user that some files will be skipped because are not PDB files: *Skipping file. Input must contain PDB file/s.*

As a result of this, output files will just correspond to PDB files:

- 5fck_binding_residues.cmd
- 5fck_binding_residues.txt
- 6upj_binding_residues.cmd
- 6upj_binding_residues.txt
- 7std_binding_residues.cmd
- 7std_binding_residues.txt

#### Single PDB file without Psiblast output, using -k option

Now, we will use a PDB file with no Psiblast result. To do so, we will use ‘5d1d.pdb’:

```bash
python predict-rf.py -i 5d1d.pdb -db ~/Documents/swissprot -k
```
In this situation, ‘5d1d.pdb’ contains four chains (A, B, C, D) but only chains A and B obtain a PSSM matrix when running Psiblast. As a result, the terminal will show a message and chains will be skipped: *Empty PSSM... Skipping chain...*

Thus, the output files (both .txt and .cmd) will just contain chains A and B, as chains C and D are not passed to prediction functions. In this case, however, as -k option has been used, there will also appear other files:

- Fasta file for each chain: 5d1dA.fa, 5d1dB.fa, 5d1dC.fa and 5d1dD.fa
- DSSP file per PDB file: 5d1d_protein.dssp
- Psiblast PSSM file for each successful chain: 5d1dA_sprot.pssm and 5d1dB_sprot.pssm

*Note:* PSSM missing output can be due to database or threshold selection, try to modify -db and/or -e options.

#### Directory with PDB files, using -v and -e options

In this example, we are using a folder called ‘pdb’ that just contains PDB files. We will also use -v and -e options with a threshold of 0.0001, so we will execute:

```bash
python predict-rf.py -i pdb -db ~/Documents/swissprot -e 0.0001 -v
```
Using -v option will inform us on what the program is doing, and the terminal will show the following information:

*Starting with: 1sbi*

*Getting chains...*

*Done*

*Getting sequence positions...*

*Done*

*Getting sequence...*

*Done*

*Getting residues hydrophobicity...*

*Done*

*Getting residues polarity...*

*Done*

*Getting residues positive charge...*

*Done*

*Getting residues negative charge...*

*Done*

*Getting residues isoelectric point...*

*Done*

*Getting residues secondary structure...*

*Done*

*Getting residues solvent accesible surface area...*

*Done*

*Getting PSSM...*

*Warning: [psiblast] Query_1 1sbi: Composition-based score adjustment conditioned on sequence properties and unconditional composition-based score adjustment is not supported with PSSMs, resetting to default value of standard composition-based statistics*

*Done*

*Note:*  Psi blast warning will always appear (when using -v option and without -v option). Other warnings can also appear.

## Documentation

You can find a background and a scientific explanation concerning the importance of accurately identifying the ligand binding sites of proteins, as well as an extended explanation of the program creation (i.e. Material and methods section), the results obtained (i.e. Results section) and the conclusions extracted (i.e. Discussion and Conclusions sections) in the following link:

[Link Text](./MLbased_algorithm_for_LBSs_prediction.pdf)

## Authors

- [@pmayog](https://www.github.com/octokatherine)
- [@codernona](https://www.github.com/octokatherine)

