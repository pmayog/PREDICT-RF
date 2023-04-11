import sys
import os
import requests
import math
import pandas as pd
import numpy as np
import urllib.request
from tables import dssp_structures
from tables import amino_acid_name
from Bio.PDB import PDBParser
from SASA import ShrakeRupley
from aaindex import aaindex1


def parseArgs():
    '''
    Parse command line arguments
    '''
    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Compute ligand binding site residues within a PDB.')

        parser.add_argument('-i',
                            '--input',
                            action='store',
                            required=True,
                            help='File must be a PDB')

        parser.add_argument('-db',
                            '--database',
                            action='store',
                            required=True,
                            help='Specify the database used when running blastp and psiblast')

    except:
            print ("An exception occurred with argument parsing. Check your provided options.")

    return parser.parse_args()

def get_pdb (prefix):
  '''
  Function that gets PDB, taking a prefix as an input.
  '''
  print('Getting PDB...')
  pdb_file=prefix + ".pdb"
  urllib.request.urlretrieve(
    f'http://files.rcsb.org/download/{prefix}.pdb', pdb_file)
  print('Done')



def merge_residue_and_ligand (prefix, pocket, ligand, protein):
  """
  Function that merges in a file pocket and ligand lines containing coordinates
  """

  all_results = []
  ligand_path="Done/"+ligand
  pocket_path="Done/"+pocket
  protein_path="Done/"+protein
  output_file= prefix+'_merge.pdb'
  merged_file = open(output_file, 'w')

  with open(ligand_path, 'r') as lig:
    for ligand_line in lig:
      ligand_split = ligand_line.split()
      if len(ligand_split) > 8 and '###' not in ligand_split:
        merged_file.write(ligand_line)

  with open(pocket_path, 'r') as pock:
    for pocket_line in pock:
      if pocket_line.startswith('ATOM'):
        merged_file.write(pocket_line)

  with open(protein_path, 'r') as prot:
    for prot_line in prot:
      specific_residue = "%s %s" % (prot_line[17:20], prot_line[23:26])
      if prot_line.startswith('ATOM') and specific_residue not in all_results:
        all_results.append(specific_residue)


  return (all_results)

def get_sequence_position (pdb):
  '''
  Function that gets the sequence from pdb input file.
  '''
  print ('Getting sequence positions...')
  all_results = []
  with open(pdb, 'r') as prot:
    for prot_line in prot:
      specific_residue = "%s %s" % (prot_line[17:20], prot_line[23:26])
      if prot_line.startswith('ATOM') and specific_residue not in all_results:
        all_results.append(specific_residue)

  print('Done')
  return (all_results)

def get_sequence (sequence_position, prefix):
  '''
  Function that gets a list from get_sequence_position function as input.
  '''
  print ('Getting sequence...')
  count=0
  fasta=prefix + ".fa"
  title=">"+ prefix+"\n"
  fasta_file = open(fasta, 'w')
  inp=title
  sequence=""
  for residue in sequence_position:
    res=residue[:3]
    for key, value in amino_acid_name.items():
      if res == key:
        sequence=sequence+value  


  ## Fasta file creation
  for res in sequence:
    count+=1
    if count%60 !=0:
      inp=inp+(res)
    else:
      inp=inp+res+"\n"
  inp=inp+"\n"
  fasta_file.write(inp)
  fasta_file.close()


  print('Done')
  return (sequence)


def residue_hydrophobicity (sequence):
  '''
  Function that retrieves residues hydrophobicity from aaindex1.
  '''
  print ('Getting residues hydrophobicity...')
  hydrophobicity=[]
  hydrophobicity_values = aaindex1['PRAM900101'].values
  for aa in sequence:
    for residue, values in hydrophobicity_values.items():
      if aa == residue:
        hydrophobicity.append(values)
  
  print('Done')
  return (hydrophobicity)


def residue_polarity (sequence):
  '''
  Function that retrieves residues polarity from aaindex1.
  '''
  print ('Getting residues polarity...')
  polarity=[]
  mean_polarity_values = aaindex1['RADA880108'].values
  for aa in sequence:
    for residue, values in mean_polarity_values.items():
      if aa == residue:
        polarity.append(values)

  print('Done')
  return (polarity)


def residue_negative (sequence):
  '''
  Function that retrieves residues negativity from aaindex1.
  '''
  print ('Getting residues negative charge...')
  negative=[]
  negative_charge_values = aaindex1['FAUJ880112'].values
  for aa in sequence:
    for residue, values in negative_charge_values.items():
      if aa == residue:
        negative.append(values)

  print('Done')
  return(negative)


def residue_positive (sequence):
  '''
  Function that retrieves residues positive charge from aaindex1.
  '''
  print ('Getting residues positive charge...')
  positive=[]
  positive_charge_values = aaindex1['FAUJ880111'].values
  for aa in sequence:
    for residue, values in positive_charge_values.items():
      if aa == residue:
        positive.append(values)

  print('Done')
  return (positive)
  
def residue_isoelectric_point (sequence):
  '''
  Function that retrieves residues isoelectric point from aaindex1.
  '''
  print ('Getting residues isoelectric point...')
  isoelectric_point=[]
  isoelectric_point_values = aaindex1['ZIMJ680104'].values
  for aa in sequence:
    for residue, values in isoelectric_point_values.items():
      if aa == residue:
        isoelectric_point.append(values)

  print('Done')
  return (isoelectric_point)

def distance_residue_and_ligand (merged, sequence_residues):
  '''
  Function that takes residues with distance < 4 with respect to the ligand 
  '''

  action=[]
  res = 0
  ligand_coords = np.empty([0, 0])
  close_res = []
  already_known = set()
  with open(merged, 'r') as merg:
    for merged_line in merg:
      if 'ATOM' not in merged_line:
        x1 = float(merged_line[18:27])
        y1 = float(merged_line[27:36])
        z1 = float(merged_line[37:46])
        ligand_coords = np.append(ligand_coords, [x1, y1, z1])

      else:
        x2 = float(merged_line[31:38])
        y2 = float(merged_line[39:46])
        z2 = float(merged_line[47:54])
        coords = np.array((x2, y2, z2))

        if merged_line[23:26] == res:
          for i in range(0, len(ligand_coords) - 3, 3):
            if res not in already_known:
              distance = np.linalg.norm(
                [ligand_coords[i], ligand_coords[i + 1], ligand_coords[i +
                                                                       2]] -
                coords)
              if distance < 4:
                specific_residue = "%s %s" % (merged_line[17:20],
                                              merged_line[23:26])
                close_res.append(specific_residue)
                already_known.add(merged_line[23:26])

        else:
          res = merged_line[23:26]
          already_known = set()
          for i in range(0, len(ligand_coords) - 3, 3):
            if res not in already_known:
              distance = np.linalg.norm(
                [ligand_coords[i], ligand_coords[i + 1], ligand_coords[i +
                                                                       2]] -
                coords)
              if distance < 4:
                specific_residue = "%s %s" % (merged_line[17:20],
                                              merged_line[23:26])
                close_res.append(specific_residue)
                already_known.add(merged_line[23:26])


  for res in sequence_residues:
    if res in close_res:
      #print('Binding',res)
      action.append("Binding")
    else:
      #print('Non-inding',res)
      action.append("Non-binding")
  
  #print('All residues:',len(all_results))
  #print('Activity:',len(action) ,action)

  return (action)



def extract_secondary_structure(prefix, seq_len):
  '''
  Function that runs dssp and extracts secondary structure
  '''
  print ('Getting secondary structures...')
  secondary = []
  seq_started = False
  dssp_file_name = prefix + ".dssp"
  pdb_file_name = prefix +".pdb"

  urllib.request.urlretrieve(
    f'http://files.rcsb.org/download/{pdb_file_name}')
  
  os.system(f'dssp -i {pdb_file_name} -o {dssp_file_name}')

  count=0
  with open(dssp_file_name, 'r') as out_dssp:
    for dssp_line in out_dssp:
      split_line = dssp_line.split()
      if seq_started is True:
        count+=1       
        if count <= seq_len:
          current_structure = dssp_line[16]
          if dssp_line[16] not in dssp_structures:
            current_structure = "-"
          secondary.append(current_structure)
      elif split_line[0] == "#":
        seq_started = True

  for i in range(len(secondary)):
      secondary[i] = secondary[i].replace('-', 'C')
      secondary[i] = secondary[i].replace('I', 'C')
      secondary[i] = secondary[i].replace('T', 'C')
      secondary[i] = secondary[i].replace('S', 'C')
      secondary[i] = secondary[i].replace('G', 'H')
      secondary[i] = secondary[i].replace('B', 'E')

  print('Done')
  return secondary


# The DSSP codes for secondary structure used here are:
# =====     ====
# Code      Structure
# =====     ====
# H         Alpha helix (4-12)
# B         Isolated beta-bridge residue
# E         Strand
# G         3-10 helix
# I         Pi helix
# T         Turn
# S         Bend
# -         None
# =====     ====
#http://www.compbio.dundee.ac.uk/jpred/references/prot_html/node17.html
# S-Scheme 2: H,G->H ; E,B->E ; I,T,S->C



def SASA(pdb, seq_len):
  '''
  #Function that calculates solvent accessible surface areas
  '''
  print ('Getting SASA...')
  sasa = []
  #protein_id = pdb_file_name.split("_")[0]

  p = PDBParser(QUIET=1)
  struct = p.get_structure(pdb, pdb)
  sr = ShrakeRupley()
  for model in struct:
    for chain in model:
      for residue in chain:
        sr.compute(residue, level="R")
        sasa.append(round(residue.sasa, 2))
  
  print('Done')
  return (sasa[0:seq_len])


def pssm_calculation (fasta_file, path_to_db, prefix):
  ''' Creates a pssm for a given protein using a fasta file, a given database and a prefix as input and creates a pssm matrix file as output.
  '''
  print('Getting PSSM...')

  A=[]
  R=[]
  N=[]
  D=[]
  C=[]
  Q=[]
  E=[]
  G=[]
  H=[]
  I=[]
  L=[]
  K=[]
  M=[]
  F=[]
  P=[]
  S=[]
  T=[]
  W=[]
  Y=[]
  V=[]
  entropy=[]

  out_file = prefix + "_sprot.pssm"
  print(fasta_file)
  os.system("psiblast -query " + fasta_file + " -evalue 0.001 -num_iterations 3 -out_ascii_pssm " +
            out_file + " -db " + path_to_db+" -out pssm.out")
  
  if os.path.isfile(out_file):
    with open (out_file, 'r') as pssm_file:
      for row in pssm_file:
        row.strip()
        split_row = row.split()
        single_entropy=0
        if 'Last' not in split_row and 'Lambda' not in split_row and 'Standard' not in split_row and 'PSI' not in split_row and len(split_row)>40:
          A.append(split_row[2])
          R.append(split_row[3])
          N.append(split_row[4])
          D.append(split_row[5])
          C.append(split_row[6])
          Q.append(split_row[7])
          E.append(split_row[8])
          G.append(split_row[9])
          H.append(split_row[10])
          I.append(split_row[11])
          L.append(split_row[12])
          K.append(split_row[13])
          M.append(split_row[14])
          F.append(split_row[15])
          P.append(split_row[16])
          S.append(split_row[17])
          T.append(split_row[18])
          W.append(split_row[19])
          Y.append(split_row[20])
          V.append(split_row[21])
          for i in range(22,42):
            if int(split_row[i]) >0:
              single_entropy = int(split_row[i])*math.log10(int(split_row[i])) + single_entropy
          entropy.append(single_entropy)


    print('Done')
    return(A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V, entropy)
  else:
    return(False)




def main():
  '''Compute LBS-forming residues within a PDB file'''
  import pandas as pd

  ##Parse arguments  
  args = parseArgs()
  pdb = args.input
  path_to_db = args.database
  #delete = args.delete
  no_result_pssm=[]
  files=set()


  if os.path.isdir(pdb):
    for pdb_file in os.listdir(pdb):
      prefix = pdb_file.split("_")[0]
      files.add(prefix)

    for file in list(files):
      
      print('Starting with',file)
      pdb_last= file+".pdb"
      fasta_file = file + ".fa"
      pocket = file + "_pocket.pdb"
      ligand = file + "_ligand.mol2"
      protein = file + "_protein.pdb"
      merged = file + "_merge.pdb"
      output=file+".csv"

      get_pdb (file)
      merge_residue_and_ligand (file, pocket, ligand, protein)
      sequence_residues = get_sequence_position (pdb_last)
      action = distance_residue_and_ligand (merged, sequence_residues)
      sequence = get_sequence (sequence_residues, file)
      hydrophobicity = residue_hydrophobicity (sequence)
      polarity = residue_polarity (sequence)
      positive = residue_positive (sequence)
      negative = residue_negative (sequence)
      isoelectric_point = residue_isoelectric_point (sequence)
      secondary_structure = extract_secondary_structure(file,len(sequence))
      sasa = SASA(pdb_last, len(sequence))
      pssm = pssm_calculation (fasta_file, path_to_db, file)

      if pssm == False:
        no_result_pssm.append(file)
      else:

        ##Creating the dataframe with the previous lists
        table = {
          'residue': sequence_residues,
          'hydrophobicity': hydrophobicity,
          'polarity': polarity,
          'positive': positive,
          'negative': negative,
          'isoelectric point': isoelectric_point,
          'secondary structure': secondary_structure,
          'sasa': sasa,
          'A':pssm[0],
          'R':pssm[1],
          'N':pssm[2],
          'D':pssm[3],
          'C':pssm[4],
          'Q':pssm[5],
          'E':pssm[6],
          'G':pssm[7],
          'H':pssm[8],
          'I':pssm[9],
          'L':pssm[10],
          'K':pssm[11],
          'M':pssm[12],
          'F':pssm[13],
          'P':pssm[14],
          'S':pssm[15],
          'T':pssm[16],
          'W':pssm[17],
          'Y':pssm[18],
          'V':pssm[19],
          'entropy':pssm[20],
          'activity': action
        }
          
        result=pd.DataFrame(data=table)
        result.to_csv(output, index=False)

    print(no_result_pssm)




##Ara suposo que s'hauria d'aplicar el model de ML i fer un output per dir quins són els residus que formen el LBS.


if __name__ == '__main__':
    main()