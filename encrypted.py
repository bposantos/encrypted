#!/usr/bin/python3
# -*-coding: utf-8 -*-

from __future__ import division
# Goal: Make a fragment list of peptides from proteins
# Author: Bruno Santos
# Comment: Viva la revolucion.
import os
import argparse
import csv
import sys
import math
import pandas as pd
import time
import pickle
import seaborn as sns
import shutil
import matplotlib.pyplot as plt
from classes import *
from sharedstuff import *

start = time.perf_counter()

#----------------------------------------------------------------------------------------------
# Proteolysis
#----------------------------------------------------------------------------------------------

if argstype == 'protein':

	for arquivo in sequencia:

		if not args.cleavage:

			import not_cleavage

		if args.cleavage:

			import cleavage

if argstype == 'peptide':

	import peptide_function

##----------------------------------------------------------------------------------------------
## Exporting CSV file
##----------------------------------------------------------------------------------------------

with open('file_original.csv', mode='w') as csv_file:
    fieldnames = list(dicio_frags3.keys())
    csvwriter = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for session in dicio_frags3:
        for item in dicio_frags3[session]:
            csvwriter.writerow([session, item, dicio_frags3[session][item]])

##----------------------------------------------------------------------------------------------
## Removing "fragments" title
##----------------------------------------------------------------------------------------------

bad_words = ['Fragments']

with open('file_original.csv') as oldfile, open('newfile.csv', 'w') as newfile:
    for line in oldfile:
        if not any(bad_word in line for bad_word in bad_words):
            newfile.write(line)

with open('newfile.csv') as file:
	df = pd.read_csv('newfile.csv')

df1 = df.columns
#print(df)
df.columns = ['Peptides','Features', 'Values']
df.loc[-1] = df1
df.index = df.index +1
df = df.sort_index()
df_reorder = df[['Features','Values', 'Peptides']] # rearrange column here
df_sort = df_reorder.sort_values('Peptides')
df2=pd.DataFrame(df_sort.groupby('Features')['Peptides'].apply(list).to_dict())
df3=pd.DataFrame(df_sort.groupby('Features')['Values'].apply(list).to_dict())
df3['Peptides'] = df2[u'Amidated_Mass']

########################### Naming the peptides, in case the --type = peptide ##########################################

lista_column_names_df = []
if args.type == 'peptide':
	for idx, row in df3['Peptides'].items():
		if row in d.values():
			lista_column_names_df.append(list(d.keys())[list(d.values()).index(row)])
		else:
			lista_column_names_df.append(None)
	#print('lista',lista_column_names_df)

	df3['Name'] = lista_column_names_df

########################################################################################################################

df3.to_csv('dataframe'+args.number+'.csv', index=True)

##----------------------------------------------------------------------------------------------
## Removing unecessary files
##----------------------------------------------------------------------------------------------

os.remove('file_original.csv')
os.remove('newfile.csv')

##----------------------------------------------------------------------------------------------
## MS/MS table
##----------------------------------------------------------------------------------------------

msms.b_series(df['Peptides'])

msms.y_series(df['Peptides'])

msms.a_series(df['Peptides'])

with open('ms_ms_series'+args.number+'.csv', mode='w') as csv_file:
    fieldnames = list(dic3.keys())
    csvwriter = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for session in dic3:
        for item in dic3[session]:
            csvwriter.writerow([session, item, dic3[session][item]])

##----------------------------------------------------------------------------------------------
## Classifier
##----------------------------------------------------------------------------------------------

df4 = pd.read_csv('dataframe'+args.number+'.csv', index_col=0)

df4['Amidated_Mass']=round(df4['Amidated_Mass']/df4['Amidated_Mass'].max(),2)
df4['Acidic_Residues_(%)']=0.00
df4['Basic_Residues_(%)']=round(df4['Basic_Residues_(%)']/df4['Basic_Residues_(%)'].max(),2)
df4['Charged_Residues_(%)']=round(df4['Charged_Residues_(%)']/df4['Charged_Residues_(%)'].max(),2)
df4['Hydrophobic_Moment']=round(df4['Hydrophobic_Moment']/df4['Hydrophobic_Moment'].max(),2)
df4['Hydrophobicity']=round(df4['Hydrophobicity']/df4['Hydrophobicity'].max(),2)
df4['Isoeletric_Point']=round(df4['Isoeletric_Point']/df4['Isoeletric_Point'].max(),2)
df4['Nonpolar_Residues_(%)']=round(df4['Nonpolar_Residues_(%)']/df4['Nonpolar_Residues_(%)'].max(),2)
df4['Polar_Residues_(%)']=round(df4['Polar_Residues_(%)']/df4['Polar_Residues_(%)'].max(),2)
df4['Polar_Residues_(%)']=round(df4['Polar_Residues_(%)']/df4['Polar_Residues_(%)'].max(),2)
df4['aromatic_rsn'] = round(df4['aromatic_rsn']/df4['Amidated_Mass'].max(),2)
df4['cf_alpha'] = round(df4['cf_alpha']/df4['cf_alpha'].max(),2)
df4['cf_beta'] = round(df4['cf_beta']/df4['cf_beta'].max(),2)
df4['cf_turn'] = round(df4['cf_turn']/df4['cf_turn'].max(),2)
df4['volume'] = round(df4['volume']/df4['volume'].max(),2)
df4['pct_buried'] = round(df4['pct_buried']/df4['pct_buried'].max(),2)
df4['solvent_accessibility_main_chain'] = round(df4['solvent_accessibility_main_chain']/df4['solvent_accessibility_main_chain'].max(),2)
df4['solvent_accessibility_nonpolar'] = round(df4['solvent_accessibility_nonpolar']/df4['solvent_accessibility_nonpolar'].max(),2)
df4['solvent_accessibility_polar'] = round(df4['solvent_accessibility_polar']/df4['solvent_accessibility_polar'].max(),2)

modelo = ''

if args.activity == "antibacterial":
	modelo = "model_antibacterial"
	X = df4.drop(columns=["Monoisotopic_Mass","Amidated_Mass","Peptides"])
	loaded_model = pickle.load(open(modelo,"rb"))
	loaded_model.predict(X)
	#sns.catplot(x=loaded_model.predict(X), kind="count", data=df)
	#plt.show()
	df4['Predicted_Peptides'] = loaded_model.predict(X)
	df4.to_csv("predicted_peptides"+args.number+".csv")
else:
	None

##----------------------------------------------------------------------------------------------
## End
##----------------------------------------------------------------------------------------------

finish = time.perf_counter()
print(f'Finished in {round(finish-start,2)} second(s)')

# Remove the temporary directory
try:
    shutil.rmtree(os.path.join(os.getcwd(),"temp/"))
except OSError as e:
    print("Error: %s : %s" % (dir_path, e.strerror))