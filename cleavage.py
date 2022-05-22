from __future__ import division
# Goal: Make a fragment list of peptides from proteins
# Author: Bruno Santos
# Comment: Viva la revolucion.

import argparse
import csv
import sys
import math
import pandas as pd
import time
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from classes import *
from sharedstuff import *

total = []

for arquivo in sequencia:

	def argCproteinase(sequence):
		for i in sequence:
			for elemento in argC:
				total.append(i.split(elemento))
		return total

	def trypsinKproteinase(sequence):
		for i in sequence:
			for elemento in trypsinK:
				total.append(i.split(elemento))
		return total

	def trypsinGKproteinase(sequence):
		for i in sequence:
			for elemento in trypsinGK:
				total.append(i.split(elemento))
		return total

	def trypsinRproteinase(sequence):
		for i in sequence:
			for elemento in trypsinR:
				total.append(i.split(elemento))
		return total

	def chymotrypsin_high_proteinase(sequence):
		for i in sequence:
			for elemento in chymotrypsin_high_specificity:
				total.append(i.split(elemento))
		return total

	def chymotrypsin_low_proteinase(sequence):
		for i in sequence:
			for elemento in chymotrypsin_low_specificity:
				total.append(i.split(elemento))
		return total

	def pepsin_ph_2_proteinase(sequence):
		for i in sequence:
			for elemento in pepsin_ph_higher_2:
				total.append(i.split(elemento))
		return total

	if args.cleavage == '1':
		total = argCproteinase(sequencia)
	elif args.cleavage == '2':
		total = trypsinKproteinase(sequencia)
	elif args.cleavage == '3':
		total = trypsinRproteinase(sequencia)
	elif args.cleavage == '4':
		total = chymotrypsin_high_proteinase(sequencia)
	elif args.cleavage == '5':
		total = chymotrypsin_low_proteinase(sequencia)
	elif args.cleavage == '6':
		total = pepsin_ph_2_proteinase(sequencia)
	elif args.cleavage == '7':
		total = trypsinGKproteinase(sequencia)
	else:
		print("No Enzyme Choosen")

	seq_string = []
	for i in list(total):
		seq_string += i

	while '' in seq_string:
		seq_string.remove('')

	#----------------------------------------------------------------------------------------------
	# Functions
	#----------------------------------------------------------------------------------------------

	setting_protein.nome_proteina(arquivo) #cada fragmento é o valor da chave que é o nome do fragmento. ex: {fragmento1: KKKIKGED, fragmento2: KKKIKG}

	setting_protein.nested_dic(seq_string) #um outro dicionário (só que dicionário de dicionários) em que a sequencia é a chave e o valor será preenchido com os parâmentros. ex: {KKKIKGED: {}, KKKIKG: {}}

	for i in range(0, len(seq_string)):
		exec("dic_frag%d = %s" % (i + 1, {}));
	#----------------------------------------------------------------------------------------------
	# Dicionários que atribuirão as funções
	#----------------------------------------------------------------------------------------------

	##----------------------------------------------------------------------------------------------
	## Calcular ponto isoelétrico
	##----------------------------------------------------------------------------------------------

	dic_fragmentos_gk ={}

	setting_protein_enzymes.pka_residues(seq_string)

	setting_protein_enzymes.resulting_charge_ph0(seq_string)

	setting_protein_enzymes.charged_residues_quantity(seq_string)

	setting_protein_enzymes.CR_list(seq_string)

	setting_protein_enzymes.pI(seq_string) #neste ponto o dicionário de dicionários contém os parametros do pI para cada fragmento


	##----------------------------------------------------------------------------------------------
	## Calcular massa monoisotópica
	##----------------------------------------------------------------------------------------------

	setting_protein.monoisotopic_mass(seq_string)

	setting_protein_enzymes.massa_amidada(seq_string)

	setting_protein_enzymes.m_h(seq_string)


	##----------------------------------------------------------------------------------------------
	## Calcular Hidrofobicidade - GRAVY
	##----------------------------------------------------------------------------------------------

	setting_protein_enzymes.hydropathicity(seq_string) #neste ponto o dicionário de dicionários contém os parametros do pI, massa ou hidrofobicidade para cada fragmento

	setting_protein.hydrophobic_moment_2(seq_string)

	setting_protein_enzymes.aa_percentage(seq_string)

	##----------------------------------------------------------------------------------------------
	## Parte 2 - Redução dos fragmentos
	##----------------------------------------------------------------------------------------------

	setting_protein_enzymes.sistematic_reduction_deep(seq_string)

	setting_protein.nested_dic_selection(seq_string)

	setting_protein.monoisotopic_mass_selection(seq_string)

	LOWEST_MASS, HIGHEST_MASS = args.length

	setting_protein.nested_dic2(seq_string) #outro nested dictionary com os nested fragmentos. ex: {KKKIKGED: {}, KKKIKGE: {}, KKKIKG: {}, KKKIK: {},... KKKIKG: {}}

	##----------------------------------------------------------------------------------------------
	## Ponto Isoeletrico dos fragmentos
	##----------------------------------------------------------------------------------------------

	setting_protein.pka_residues2(seq_string)

	setting_protein.resulting_charge_ph0_2(seq_string)

	setting_protein.charged_residues_quantity_2(seq_string)

	setting_protein.CR_list_2(seq_string)

	setting_protein.pI_2(seq_string)


	##----------------------------------------------------------------------------------------------
	## Calcular massa monoisotópica dos fragmentos
	##----------------------------------------------------------------------------------------------

	setting_protein.massa_amidada_2(seq_string)

	setting_protein.m_h_2(seq_string)

	##----------------------------------------------------------------------------------------------
	## Calcular Hidrofobicidade - GRAVY
	##----------------------------------------------------------------------------------------------

	setting_protein.hydropathicity_2(seq_string)

	setting_protein.hydrophobic_moment_2(seq_string)

	setting_protein.aa_percentage_2(seq_string)

	##----------------------------------------------------------------------------------------------
	## Others
	##----------------------------------------------------------------------------------------------

	setting_protein.aromatic_2(seq_string)

	setting_protein.asa_mainchain_2(seq_string)

	setting_protein.asa_mainchain_nonpolar_2(seq_string)

	setting_protein.asa_mainchain_polar_2(seq_string)

	setting_protein.pct_buried_2(seq_string)

	setting_protein.cf_alpha_2(seq_string)

	setting_protein.cf_beta_2(seq_string)

	setting_protein.cf_turn_2(seq_string)

	setting_protein.charge_2(seq_string)

	setting_protein.volume_2(seq_string)

	##----------------------------------------------------------------------------------------------
	## Novo dicionário com os atributos que serão exportados
	##----------------------------------------------------------------------------------------------

	setting_protein_enzymes.nested_dic_4(seq_string)

	setting_protein_enzymes.attributes_exportation2(seq_string)