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

#----------------------------------------------------------------------------------------------
# Dicionários que atribuirão as funções
#----------------------------------------------------------------------------------------------

for arquivo in sequencia:

	seq_string = sequencia

	setting_protein.nome_proteina(arquivo) #cada fragmento é o valor da chave que é o nome do fragmento. ex: {fragmento1: KKKIKGED, fragmento2: KKKIKG}

	setting_protein.nested_dic(seq_string) #um outro dicionário (só que dicionário de dicionários) em que a sequencia é a chave e o valor será preenchido com os parâmentros. ex: {KKKIKGED: {}, KKKIKG: {}}

	for i in range(0, len(arquivo)):
		exec("dic_frag%d = %s" % (i + 1, {}));

	##----------------------------------------------------------------------------------------------
	## Calcular massa monoisotópica
	##----------------------------------------------------------------------------------------------
	#Dicionário contém lista de resíduos (aa - 18 Da)

	setting_protein.monoisotopic_mass(seq_string)

	##----------------------------------------------------------------------------------------------
	## Parte 2 - Redução dos fragmentos
	##----------------------------------------------------------------------------------------------

	setting_protein.sistematic_reduction(seq_string)

	setting_protein.nested_dic_selection(seq_string)

	setting_protein.monoisotopic_mass_selection(seq_string)

	LOWEST_MASS, HIGHEST_MASS = args.length

	setting_protein.nested_dic2(seq_string) #outro nested dictionary com os nested fragmentos. ex: {KKKIKGED: {}, KKKIKGE: {}, KKKIKG: {}, KKKIK: {},... KKKIKG: {}}

	##----------------------------------------------------------------------------------------------
	## Isoeletric point
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
	## Hydrophobicity Calculation - GRAVY
	##----------------------------------------------------------------------------------------------

	setting_protein.hydropathicity_2(seq_string)

	setting_protein.hydrophobic_moment_2(seq_string)

	##----------------------------------------------------------------------------------------------
	## Amino acids residues percentage
	##----------------------------------------------------------------------------------------------

	setting_protein.aa_percentage_2(seq_string)

	setting_protein.aromatic_2(seq_string)

	setting_protein.cf_alpha_2(seq_string)

	setting_protein.cf_beta_2(seq_string)

	setting_protein.cf_turn_2(seq_string)

	setting_protein.pct_buried_2(seq_string)

	setting_protein.charge_2(seq_string)

	setting_protein.volume_2(seq_string)

	setting_protein.asa_mainchain_2(seq_string)

	setting_protein.asa_mainchain_nonpolar_2(seq_string)

	setting_protein.asa_mainchain_polar_2(seq_string)

	##----------------------------------------------------------------------------------------------
	## New dictionary with the attributes that shall be exported
	##----------------------------------------------------------------------------------------------

	setting_protein.nested_dic3(seq_string)

	setting_protein.attributes_exportation(seq_string)