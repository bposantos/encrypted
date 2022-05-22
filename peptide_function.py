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

for arquivo in sequencia:

	seq_string = sequencia

	while '' in seq_string:
		seq_string.remove('')

	setting_protein.nome_proteina(arquivo)

	setting_protein.nested_dic(seq_string)

	for i in range(0, len(arquivo)):
		exec("dic_frag%d = %s" % (i + 1, {}));

	##----------------------------------------------------------------------------------------------
	## Calcular massa monoisotópica
	##----------------------------------------------------------------------------------------------
	#Dicionário contém lista de resíduos (aa - 18 Da)

	setting_protein.monoisotopic_mass(seq_string)

	##----------------------------------------------------------------------------------------------
	## Isoeletric point
	##----------------------------------------------------------------------------------------------

	setting_peptide.pka_residues3(seq_string)	

	setting_peptide.resulting_charge_ph0_3(seq_string)

	setting_peptide.charged_residues_quantity_3(seq_string)	

	setting_peptide.CR_list_3(seq_string)

	setting_peptide.pI_3(seq_string)

	##----------------------------------------------------------------------------------------------
	## Calcular massa monoisotópica dos peptideos
	##----------------------------------------------------------------------------------------------

	setting_peptide.massa_amidada_3(seq_string)

	setting_peptide.m_h_3(seq_string)

	##----------------------------------------------------------------------------------------------
	## Hydrophobicity Calculation - GRAVY
	##----------------------------------------------------------------------------------------------

	setting_peptide.hydropathicity_3(seq_string)

	setting_peptide.hydrophobic_moment_3(seq_string)

	##----------------------------------------------------------------------------------------------
	## More
	##----------------------------------------------------------------------------------------------

	setting_peptide.aromatic_3(seq_string)

	setting_peptide.cf_alpha_3(seq_string)

	setting_peptide.cf_beta_3(seq_string)

	setting_peptide.cf_turn_3(seq_string)

	setting_peptide.volume_3(seq_string)				

	setting_peptide.charge_3(seq_string)

	setting_peptide.pct_buried_3(seq_string)

	setting_peptide.asa_mainchain_3(seq_string)

	setting_peptide.asa_mainchain_nonpolar_3(seq_string)

	setting_peptide.asa_mainchain_polar_3(seq_string)

	##----------------------------------------------------------------------------------------------
	## Amino acids residues percentage
	##----------------------------------------------------------------------------------------------

	setting_peptide.aa_percentage_3(seq_string)	

	##----------------------------------------------------------------------------------------------
	## New dictionary with the attributes that shall be exported
	##----------------------------------------------------------------------------------------------

	setting_peptide.nested_dic3_2(seq_string)

	setting_peptide.attributes_exportation_3(seq_string)