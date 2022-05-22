from sharedstuff import *
import math

class setting_protein:

	def nome_proteina(lista_fragmentos): #cria dicionário ex: {fragmento1: KKKIKGED, fragmento2: KKKIKG}
		for elemento in lista_fragmentos:
			indice = "fragmento_" + str(lista_fragmentos.index(elemento)+1)
			dic[indice] = elemento

	def nested_dic(lista_peptides):
		for elemen in lista_peptides:
			dicio_frags1[elemen] = {}

	def monoisotopic_mass(fragments_list):
		for fragment in fragments_list:
			n =sum([massa_monoisotopica[residue] for residue in fragment if residue in massa_monoisotopica])
			dicio_frags1[fragment]['monoisotopic_mass'] = round(n+18,2)

	def sistematic_reduction(fragments_list):
		for fragment in fragments_list:
			alist = []
			alist.extend([fragment[i:j] for i in range(len(fragment)) for j in range(i+1,len(fragment)+1)])
			dicio_frags1[fragment]['reduction_list'] = alist

	def nested_dic_selection(fragments_list):
		for fragment in fragments_list:
			for elemen in dicio_frags1[fragment]['reduction_list']:
				dicio_frags_sel[elemen] = {}

	def monoisotopic_mass_selection(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1[fragment]['reduction_list']:
				n =sum([massa_monoisotopica[residue] for residue in element if residue in massa_monoisotopica])
				dicio_frags_sel[element]['monoisotopic_mass'] = round(n+18,2)

	def nested_dic2(fragments_list):
		for elemen in dicio_frags_sel:
			if args.length[0] < dicio_frags_sel[elemen]['monoisotopic_mass'] < args.length[1]:
				dicio_frags2[elemen] = {}
				dicio_frags2[elemen]['monoisotopic_mass'] = dicio_frags_sel[elemen]['monoisotopic_mass']

	def pka_residues2(fragments_list):
		for fragment in fragments_list:
			for elemen in dicio_frags2:
				list_pka = [3.55,7.5]
				list_pka.extend([positive_pKs[residue] for residue in elemen if residue in positive_pKs])
				list_pka.extend([negative_pKs[residue] for residue in elemen if residue in negative_pKs])
				sorted_pka_list = sorted(list_pka)
				dicio_frags2[elemen]['list_pka'] = sorted_pka_list
				unique_list_pka = sorted(list(set(list_pka)))

	def resulting_charge_ph0_2(fragments_list):
		for fragment in fragments_list:
			for elemen in dicio_frags2:
				CR_sum = []
				CR_sum = sum([dic_positive_charge_pKs[element] for element in dicio_frags2[elemen]['list_pka'] if element in dic_positive_charge_pKs])
				dicio_frags2[elemen]['CR_sum'] = CR_sum

	def charged_residues_quantity_2(fragments_list):
		for fragment in fragments_list:
			for elemen in dicio_frags2:
				qtty_charge = sum([1 for residue in elemen if residue in charged_pos_aas or residue in charged_neg_aas])
				qtty_charge += 2
				dicio_frags2[elemen]['quantity_charged_groups'] = qtty_charge

	def CR_list_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				diff = dicio_frags2[element]['CR_sum'] - dicio_frags2[element]['quantity_charged_groups']
				diff += -1
				list_CR = list(range(dicio_frags2[element]['CR_sum'], diff, -1))
				dicio_frags2[element]['list_CR'] = list_CR
				list_CR = 0

	def pI_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				p1 = dicio_frags2[element]['list_CR'].index(0)
				p1 += -1
				p2 = dicio_frags2[element]['list_CR'].index(0)
				p2 += -2
				pI = (dicio_frags2[element]['list_pka'][p1] + dicio_frags2[element]['list_pka'][p2])/2
				if p2 < 0:
					pI = (dicio_frags2[element]['list_pka'][p1])/2
				dicio_frags2[element]['Isoeletric_point'] = pI

	def massa_amidada_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([massa_monoisotopica[residue] for residue in element if residue in massa_monoisotopica])
				dicio_frags2[element]['amidated_mass'] = round(n+17,2)

	def m_h_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([massa_monoisotopica[residue] for residue in element if residue in massa_monoisotopica])
				dicio_frags2[element]['[M+H+]'] = round(n+18,2)

	def hydropathicity_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([hydropathicity_aa[residue] for residue in element if residue in hydropathicity_aa])
				qtdd_rsn = len(element)
				hydropathicity = round(n/qtdd_rsn, 2)
				dicio_frags2[element]['Hydrophobicity'] = hydropathicity

	def hydrophobic_moment_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:		
				vcos = sum([scale_fauchere_pliska[ins]*math.cos(((position*100)*math.pi)/180) for position, ins in enumerate(element)])
				vsin = sum([scale_fauchere_pliska[ins]*math.sin(((position*100)*math.pi)/180) for position, ins in enumerate(element)])
				mH = round(math.sqrt(vsin**2+vcos**2)/len(element),2)
				dicio_frags2[element]['Hydrophobic_moment'] = mH

	def aa_percentage_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				rsn_list = []
				qtty = sum([1 for residue in element if residue in charged_aas])
				rsn_list.extend(residue for residue in element if residue in charged_aas)		
				dicio_frags2[element]['Charged Residues'] = {}
				dicio_frags2[element]['Charged Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []		
				qtty = sum([1 for residue in element if residue in aliphatic_aas])
				rsn_list.extend(residue for residue in element if residue in aliphatic_aas)		
				dicio_frags2[element]['Aliphatic Residues'] = {}
				dicio_frags2[element]['Aliphatic Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []		
				qtty = sum([1 for residue in element if residue in aromatic_aas])
				rsn_list.extend(residue for residue in element if residue in aromatic_aas)
				dicio_frags2[element]['Aromatic Residues'] = {}
				dicio_frags2[element]['Aromatic Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in acidic_aas])
				rsn_list.extend(residue for residue in element if residue in acidic_aas)		
				dicio_frags2[element]['Acidic Residues'] = {}
				dicio_frags2[element]['Acidic Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in basic_aas])
				rsn_list.extend(residue for residue in element if residue in basic_aas)		
				dicio_frags2[element]['Basic Residues'] = {}
				dicio_frags2[element]['Basic Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in hydroxilic_aa])
				rsn_list.extend(residue for residue in element if residue in hydroxilic_aa)		
				dicio_frags2[element]['Hydrophilic Residues'] = {}
				dicio_frags2[element]['Hydrophilic Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in polar_aas])
				rsn_list.extend(residue for residue in element if residue in polar_aas)		
				dicio_frags2[element]['Polar Residues'] = {}
				dicio_frags2[element]['Polar Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in nonpolar_aas])
				rsn_list.extend(residue for residue in element if residue in nonpolar_aas)	
				dicio_frags2[element]['Nonpolar Residues'] = {}
				dicio_frags2[element]['Nonpolar Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []

	def aromatic_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([aromatic[residue] for residue in element if residue in aromatic])
				dicio_frags2[element]['aromatic_rsn'] = round(n+18,2)

	def asa_mainchain_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([asa_mainchain[residue] for residue in element if residue in asa_mainchain])
				dicio_frags2[element]['solvent_accessibility_main_chain'] = round(n+18,2)

	def asa_mainchain_nonpolar_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([asa_sidechain_nonpolar[residue] for residue in element if residue in asa_sidechain_nonpolar])
				dicio_frags2[element]['solvent_accessibility_nonpolar'] = round(n+18,2)

	def asa_mainchain_polar_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([asa_sidechain_polar[residue] for residue in element if residue in asa_sidechain_polar])
				dicio_frags2[element]['solvent_accessibility_polar'] = round(n+18,2)

	def cf_alpha_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([cf_alpha[residue] for residue in element if residue in cf_alpha])
				dicio_frags2[element]['cf_alpha'] = round(n+18,2)

	def pct_buried_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([pct_buried[residue] for residue in element if residue in pct_buried])
				dicio_frags2[element]['pct_buried'] = round(n+18,2)	

	def cf_beta_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([cf_beta[residue] for residue in element if residue in cf_beta])
				dicio_frags2[element]['cf_beta'] = round(n+18,2)

	def cf_turn_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([cf_turn[residue] for residue in element if residue in cf_turn])
				dicio_frags2[element]['cf_turn'] = round(n+18,2)

	def charge_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([charge[residue] for residue in element if residue in charge])
				dicio_frags2[element]['charge'] = round(n+18,2)

	def volume_2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:
				n = sum([volume[residue] for residue in element if residue in volume])
				dicio_frags2[element]['volume'] = round(n+18,2)

	def nested_dic3(fragments_list):
		for fragments_list in dicio_frags2:
				dicio_frags3[fragments_list] = {}

	def attributes_exportation(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:				
				dicio_frags3[element]['[M+H+]'] = dicio_frags2[element]['[M+H+]']
				dicio_frags3[element]['Amidated_Mass'] = dicio_frags2[element]['amidated_mass']
				dicio_frags3[element]['Monoisotopic_Mass'] = dicio_frags2[element]['monoisotopic_mass']				
				dicio_frags3[element]['Hydrophobicity'] = dicio_frags2[element]['Hydrophobicity']
				dicio_frags3[element]['Hydrophobic_Moment'] = dicio_frags2[element]['Hydrophobic_moment']				
				dicio_frags3[element]['Isoeletric_Point'] = dicio_frags2[element]['Isoeletric_point']
				dicio_frags3[element]['Nonpolar_Residues_(%)'] = dicio_frags2[element]['Nonpolar Residues']['Percentage (%)']				
				dicio_frags3[element]['Polar_Residues_(%)'] = dicio_frags2[element]['Polar Residues']['Percentage (%)']				
				dicio_frags3[element]['Acidic_Residues_(%)'] = dicio_frags2[element]['Acidic Residues']['Percentage (%)']
				dicio_frags3[element]['Basic_Residues_(%)'] = dicio_frags2[element]['Basic Residues']['Percentage (%)']
				dicio_frags3[element]['Charged_Residues_(%)'] = dicio_frags2[element]['Charged Residues']['Percentage (%)']
				dicio_frags3[element]['aromatic_rsn'] = dicio_frags2[element]['aromatic_rsn']
				dicio_frags3[element]['cf_alpha'] = dicio_frags2[element]['cf_alpha']
				dicio_frags3[element]['cf_beta'] = dicio_frags2[element]['cf_beta']
				dicio_frags3[element]['cf_turn'] = dicio_frags2[element]['cf_turn']
				dicio_frags3[element]['volume'] = dicio_frags2[element]['volume']
				dicio_frags3[element]['charge'] = dicio_frags2[element]['charge']
				dicio_frags3[element]['pct_buried'] = dicio_frags2[element]['pct_buried']
				dicio_frags3[element]['solvent_accessibility_main_chain'] = dicio_frags2[element]['solvent_accessibility_main_chain']
				dicio_frags3[element]['solvent_accessibility_nonpolar'] = dicio_frags2[element]['solvent_accessibility_nonpolar']
				dicio_frags3[element]['solvent_accessibility_polar'] = dicio_frags2[element]['solvent_accessibility_polar']
			indice = "fragmento_" + str(fragments_list.index(fragment)+1)
			dic2[indice + cabecalho4] = fragments_list
			dicio_frags3['Fragments'] =  dic2

class setting_protein_enzymes:

	def pka_residues(fragments_list):
		for fragment in fragments_list:
			list_pka = [3.55,7.5]
			list_pka.extend([positive_pKs[residue] for residue in fragment if residue in positive_pKs])
			list_pka.extend([negative_pKs[residue] for residue in fragment if residue in negative_pKs])
			sorted_pka_list = sorted(list_pka)
			dicio_frags1[fragment]['list_pka'] = sorted_pka_list
			unique_list_pka = sorted(list(set(list_pka)))
			dicio_frags1[fragment]['unique_list_pka'] = unique_list_pka

	def resulting_charge_ph0(fragments_list):
		for fragment in fragments_list:
			CR_sum = []
			CR_sum = sum([dic_positive_charge_pKs[element] for element in dicio_frags1[fragment]['list_pka'] if element in dic_positive_charge_pKs])
			dicio_frags1[fragment]['CR_sum'] = CR_sum

	def charged_residues_quantity(fragments_list):
		for fragment in fragments_list:
			qtty_charge = sum([1 for residue in fragment if residue in charged_pos_aas or residue in charged_neg_aas])
			qtty_charge += 2
			dicio_frags1[fragment]['quantity_charged_groups'] = qtty_charge

	def CR_list(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1[fragment]:
				diff = dicio_frags1[fragment]['CR_sum'] - dicio_frags1[fragment]['quantity_charged_groups']
				diff += -1
				list_CR = list(range(dicio_frags1[fragment]['CR_sum'], diff, -1))
			dicio_frags1[fragment]['list_CR'] = list_CR
			list_CR = 0

	def pI(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1[fragment]['list_CR']:
				position_1 = dicio_frags1[fragment]['list_CR'].index(0)
				position_1 += -1
				position_2 = dicio_frags1[fragment]['list_CR'].index(0)
				position_2 += -2
				pI = (dicio_frags1[fragment]['list_pka'][position_1] + dicio_frags1[fragment]['list_pka'][position_2])/2
				if position_2 == 0:
					pI = (dicio_frags1[fragment]['list_pka'][position_1] + 0)/2
				elif position_2 < 0:
					pI = (dicio_frags1[fragment]['list_pka'][position_1] + 0)/2
			dicio_frags1[fragment]['Isoeletric_point'] = pI

	def massa_amidada(fragments_list):
		for fragment in fragments_list:
			n = sum([massa_monoisotopica[residue] for residue in fragment if residue in massa_monoisotopica])
			dicio_frags1[fragment]['amidated_mass'] = round(n+17,2)

	def m_h(fragments_list):
		for fragment in fragments_list:
			n = sum([massa_monoisotopica[residue] for residue in fragment if residue in massa_monoisotopica])
			dicio_frags1[fragment]['[M+H+]'] = round(n+18,2)

	def aromatic(fragments_list):
		for fragment in fragments_list:
			n = sum([aromatic[residue] for residue in fragment if residue in aromatic])
			dicio_frags1[fragment]['aromatic_rsn'] = round(n+18,2)

	def asa_mainchain(fragments_list):
		for fragment in fragments_list:
			n = sum([asa_mainchain[residue] for residue in fragment if residue in asa_mainchain])
			dicio_frags1[fragment]['solvent_accessibility_main_chain'] = round(n+18,2)

	def asa_mainchain_nonpolar(fragments_list):
		for fragment in fragments_list:
			n = sum([asa_sidechain_nonpolar[residue] for residue in fragment if residue in asa_sidechain_nonpolar])
			dicio_frags1[fragment]['solvent_accessibility_nonpolar'] = round(n+18,2)

	def asa_mainchain_polar(fragments_list):
		for fragment in fragments_list:
			n = sum([asa_sidechain_polar[residue] for residue in fragment if residue in asa_sidechain_polar])
			dicio_frags1[fragment]['solvent_accessibility_polar'] = round(n+18,2)

	def cf_alpha(fragments_list):
		for fragment in fragments_list:
			n = sum([cf_alpha[residue] for residue in fragment if residue in cf_alpha])
			dicio_frags1[fragment]['cf_alpha'] = round(n+18,2)

	def pct_buried(fragments_list):
		for fragment in fragments_list:
			n = sum([pct_buried[residue] for residue in fragment if residue in pct_buried])
			dicio_frags1[fragment]['pct_buried'] = round(n+18,2)	

	def cf_beta(fragments_list):
		for fragment in fragments_list:
			n = sum([cf_beta[residue] for residue in fragment if residue in cf_beta])
			dicio_frags1[fragment]['cf_beta'] = round(n+18,2)

	def cf_turn(fragments_list):
		for fragment in fragments_list:
			n = sum([cf_turn[residue] for residue in fragment if residue in cf_turn])
			dicio_frags1[fragment]['cf_turn'] = round(n+18,2)

	def charge(fragments_list):
		for fragment in fragments_list:
			n = sum([charge[residue] for residue in fragment if residue in charge])
			dicio_frags1[fragment]['charge'] = round(n+18,2)

	def volume(fragments_list):
		for fragment in fragments_list:
			n = sum([volume[residue] for residue in fragment if residue in volume])
			dicio_frags2[fragment]['volume'] = round(n+18,2)

	def hydropathicity(fragments_list):
		for fragment in fragments_list:
			n = sum([hydropathicity_aa[residue] for residue in fragment if residue in hydropathicity_aa])
			qtdd_rsn = len(fragment)
			hydropathicity = round(n/qtdd_rsn,2)			
			dicio_frags1[fragment]['Hydrophobicity'] = hydropathicity

	def aa_percentage(fragments_list):
		qtty = 0
		for fragment in fragments_list:
			rsn_list = []
			qtty = sum([1 for residue in fragment if residue in charged_aas])
			rsn_list.extend(residue for residue in fragment if residue in charged_aas)		
			dicio_frags1[fragment]['Charged Residues'] = {}
			dicio_frags1[fragment]['Charged Residues']['Percentage (%)'] = round((qtty/len(fragment)*100),2)	
			qtty = 0
			rsn_list = []		
			qtty = sum([1 for residue in fragment if residue in aliphatic_aas])
			rsn_list.extend(residue for residue in fragment if residue in aliphatic_aas)		
			dicio_frags1[fragment]['Aliphatic Residues'] = {}
			dicio_frags1[fragment]['Aliphatic Residues']['Percentage (%)'] = round((qtty/len(fragment)*100),2)		
			qtty = 0
			rsn_list = []		
			qtty = sum([1 for residue in fragment if residue in aromatic_aas])
			rsn_list.extend(residue for residue in fragment if residue in aromatic_aas)		
			dicio_frags1[fragment]['Aromatic Residues'] = {}
			dicio_frags1[fragment]['Aromatic Residues']['Percentage (%)'] = round((qtty/len(fragment)*100),2)		
			qtty = 0
			rsn_list = []
			qtty = sum([1 for residue in fragment if residue in acidic_aas])
			rsn_list.extend(residue for residue in fragment if residue in acidic_aas)		
			dicio_frags1[fragment]['Acidic Residues'] = {}
			dicio_frags1[fragment]['Acidic Residues']['Percentage (%)'] = round((qtty/len(fragment)*100),2)		
			qtty = 0
			rsn_list = []
			qtty = sum([1 for residue in fragment if residue in basic_aas])
			rsn_list.extend(residue for residue in fragment if residue in basic_aas)		
			dicio_frags1[fragment]['Basic Residues'] = {}
			dicio_frags1[fragment]['Basic Residues']['Percentage (%)'] = round((qtty/len(fragment)*100),2)		
			qtty = 0
			rsn_list = []
			qtty = sum([1 for residue in fragment if residue in hydroxilic_aa])
			rsn_list.extend(residue for residue in fragment if residue in hydroxilic_aa)		
			dicio_frags1[fragment]['Hydroxilic Residues'] = {}
			dicio_frags1[fragment]['Hydroxilic Residues']['Percentage (%)'] = round((qtty/len(fragment)*100),2)		
			qtty = 0
			rsn_list = []
			qtty = sum([1 for residue in fragment if residue in polar_aas])
			rsn_list.extend(residue for residue in fragment if residue in polar_aas)	
			dicio_frags1[fragment]['Polar Residues'] = {}
			dicio_frags1[fragment]['Polar Residues']['Percentage (%)'] = round((qtty/len(fragment)*100),2)		
			qtty = 0
			rsn_list = []
			qtty = sum([1 for residue in fragment if residue in nonpolar_aas])
			rsn_list.extend(residue for residue in fragment if residue in nonpolar_aas)		
			dicio_frags1[fragment]['Nonpolar Residues'] = {}
			dicio_frags1[fragment]['Nonpolar Residues']['Percentage (%)'] = round((qtty/len(fragment)*100),2)		
			qtty = 0
			rsn_list = []

	def sistematic_reduction_deep(fragments_list):
		for fragment in fragments_list:
			length = len(fragment)
			alist = []
			alist.extend([fragment[i:] for i in range(length)])
			dicio_frags1[fragment]['reduction_list'] = alist

	def nested_dic_4(fragments_list):
		for fragment in fragments_list:
			for elemen in dicio_frags2:
					dicio_frags3[elemen] = {}

	def attributes_exportation2(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags2:				
				dicio_frags3[element]['[M+H+]'] = dicio_frags2[element]['[M+H+]']
				dicio_frags3[element]['Amidated_Mass'] = dicio_frags2[element]['amidated_mass']
				dicio_frags3[element]['Monoisotopic_Mass'] = dicio_frags2[element]['monoisotopic_mass']				
				dicio_frags3[element]['Hydrophobicity'] = dicio_frags2[element]['Hydrophobicity']
				dicio_frags3[element]['Hydrophobic_Moment'] = dicio_frags2[element]['Hydrophobic_moment']				
				dicio_frags3[element]['Isoeletric_Point'] = dicio_frags2[element]['Isoeletric_point']
				dicio_frags3[element]['Nonpolar_Residues_(%)'] = dicio_frags2[element]['Nonpolar Residues']['Percentage (%)']				
				dicio_frags3[element]['Polar_Residues_(%)'] = dicio_frags2[element]['Polar Residues']['Percentage (%)']				
				dicio_frags3[element]['Acidic_Residues_(%)'] = dicio_frags2[element]['Acidic Residues']['Percentage (%)']
				dicio_frags3[element]['Basic_Residues_(%)'] = dicio_frags2[element]['Basic Residues']['Percentage (%)']
				dicio_frags3[element]['Charged_Residues_(%)'] = dicio_frags2[element]['Charged Residues']['Percentage (%)']
				dicio_frags3[element]['aromatic_rsn'] = dicio_frags2[element]['aromatic_rsn']
				dicio_frags3[element]['cf_alpha'] = dicio_frags2[element]['cf_alpha']
				dicio_frags3[element]['cf_beta'] = dicio_frags2[element]['cf_beta']
				dicio_frags3[element]['cf_turn'] = dicio_frags2[element]['cf_turn']
				dicio_frags3[element]['volume'] = dicio_frags2[element]['volume']
				dicio_frags3[element]['charge'] = dicio_frags2[element]['charge']
				dicio_frags3[element]['pct_buried'] = dicio_frags2[element]['pct_buried']
				dicio_frags3[element]['solvent_accessibility_main_chain'] = dicio_frags2[element]['solvent_accessibility_main_chain']
				dicio_frags3[element]['solvent_accessibility_nonpolar'] = dicio_frags2[element]['solvent_accessibility_nonpolar']
				dicio_frags3[element]['solvent_accessibility_polar'] = dicio_frags2[element]['solvent_accessibility_polar']
				dic = {}
			indice = "fragmento_" + str(fragments_list.index(fragment)+1)
			dic2[indice + cabecalho4] = fragments_list
			dicio_frags3['Fragments'] =  dic2

class setting_peptide:

	def pka_residues3(fragments_list):
		for fragment in fragments_list:
			for elemen in dicio_frags1:
				list_pka = [3.55,7.5]
				list_pka.extend([positive_pKs[residue] for residue in elemen if residue in positive_pKs])
				list_pka.extend([negative_pKs[residue] for residue in elemen if residue in negative_pKs])
				sorted_pka_list = sorted(list_pka)
				dicio_frags1[elemen]['list_pka'] = sorted_pka_list
				unique_list_pka = sorted(list(set(list_pka)))

	def resulting_charge_ph0_3(fragments_list):
		for fragment in fragments_list:
			for elemen in dicio_frags1:
				CR_sum = []
				CR_sum = sum([dic_positive_charge_pKs[element] for element in dicio_frags1[elemen]['list_pka'] if element in dic_positive_charge_pKs])
				dicio_frags1[elemen]['CR_sum'] = CR_sum

	def charged_residues_quantity_3(fragments_list):
		for fragment in fragments_list:
			for elemen in dicio_frags1:
				qtty_charge = sum([1 for residue in elemen if residue in charged_pos_aas or residue in charged_neg_aas])
				qtty_charge += 2
				dicio_frags1[elemen]['quantity_charged_groups'] = qtty_charge

	def CR_list_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				diff = dicio_frags1[element]['CR_sum'] - dicio_frags1[element]['quantity_charged_groups']
				diff += -1
				list_CR = list(range(dicio_frags1[element]['CR_sum'], diff, -1))
				dicio_frags1[element]['list_CR'] = list_CR
				list_CR = 0

	def pI_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				p1 = dicio_frags1[element]['list_CR'].index(0)
				p1 += -1
				p2 = dicio_frags1[element]['list_CR'].index(0)
				p2 += -2
				pI = (dicio_frags1[element]['list_pka'][p1] + dicio_frags1[element]['list_pka'][p2])/2
				if p2 < 0:
					pI = (dicio_frags1[element]['list_pka'][p1])/2
				dicio_frags1[element]['Isoeletric_point'] = pI

	def massa_amidada_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([massa_monoisotopica[residue] for residue in element if residue in massa_monoisotopica])
				dicio_frags1[element]['amidated_mass'] = round(n+17,2)

	def m_h_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([massa_monoisotopica[residue] for residue in element if residue in massa_monoisotopica])
				dicio_frags1[element]['[M+H+]'] = round(n+18,2)

	def aromatic_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([aromatic[residue] for residue in element if residue in aromatic])
				dicio_frags1[element]['aromatic_rsn'] = round(n+18,2)

	def asa_mainchain_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([asa_mainchain[residue] for residue in element if residue in asa_mainchain])
				dicio_frags1[element]['solvent_accessibility_main_chain'] = round(n+18,2)

	def asa_mainchain_nonpolar_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([asa_sidechain_nonpolar[residue] for residue in element if residue in asa_sidechain_nonpolar])
				dicio_frags1[element]['solvent_accessibility_nonpolar'] = round(n+18,2)

	def asa_mainchain_polar_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([asa_sidechain_polar[residue] for residue in element if residue in asa_sidechain_polar])
				dicio_frags1[element]['solvent_accessibility_polar'] = round(n+18,2)

	def cf_alpha_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([cf_alpha[residue] for residue in element if residue in cf_alpha])
				dicio_frags1[element]['cf_alpha'] = round(n+18,2)

	def pct_buried_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([pct_buried[residue] for residue in element if residue in pct_buried])
				dicio_frags1[element]['pct_buried'] = round(n+18,2)	

	def cf_beta_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([cf_beta[residue] for residue in element if residue in cf_beta])
				dicio_frags1[element]['cf_beta'] = round(n+18,2)

	def cf_turn_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([cf_turn[residue] for residue in element if residue in cf_turn])
				dicio_frags1[element]['cf_turn'] = round(n+18,2)

	def charge_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([charge[residue] for residue in element if residue in charge])
				dicio_frags1[element]['charge'] = round(n+18,2)

	def volume_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([volume[residue] for residue in element if residue in volume])
				dicio_frags1[element]['volume'] = round(n+18,2)

	def hydropathicity_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				n = sum([hydropathicity_aa[residue] for residue in element if residue in hydropathicity_aa])
				qtdd_rsn = len(element)
				hydropathicity = round(n/qtdd_rsn, 2)
				dicio_frags1[element]['Hydrophobicity'] = hydropathicity

	def hydrophobic_moment_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:		
				vcos = sum([scale_fauchere_pliska[ins]*math.cos(((position*100)*math.pi)/180) for position, ins in enumerate(element)])
				vsin = sum([scale_fauchere_pliska[ins]*math.sin(((position*100)*math.pi)/180) for position, ins in enumerate(element)])
				mH = round(math.sqrt(vsin**2+vcos**2)/len(element),2)
				dicio_frags1[element]['Hydrophobic_moment'] = mH

	def aa_percentage_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:
				rsn_list = []
				qtty = sum([1 for residue in element if residue in charged_aas])
				rsn_list.extend(residue for residue in element if residue in charged_aas)		
				dicio_frags1[element]['Charged Residues'] = {}
				dicio_frags1[element]['Charged Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []		
				qtty = sum([1 for residue in element if residue in aliphatic_aas])
				rsn_list.extend(residue for residue in element if residue in aliphatic_aas)		
				dicio_frags1[element]['Aliphatic Residues'] = {}
				dicio_frags1[element]['Aliphatic Residues']['Percentage (%)'] = round((qtty/len(element)*100),2)		
				qtty = 0
				rsn_list = []		
				qtty = sum([1 for residue in element if residue in aromatic_aas])
				rsn_list.extend(residue for residue in element if residue in aromatic_aas)
				dicio_frags1[element]['Aromatic Residues'] = {}
				dicio_frags1[element]['Aromatic Residues']['Percentage (%)'] = round((qtty/len(element)*100),3)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in acidic_aas])
				rsn_list.extend(residue for residue in element if residue in acidic_aas)		
				dicio_frags1[element]['Acidic Residues'] = {}
				dicio_frags1[element]['Acidic Residues']['Percentage (%)'] = round((qtty/len(element)*100),3)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in basic_aas])
				rsn_list.extend(residue for residue in element if residue in basic_aas)		
				dicio_frags1[element]['Basic Residues'] = {}
				dicio_frags1[element]['Basic Residues']['Percentage (%)'] = round((qtty/len(element)*100),3)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in hydroxilic_aa])
				rsn_list.extend(residue for residue in element if residue in hydroxilic_aa)		
				dicio_frags1[element]['Hydrophilic Residues'] = {}
				dicio_frags1[element]['Hydrophilic Residues']['Percentage (%)'] = round((qtty/len(element)*100),3)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in polar_aas])
				rsn_list.extend(residue for residue in element if residue in polar_aas)		
				dicio_frags1[element]['Polar Residues'] = {}
				dicio_frags1[element]['Polar Residues']['Percentage (%)'] = round((qtty/len(element)*100),3)		
				qtty = 0
				rsn_list = []
				qtty = sum([1 for residue in element if residue in nonpolar_aas])
				rsn_list.extend(residue for residue in element if residue in nonpolar_aas)	
				dicio_frags1[element]['Nonpolar Residues'] = {}
				dicio_frags1[element]['Nonpolar Residues']['Percentage (%)'] = round((qtty/len(element)*100),3)		
				qtty = 0
				rsn_list = []

	def nested_dic3_2(fragments_list):
		for fragment in fragments_list:
			for elemen in dicio_frags1:
					dicio_frags3[elemen] = {}

	def attributes_exportation_3(fragments_list):
		for fragment in fragments_list:
			for element in dicio_frags1:				
				dicio_frags3[element]['[M+H+]'] = dicio_frags1[element]['[M+H+]']
				dicio_frags3[element]['Amidated_Mass'] = dicio_frags1[element]['amidated_mass']
				dicio_frags3[element]['Monoisotopic_Mass'] = dicio_frags1[element]['monoisotopic_mass']				
				dicio_frags3[element]['Hydrophobicity'] = dicio_frags1[element]['Hydrophobicity']
				dicio_frags3[element]['Hydrophobic_Moment'] = dicio_frags1[element]['Hydrophobic_moment']				
				dicio_frags3[element]['Isoeletric_Point'] = dicio_frags1[element]['Isoeletric_point']
				dicio_frags3[element]['Nonpolar_Residues_(%)'] = dicio_frags1[element]['Nonpolar Residues']['Percentage (%)']				
				dicio_frags3[element]['Polar_Residues_(%)'] = dicio_frags1[element]['Polar Residues']['Percentage (%)']				
				dicio_frags3[element]['Acidic_Residues_(%)'] = dicio_frags1[element]['Acidic Residues']['Percentage (%)']
				dicio_frags3[element]['Basic_Residues_(%)'] = dicio_frags1[element]['Basic Residues']['Percentage (%)']
				dicio_frags3[element]['Charged_Residues_(%)'] = dicio_frags1[element]['Charged Residues']['Percentage (%)']
				dicio_frags3[element]['aromatic_rsn'] = dicio_frags1[element]['aromatic_rsn']
				dicio_frags3[element]['cf_alpha'] = dicio_frags1[element]['cf_alpha']
				dicio_frags3[element]['cf_beta'] = dicio_frags1[element]['cf_beta']
				dicio_frags3[element]['cf_turn'] = dicio_frags1[element]['cf_turn']
				dicio_frags3[element]['volume'] = dicio_frags1[element]['volume']
				dicio_frags3[element]['charge'] = dicio_frags1[element]['charge']
				dicio_frags3[element]['pct_buried'] = dicio_frags1[element]['pct_buried']
				dicio_frags3[element]['solvent_accessibility_main_chain'] = dicio_frags1[element]['solvent_accessibility_main_chain']
				dicio_frags3[element]['solvent_accessibility_nonpolar'] = dicio_frags1[element]['solvent_accessibility_nonpolar']
				dicio_frags3[element]['solvent_accessibility_polar'] = dicio_frags1[element]['solvent_accessibility_polar']
				dic = {}
			indice = "fragmento_" + str(fragments_list.index(fragment)+1)
			dic2[indice + cabecalho4] = fragment
			dicio_frags3['Fragments'] =  dic2

class msms:

	def massa(proteina):
		n = 0
		for aa in proteina:
			n = round((n + massa_monoisotopica[aa]),2)
		return n

	def massa_amidada_orig(proteina):
		n = 0
		for aa in proteina:
			n = n + massa_monoisotopica[aa]
			m = round((n + 17),2)
		return m

	def y_series(peptides_column):
		for element in peptides_column:	
			other_fragment = element
			serie_y = []
			while other_fragment:
				serie_y.append(msms.massa_amidada_orig(other_fragment)+2)
				other_fragment = other_fragment[1:]
			dic3[element]['y_series'] = serie_y

	def a_series(peptides_column):
		for element in peptides_column:
			other_fragment = element
			serie_a = []
			while other_fragment:
				serie_a.append(msms.massa_amidada_orig(other_fragment)-44)
				other_fragment = other_fragment[:-1]
			dic3[element]['a_series'] = serie_a[::-1]

	def b_series(peptides_column):
		for element in peptides_column:
			dic3[element] = {}
			dic3[element]['massa_monoisotopica'] = round((msms.massa_amidada_orig(element)+1),2)
			other_fragment = element
			serie_b = []
			while other_fragment:
				serie_b.append(msms.massa_amidada_orig(other_fragment)-16)
				other_fragment = other_fragment[:-1]
			dic3[element]['b_series'] = serie_b[::-1]