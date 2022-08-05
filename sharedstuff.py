from itertools import chain, islice
import glob
import argparse
import os


#----------------------------------------------------------------------------------------------
# Read FASTA file
#----------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Proteolysis/Physical-Chemical Features/Antimicrobial Peptides Predictor')
parser.add_argument("-f", '--file', help="sequence in fasta format", action="store")
parser.add_argument("-t", '--type', help="peptide or protein? Type the first or the second one.", action="store")
parser.add_argument("-cl","--cleavage", help="(ONLY FOR protein) optional: breaks the sequence at the typed enzymes:\n 1 = argCproteinase,\n / 2 = trypsinKproteinase, 3 - trypsinRproteinase, 4 - chymotrypsin_high_proteinase, 5 - chymotrypsin_low_proteinase, 6 - pepsin_ph_2_proteinase, 7 = trypsinGKproteinase", action="store")
parser.add_argument("-l", '--length', help="(ONLY FOR protein) sets the peptides range in Da. Suggested> -l 699 3001", action="store", type=int, nargs=2, metavar=('LOWEST_MASS', 'HIGHEST_MASS'))
parser.add_argument("-act", '--activity', help="select the activity", action="store")
parser.add_argument("-n", '--number', help="number on files to differentiation - optional", action="store")
parser.add_argument("-m", '--multiprocessing', help="chopse the number of CPU you want to use. This helps in acceleration. Optional", action="store")
args = parser.parse_args()

#----------------------------------------------------------------------------------------------
# Args Type
#----------------------------------------------------------------------------------------------

argstype = str(args.type)
nome_arquivo = args.file

#----------------------------------------------------------------------------------------------
# Identify header and sequences from FASTA file
#----------------------------------------------------------------------------------------------

print("\nHello. Your Sequence is being fragmented in the selected range and specification. Wait some minutes until is done.\n")

cabecalho =[]
sequencia = []
aux = ""

#############################################################################
# Função para dividir o texto em linhas
#############################################################################

def chunks(iterable, n):
	iterable = iter(iterable)
	while True:
		try:
			yield chain([next(iterable)], islice(iterable, n-1))
		except StopIteration:
			return

#############################################################################
# Abrir o arquivo original, contar a quantidade de linhas e dividir no
# inteiro par mais próximo
#############################################################################

with open(args.file) as f:
    for i, l in enumerate(f):
        pass
        qtdd_linhas = i + 1

if args.multiprocessing:
	mp = int(args.multiprocessing)

if not args.multiprocessing:
	mp = 1

num = qtdd_linhas/mp
if (round(num) % 2) == 0:
  line_file = round(num)
else:
  line_file = round(num)+1

#############################################################################
# Criar arquivos temporários
#############################################################################

dir = './temp'       
os.makedirs(dir)

with open(args.file) as bigfile:
    for i, lines in enumerate(chunks(bigfile, line_file)):
        file_split = '{}.{}'.format(args.file, i)
        with open("temp/"+file_split, 'w') as f:
            f.writelines(lines)

#############################################################################
# Multiprocessamento
#############################################################################

#os("python3 encrypted.py -f teste2.fasta.0 -t peptide -n 3 -m 4 -act antibacterial | python3 encrypted.py -f teste2.fasta.0 -t peptide -n 3 -m 4 -act antibacterial")

#############################################################################
# Ler cada arquivo temporário
#############################################################################

path = './temp'
for filename in glob.glob(os.path.join(path, '*.fasta.*')):
  with open(os.path.join(os.getcwd(), filename), 'r') as files:
    for linha in files:
      if linha[0] == ">":
        cabecalho.append(linha)
        if aux != "":
          sequencia.append(aux)
        aux = ""
      else:
        aux += linha.strip() #aux = aux + linha.strip (same thing)

  sequencia.append(aux)
  cabecalho2 = str(cabecalho).strip('['']')
  cabecalho3 = cabecalho2.split(" ")
  cabecalho4 = str(cabecalho3[0]).strip('['']')

#----------------------------------------------------------------------------------------------
# Parameters
#----------------------------------------------------------------------------------------------

hydropathicity_aa = {
	'A': 1.800,
	'C': 2.500,
	'D': -3.500,
	'E': -3.500,
	'F': 2.800,
	'G': -0.400,
	'H': -3.200,
	'I': 4.500,
	'K': -3.900,
	'L': 3.800,
	'M': 1.900,
	'N': -3.500,
	'P': -1.600,
	'Q': -3.500,
	'R': -4.500,
	'S': -0.800,
	'T': -0.700,
	'V': 4.200,
	'W': -0.900,
	'Y': -1.300,
}

#Dicionário contém lista de resíduos (aa - 18 Da)
massa_monoisotopica = {
	'A': 71.03711,
	'C': 103.00919,
	'D': 115.02694,
	'E': 129.04259,
	'F': 147.06841,
	'G': 57.02146,
	'H': 137.05891,
	'I': 113.08406,
	'K': 128.09496,
	'L': 113.08406,
	'M': 131.04049,
	'N': 114.04293,
	'P': 97.05276,
	'Q': 128.05858,
	'R': 156.10111,
	'S': 87.03203,
	'T': 101.04768,
	'V': 99.09841,
	'W': 186.07931,
	'Y': 163.06333,
}

charged_pos_aas = ('K', 'R', 'H') 
charged_neg_aas = ('D', 'E', 'C', 'Y') 
positive_pKs = {'Nterm': 7.5, 'K': 10.0, 'R': 12.0, 'H': 5.98} 
negative_pKs = {'Cterm': 3.55, 'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0}
dic_positive_charge_pKs = {7.5: +1, 10.0: +1, 12.0: +1, 5.98: +1} 
dic_negative_charge_pKs = {3.55: -1, 4.05: -1, 4.45: -1, 9.0: -1, 10.0: -1}

scale_fauchere_pliska = {
	'A':  0.31, 
	'R': -1.01, 
	'N': -0.60,
	'D': -0.77, 
	'C':  1.54, 
	'Q': -0.22,
	'E': -0.64, 
	'G':  0.00, 
	'H':  0.13,
	'I':  1.80, 
	'L':  1.70, 
	'K': -0.99,
	'M':  1.23, 
	'F':  1.79, 
	'P':  0.72,
	'S': -0.04, 
	'T':  0.26, 
	'W':  2.25,
	'Y':  0.96, 
	'V':  1.22,
}

polar_aas = ('R','Q','N','E','K','H','D','S','T','Y','C','Y','M','G')
nonpolar_aas = ('A','I','L','M','F','V','P','W')
charged_aas = ('K', 'R', 'H', 'D', 'E')
aliphatic_aas = ('A','G','I', 'L','P','V')
aromatic_aas = ('F','W','Y')
acidic_aas = ('D','E')
basic_aas = ('R','H','K')
hydroxilic_aa = ('S','T')

argC = ['R']
chymotrypsin_high_specificity = ['F','Y','W']
chymotrypsin_low_specificity = ['F','Y','M','L']
pepsin_ph_higher_2 = ['F','L','W','Y']
trypsinK = ['K']
trypsinR = ['R']
trypsinGK = ['GK']


#description": "is amino acid aromatic?
aromatic = {
	"A": 0.0,
	"C": 0.0,
	"D": 0.0,
	"E": 0.0,
	"F": 1.0,
	"G": 0.0,
	"H": 0.0,
	"I": 0.0,
	"K": 0.0,
	"L": 0.0,
	"M": 0.0,
	"N": 0.0,
	"P": 0.0,
	"Q": 0.0,
	"R": 0.0,
	"S": 0.0,
	"T": 0.0,
	"V": 0.0,
	"W": 1.0,
	"Y": 1.0
	}

#"description": "solvent accessibility of main chain",
#"notes": "calculated using NACCESS on extended ALA-X-ALA tripeptide, water radius of 1.4"

asa_mainchain = {
	"A": 38.54,
	"C": 37.53,
	"D": 37.7,
	"E": 37.51,
	"F": 35.37,
	"G": 47.77,
	"H": 35.8,
	"I": 37.16,
	"K": 37.51,
	"L": 37.51,
	"M": 37.51,
	"N": 37.7,
	"P": 16.23,
	"Q": 37.51,
	"R": 37.51,
	"S": 38.4,
	"T": 37.57,
	"V": 37.16,
	"W": 38.1,
	"Y": 35.38
	}

#"description": "solvent accessibility of nonpolar atoms in sidechain",
#"notes": "calculated using NACCESS on extended ALA-X-ALA tripeptide, water radius of 1.4",
#"refs": "Hubbard, Simon J., and Janet M. Thornton. 'Naccess.' Computer Program, Department of Biochemistry and Molecular Biology, University College London 2.1 (1993).",
asa_sidechain_nonpolar = {
	"A": 71.38,
	"C": 97.93,
	"D": 49.24,
	"E": 60.29,
	"F": 165.25,
	"G": 37.55,
	"H": 97.15,
	"I": 139.14,
	"K": 116.57,
	"L": 142.31,
	"M": 157.84,
	"N": 46.23,
	"P": 120.95,
	"Q": 52.22,
	"R": 77.8,
	"S": 48.55,
	"T": 75.72,
	"V": 115.47,
	"W": 189.67,
	"Y": 136.5
	}

#"description": "solvent accessibility of polar atoms in sidechain",
#"notes": "calculated using NACCESS on extended ALA-X-ALA tripeptide, water radius of 1.4",
#"refs": "Hubbard, Simon J., and Janet M. Thornton. 'Naccess.' Computer Program, Department of Biochemistry and Molecular Biology, University College London 2.1 (1993).",
asa_sidechain_polar = {
	"A": 36.58,
	"C": 36.35,
	"D": 91.15,
	"E": 111.96,
	"F": 34.23,
	"G": 42.55,
	"H": 85.73,
	"I": 35.98,
	"K": 84.24,
	"L": 36.32,
	"M": 36.32,
	"N": 97.72,
	"P": 15.19,
	"Q": 126.28,
	"R": 160.97,
	"S": 67.95,
	"T": 63.55,
	"V": 35.97,
	"W": 59.69,
	"Y": 76.26
    }

#"description": "Chou-Fasman alpha-helix propensity",
#"notes": "",
#"refs": "Chou, Peter Y., and Gerald D. Fasman. 'Empirical predictions of protein conformation.' Annual review of biochemistry 47.1 (1978): 251-276.",
cf_alpha = {
	"A": 1.42,
	"C": 0.7,
	"D": 1.01,
	"E": 1.51,
	"F": 1.13,
	"G": 0.57,
	"H": 1.0,
	"I": 1.08,
	"K": 1.14,
	"L": 1.21,
	"M": 1.45,
	"N": 0.67,
	"P": 0.57,
	"Q": 1.11,
	"R": 0.98,
	"S": 0.77,
	"T": 0.83,
	"V": 1.06,
	"W": 0.81,
	"Y": 0.68
	}


#"description": "fraction of time this amino acid is buried in structures",
#"notes": "",
#"refs": "Schein, Catherine H. 'Solubility as a function of protein structure and solvent components.' Nature Biotechnology 8.4 (1990): 308-317."
pct_buried = {
	"A": 0.38,
	"C": 0.47,
	"D": 0.145,
	"E": 0.2,
	"F": 0.48,
	"G": 0.37,
	"H": 0.19,
	"I": 0.63,
	"K": 0.042,
	"L": 0.41,
	"M": 0.5,
	"N": 0.1,
	"P": 0.24,
	"Q": 0.63,
	"R": 0.0,
	"S": 0.24,
	"T": 0.25,
	"V": 0.56,
	"W": 0.23,
	"Y": 0.13
    }

#"description": "Chou-Fasman beta-sheet propensity",
#"refs": "Chou, Peter Y., and Gerald D. Fasman. 'Empirical predictions of protein conformation.' Annual review of biochemistry 47.1 (1978): 251-276."
cf_beta = {
	"A": 0.83,
	"C": 1.19,
	"D": 0.54,
	"E": 0.37,
	"F": 1.38,
	"G": 0.75,
	"H": 0.87,
	"I": 1.6,
	"K": 0.74,
	"L": 1.3,
	"M": 1.05,
	"N": 0.89,
	"P": 0.55,
	"Q": 1.11,
	"R": 0.93,
	"S": 0.75,
	"T": 1.19,
	"V": 1.7,
	"W": 1.19,
	"Y": 1.47
	}

#"description": "Chou-Fasman turn propensity",
#"refs": "Chou, Peter Y., and Gerald D. Fasman. 'Empirical predictions of protein conformation.' Annual review of biochemistry 47.1 (1978): 251-276.",
cf_turn = {
	"A": 0.66,
	"C": 1.19,
	"D": 1.46,
	"E": 0.74,
	"F": 0.6,
	"G": 1.56,
	"H": 0.95,
	"I": 0.47,
	"K": 1.01,
	"L": 0.59,
	"M": 0.6,
	"N": 1.56,
	"P": 1.52,
	"Q": 0.98,
	"R": 0.95,
	"S": 1.43,
	"T": 0.96,
	"V": 0.5,
	"W": 0.96,
	"Y": 1.14
	}
#"description": "side chain charge (when ionized)"

charge = {
	"A": 0.0,
	"C": -1.0,
	"D": -1.0,
	"E": -1.0,
	"F": 0.0,
	"G": 0.0,
	"H": 1.0,
	"I": 0.0,
	"K": 1.0,
	"L": 0.0,
	"M": 0.0,
	"N": 0.0,
	"P": 0.0,
	"Q": 0.0,
	"R": 1.0,
	"S": 0.0,
	"T": 0.0,
	"V": 0.0,
	"W": 0.0,
	"Y": -1.0
	}

#"description": "sidechain volume",
#"notes": "",
#"refs": "Richards, Frederic M. 'Areas, volumes, packing, and protein structure.' Annual review of biophysics and bioengineering 6.1 (1977): 151-176."
volume = {
	"A": 67.0,
	"C": 86.0,
	"D": 91.0,
	"E": 109.0,
	"F": 135.0,
	"G": 48.0,
	"H": 118.0,
	"I": 124.0,
	"K": 135.0,
	"L": 124.0,
	"M": 124.0,
	"N": 96.0,
	"P": 90.0,
	"Q": 114.0,
	"R": 148.0,
	"S": 73.0,
	"T": 93.0,
	"V": 105.0,
	"W": 163.0,
	"Y": 141.0
    }
#----------------------------------------------------------------------------------------------
# Dictionaries
#----------------------------------------------------------------------------------------------

dic = {}

dicio_frags2 = {}

dicio_frags3 = {}

dic2 = {}

dicio_frags_sel = {}

dicio_frags1 = {}

dic3 = {}
