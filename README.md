# encrypted
![Encrypted Logo](https://github.com/bposantos/encrypted/blob/main/encrypted-logo.png?raw=true)
Code for cryptide's sequences identification.

Quote:
>Santos BPO, Alves ESF, Ferreira CS, Ferreira-Silva A, Góes-Neto A, Verly RM, Lião LM, Oliveira SC, de Magalhães MTQ. Schistocins: Novel antimicrobial peptides encrypted in the Schistosoma mansoni Kunitz Inhibitor SmKI-1. Biochim Biophys Acta Gen Subj. 2021 Nov;1865(11):129989. [doi: 10.1016/j.bbagen.2021.129989](https://www.sciencedirect.com/science/article/pii/S0304416521001483). Epub 2021 Aug 10. PMID: 34389467.

For usage, python3 is necessary.

### Example: 
```
python3 encrypted.py -f fastafile.fasta -t protein -cl 2 -l 500 3000 -act antibacterial -n 1
```

## Decoding the arguments:
- **File (-f):** a fasta or multifasta file.

- **Activity types (-act):** in the moment, only antibacterial is available.

- **Cleavage types (-cl):** 1 for argCproteinase; 2 for trypsinKproteinase; 3 for trypsinRproteinase, 4 and 5 for chymotrypsin, 6 pepsin and 7 for trypsinGKproteinase.

- **Length (-l):** two enters are necessary, being the first the lower peptides length and the second, the higher length (in Daltons).

- **Type (-t):** 'protein' or 'peptide'. The first option cleavages the sequence(s) while the second jumps this step.

- **Number (-n):** a way to separate different predictions. You can type any number.

- **Help (-help):** print in the screen the information above.

## Exit files:

CSV table with predicted peptides and their physicochemical properties;

CSV table with MS/MS peptides fragmentation (mass spectrometry data).

### Plotting

If the user wants to visualize the data, a three dimensional plot is additionally available.

- **File (-f):** a csv file originated by encrypted.py.

- **Number of clusters (-cl):** kmeans clusterization. User defined number of clusters.

- **Output (-o):** name of the output file.


### Example: 
```
python3 pca_v2.py -f dataframe.csv -cl 3 -o clusters01
```
