# encrypted
Code for cryptide's sequences identification.

For usage, python3 is necessary.

### Example: 
```
python3 encrypted.py -f fastafile.fasta -t protein -cl 2 -l 500 3000 -act antibacterial -n 1
```

## Decoding the arguments:
**File (-f):** a fasta or multifasta file.

**Activity types (-act):** in the moment, only antibacterial is available.

**Cleavage types (-cl):** 1 for argCproteinase; 2 for trypsinKproteinase; 3 for trypsinRproteinase, 4 and 5 for chymotrypsin, 6 pepsin and 7 for trypsinGKproteinase.

**Length (-l):** two enters are necessary, being the first the lower peptides length and the second, the higher length (in Daltons).

**Type (-t):** 'protein' or 'peptide'. The first option cleavages the sequence(s) while the second jumps this step.

**Number (-n):** a way to separate different predictions. You can type any number.

**Help (-help):** print in the screen the information above.

## Exit files:

CSV table with predicted peptides and their physicochemical properties;

CSV table with MS/MS peptides fragmentation (mass spectrometry data).

### Plotting

If the user wants to visualize the data, a three dimensional plot is additionally available. The csv file name needs to be added to the code (pca_01.py) line 9, under '    '.
