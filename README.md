# SauerEnrichmentAnalysis

* This is a refactor of the original Sauer enrichment analysis tool, migrated to Python3.
* Implementation of the metabolic flux analysis procedures described by Sauer and co-workers (https://www.ncbi.nlm.nih.gov/pubmed/17035687).

## Usage
Two different input files are required.

1. **file1**: A file defining the name and chemical formula of the analyte molecule(s) as well as the number of atoms that can possibly be exchanged.
2. **file2**: A file providing the intensities for each isotopic peak for each analyte defined in file 1 in each sample.

Sample input and output files are provided in the 'sample_files' directory.

```
python3 proc.py
```

**Note the following points:**
1. The number of possible enrichment sites for an analyte defined in file 1 **must** match the number of isotope data series' provided in file 2. Note that the M0 isotope is not included in this number.
   - For example, if M0, M+1 and M+2 are defined for a molecule, the number of enrichment sites specified in file 1 must be 2.
2. The names of metabolites defined in file 1 must exactly match those in file 2.
3. In file 2, the headers for the isotope data series must be given exactly as follows:

```
deoxyribose M0 Results, deoxyribose M+1 Results, deoxyribose M+2 Results etc...
```

   Where 'M0' represents the monoisotopic peak, 'M+1' is the first isotope peak, 'M+2' is the second isotope peak and so on. The current script is very sensitive to the formatting of these headers so be sure to change the example files carefully.
4. When specifying chemical formulae in file 1:
   - Place spaces between entries for different elements. Eg. C6 H12 O6 rather than C6H12O6.
   - If there is only one atom of a given element X, it must be entered as X1. Eg. C2 O2 H5 N1 rather than C2 O2 H5 N.
   - An internal atom type dictionary seems to go looking for Si atoms. The script crashes if Si0 is not added for a compound with no silicione. Fix this in a future update.
5. No cells in the file 2 input can be blank. Must replace with 0.

### To-do
1. Add command-line arg parsing
2. Try to increase flexibility of input file type definitions
  * change reading logic to regex, instead of looking for ` M`.
3. Sort out bug for Si0 compound formulae definitions (pt. 4 above)
4. Add Python Dash front end?
