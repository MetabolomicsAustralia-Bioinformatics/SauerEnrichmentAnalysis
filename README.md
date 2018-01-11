# SauerEnrichmentAnalysis

Implementation of the metabolic flux analysis procedures described by Sauer and co-workers (https://www.ncbi.nlm.nih.gov/pubmed/17035687).

## Usage
Two different input files are required.

1. A file defining the name and chemical formula of the analyte molecule(s) as well as the number of atoms that can possibly be exchanged.
2. A file providing the intensities for each isotopic peak for each analyte defined in file 1 in each sample.

Sample input and output files are provided in the 'sample_files' directory.

```
python proc.py
```

**Note the following points:**
1. The number of possible enrichment sites for an analyte defined in file 1 **must** match the number of isotope data series' provided in file 2.
2. The names of metabolites defined in file 1 must exactly match those in file 2.
3. In file 2, the headers for the isotope data series must be given exactly as follows:
~~~
deoxyribose M0 Results, deoxyribose M+1 Results,	deoxyribose M+2 Results etc...
~~~
   Where 'M9' represents the monoisotopic peak, 'M+1' is the first isotope peak, 'M+2' is the second isotope peak and so on. The current script is very sensitive to the formatting of these headers so be sure to adapt the example files carefully.
4. When specifying chemical formulae in file 1:
   - Place spaces between entries for different elements. Eg. C6 H12 O6 rather than C6H12O6.
   - If there is only one atom of a given element, it must be entered as X1. Eg. C2 O2 H5 N1 rather than C2 O2 H5 N.
   - An internal atom type dictionary seems to go looking for Si atoms. The script crashes if Si0 is not added for a compound with no silicione. Fix this in a future update.
