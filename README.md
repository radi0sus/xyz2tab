# xyz2tab
A Python 3 script for printing tables of bond lengths and angles from xyz files to the console. The script furthermore calculates average values, including a variety of statistical parameters, and is able to group bonding parameters. Named atoms or elements can be excluded from bond or angle tables. Contacts of two or more atoms can be included. The output should result in nicely rendered mark down tables. 

## External modules
 `gemmi`,  `pandas`, `numpy`, `scipy`, `tabulate`
 
## Quick start
 Start the script with:
```console
python3 xyz2tab.py filename.xyz
```
to open the XYZ. It gives the following output:

Tables with general information. Please note that the covalent radius (column `Cov. radius`) has been increased by `Covalent Radius +  = 11.15 %` and bonds have been calculated by the sum of the radii given in the column `Cov. radius +`.
```
-------------------  ------------
Filename          :  asa.xyz
Number of atoms   :  21
Sum formula       :  C₉H₈O₄
Formula weight    :  180.16 g/mol
Excluded atoms    :  None
Excluded elements :  None
Included contacts :  None
Covalent Radius + :  11.15 %
-------------------  ------------

| Element   |   Atom count |   Mass fraction /% |   Cov. radius /Å |   Cov. radius + /Å |
|-----------|--------------|--------------------|------------------|--------------------|
| C         |            9 |              60.00 |             0.73 |               0.81 |
| H         |            8 |               4.48 |             0.31 |               0.34 |
| O         |            4 |              35.52 |             0.66 |               0.73 |
```
A table of bond lengths. Please note that the number after the element indicates the position in the xyz file. The numbering starts with zero due to [ORCA](https://orcaforum.kofo.mpg.de) conventions.
```
| Atoms   |   Bond length /Å |
|---------|------------------|
| C0–C2   |           1.3983 |
| C0–C5   |           1.4051 |
| ...     | ....             |
```
A table with summarized general bond lengths.
```
| Atoms   | Bond lengths /Å   |
|---------|-------------------|
| C–C     | 1.3933 - 1.4993   |
| C–H     | 1.0865 - 1.0945   |
| C–O     | 1.2189 - 1.3969   |
| O–H     | 0.9806            |
```
A table with statistical parameters.
```
| Atoms   |   Count |   Mean /Å |   Median /Å |   Sam. std. dev. |   Pop. std. dev. |   Std. error |   Skewness |
|---------|---------|-----------|-------------|------------------|------------------|--------------|------------|
| C–C     |       8 |    1.4221 |      1.4009 |           0.0436 |           0.0408 |       0.0154 |     1.4385 |
| C–H     |       7 |    1.09   |      1.0893 |           0.0032 |           0.003  |       0.0012 |     0.289  |
| C–O     |       5 |    1.3147 |      1.3438 |           0.0876 |           0.0783 |       0.0392 |    -0.3848 |
| O–H     |       1 |    0.9806 |      0.9806 |         nan      |           0      |     nan      |   nan      |
```
Four tables with angles in a likewise manner.

Start the script with:
```console
python3 cifpal.py filename.xyz > filename.md
```
will save the output in markdown format.

Convert markdown to docx (install [PANDOC](https://pandoc.org) first):
```console
pandoc filename.md -o filename.docx
```
This will convert the markdown file to a docx file. Open it with your favorite
word processor. Convert the file to even more formats such as HTML, PDF or TeX with PANDOC.

## Command-line options
- `filename` , required: filename, e.g. `my_xyz.xyz`, first two lines will be ignored, file format must be `element x y z`, cartesian coordinates, units in Å
- `-ea` `atom(s)`, optional: exclude atoms, e.g. `-ea H18` exclude bonds to H18, `-ea H18 H19` exclude bonds to H18 and H19
- `-ee` `elements(s)`,  optional: exlude elements,  e.g. `-ee H` exclude bonds to hydrogen atoms, `-ea H O` exclude bonds to hydrogen and oxygen atoms
- `-sa`, optional: sort values for bond lengths and angles ascending
- `-sd`, optional: sort values for bond lengths and angles descending
- `-sae`, optional: ascending alphabetical sort of elements
- `-sde`, optional:  descending alphabetical sort of elements
- `-ic` `atoms`, optional: include contacts of named atoms, e.g. `-ic O10 O11`, include the distance O10-O11, also include the angles X-O10-O11 and X-O11-O10. Input of more than two atoms is possible, e.g. `-ic O10 O11 O12`
- `-r` `N`, increase the covalent radii by `N` %, e.g.  `-r 20.1`, increase the covalent radii by 20.1 %. The default `N` is `11.15` %. The covalent radii used for the calculation of the bond length (bond length = rA + rB) is given in the last column of the summary table.
