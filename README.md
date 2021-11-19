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
- `-r` `N`, increase the covalent radii by `N` %, e.g.  `-r 20.1`, increase the covalent radii by 20.1 %. The default `N` is `11.15` %. The covalent radii used for the calculation of the bond length (bond length of A-B = rA + rB) is given in the last column of the summary table (`Cov. radius +`).
- `-v`, optional:  include two more tables (tables with general bond lengths and angles)

## Statistics
Statistics are derived from the values of the bonding parameters. 

Sam. std. dev. = Sample standard deviation, Pop. std. dev. = Population standard deviation, Std. error = Standard error or standard error of mean. Please refer to literature or Wikipedia for the meaning of these terms. The population standard deviation is probably the value you are looking for.

## Remarks
- The format of the tabular output can be easily changed in the script using another formatting option of the `tabulate` module.
- The `gemmi` library is only needed for covalent radii and molecular weigths.

## Known Issues
- The script makes extensive use of Unicode characters, which can cause problems with output or conversion.
- Verbose output (`-v` option) can be very large and confusing (look nicer after formatting).

## Examples

### Example 1:
```console
python3 cifpal.py asa.xyz
```
Open `asa.xyz` and show tables.


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

| Atoms   |   Bond length /Å |
|---------|------------------|
| C0–C2   |           1.3983 |
| C0–C5   |           1.4051 |
| C0–H14  |           1.0893 |
| C1–C3   |           1.3936 |
| C1–C6   |           1.4006 |
| C1–H15  |           1.0865 |
| C2–C3   |           1.3933 |
| C2–H16  |           1.0875 |
| C3–H17  |           1.0871 |
| C4–C8   |           1.4993 |
| C4–H18  |           1.0929 |
| C4–H19  |           1.0945 |
| C4–H20  |           1.0924 |
| C5–C6   |           1.4013 |
| C5–C7   |           1.4850 |
| C6–O12  |           1.3969 |
| C7–O9   |           1.2189 |
| C7–O11  |           1.3438 |
| C8–O10  |           1.2241 |
| C8–O12  |           1.3900 |
| O11–H13 |           0.9806 |

| Atoms   | Bond lengths /Å   |
|---------|-------------------|
| C–C     | 1.3933 - 1.4993   |
| C–H     | 1.0865 - 1.0945   |
| C–O     | 1.2189 - 1.3969   |
| O–H     | 0.9806            |

| Atoms   |   Count |   Mean /Å |   Median /Å |   Sam. std. dev. |   Pop. std. dev. |   Std. error |   Skewness |
|---------|---------|-----------|-------------|------------------|------------------|--------------|------------|
| C–C     |       8 |    1.4221 |      1.4009 |           0.0436 |           0.0408 |       0.0154 |     1.4385 |
| C–H     |       7 |    1.09   |      1.0893 |           0.0032 |           0.003  |       0.0012 |     0.289  |
| C–O     |       5 |    1.3147 |      1.3438 |           0.0876 |           0.0783 |       0.0392 |    -0.3848 |
| O–H     |       1 |    0.9806 |      0.9806 |         nan      |           0      |     nan      |   nan      |

| Atoms      |   Angle /° |
|------------|------------|
| C2–C0–C5   |     120.83 |
| C2–C0–H14  |     118.97 |
| C5–C0–H14  |     120.20 |
| C3–C1–C6   |     119.91 |
| C3–C1–H15  |     119.83 |
| C6–C1–H15  |     120.26 |
| C0–C2–C3   |     120.08 |
| C0–C2–H16  |     119.86 |
| C3–C2–H16  |     120.06 |
| C1–C3–C2   |     119.90 |
| C1–C3–H17  |     120.06 |
| C2–C3–H17  |     120.04 |
| C8–C4–H18  |     109.89 |
| C8–C4–H19  |     109.30 |
| C8–C4–H20  |     109.85 |
| H18–C4–H19 |     108.47 |
| H18–C4–H20 |     110.68 |
| H19–C4–H20 |     108.62 |
| C0–C5–C6   |     118.30 |
| C0–C5–C7   |     116.52 |
| C6–C5–C7   |     125.18 |
| C1–C6–C5   |     120.97 |
| C1–C6–O12  |     116.28 |
| C5–C6–O12  |     122.70 |
| C5–C7–O9   |     125.02 |
| C5–C7–O11  |     114.98 |
| O9–C7–O11  |     120.00 |
| C4–C8–O10  |     124.34 |
| C4–C8–O12  |     108.76 |
| O10–C8–O12 |     126.89 |
| C7–O11–H13 |     103.14 |
| C6–O12–C8  |     112.52 |

| Atoms   | Angle /°        |
|---------|-----------------|
| C–C–C   | 116.52 - 125.18 |
| C–C–H   | 109.30 - 120.26 |
| H–C–H   | 108.47 - 110.68 |
| C–C–O   | 108.76 - 125.02 |
| O–C–O   | 120.00 / 126.89 |
| C–O–H   | 103.14          |
| C–O–C   | 112.52          |

| Atoms   |   Count |   Mean /° |   Median /° |   Sam. std. dev. |   Pop. std. dev. |   Std. error |   Skewness |
|---------|---------|-----------|-------------|------------------|------------------|--------------|------------|
| C–C–C   |       8 |    120.21 |      119.99 |             2.48 |             2.32 |         0.88 |       0.82 |
| C–C–H   |      11 |    117.12 |      119.86 |             4.79 |             4.57 |         1.45 |      -1.17 |
| H–C–H   |       3 |    109.26 |      108.62 |             1.24 |             1.01 |         0.72 |       1.7  |
| C–C–O   |       6 |    118.68 |      119.49 |             6.42 |             5.86 |         2.62 |      -0.61 |
| O–C–O   |       2 |    123.44 |      123.44 |             4.88 |             3.45 |         3.45 |     nan    |
| C–O–H   |       1 |    103.14 |      103.14 |           nan    |             0    |       nan    |     nan    |
| C–O–C   |       1 |    112.52 |      112.52 |           nan    |             0    |       nan    |     nan    |

### Example 2:
```console
python3 cifpal.py 2103396.cif Cu1 Cu2 -ee K -f 12
