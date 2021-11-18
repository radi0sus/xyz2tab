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

A tables with general information.
```
-------------------  ------------
Filename          :  aspirine.xyz
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
A table of bond lengths. Please note that the number after the element indicates the position in the xyz file. The numbering starts with zero due to ORCA conventions.
```
| Atoms   |   Bond length /Å |
|---------|------------------|
| C0–C2   |           1.3983 |
| C0–C5   |           1.4051 |
| C0–H14  |           1.0893 |
| ...     | ....             |
```

