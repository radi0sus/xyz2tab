# xyz2tab
A Python 3 script for printing tables of bond lengths and angles from xyz files to the console. The script furthermore calculates average values, including a variety of statistical parameters, and is able to group bonding parameters. Named atoms or elements can be excluded from bond or angle table. Contacts of two or more atoms can be included. The output should result in nicely rendered mark down tables. 

## External modules
 `gemmi`,  `pandas`, `numpy`, `scipy`, `tabulate`
 
## Quick start
 Start the script with:
```console
python3 xyz2tab.py filename.xyz
```
to open the XYZ. It gives the following output:
