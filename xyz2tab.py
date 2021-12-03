#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#xyz2tab

import sys                                                      #stdout
import argparse                                                 #argument parser
import re                                                       #regex
import itertools                                                #for r-length tuples, in sorted order, no repeated elements
import pandas as pd                                             #pandas tables
import numpy as np                                              #for calculations
from scipy.spatial.distance import pdist, squareform, cosine    #for the calculations of the distance matrix and angles (cosine)
from tabulate import tabulate                                   #nice table output

#for windows console
sys.stdout.reconfigure(encoding='utf-8')  

#pd.set_option("display.max_rows", None, "display.max_columns", None)

#+x % to covalence radius
#radii_ext = 8 / 100 

#covalent radii from Alvarez (2008)
#DOI: 10.1039/b801115j
covalent_radii = {
	'H': 0.31, 'He': 0.28, 'Li': 1.28,
	'Be': 0.96, 'B': 0.84, 'C': 0.76, 
	'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
	'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
	'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
	'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
	'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
	'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
	'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
	'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
	'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
	'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
	'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
	'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
	'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
	'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
	'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
	'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
	'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
	'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
	'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
	'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
	'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
	'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
	'Am': 1.80, 'Cm': 1.69
}

#atomic weights
atomic_weights = {
	'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,
	'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
	'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,
	'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,
	'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,
	'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,
	'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,
	'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,
	'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,
	'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,
	'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,
	'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,
	'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,
	'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,
	'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,
	'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,
	'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,
	'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,
	'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,
	'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,
	'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,
	'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,
	'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,
	'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247
}

#dict for numbers to subscript numbers
utf_sub_dict = {
	"0" : "₀",
	"1" : "₁",
	"2" : "₂",
	"3" : "₃",
	"4" : "₄",
	"5" : "₅",
	"6" : "₆",
	"7" : "₇",
	"8" : "₈",
	"9" : "₉",
}

#numbers to subscript (utf8) numbers
def num_to_subnum(number):
	utf_number=''
	for letter in str(number):
		utf_letter=utf_sub_dict[letter]
		utf_number=utf_number+utf_letter
	return(utf_number)

#removes the upper triangle of the distance matrix and zeros
#e.g., from d(A-B) = 1.234 Å = d(B-A) =1.234 Å, d(B-A) will be removed 
#d(A-B) = 0 Å will be removed as well
def dm_to_series1(df):
	df = df.astype(float) # do not comment this, angle list will be incomplete
	df.values[np.triu_indices_from(df, k=1)] = np.nan
	#replace zeros with nan
	df = df.replace(0, np.nan)
	#return and drop all nan
	return df.unstack().dropna()

#calculate angle from 3 vectors / atomic coordinates: i(x,y,z); j(x,y,z); k(x,y,z) 
#xyzarr is the array of all atomic coordinates
def calc_angle(xyzarr, i, j, k):
	rij = xyzarr[i] - xyzarr[j]
	rkj = xyzarr[k] - xyzarr[j]
	#remove if cosine fails
	#cos_theta = np.dot(rij, rkj)
	#sin_theta = np.linalg.norm(np.cross(rij, rkj))
	#theta = np.arctan2(sin_theta, cos_theta)
	#scipy pdist cosine instead of the 3 lines above
	theta = cosine(rij,rkj)
	theta = np.arccos(1-theta)
	return 	np.degrees(theta)

#calculate the dihedral angle from 4 vectors / atomic coordinates: i(x,y,z); j(x,y,z); k(x,y,z); l(x,y,z)
def calc_d_angle(xyzarr, i, j, k, l):
	#no warning if division by zero
	np.seterr(invalid='ignore')
	rji = -1*(xyzarr[j] - xyzarr[i])
	rkj = xyzarr[k] - xyzarr[j]
	rlk = xyzarr[l] - xyzarr[k]
	rkj /= np.linalg.norm(rkj)
	v = rji - np.dot(rji, rkj)*rkj
	w = rlk - np.dot(rlk, rkj)*rkj
	x = np.dot(v, w)
	y = np.dot(np.cross(rkj, v), w)
	return 	np.degrees(np.arctan2(y,x))

#calculation of the best-fit plane
#https://gist.github.com/bdrown/a2bc1da0123b142916c2f343a20784b4
def svd_fit(X):
	C = np.average(X, axis=0)
	# Create CX vector (centroid to point) matrix
	CX = X - C
	# Singular value decomposition
	U, S, V = np.linalg.svd(CX)
	# The last row of V matrix indicate the eigenvectors of
	# smallest eigenvalues (singular values).
	N = V[-1]
	return C, N

#argument parser
parser = argparse.ArgumentParser(
		prog='xyz2tab', 
		description = "Print bond, lengths angles and more from xyz files.")

#filename is required
parser.add_argument("filename", 
	help = "filename, xyz; e.g. mymolecule.xyz")

#exclude atoms
parser.add_argument('-ea','--excludeAt',
	nargs="+",
	type=str,
	help='exclude bonds and angles to specified atoms; e.g. -ea N1 or -ea N1 N2')

#exclude elements
parser.add_argument('-ee','--excludeEl',
	nargs="+",
	type=str,
	help='exclude bonds and angles to specified elements; e.g. -ee C or -ee C N')

#sort by value
parser.add_argument('-sa','--sortasc',
	default=0,
	action='store_true',
	help='sort values for bond lengths and angles ascending')

#sort by value
parser.add_argument('-sd','--sortdes',
	default=0, 
	action='store_true',
	help='sort values for bond lengths and angles descending')

#sort by name
parser.add_argument('-sae','--sortascEl',
	default=0, 
	action='store_true',
	help='ascending alphabetical sort of elements')

#sort by name
parser.add_argument('-sde','--sortdesEl',
	default=0, 
	action='store_true',
	help='descending alphabetical sort of elements')

#include contacts
parser.add_argument('-ic','--includeCon',
	nargs="+",
	type=str,
	help='include contacts, e.g. -ic C0 C11 or -ic C0 C11 N14')

#calculate dihedral angle of selected atoms
parser.add_argument('-d','--dihedral',
	nargs=4,
	type=str,
	help='calculate the dihedral angle of 4 atoms, e.g. -d C0 C11 C12 N13')
	
#calculate the best plane 1
parser.add_argument('-p1','--plane1',
	nargs='+',
	type=str,
	help='calculate the best plane through selected (A B C), a range of (A : C) or all atoms (:)\
		and show the distances to plane 1, \
		e.g. -p1 C0 C11 C12 N13 or -p1 C0 : N13 or -p1 :')

#calculate the best plane 2
parser.add_argument('-p2','--plane2',
	nargs='+',
	type=str,
	help='calculate the best plane through selected (A B C), a range of (A : C) or all atoms (:)\
		and show the distances to plane 1 and 2 and the angle to plane 1, \
		e.g. -p2 C21 C22 C23 N23 or -p2 C21 : N23 or -p2 :')

#add +x% to radius
parser.add_argument('-r','--radius',
	default=8,
	type=float,
	help='enlarge atomic radii by x %%, e.g. -r 15.2, default is 8 %%')

#verbose
parser.add_argument('-v','--verbose',
	default=0, action='store_true',
	help='verbose print, includes 2 more tables')

#parse arguments
args = parser.parse_args()

#read xyz into data frame
#skip first two rows of the xyz file
#only XMol xyz is supportet, atom(as element) x y z, e.g. C 1.58890 -1.44870 -0.47000
try:
	xyz_df = pd.read_csv(args.filename, 
			 delim_whitespace=True, 
			 skiprows=2, 
			 names=["element", "x", "y", "z"])
#file not found
except IOError:
	print(f"'{args.filename}'" + " not found")
	sys.exit(1)

#xx.x% --> 0.xxx
radii_ext = args.radius / 100

#element + position in xyz = atom + number, e.g. C --> C0; change to C(0) can be arranged here
xyz_df['atom1_idx'] = ["{}{}".format(atm, idx) for atm, idx in zip(xyz_df.element, xyz_df.index.array)]
#atom 1 is atom 2
xyz_df['atom2_idx'] = xyz_df['atom1_idx']
#atomic weight gemmi
#xyz_df['weight'] = xyz_df['element'].apply(lambda x: gemmi.Element(x).weight)
#atomic weight dict
xyz_df['weight'] = xyz_df['element'].apply(lambda x: atomic_weights[x])
#covalent radius gemmi
#xyz_df['cov_radius'] = xyz_df['element'].apply(lambda x: gemmi.Element(x).covalent_r)
#covalent radius dict
xyz_df['cov_radius'] = xyz_df['element'].apply(lambda x: covalent_radii[x])
#reorder data frame
xyz_df=xyz_df[['atom1_idx','atom2_idx','element','x','y','z','weight','cov_radius']]
#formula weight from sum of atomic weights
fw = xyz_df['weight'].sum()

#get total number of each element by counting elements and grouping
grouped=xyz_df.groupby('element',sort=False)['atom1_idx'].count()
#results of grouping into new data frame
info_df=grouped.reset_index()
#numbers to utf8 subscript number, e.g. 2 --> ₂ for sum formula
info_df['utf_num']=info_df['atom1_idx'].apply(lambda x: num_to_subnum(x))
#replace '1' with '', e.g. C1O2 --> CO2
info_df.loc[info_df['utf_num'] == '₁', 'utf_num'] = ''
#sort elements alphabetically
info_df.sort_values(by=['element'], inplace=True)
#generate the sum formula out of elements and the number of each element, e.g. C 5 H 12 --> C5H12
sum_formula_list=list("{}{}".format(element, number) for element, number in zip(info_df.element,info_df.utf_num))
#drop the utf8 subscript number column
info_df.drop(columns=['utf_num'], inplace=True)
#calculate mass fraction of each element
#gemmi
#info_df['atom%']=info_df.apply(lambda x: gemmi.Element(x.element).weight*x.atom1_idx/fw*100, axis=1)
#dict
info_df['atom%']=info_df.apply(lambda x: atomic_weights[x.element]*x.atom1_idx/fw*100, axis=1)
#atomic radii from the gemmi table
#gemmi
#info_df['radii']=info_df.apply(lambda x: gemmi.Element(x.element).covalent_r, axis=1)
#dict
info_df['radii']=info_df.apply(lambda x: covalent_radii[x.element], axis=1)
#atomic radii + x% used for calculation of bonds
info_df['radii_plus']=info_df.apply(lambda x: x.radii + x.radii*radii_ext, axis=1)

#print summary table
print('')
print(tabulate([['Filename          :', args.filename],
				['Number of atoms   :', xyz_df.shape[0]],
	            ['Sum formula       :', ''.join(sum_formula_list)],
				['Formula weight    :', '{:.2f} g/mol'.format(fw)],
				['Excluded atoms    :', re.sub(r'[^a-zA-Z0-9,]','',str(args.excludeAt))],
				['Excluded elements :', re.sub(r'[^a-zA-Z0-9,]','',str(args.excludeEl))],
				['Included contacts :', re.sub(r'[^a-zA-Z0-9,]','',str(args.includeCon))],
				['Covalent Radius + :', '{:.2f} %'.format(args.radius)]],
				tablefmt='simple'))

#print info table
print('')
print(tabulate(info_df,
			headers=['Element','Atom count','Mass fraction /%', 'Cov. radius /Å', 'Cov. radius + /Å'],
			tablefmt='github',
			floatfmt=(".2f"),
			showindex=False))

#calculate the full distance matrix & put to square form, e.g.:
#
#   C0  C1  C2
#C0 0.0 1.1 2.3
#C1 1.1 0.0 1.5
#C2 2.3 1.5 0.0
#iloc [:,3:6] contains xyz coordinates
dist_mat_full=pd.DataFrame(squareform(pdist(xyz_df.iloc[:,3:6],'euclid')),
			  columns = xyz_df[['atom1_idx','element','cov_radius']],
			  index = xyz_df[['atom2_idx','element','cov_radius']])

#remove the upper triangle and zeros, e.g.:
#   C0  C1  C2
#C0 NaN NaN NaN
#C1 1.1 NaN NaN
#C2 2.3 1.5 NaN
dist_mat_red = dm_to_series1(dist_mat_full)
#bring it to a "normal" form
dist_mat_red = dist_mat_red.reset_index(level=[1])
#bring it to a "normal" form, distance matrix --> disctance data frame
dist_df=pd.DataFrame(dist_mat_red.to_records())
#bring it to a "normal" form ...
dist_df[['atom1_idx','element1','cov_radius1']]=pd.DataFrame(dist_df['index'].tolist(), index=dist_df.index)
dist_df[['atom2_idx','element2','cov_radius2']]=pd.DataFrame(dist_df['level_1'].tolist(), index=dist_df.index)
dist_df.drop(['index', 'level_1'], axis=1,inplace=True)
dist_df.rename(columns={'0':'distance_calc'},inplace=True)
#reorder data frame
dist_df=dist_df[['atom1_idx','element1','cov_radius1','atom2_idx','element2','cov_radius2','distance_calc']]

#column with the sum of the atomic radii from two elements /atoms + x%
dist_df['distance_radii'] = (dist_df['cov_radius1'] + 
							dist_df['cov_radius2'])+(dist_df['cov_radius1'] + 
							dist_df['cov_radius2'])*radii_ext

#distance is considered as bond if the calculated distance is smaller than the sum of the atomic radii
#set to 'True' if this is a bond
dist_df['is_bond']=(dist_df['distance_calc'] < dist_df['distance_radii'])

#include distances from selected atom (pairs) from args include connections
#sets 'is_bond' to 'True' if 'False'
if args.includeCon:
	#must have atoms 1 and 2 and 'is_bond' must be 'False'
	dist_df.loc[(dist_df.atom1_idx.isin(args.includeCon) & 
				dist_df.atom2_idx.isin(args.includeCon) & 
			   ~dist_df.is_bond), 'is_bond'] = True
	#np variant of the above, slightly slower
	#dist_df['is_bond'] = np.where((dist_df.atom1_idx.isin(args.includeCon)) &
	#(dist_df.atom2_idx.isin(args.includeCon) & (~dist_df.is_bond)), True, dist_df['is_bond'])
	
#fusion char, A B --> A-B
dist_df['fusion_char']='–'

#A0 B1 - --> A0-B1
dist_df['A-B'] = dist_df['atom1_idx']+dist_df['fusion_char']+dist_df['atom2_idx']

#A B - --> A-B
dist_df['El1-El2'] = dist_df['element1']+dist_df['fusion_char']+dist_df['element2']

#A B --> AB
dist_df['ElEl'] = dist_df['El1-El2'].apply(lambda x: ''.join(sorted(x)))

#dist_df data frame --> all_dist data frame if 'is_bond' is True
all_dist = pd.DataFrame(dist_df[(dist_df.is_bond == True)])
#all_dist data frame --> sel_dist data frame
sel_dist=all_dist
#sel_dist.reset_index(drop=True,inplace=True)

############ Exclude

#exclude named atoms from input (in Atom1 and Atom2)
if args.excludeAt:
	sel_dist=all_dist[~all_dist.atom1_idx.isin(args.excludeAt) & ~all_dist.atom2_idx.isin(args.excludeAt)]
	
#exclude named elements from input (in Element1 and Element2)
if args.excludeEl:
	sel_dist=all_dist[~all_dist.element1.isin(args.excludeEl) & ~all_dist.element2.isin(args.excludeEl)]

#exit if selected bonds data frame is empty 
if len(sel_dist) == 0:
	print("No bonds found. Include more atoms or elements. Exit.")
	sys.exit(1)

############ Sort

#sort bond length values ascending
if args.sortasc:
	sel_dist=sel_dist.sort_values(by=['distance_calc'])
	
#sort bond length values descending
if args.sortdes:
	sel_dist=sel_dist.sort_values(by=['distance_calc'],ascending=False)
	
#sort by elements ascending, A --> Z (not PSE like)
if args.sortascEl:
	sel_dist=sel_dist.sort_values(by=['element1','element2','distance_calc','atom1_idx','atom2_idx'])
	
#sort by elements descending, A --> Z (not PSE like)
if args.sortdesEl:
	sel_dist=sel_dist.sort_values(by=['element1','element2','distance_calc','atom1_idx','atom2_idx'],ascending=False)

############ Print

#table with all selected distances, A-B | 1.234 Å
pr_sel_dist=sel_dist[['A-B','distance_calc']]
print('')
print(tabulate(pr_sel_dist,
	  headers=['Atoms','Bond length /Å'],
	  tablefmt='github',
	  floatfmt=(".4f"),
	  showindex=False))

#lists for printed tables
summary_bond_table_1 = list()
summary_bond_table_2 = list()
summary_bond_table_3 = list()

#group El-El and distances by ElEl, e.g. C-N 1.234 Å by CN
grouped = sel_dist[['El1-El2','distance_calc']].groupby(sel_dist['ElEl'],sort=False)

#verbose table El1-El2 | bond length, e.g. C-C 1.223, 1.456, 1.511
for groups in grouped:
	summary_bond_table_1.append([groups[1].iloc[0].tolist()[0], 
		', '.join(groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist())])

#short table El1-El2 | bond length, e.g. C-C 1.223 - 1.511 (for >2), or C-C 1.223 / 1.511 (for 2), C-C 1.223 (for one)
# float to 4 decimals, e.g. 0.1234 
for groups in grouped:
	if len(groups[1]) == 1:
		summary_bond_table_2.append([groups[1].iloc[0].tolist()[0], 
			groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[0]])
	elif len(groups[1]) == 2:
		summary_bond_table_2.append([groups[1].iloc[0].tolist()[0], 
			groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[0] + 
			' / ' + groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[-1]])
	else:
		summary_bond_table_2.append([groups[1].iloc[0].tolist()[0], 
			groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[0] + 
			' - ' + groups[1].sort_values(by=['distance_calc']).distance_calc.apply(lambda x: '{:.4f}'.format(x)).tolist()[-1]])

#grouped = sel_dist[['El1-El2','distance_calc']].groupby(sel_dist['ElEl'],sort=False)

#generate table with statistics, | El1-El2 | Count | Mean | Median | Sam. std. dev. | Pop. std. dev. | Std. error |
for groups in grouped:
	summary_bond_table_3.append([groups[1].iloc[0].tolist()[0], groups[1].distance_calc.count(),f'{groups[1].distance_calc.mean():.4f}', \
		f'{groups[1].distance_calc.median():.4f}', f'{groups[1].distance_calc.std():.4f}', f'{groups[1].distance_calc.std(ddof=0):.4f}', \
		f'{groups[1].distance_calc.sem():.4f}',f'{groups[1].distance_calc.skew():.4f}'])

#print verbose table El1-El2 | bond length only by request
if args.verbose:
	print('')
	print(tabulate(summary_bond_table_1,
		  headers=['Atoms','Bond lengths /Å'], 
		  tablefmt='github',
		  floatfmt=(".4f"),
		  showindex=False))

#print short table
print('')
print(tabulate(summary_bond_table_2,
	  headers=['Atoms','Bond lengths /Å'],
	  tablefmt='github',
	  floatfmt=(".4f"),
	  showindex=False))

#print statistics table
print('')
print(tabulate(summary_bond_table_3,
	  headers=['Atoms','Count','Mean /Å', 'Median /Å','Sam. std. dev.', \
	  'Pop. std. dev.','Std. error','Skewness'], 
	  tablefmt='github',
	  showindex=False))

############ Angles

#in case of problems (lost angles) remove comment from next line ⇩
#dist_mat2=pd.DataFrame(squareform(pdist(xyz_df.iloc[:,3:6],'euclid')),columns = xyz_df[['atom1_idx','element','cov_radius']],index = xyz_df[['atom2_idx','element','cov_radius']])
# and comment next line here ⇩
#full distance matrix is needed for the angle calculation
dist_mat2 = dist_mat_full
dist_mat2 = dist_mat2.unstack()
dist_mat2 = dist_mat2.reset_index()
#copy to dist2_df data frame
dist2_df=dist_mat2

#recover 'atom1_idx'... etc from distance matrix
dist2_df[['atom1_idx','element1','cov_radius1']]=pd.DataFrame(dist2_df['level_0'].tolist(), index=dist2_df.index)
dist2_df[['atom2_idx','element2','cov_radius2']]=pd.DataFrame(dist2_df['level_1'].tolist(), index=dist2_df.index)
dist2_df.rename(columns={0:'distance_calc'},inplace=True)
dist2_df.drop(['level_0'], axis=1,inplace=True)
dist2_df.drop(['level_1'], axis=1,inplace=True)
#column with the sum of the atomic radii from two elements /atoms + x%
dist2_df['distance_radii'] = (dist2_df['cov_radius1'] + 
							  dist2_df['cov_radius2'])+(dist2_df['cov_radius1'] + 
							  dist2_df['cov_radius2'])*radii_ext
#distance is considered as bond if the calculated distance is smaller than the sum of the atomic radii
#set to 'True' if this is a bond
dist2_df['is_bond']=((dist2_df['distance_calc'] > 0) & (dist2_df['distance_calc'] < dist2_df['distance_radii']))

#include distances from selected atom (pairs) from args include connections
#sets 'is_bond' to 'True' if 'False'
if args.includeCon:
	#must have Atoms 1 and 2 and is_bond must be false
	dist2_df.loc[((dist2_df.distance_calc > 0) & 
				   dist2_df.atom1_idx.isin(args.includeCon) & 
				   dist2_df.atom2_idx.isin(args.includeCon) & 
				  ~dist2_df.is_bond), 'is_bond'] = True
	#np variant of the above, slightly slower
	#dist2_df['is_bond'] = np.where((dist2_df.distance_calc > 0) & 
	#(dist2_df.atom1_idx.isin(args.includeCon)) & 
	#(dist2_df.atom2_idx.isin(args.includeCon) & 
	#(~dist2_df.is_bond)), True, dist2_df['is_bond'])

#reoder data frame - maybe not necessary
dist2_df=dist2_df[['atom1_idx','cov_radius1','atom2_idx','cov_radius2','distance_calc','distance_radii','is_bond']]

#dist2_df data frame --> all_dist2 data frame if 'is_bond' is True
all_dist2 = pd.DataFrame(dist2_df[(dist2_df.is_bond == True)])
#all_dist2 data frame --> sel_dist2 data frame
sel_dist2 = all_dist2

#group atom2 by atom1,e.g. C0 C1, C0 C2, C0 C3...
group1 = sel_dist2.groupby('atom1_idx',sort=False)['atom2_idx']
#group2=sel_dist2.groupby('atom2_idx',sort=False)['atom1_idx']

#xyz array with coordinates of all atoms is needed
#needed for (dihedral) angle calculation
xyzarr = xyz_df.iloc[:,3:6].to_numpy()
#angle calculation is not as fast as bond length calculation

#make 4 empty lists
atom1 = list()
atom2 = list()
atom3 = list()
anglelist = list()

#middle atom 'B' (for angle A-B-C) is in name of the group
#get x,y,z coordinates (a2) of 'B' from xyz data frame
for name, group in group1:
	a2=xyz_df.index[xyz_df['atom1_idx'] == name].tolist()
	#'A' and 'C' atoms are in the group
	#get x,y,z coordinates of 'A' (a1) and 'C' (a3) from xyz data frame 
	for s in itertools.combinations(group,2):
		#very nice itertool, e.g.:
		#a1 (central atom) binds to a2, a3, a4, a5
		#angles will be a2-a1-a3, a2-a1-a4, a2-a1-a5, a3-a1-a4, a3-a1-a5, a4-a1-a5
		#exludes double entries like a3-a1-a2, a4-a1-a2,....
		a1=xyz_df.index[xyz_df['atom1_idx'] == s[0]].tolist()
		a3=xyz_df.index[xyz_df['atom1_idx'] == s[1]].tolist()
		#calculate the angle
		angle = calc_angle(xyzarr, *a1, *a2, *a3)
		#name of atom1 ('A') --> atom1 list
		atom1.append(s[0])
		#name of atom2 ('B') --> atom2 list
		atom2.append(name)
		#name of atom3 ('C') --> atom3 list
		atom3.append(s[1])
		#calculated angle to list of angles
		anglelist.append(angle)

#all 4 lists in angles_df data frame
angles_df=pd.DataFrame(({'atom1_idx': atom1, 'atom2_idx': atom2, 'atom3_idx': atom3, 'angle_calc': anglelist}))
#construct elements from atom names, e.g. C1 --> C, Fe13 --> Fe
angles_df['element1']=angles_df['atom1_idx'].apply(lambda x: re.sub(r'\d+','', x))
angles_df['element2']=angles_df['atom2_idx'].apply(lambda x: re.sub(r'\d+','', x))
angles_df['element3']=angles_df['atom3_idx'].apply(lambda x: re.sub(r'\d+','', x))
#fuse atom names A B C by '-' --> A-B-C
angles_df['fusion_char'] = '–'
#A0 B1 C2- --> A0-B1-C2
angles_df['A-B-C'] = angles_df['atom1_idx']+angles_df['fusion_char']+angles_df['atom2_idx']+angles_df['fusion_char']+angles_df['atom3_idx']
#A B C--> A-B-C
angles_df['El1-El2-El3'] = angles_df['element1']+angles_df['fusion_char']+angles_df['element2']+angles_df['fusion_char']+angles_df['element3']
#A (B) C--> A-C, for grouping
angles_df['El1-El3'] = angles_df['element1']+angles_df['fusion_char']+angles_df['element3']
#A-B-C--> ABC, for grouping
angles_df['ElElEl'] = angles_df['El1-El2-El3'].apply(lambda x: ''.join(sorted(x)))
#A-C--> AC, for grouping
angles_df['ElEl'] = angles_df['El1-El3'].apply(lambda x: ''.join(sorted(x)))
#AC + ABC --> ACABC, for grouping
angles_df['ElEl_ElElEl'] = angles_df['ElEl'] + angles_df['ElElEl']
#angles_df data frame --> sel_angles data frame
sel_angles=angles_df

############ Exclude

#exclude named atoms from input (in Atom1 and Atom2 and Atom3)
if args.excludeAt:
	sel_angles=angles_df[~angles_df.atom1_idx.isin(args.excludeAt) & ~angles_df.atom2_idx.isin(args.excludeAt) & ~angles_df.atom3_idx.isin(args.excludeAt)] 
	
#exclude named elements from input (in Element1 and Element2 and Element3)
if args.excludeEl:
	sel_angles=angles_df[~angles_df.element1.isin(args.excludeEl) & ~angles_df.element2.isin(args.excludeEl) & ~angles_df.element3.isin(args.excludeEl)]
	
#exit if selected angles data frame is empty 
if len(sel_angles) == 0:
	print("No angles found. Include more atoms or elements. Exit.")
	sys.exit(1)

############ Sort

#sort angle values ascending
if args.sortasc:
	sel_angles=sel_angles.sort_values(by=['angle_calc'])
	
#sort angle values descending
if args.sortdes:
	sel_angles=sel_angles.sort_values(by=['angle_calc'],ascending=False)
	
#sort by elements ascending, A --> Z (not PSE like)
if args.sortascEl:
	sel_angles=sel_angles.sort_values(by=['element1','element2','element3','angle_calc','atom1_idx','atom2_idx','atom3_idx'])
	
#sort by elements descending, A --> Z (not PSE like)
if args.sortdesEl:
	sel_angles=sel_angles.sort_values(by=['element1','element2','element3','angle_calc','atom1_idx','atom2_idx','atom3_idx'],ascending=False)

############ Print

#table with all selected angles A-B-C | 123.45°
pr_sel_angles=sel_angles[['A-B-C','angle_calc']]
print('')
print(tabulate(pr_sel_angles,
	  headers=['Atoms','Angle /°'],
	  tablefmt='github',
	  floatfmt=(".2f"),
	  showindex=False))

#lists for printed tables
summary_angle_table_1 = list()
summary_angle_table_2 = list()
summary_angle_table_3 = list()

#group El1-El2-El3 and angles by ElElEl, e.g. O-C-N 1.234 Å by OCOCN sorted CCNOO
grouped = sel_angles[['El1-El2-El3','angle_calc']].groupby(sel_angles['ElEl_ElElEl'],sort=False)

#verbose table El1-El2-El3 | angle, e.g. C-C-C 122.31, 145.61, 151.11
for groups in grouped:
	summary_angle_table_1.append([groups[1].iloc[0].tolist()[0], 
		', '.join(groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist())])

#short table El1-El2-El3 | angle, e.g.C-C-C 122.32 - 151.11 (for >2), 
#or C-C-C 122.32 / 151.11 (for 2), C-C-C 122.32 (for one)
for groups in grouped:
	if len(groups[1]) == 1:
		summary_angle_table_2.append([groups[1].iloc[0].tolist()[0], 
			groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[0]])
	elif len(groups[1]) == 2:
		summary_angle_table_2.append([groups[1].iloc[0].tolist()[0], 
			groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[0] + 
			' / ' + groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[-1]])
	else:
		summary_angle_table_2.append([groups[1].iloc[0].tolist()[0], 
			groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[0] + 
			' - ' + groups[1].sort_values(by=['angle_calc']).angle_calc.apply(lambda x: '{:.2f}'.format(x)).tolist()[-1]])
		
#grouped = sel_angles[['El1-El2-El3','angle_calc']].groupby(sel_angles['ElEl_ElElEl'],sort=False)

#generate table with statistics, | El1-El2-El3 | Count | Mean | Median | Sam. std. dev. | Pop. std. dev. | Std. error |
for groups in grouped:
	summary_angle_table_3.append([groups[1].iloc[0].tolist()[0], groups[1].angle_calc.count(),f'{groups[1].angle_calc.mean():.2f}', \
		f'{groups[1].angle_calc.median():.2f}', f'{groups[1].angle_calc.std():.2f}', f'{groups[1].angle_calc.std(ddof=0):.2f}', \
		f'{groups[1].angle_calc.sem():.2f}',f'{groups[1].angle_calc.skew():.2f}'])

#print verbose table El1-El2-El3 | angle only by request
if args.verbose:
	print('')
	print(tabulate(summary_angle_table_1,
		  headers=['Atoms','Angle /°'], 
		  tablefmt='github',
		  floatfmt=(".2f"),
		  showindex=False))

#print short table
print('')
print(tabulate(summary_angle_table_2,
	  headers=['Atoms','Angle /°'], 
	  tablefmt='github',
	  floatfmt=(".2f"),
	  showindex=False))

#print statistics table
print('')
print(tabulate(summary_angle_table_3,
	  headers=['Atoms','Count','Mean /°', 'Median /°','Sam. std. dev.', \
	 'Pop. std. dev.','Std. error','Skewness'], 
	  tablefmt='github',
	  showindex=False))
	
#print the dihedral angle on request
if args.dihedral:
	a1=xyz_df.index[xyz_df['atom1_idx'] == args.dihedral[0]].tolist()
	a2=xyz_df.index[xyz_df['atom1_idx'] == args.dihedral[1]].tolist()
	a3=xyz_df.index[xyz_df['atom1_idx'] == args.dihedral[2]].tolist()
	a4=xyz_df.index[xyz_df['atom1_idx'] == args.dihedral[3]].tolist()
	try:
		d_angle = calc_d_angle(xyzarr, *a1, *a2, *a3, *a4)
	except TypeError:
		print('')
		print('Warning! Dihedral angle: One or more atoms could not be found in the input file.')
		sys.exit(1)
	print('')
	print('Dihedral angle ' + args.dihedral[0] + '-' +args.dihedral[1] + '-' + \
		    args.dihedral[2] + '-' + args.dihedral[3] + ':', f'{d_angle:.2f}°')

#Plane No. 1 through selected or all atoms on request
if args.plane1:
	#get the index of the selected atoms
	a1 = xyz_df.index[xyz_df['atom1_idx'].isin(args.plane1)].tolist()
	#check for range notation symbol ':'
	if ':' in args.plane1:
		#get indices for ':'
		sep_indices = [i for i, x in enumerate(args.plane1) if x == ":"]
		#avoid error message for strange input
		start = xyz_df.index[0]
		#loop over indices for ':'
		for sep_idx in sep_indices: 
			#start is at atom index 0 if 0
			if sep_idx == 0:
				start = xyz_df.index[0]
			else:
				try:
					#start index is the atom index of the atom in input, left from ':'
					start = xyz_df.index[xyz_df['atom1_idx'] == args.plane1[sep_idx-1]].values.astype(int)[0]
				except IndexError:
					#catch malformed inputs
					print('')
					print('Warning! Malformed input.')
			try: 
				#stop index is the atom index of the atom in input, right from ':'
				stop = xyz_df.index[xyz_df['atom1_idx'] == args.plane1[sep_idx+1]].values.astype(int)[0]
			except IndexError:
				#if there is no atom right from ':' take the last atom
				stop = xyz_df.index[-1]	
			#atom index range, stop+1 because of numpy 
			atm_range = np.arange(start, stop+1)
			#get unique index, avoid duplicates from range input and singel ato input 
			#e.g. 1 2 4 + range 3 4 5 --> 1 2 3 4 4 5 --> 1 2 3 4 5
			a1 = np.unique(list(a1) + list(atm_range))
		
	#get atom names from the selected atoms
	atom_names1 = xyz_df.iloc[a1,1].tolist()
	#get the x,y,z coordinates from the selected atoms
	xyz_pl1_arr = xyzarr[a1]
	
	#check if atom names from input are in the atom list
	#if not, warning
	atom_names_in_arg = [x for x in args.plane1 if x != ':']
	if not all(elem in atom_names1 for elem in atom_names_in_arg):
		print('')
		print('Warning! One or more atoms could not be found in the input file.')
	
	#no atoms / coordinates for calculation --> exit
	if len(xyz_pl1_arr) == 0:
		#if the list is empty --> exit
		print('')
		print('Warning! No atoms for Plane 1. Exit.')
		sys.exit(1)
	
	#calculate the best fit plane
	c1, n1 = svd_fit(np.asarray(xyz_pl1_arr))
	
	#create data frame for output
	plane1_df = pd.DataFrame()
	plane1_df['Atom'] = atom_names1
	#calculate the distance of each atom to the plane 
	plane1_df['Distance'] = np.dot(np.asarray(xyz_pl1_arr)-c1, n1)
	#sum of squares error
	#sosqf_p1  = plane1_df['Distance'].apply(lambda x: abs(x)**2)
	
	print('')
	print('Best-fit Plane 1 through', len(atom_names1), 'atoms.')
	
	#print some plane related parameters
	#print('')
	#print('Centroid: ', *c1, 'Å')
	#print('Plane normal: ', *n1, 'Å')
	#print('Sum-of-squares error:', f'{sosqf_p1.sum():.4f} Å²')
	
	#print the table with atom names and distances to the plane
	print('')
	print(tabulate(plane1_df,
	headers=['Atoms','Distances to Plane 1 /Å'], 
	tablefmt='github',
	floatfmt=(".4f"),
	showindex=False))
	
#Plane No. 2 through selected atoms on request
#check comments for plane 1
if args.plane2:
	#get the index of the selected atoms
	a2 = xyz_df.index[xyz_df['atom1_idx'].isin(args.plane2)].tolist()
	
	if ':' in args.plane2:
		
		sep_indices2 = [i for i, x in enumerate(args.plane2) if x == ":"]
		start2 = xyz_df.index[0]
		for sep_idx2 in sep_indices2: 
			
			if sep_idx2 == 0:
				start2 = xyz_df.index[0]
			else:
				try:
					start2 = xyz_df.index[xyz_df['atom1_idx'] == args.plane2[sep_idx2-1]].values.astype(int)[0]
				except IndexError:
					print('')
					print('Warning! Malformed input.')
			try: 
				stop2 = xyz_df.index[xyz_df['atom1_idx'] == args.plane2[sep_idx2+1]].values.astype(int)[0]
			except IndexError:
				stop2 = xyz_df.index[-1]	
				
			atm_range2 = np.arange(start2, stop2+1)
			a2 = np.unique(list(a2) + list(atm_range2))
			
	#get atom names from the selected atoms
	atom_names2 = xyz_df.iloc[a2,1].tolist()
	#get the x,y,z coordinates from the selected atoms
	xyz_pl2_arr = xyzarr[a2]
	
	atom_names_in_arg2 = [x for x in args.plane2 if x != ':']
	if not all(elem in atom_names2 for elem in atom_names_in_arg2):
		print('')
		print('Warning! One or more atoms could not be found in the input file.')

	if len(xyz_pl2_arr) == 0:
		#if the list is empty --> exit
		print('')
		print('Warning! No atoms for Plane 2. Exit.')
		sys.exit(1)
		
	#calculate the best fit plane
	c2, n2 = svd_fit(np.asarray(xyz_pl2_arr))
	
	#create data frame for output
	plane2_df = pd.DataFrame()
	plane2_df['Atom'] = atom_names2
	#calculate the distance of each atom to plane 2 
	plane2_df['DistanceP2'] = np.dot(np.asarray(xyz_pl2_arr)-c2, n2)
	#calculate the distance of each atom to plane 1 
	plane2_df['DistanceP1'] = np.dot(np.asarray(xyz_pl2_arr)-c1, n1)
	#sum of squares error
	#sosqf_p2 = plane2_df['DistanceP2'].apply(lambda x: abs(x)**2)
	
	print('')
	print('Best-fit Plane 2 through', len(atom_names2), 'atoms.')
	
	#print some plane related parameters
	#print('')
	#print('Centroid: ', *c2, 'Å')
	#print('Plane normal: ', *n2, 'Å')
	#print('Sum-of-squares error:', f'{sosqf_p2.sum():.4f} Å²')
	
	#print the table with atom names and distances to the plane 2 and plane 1 
	print('')
	print(tabulate(plane2_df,
	headers=['Atoms','Distances to Plane 2 /Å','Distances to Plane 1 /Å'], 
	tablefmt='github',
	floatfmt=(".4f"),
	showindex=False))
	
	#calculate the angle between plane 1 and plane 2
	#no warning if senseless result, e.g. plane 1 = plane 2
	np.seterr(invalid='ignore')
	phi = np.arccos(np.dot(n1,n2))
	print('')
	print('Angle between Plane 1 and Plane 2:', f'{np.degrees(phi):.2f}°')
