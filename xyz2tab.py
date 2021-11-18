#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#xyz2tab

import sys                                                  #stdout
import argparse                                             #argument parser
import re                                                   #regex
import itertools	                                    #for r-length tuples, in sorted order, no repeated elements
import pandas as pd                                         #pandas tables
import numpy as np                                          #for calculations
from scipy.spatial.distance import pdist, squareform	    #for the calculations of the distance matrix
from tabulate import tabulate                               #nice table output
import gemmi	                                            #covalence radii, molecular weight; can be replaced by own or other
                                                            #covalent radii and mol. weight lists

#for windows console
sys.stdout.reconfigure(encoding='utf-8')  

#pd.set_option("display.max_rows", None, "display.max_columns", None)

#+x % to covalence radius
radii_ext = 11.15 / 100 

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

#calculate angle from 3 vectors / atomic coordinates: j(x,y,z); k(x,y,z); l(x,y,z) 
#xyzarr is the array of all atomic coordinates
def calc_angle(xyzarr, i, j, k):
	rij = xyzarr[i] - xyzarr[j]
	rkj = xyzarr[k] - xyzarr[j]
	cos_theta = np.dot(rij, rkj)
	sin_theta = np.linalg.norm(np.cross(rij, rkj))
	theta = np.arctan2(sin_theta, cos_theta)
	return 	np.degrees(theta)

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
	help='include contacts, e.g. -ee C0 C11 or -ee C0 C11 N14')

#include contacts
parser.add_argument('-r','--radius',
	default=11.15,
	type=float,
	help='enlarge atomic radii by x %%, e.g. -r 15.2, default is 11.15 %%')

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
#atomic weight
xyz_df['weight'] = xyz_df['element'].apply(lambda x: gemmi.Element(x).weight)
#covalent radius
xyz_df['cov_radius'] = xyz_df['element'].apply(lambda x: gemmi.Element(x).covalent_r)
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
info_df['atom%']=info_df.apply(lambda x: gemmi.Element(x.element).weight*x.atom1_idx/fw*100, axis=1)
#atomic radii from the gemmi table
info_df['radii']=info_df.apply(lambda x: gemmi.Element(x.element).covalent_r, axis=1)
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
#angle calculation is not as fast as bond length calculation
xyzarr = xyz_df.iloc[:,3:6].to_numpy()

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
