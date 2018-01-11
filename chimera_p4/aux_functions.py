#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import chimera
from chimera import Vector

import numpy as np

from rdkit import Chem, Geometry
from rdkit.Chem import SanitizeMol
from rdkit.Chem.FeatMaps import FeatMapUtils as FMU
from rdkit.Chem.Features import FeatDirUtilsRD as FeatDirUtils

# The RDKit_BondType dictionary is defined to convert float bond orders into RDKit bond types
RDKit_BondType = {
	1.0 : Chem.BondType.SINGLE,
	2.0 : Chem.BondType.DOUBLE,
	3.0 : Chem.BondType.TRIPLE,
	4.0 : Chem.BondType.QUADRUPLE,
	5.0 : Chem.BondType.QUINTUPLE,
	6.0 : Chem.BondType.HEXTUPLE,
	1.5 : Chem.BondType.ONEANDAHALF,
	2.5 : Chem.BondType.TWOANDAHALF,
	3.5 : Chem.BondType.THREEANDAHALF,
	4.5 : Chem.BondType.FOURANDAHALF,
	5.5 : Chem.BondType.FIVEANDAHALF
}

def feq(v1, v2, tol=1e-4): 
	return abs(v1 - v2) < tol 

def _ArbAxisRotation(theta, ax, pt): 
	theta = math.pi * theta / 180 
	c = math.cos(theta) 
	s = math.sin(theta) 
	t = 1 - c 
	X = ax.x 
	Y = ax.y 
	Z = ax.z 
	mat = [[t * X * X + c, t * X * Y + s * Z, t * X * Z - s * Y], 
		   [t * X * Y - s * Z, t * Y * Y + c, t * Y * Z + s * X], 
		   [t * X * Z + s * Y, t * Y * Z - s * X, t * Z * Z + c]] 
	mat = np.array(mat) 
	if isinstance(pt, Geometry.Point3D): 
		pt = np.array((pt.x, pt.y, pt.z)) 
		tmp = np.dot(mat, pt) 
		res = Geometry.Point3D(tmp[0], tmp[1], tmp[2]) 
	elif isinstance(pt, list) or isinstance(pt, tuple): 
		pts = pt 
		res = [] 
		for pt in pts: 
			pt = np.array((pt.x, pt.y, pt.z)) 
			tmp = np.dot(mat, pt) 
			res.append(Geometry.Point3D(tmp[0], tmp[1], tmp[2])) 
	else: 
		res = None 
	return res 

def _findAvgVec(conf, center, nbrs): 
	# find the average of the normalized vectors going from the center atoms to the 
	# neighbors 
	# the average vector is also normalized 
	avgVec = 0 
	for nbr in nbrs: 
		nid = nbr.GetIdx() 
		pt = conf.GetAtomPosition(nid) 
		pt -= center 
		pt.Normalize() 
		if (avgVec == 0): 
			avgVec = pt 
		else: 
			avgVec += pt 
   
	avgVec.Normalize() 
	return avgVec 

def _MergeFeatPoints(fm, mergeMetric=FMU.MergeMetric.NoMerge, mergeTol=1.5, 
					  dirMergeMode=FMU.DirMergeMode.NoMerge, mergeMethod=FMU.MergeMethod.WeightedAverage, 
					  compatFunc=FMU.familiesMatch): 
	""" 
   
	  NOTE that mergeTol is a max value for merging when using distance-based 
	  merging and a min value when using score-based merging. 
   
	  returns whether or not any points were actually merged 
   
	""" 
	FMU.MergeMetric.valid(mergeMetric) 
	FMU.MergeMethod.valid(mergeMethod) 
	FMU.DirMergeMode.valid(dirMergeMode) 
   
	res = False 
	if mergeMetric == FMU.MergeMetric.NoMerge: 
		return res 
	dists = FMU.GetFeatFeatDistMatrix(fm, mergeMetric, mergeTol, dirMergeMode, compatFunc) 
	distOrders = [None] * len(dists) 
	for i in range(len(dists)): 
		distV = dists[i] 
		distOrders[i] = [] 
		for j, dist in enumerate(distV): 
			if dist < mergeTol: 
				distOrders[i].append((dist, j)) 
		distOrders[i].sort() 
  
	# we now know the "distances" and have rank-ordered list of 
	# each point's neighbors. Work with that. 
   
	# progressively merge nearest neighbors until there 
	# are no more points left to merge 
	featsInPlay = list(range(fm.GetNumFeatures())) 
	featsToRemove = [] 
	# print '--------------------------------' 
	while featsInPlay: 
		# find two features who are mutual nearest neighbors: 
		fipCopy = featsInPlay[:] 
		for fi in fipCopy: 
	
			mergeThem = False 
			if not distOrders[fi]: 
				featsInPlay.remove(fi) 
				continue 
			dist, nbr = distOrders[fi][0] 
			if nbr not in featsInPlay: 
				continue 
			if distOrders[nbr][0][1] == fi: 
				# print 'direct:',fi,nbr 
				mergeThem = True 
			else: 
				# it may be that there are several points at about the same distance, 
				# check for that now 
				if (feq(distOrders[nbr][0][0], dist)): 
					for distJ, nbrJ in distOrders[nbr][1:]: 
						if feq(dist, distJ): 
							if nbrJ == fi: 
								# print 'indirect: ',fi,nbr 
								mergeThem = True 
								break 
						else: 
							break 
			# print '    bottom:',mergeThem 
			if mergeThem: 
				break 
		if mergeThem: 
			res = True 
			featI = fm.GetFeature(fi) 
			nbrFeat = fm.GetFeature(nbr) 
   
			if mergeMethod == FMU.MergeMethod.WeightedAverage: 
				newPos = featI.GetPos() * featI.weight + nbrFeat.GetPos() * nbrFeat.weight 
				newPos /= (featI.weight + nbrFeat.weight) 
				newWeight = (featI.weight + nbrFeat.weight) / 2 
			elif mergeMethod == FMU.MergeMethod.Average: 
				newPos = featI.GetPos() + nbrFeat.GetPos() 
				newPos /= 2 
				newWeight = (featI.weight + nbrFeat.weight) / 2 
			elif mergeMethod == FMU.MergeMethod.UseLarger: 
				if featI.weight > nbrFeat.weight: 
					newPos = featI.GetPos() 
					newWeight = featI.weight 
				else: 
					newPos = nbrFeat.GetPos() 
					newWeight = nbrFeat.weight 
   
			featI.SetPos(newPos) 
			featI.weight = newWeight 
			
			if dirMergeMode == FMU.DirMergeMode.Sum:
				if featI.featDirs and nbrFeat.featDirs:
					ps1, fType1 = featI.featDirs
					ps2, fType2 = nbrFeat.featDirs
				else:
					ps1 = ps2 = None
				if ps1 and ps2:
					with open('test_file.txt', 'a') as f:
						print('ps1', file=f)
						print(ps1, file=f)
						print('ps2', file=f)
						print(ps2, file=f)
						if fType1 == fType2:
							sumVec1 = 0
							print('ps1 enumerate', file=f)
							for i, pt in enumerate(ps1):
								print(pt, file=f)
								if (sumVec1 == 0):
									sumVec1 = pt[1] - pt[0]
								else:
									sumVec1 += pt[1] - pt[0]
							sumVec1 = (sumVec1/(i+1))
							sumVec2 = 0
							print('ps2 enumerate', file=f)
							for i, pt in enumerate(ps2):
								print(pt, file=f)
								if (sumVec2 == 0):
									sumVec2 = pt[1] - pt[0]
								else:
									sumVec2 += pt[1] - pt[0]
							sumVec2 = (sumVec2/(i+1))
							
							sumVec = (sumVec1+sumVec2)/2
							sumVec = Geometry.Point3D(sumVec[0], sumVec[1], sumVec[2])
							sumVec.Normalize() 
							if fType1 == 'linear':
								sumVec *= 1.5
							elif fType1 == 'cone':
								sumVec *= 0.5
							sumVec += newPos
							featI.featDirs = ((newPos, sumVec), ), fType1
						else:
							del featI.featDirs
				else:
					try:
						del featI.featDirs
					except:
						pass
			
			# nbr and fi are no longer valid targets: 
			# print 'nbr done:',nbr,featsToRemove,featsInPlay 
			featsToRemove.append(nbr) 
			featsInPlay.remove(fi) 
			featsInPlay.remove(nbr) 
			for nbrList in distOrders: 
				try: 
					nbrList.remove(fi) 
				except ValueError: 
					pass 
				try: 
					nbrList.remove(nbr) 
				except ValueError: 
					pass 
		else: 
			# print ">>>> Nothing found, abort" 
			break 
	featsToRemove.sort() 
	for i, fIdx in enumerate(featsToRemove): 
		fm.DropFeature(fIdx - i) 
	return res 

def _GetAcceptor1FeatVects(conf, featAtoms, scale=1.5): 
	""" 
	Get the direction vectors for Acceptor of type 1 
   
	This is a acceptor with one heavy atom neighbor. There are two possibilities we will 
	consider here 
	1. The bond to the heavy atom is a single bond e.g. CO 
	   In this case we don't know the exact direction and we just use the inversion of this bond direction 
	   and mark this direction as a 'cone' 
	2. The bond to the heavy atom is a double bond e.g. C=O 
	   In this case the we have two possible direction except in some special cases e.g. SO2 
	   where again we will use bond direction 
		
	ARGUMENTS: 
	  featAtoms - list of atoms that are part of the feature 
	  scale - length of the direction vector 
	""" 
	assert len(featAtoms) == 1 
	aid = featAtoms[0] 
	mol = conf.GetOwningMol() 
	nbrs = mol.GetAtomWithIdx(aid).GetNeighbors() 
   
	cpt = conf.GetAtomPosition(aid) 
   
	# find the adjacent heavy atom 
	heavyAt = -1 
	for nbr in nbrs: 
		if nbr.GetAtomicNum() != 1: 
			heavyAt = nbr 
			break 
   
	singleBnd = mol.GetBondBetweenAtoms(aid, heavyAt.GetIdx()).GetBondType() == Chem.BondType.SINGLE 
   
	# special scale - if the heavy atom is a sulfur (we should proabably check phosphorous as well) 
	sulfur = heavyAt.GetAtomicNum() == 16 
   
	if singleBnd or sulfur: 
		v1 = conf.GetAtomPosition(heavyAt.GetIdx()) 
		v1 -= cpt 
		v1.Normalize() 
		v1 *= (-1.0 * scale) 
		v1 += cpt 
		return ((cpt, v1), ), 'cone' 
	else: 
		# ok in this case we will assume that 
		# heavy atom is sp2 hybridized and the direction vectors (two of them) 
		# are in the same plane, we will find this plane by looking for one 
		# of the neighbors of the heavy atom 
		hvNbrs = heavyAt.GetNeighbors() 
		hvNbr = -1 
		for nbr in hvNbrs: 
			if nbr.GetIdx() != aid: 
				hvNbr = nbr 
				break 
   
		pt1 = conf.GetAtomPosition(hvNbr.GetIdx()) 
		v1 = conf.GetAtomPosition(heavyAt.GetIdx()) 
		pt1 -= v1 
		v1 -= cpt 
		rotAxis = v1.CrossProduct(pt1) 
		rotAxis.Normalize() 
		bv1 = _ArbAxisRotation(120, rotAxis, v1) 
		bv1.Normalize() 
		bv1 *= scale 
		bv1 += cpt 
		bv2 = _ArbAxisRotation(-120, rotAxis, v1) 
		bv2.Normalize() 
		bv2 *= scale 
		bv2 += cpt 
		return ((cpt, bv1), (cpt, bv2), ), 'linear' 


def _GetDonor2FeatVects(conf, featAtoms, scale=1.5): 
	""" 
	Get the direction vectors for Donor of type 2 
   
	This is a donor with two heavy atoms as neighbors. The atom may are may not have 
	hydrogen on it. Here are the situations with the neighbors that will be considered here 
	  1. two heavy atoms and two hydrogens: we will assume a sp3 arrangement here 
	  2. two heavy atoms and one hydrogen: this can either be sp2 or sp3 
	  3. two heavy atoms and no hydrogens 
	   
	ARGUMENTS: 
	  featAtoms - list of atoms that are part of the feature 
	  scale - length of the direction vector 
	""" 
	assert len(featAtoms) == 1 
	aid = featAtoms[0] 
	mol = conf.GetOwningMol() 
   
	cpt = conf.GetAtomPosition(aid) 
   
	# find the two atoms that are neighbors of this atoms 
	nbrs = list(mol.GetAtomWithIdx(aid).GetNeighbors()) 
	assert len(nbrs) >= 2 
   
	hydrogens = []
	heavy=[] 
	for nbr in nbrs: 
		if nbr.GetAtomicNum() == 1:
			hydrogens.append(nbr) 
		else:
			heavy.append(nbr)
   
	if len(nbrs) == 2: 
		# there should be no hydrogens in this case 
		assert len(hydrogens) == 0 
		# in this case the direction is the opposite of the average vector of the two neighbors 
		bvec = _findAvgVec(conf, cpt, heavy) 
		bvec *= (-1.0 * scale) 
		bvec += cpt 
		return ((cpt, bvec), ), 'linear' 
	elif len(nbrs) == 3: 
		assert len(hydrogens) == 1 
		# this is a little more tricky we have to check if the hydrogen is in the plane of the 
		# two heavy atoms (i.e. sp2 arrangement) or out of plane (sp3 arrangement) 
   
		# one of the directions will be from hydrogen atom to the heavy atom 
		hid = hydrogens[0].GetIdx() 
		bvec = conf.GetAtomPosition(hid) 
		bvec -= cpt 
		bvec.Normalize() 
		bvec *= scale 
		bvec += cpt 
		#return ((cpt, bvec), ), 'linear'
		if FeatDirUtils._checkPlanarity(conf, cpt, nbrs): 
			# only the hydrogen atom direction needs to be used 
			return ((cpt, bvec), ), 'linear' 
		else: 
			# we have a non-planar configuration - we will assume sp3 and compute a second direction vector 
			ovec = _findAvgVec(conf, cpt, heavy) 
			ovec *= (-1.0 * scale) 
			ovec += cpt 
			return ((cpt, bvec), (cpt, ovec), ), 'linear' 
	elif len(nbrs) >= 4: 
		# in this case we should have two or more hydrogens we will simple use there directions 
		res = [] 
		for hid in hydrogens:
			hid = hid.GetIdx() 
			bvec = conf.GetAtomPosition(hid) 
			bvec -= cpt 
			bvec.Normalize() 
			bvec *= scale 
			bvec += cpt 
			res.append((cpt, bvec)) 
		return tuple(res), 'linear' 

#Function to fix Nitro groups from chimera molecules without explicit formal charges
def _fix_nitro(molecule):
	for atom in molecule.GetAtoms():
		valence = 0.0
		for bond in atom.GetBonds():
			valence += bond.GetBondTypeAsDouble()

		if atom.GetAtomicNum() == 7:
			if atom.GetFormalCharge():
				continue
			aromHolder = atom.GetIsAromatic()
			atom.SetIsAromatic(0)
			if valence == 4:
				for neighbor in atom.GetNeighbors():
					print (neighbor.GetAtomicNum())
					print(molecule.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType())
					if (neighbor.GetAtomicNum() == 8) and (molecule.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE):
						print('entra')
						atom.SetFormalCharge(1)
						neighbor.SetFormalCharge(-1)
						break
				atom.SetIsAromatic(aromHolder)

#Function to convert from Chimera to RDKit
def _chimera_to_rdkit(molecule, sanitize=True, useXformCoord=True):
	from rdkit.Chem import Mol, RWMol, Atom, Conformer
	mol = Mol()
	emol = RWMol(mol)
	emol.AddConformer(Conformer())
	atom_map = {}
	
	for atom in molecule.atoms:
		a = Atom(atom.element.number)
		atom_map[atom] = i = emol.AddAtom(a)
		if useXformCoord:
			emol.GetConformer().SetAtomPosition(i, atom.xformCoord())
		else:
			emol.GetConformer().SetAtomPosition(i, atom.coord())
	for bond in molecule.bonds:
		a1, a2 = bond.atoms
		if hasattr(bond, 'order') and bond.order:
			bond_order = RDKit_BondType[float(bond.order)]
		else:
			bond_order = Chem.BondType.SINGLE		
		emol.AddBond(atom_map[a1], atom_map[a2], bond_order)
		
	mol = Mol(emol)
	_fix_nitro(mol)

	if sanitize:
		SanitizeMol(mol)

	return mol, atom_map

#### Open3Align code
def _return_atom_positions(rdkit_mol, chimera_mol, atom_map, rdkit_confId=0):
	from chimera import Coord, Point

	atom_positions = {}
	for rdkit_atom in rdkit_mol.GetAtoms():
		chimera_atom = list(atom_map.keys())[list(atom_map.values()).index(rdkit_atom.GetIdx())]

		atom_position = Point()
		atom_position[0] = list(rdkit_mol.GetConformer(id=rdkit_confId).GetAtomPosition(rdkit_atom.GetIdx()))[0]
		atom_position[1] = list(rdkit_mol.GetConformer(id=rdkit_confId).GetAtomPosition(rdkit_atom.GetIdx()))[1]
		atom_position[2] = list(rdkit_mol.GetConformer(id=rdkit_confId).GetAtomPosition(rdkit_atom.GetIdx()))[2]

		atom_positions[chimera_atom] = Coord(atom_position)

	return atom_positions

def _apply_atom_positions(chimera_mol, new_pos_dict):
	for atom in new_pos_dict.keys():
		atom.setCoord(new_pos_dict[atom])
	chimera_mol.computeIdatmTypes()
	#Undo possible previous rotations/translations of the molecule
	chimera_mol.openState.xform = chimera.Xform.identity()

def _del_chimeraHs(molecule):
	for atom in molecule.atoms:
		if atom.element == chimera.Element("H"):
			molecule.deleteAtom(atom)
