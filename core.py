#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chimera
from AddCharge import estimateFormalCharge
from SplitMolecule.split import molecule_from_atoms
from Matrix import chimera_xform

try:
	from cStringIO import StringIO
except ImportError:
	from StringIO import StringIO
from textwrap import dedent
from Bld2VRML import openFileObject as openBildFileObject

import os
import numpy as np
import math

from rdkit import Chem, RDConfig, Geometry
from rdkit.Geometry.rdGeometry import Point3D
from rdkit.Chem import SanitizeMol, AllChem, rdmolops, MolFromPDBBlock, FastFindRings
from rdkit.Chem.rdMolAlign import GetAlignmentTransform, GetO3A, AlignMol, GetO3AForProbeConfs
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem.FeatMaps import FeatMapUtils as FMU
from rdkit.Chem.Features import FeatDirUtilsRD as FeatDirUtils
from rdkit.Chem.AllChem import EmbedMolecule, ETKDG, MMFFOptimizeMolecule, MMFFGetMoleculeProperties, EmbedMultipleConfs, AddHs, RemoveHs

FEATURES_FILE = os.path.join(os.path.dirname(__file__), 'BaseFeatures.fdef')

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

_featColors = { 
	'Donor': '0 1 1', 
	'Acceptor': '1 0 1', 
	'NegIonizable': '1 0 0', 
	'PosIonizable': '0 0 1', 
	'ZnBinder': '1 .5 .5', 
	'Aromatic': '1 .8 .2', 
	'LumpedHydrophobe': '.5 .25 0', 
	'Hydrophobe': '.5 .25 0', 
}

#Find if there is a method in chimera to do this (new function)
def _find_by_sn(molecule, serialN):
	for atom in molecule.atoms:
		if(atom.serialNumber == serialN):
			return atom

#New function
def sulfurOxygen(atom):
	if atom.idatmType != "O3-":
		return False
	try:
		s = atom.bondsMap.keys()[0]
	except IndexError:
		return False
	if s.idatmType in ['Son', 'Sxd']:
		return True
	if s.idatmType == 'Sac':
		o3s = [a for a in s.neighbors if a.idatmType == 'O3-']
		o3s.sort()
		return o3s.index(atom) > 1
	return False


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
					if fType1 == fType2:
						sumVec1 = 0
						for i, pt in enumerate(ps1):
							if (sumVec1 == 0):
								sumVec1 = pt[1] - pt[0]
							else:
								sumVec1 += pt[1] - pt[0]
						sumVec1 = (sumVec1/(i+1))
						sumVec2 = 0
						for i, pt in enumerate(ps2):
							if (sumVec2 == 0):
								sumVec2 = pt[1] - pt[0]
							else:
								sumVec2 += pt[1] - pt[0]
						sumVec2 = (sumVec2/(i+1))
						
						sumVec = (sumVec1+sumVec2)/2
						sumVec = Point3D(sumVec[0], sumVec[1], sumVec[2])
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
		return ((cpt, bvec), ), 'linear'
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
			bvec = conf.GetAtomPosition(hid) 
			bvec -= cpt 
			bvec.Normalize() 
			bvec *= scale 
			bvec += cpt 
			res.append((cpt, bvec)) 
		return tuple(res), 'linear' 

class p4_element(object):

	SUPPORTED_SHAPES = set('sphere arrow cone'.split())

	def __init__(self, shape, origin, color,
				 name='p4', parent_id=100, end=None, size=None):
		if shape not in self.SUPPORTED_SHAPES:
			raise ValueError('`shape` should be one of: '
							 '{}'.format(', '.join(self.SUPPORTED_SHAPES)))
		self.shape = shape
		self.name = name
		self.size = size
		self.origin = origin
		self.end = end
		self.color = color
		self._vrml_shape = None
		self._id = parent_id
		self._subid = 0

	def destroy(self):
		if self._vrml_shape is not None:
			chimera.openModels.close(self._vrml_shape)
			self._vrml_shape = None

	def draw(self):
		self._vrml_shape = getattr(self, '_draw_' + self.shape)()

	def _draw_sphere(self):
		x, y, z = self.origin
		bild = """
		.color {}
		.transparency {}
		.sphere {} {} {} {}
		""".format(self.color, 0.3, x, y, z, self.size)
		return self._build_vrml(bild)

	def _draw_arrow(self):
		x1, y1, z1 = self.origin
		x2, y2, z2 = self.end
		bild = """
		.color {}
		.arrow {} {} {} {} {} {} {} {} {}
		""".format(self.color, x1, y1, z1, x2, y2, z2, 0.05, 0.1, 0.8)
		return self._build_vrml(bild)

	def _draw_cone(self):
		x1, y1, z1 = self.origin
		x2, y2, z2 = self.end
		bild = """
		.color {}
		.transparency {}
		.cone {} {} {} {} {} {} {} {}
		""".format(self.color, 0.5, x1, y1, z1, x2, y2, z2, self.size, 'open')
		return self._build_vrml(bild)

	def _build_vrml(self, bild, name=None):
		if name is None:
			name = self.name
		f = StringIO(dedent(bild))
		try:
			vrml = openBildFileObject(f, '<string>', name)
		except chimera.NotABug:
			print(bild)
		else:
			chimera.openModels.add(vrml, baseId=self._id, subid=self._subid)
			self._subid += 1
		return vrml

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
def _chimera_to_rdkit(molecule, sanitize=True):
	from rdkit.Chem import Mol, RWMol, Atom, Conformer
	mol = Mol()
	emol = RWMol(mol)
	emol.AddConformer(Conformer())
	atom_map = {}
	
	for atom in molecule.atoms:
		a = Atom(atom.element.number)
		atom_map[atom] = i = emol.AddAtom(a)
		emol.GetConformer().SetAtomPosition(i, atom.coord().data())
	for bond in molecule.bonds:
		a1, a2 = bond.atoms
		if hasattr(bond, 'order'):
			bond_order =  RDKit_BondType[float(bond.order)]
		else:
			bond_order = Chem.BondType.SINGLE		
		emol.AddBond(atom_map[a1], atom_map[a2], bond_order)
		
	mol = Mol(emol)
	_fix_nitro(mol)

	if sanitize:
		SanitizeMol(mol)

	return mol, atom_map

def calc_p4map(molecules, families=('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic', 'LumpedHydrophobe'), mergeMetric=1, mergeTol=2.5, dirMergeMode=1, minRepeats=1, showVectors=True):
	rdkit_mols = []
	rdkit_maps = []
	for mol in molecules:
		rdkit_mol, rdkit_map = _chimera_to_rdkit(mol)
		rdmolops.Cleanup(rdkit_mol)
		rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
		rdkit_mols.append(rdkit_mol)
		rdkit_maps.append(rdkit_map)

	fdef = AllChem.BuildFeatureFactory(FEATURES_FILE)
	fmParams = {}
	for k in fdef.GetFeatureFamilies():
		fparams = FeatMaps.FeatMapParams()
		fmParams[k] = fparams

	keep = families
	global_fmap = FeatMaps.FeatMap(params=fmParams)
	for m in rdkit_mols:
		rawFeats=[]
		for f in fdef.GetFeaturesForMol(m):
			if showVectors:
				if f.GetFamily() == 'Acceptor':
					aids = f.GetAtomIds() 
					"""
					if len(aids) == 1: 
						featAtom = m.GetAtomWithIdx(aids[0]) 
						hvyNbrs = [x for x in featAtom.GetNeighbors() if x.GetAtomicNum() != 1] 
						if len(hvyNbrs) == 1: 
							f.featDirs = _GetAcceptor1FeatVects(m.GetConformer(-1), aids, scale=(mergeTol)) 
						elif len(hvyNbrs) == 2: 
							f.featDirs = FeatDirUtils.GetAcceptor2FeatVects(m.GetConformer(-1), aids, scale=(mergeTol)) 
						elif len(hvyNbrs) == 3: 
							f.featDirs = FeatDirUtils.GetAcceptor3FeatVects(m.GetConformer(-1), aids, scale=(mergeTol))
					"""	 
				elif f.GetFamily() == 'Donor':
					aids = f.GetAtomIds() 
					if len(aids) == 1: 
						featAtom = m.GetAtomWithIdx(aids[0]) 
						hvyNbrs = [x for x in featAtom.GetNeighbors() if x.GetAtomicNum() != 1] 
						if len(hvyNbrs) == 1: 
							f.featDirs = FeatDirUtils.GetDonor1FeatVects(m.GetConformer(-1), aids, scale=(0.5)) 
						elif len(hvyNbrs) == 2: 
							f.featDirs = _GetDonor2FeatVects(m.GetConformer(-1), aids, scale=(mergeTol)) 
						elif len(hvyNbrs) == 3: 
							f.featDirs = FeatDirUtils.GetDonor3FeatVects(m.GetConformer(-1), aids, scale=(mergeTol))
				elif f.GetFamily() == 'Aromatic':
					f.featDirs = FeatDirUtils.GetAromaticFeatVects(m.GetConformer(-1), f.GetAtomIds(), f.GetPos(-1), scale=(mergeTol))				
			rawFeats.append(f)
		# filter that list down to only include the ones we're intereted in 
		featList = [f for f in rawFeats if f.GetFamily() in keep]

		fmap = FeatMaps.FeatMap(feats = featList,weights=[1]*len(featList),params=fmParams)
		Merge = True
		while Merge == True:
			Merge = _MergeFeatPoints(fmap, mergeMetric=mergeMetric, mergeTol=mergeTol, dirMergeMode=dirMergeMode)
		global_fmap = FMU.CombineFeatMaps(global_fmap, fmap, mergeMetric=0)	
	
	matrix = FMU.GetFeatFeatDistMatrix(global_fmap, mergeMetric=mergeMetric, mergeTol=mergeTol, dirMergeMode=dirMergeMode, compatFunc=FMU.familiesMatch)
	p4map = FeatMaps.FeatMap(params=fmParams)
	for i, vector in enumerate(matrix):
		feat_indexs = [vector.index(x) for x in vector if x<=mergeTol]
		feat_indexs.append(i)
		if (len(feat_indexs) >= minRepeats):
			for feat_index in feat_indexs:
				p4map.AddFeature(global_fmap._feats[feat_index], weight=1)
	Merge = True
	while Merge == True:
		Merge = _MergeFeatPoints(p4map, mergeMetric=mergeMetric, mergeTol=mergeTol, dirMergeMode=dirMergeMode)
	return p4map

def chimera_p4(molecules_sel, mergeTol=2.5, minRepeats=1, showVectors=True):
	molecules = molecules_sel.molecules()

	p4map = calc_p4map(molecules, mergeTol=mergeTol, minRepeats=minRepeats, showVectors=showVectors)

	for feat in p4map._feats:
		if feat.GetFamily() != 'Donor':
			p4_elem = p4_element(shape="sphere", size=(mergeTol/2), origin=chimera.Point(feat.GetPos()[0],feat.GetPos()[1],feat.GetPos()[2]), color=_featColors[feat.GetFamily()])
			p4_elem.draw()
		
		if feat.featDirs:
			ps, fType = feat.featDirs			
			for tail, head in ps:
				if fType == 'linear':
					p4_elem = p4_element(shape="arrow", origin=tail, end=head, color=_featColors[feat.GetFamily()])
				elif fType =='cone':
					p4_elem = p4_element(shape="cone", origin=head, end=tail, color=_featColors[feat.GetFamily()], size=(1.33))
				p4_elem.draw()
		
	msg = "Chimera pharmacophore done"
	chimera.statusline.show_message(msg)
	
	return True

#### Open3Align code
def _apply_atom_positions(rdkit_mol, chimera_mol, atom_map, rdkit_confId=0):
	from chimera import Coord, Point

	for rdkit_atom in rdkit_mol.GetAtoms():
		chimera_atom = list(atom_map.keys())[list(atom_map.values()).index(rdkit_atom.GetIdx())]
		#chimera_atom = _find_by_sn(chimera_mol, rdkit_atom.GetIdx())

		atom_position = Point()
		atom_position[0] = list(rdkit_mol.GetConformer(id=rdkit_confId).GetAtomPosition(rdkit_atom.GetIdx()))[0]
		atom_position[1] = list(rdkit_mol.GetConformer(id=rdkit_confId).GetAtomPosition(rdkit_atom.GetIdx()))[1]
		atom_position[2] = list(rdkit_mol.GetConformer(id=rdkit_confId).GetAtomPosition(rdkit_atom.GetIdx()))[2]

		chimera_atom.setCoord(Coord(atom_position))
		
	chimera_mol.computeIdatmTypes()

#Function from plume
def _transform_molecule(molecule, xform):
	new_xform = molecule.openState.xform
	new_xform.premultiply(xform)
	molecule.openState.xform = new_xform


#Function from plume with small changes in sanitization
def align_o3a(reference, probe, transform=True, sanitize=True, nConformers=0, **kwargs):
	rdk_reference, ref_map = _chimera_to_rdkit(reference)
	rdk_probe, probe_map = _chimera_to_rdkit(probe)
	FastFindRings(rdk_reference)
	FastFindRings(rdk_probe)
	reference_params = MMFFGetMoleculeProperties(rdk_reference)
	probe_params = MMFFGetMoleculeProperties(rdk_probe)
	
	if nConformers > 0:
		rdk_probe = AddHs(rdk_probe, addCoords=True)
		cids = EmbedMultipleConfs(rdk_probe, nConformers, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
		rdk_probe = RemoveHs(rdk_probe)
		rdk_reference = RemoveHs(rdk_reference)

		o3as = GetO3AForProbeConfs(rdk_probe, rdk_reference, numThreads=0)
	
		highest_conf_score = 0.0
		conf_id = 0
		for o3a in o3as:
			score_conf = o3a.Score()
			if score_conf > highest_conf_score:
				highest_conf_score = score_conf
				o3a_result = o3a
				highest_conf_id = conf_id
			conf_id += 1
	else:
		rdk_probe = RemoveHs(rdk_probe)
		rdk_reference = RemoveHs(rdk_reference)
		o3a_result = GetO3A(rdk_probe, rdk_reference, probe_params, reference_params)
		highest_conf_id = 0
	
	#rmsd, xform = o3a_result.Trans()
	if transform:
		o3a_result.Align()
		_apply_atom_positions(rdk_probe, probe, probe_map, rdkit_confId=highest_conf_id)
		#_transform_molecule(probe, chimera_xform(xform[:3]))
	
	return o3a_result   #return the alignment

#New function
def open3align(molecules_sel, transform=True, nConformers=0):
	molecules = molecules_sel.molecules()

	if not len(molecules) > 1:
		raise chimera.UserError("At least 2 molecules are needed to do an alignment")
	
	#Calculating the best scored alignment
	max_score = 0.0
	for reference in molecules:

		align_score = 0.0
		for probe in molecules:
			if(molecules.index(probe) is not molecules.index(reference)):
				o3a = align_o3a(reference, probe, transform=False, nConformers=nConformers)
				align_score += o3a.Score()
		
		if align_score > max_score:
			max_score = align_score
			if transform:
				for probe in molecules:
					if(molecules.index(probe) is not molecules.index(reference)):
						o3a = align_o3a(reference, probe, transform=True, nConformers=nConformers)

	msg = "The score of the best alignment is {}".format(max_score)
	chimera.statusline.show_message(msg)
	
	return max_score

### GUI Code
class Controller(object):

	"""
	The controller manages the communication between the UI (graphic interface)
	and the data model. Actions such as clicks on buttons, enabling certain areas, 
	or running external programs, are the responsibility of the controller.
	"""
	def __init__(self, *args, **kwargs):
		return


class Model(object):

	"""
	The model controls the data we work with. Normally, it'd be a Chimera molecule
	and some input files from other programs. The role of the model is to create
	a layer around those to allow the easy access and use to the data contained in
	those files
	"""

	def __init__(self, *args, **kwargs):
		return