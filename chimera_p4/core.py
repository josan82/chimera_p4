#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chimera
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
from rdkit.Chem import SanitizeMol
from rdkit.Chem import AllChem
from rdkit.Chem.FeatMaps import FeatMaps
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
def feq(v1, v2, tol=1e-4): 
	return abs(v1 - v2) < tol 

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
			"""
			if dirMergeMode == FMU.DirMergeMode.Sum:
				if hasattr(featI, 'featDirs') and hasattr(nbrfeat, 'featDirs'):
					sumDirs = featI.featDirs + nbrfeat.featDirs
					ps, fType = sumDirs
					for tail, head in ps:
						tails
						featI.featDirs = numpy
			"""
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
	tmp = [] 
	while len(nbrs): 
		nbr = nbrs.pop() 
		if nbr.GetAtomicNum() == 1: 
			hydrogens.append(nbr) 
	#	else: 
	#		tmp.append(nbr) 
	#nbrs = tmp 
   
	if len(nbrs) == 2: 
		# there should be no hydrogens in this case 
		assert len(hydrogens) == 0 
		# in this case the direction is the opposite of the average vector of the two neighbors 
		bvec = _findAvgVec(conf, cpt, nbrs) 
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
		if FeatDirUtils._checkPlanarity(conf, cpt, nbrs): 
			# only the hydrogen atom direction needs to be used 
			return ((cpt, bvec), ), 'linear' 
		else: 
			# we have a non-planar configuration - we will assume sp3 and compute a second direction vector 
			ovec = _findAvgVec(conf, cpt, nbrs) 
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

	SUPPORTED_SHAPES = set('sphere arrow'.split())

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

#Function from plume with sanitization and bond order=1
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
	
	if sanitize:
		SanitizeMol(mol)

	return mol, atom_map

def calc_p4map(molecules, families=('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic', 'LumpedHydrophobe'), mergeMetric=1, mergeTol=2.5, dirMergeMode=1, minRepeats=1):
	rdkit_mols = []
	rdkit_maps = []
	for mol in molecules:
		rdkit_mol, rdkit_map = _chimera_to_rdkit(mol)
		rdkit_mol = Chem.AddHs(rdkit_mol)
		rdkit_mols.append(rdkit_mol)
		rdkit_maps.append(rdkit_map)

	fdef = AllChem.BuildFeatureFactory('/home/jose/Pharmacophore/Notebooks/BaseFeatures.fdef')
	fmParams = {}
	for k in fdef.GetFeatureFamilies():
		fparams = FeatMaps.FeatMapParams()
		fmParams[k] = fparams

	keep = families
	global_fmap = FeatMaps.FeatMap(params=fmParams)
	for m in rdkit_mols:
		rawFeats=[]
		for f in fdef.GetFeaturesForMol(m):
			
			if f.GetFamily() == 'Acceptor':
				aids = f.GetAtomIds() 
				if len(aids) == 1: 
					featAtom = m.GetAtomWithIdx(aids[0]) 
					hvyNbrs = [x for x in featAtom.GetNeighbors() if x.GetAtomicNum() != 1] 
					if len(hvyNbrs) == 1: 
						f.featDirs = FeatDirUtils.GetAcceptor1FeatVects(m.GetConformer(-1), aids, scale=(0.6*mergeTol)) 
					elif len(hvyNbrs) == 2: 
						f.featDirs = FeatDirUtils.GetAcceptor2FeatVects(m.GetConformer(-1), aids, scale=(0.6*mergeTol)) 
					elif len(hvyNbrs) == 3: 
						f.featDirs = FeatDirUtils.GetAcceptor3FeatVects(m.GetConformer(-1), aids, scale=(0.6*mergeTol)) 
			elif f.GetFamily() == 'Donor':
				aids = f.GetAtomIds() 
				if len(aids) == 1: 
					featAtom = m.GetAtomWithIdx(aids[0]) 
					hvyNbrs = [x for x in featAtom.GetNeighbors() if x.GetAtomicNum() != 1] 
					if len(hvyNbrs) == 1: 
						f.featDirs = FeatDirUtils.GetDonor1FeatVects(m.GetConformer(-1), aids, scale=(0.6*mergeTol)) 
					#elif len(hvyNbrs) == 2: 
					#	f.featDirs = _GetDonor2FeatVects(m.GetConformer(-1), aids, scale=(0.6*mergeTol)) 
					elif len(hvyNbrs) == 3: 
						f.featDirs = FeatDirUtils.GetDonor3FeatVects(m.GetConformer(-1), aids, scale=(0.6*mergeTol))
			elif f.GetFamily() == 'Aromatic':
				f.featDirs = FeatDirUtils.GetAromaticFeatVects(m.GetConformer(-1), f.GetAtomIds(), f.GetPos(-1), scale=(0.6*mergeTol))
			
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
		feat_indexs = [vector.index(x) for x in vector if x<mergeTol]
		feat_indexs.append(i)
		if (len(feat_indexs) >= minRepeats):
			for feat_index in feat_indexs:
				p4map.AddFeature(global_fmap._feats[feat_index], weight=1)
	Merge = True
	while Merge == True:
		Merge = _MergeFeatPoints(p4map, mergeMetric=mergeMetric, mergeTol=mergeTol, dirMergeMode=dirMergeMode)
	return p4map

def chimera_p4(molecules_sel, mergeTol=2.5, minRepeats=1):
	molecules = molecules_sel.molecules()

	p4map = calc_p4map(molecules, mergeTol=mergeTol, minRepeats=minRepeats)

	for feat in p4map._feats:
		size_sphere = mergeTol/4
		p4_elem = p4_element(shape="sphere", size=size_sphere, origin=chimera.Point(feat.GetPos()[0],feat.GetPos()[1],feat.GetPos()[2]), color=_featColors[feat.GetFamily()])
		p4_elem.draw()
		
		if feat.featDirs:
			ps, fType = feat.featDirs
			for tail, head in ps:
				p4_elem = p4_element(shape="arrow", origin=tail, end=head, color=_featColors[feat.GetFamily()])
				p4_elem.draw() 

	msg = "Chimera pharmacophore done"
	chimera.statusline.show_message(msg)
	
	return True

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