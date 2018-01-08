#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import chimera
from chimera import runCommand
from SimpleSession import registerAttribute

try:
	from cStringIO import StringIO
except ImportError:
	from StringIO import StringIO
from textwrap import dedent
from Bld2VRML import openFileObject as openBildFileObject

import os
import math

from rdkit import Chem
from rdkit.Chem import FastFindRings
from rdkit.Chem.rdMolAlign import GetO3A, AlignMol, GetO3AForProbeConfs
from rdkit.Chem.FeatMaps import FeatMaps, FeatMapUtils as FMU
from rdkit.Chem.Features import FeatDirUtilsRD as FeatDirUtils
from rdkit.Chem.AllChem import BuildFeatureFactory, MMFFGetMoleculeProperties, EmbedMultipleConfs, AddHs, RemoveHs

from aux_functions import _chimera_to_rdkit, _GetAcceptor1FeatVects, _GetDonor2FeatVects, _MergeFeatPoints, _return_atom_positions, _apply_atom_positions, _del_chimeraHs


FEATURES_FILE = os.path.join(os.path.dirname(__file__), 'BaseFeatures.fdef')

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
			chimera.openModels.add(vrml, baseId=self._id, subid=self._subid, shareXform=0)
			self._subid += 1
		return vrml

def draw_p4legend(families):
	#Delete previous labels
	runCommand("2dlabels delete *")

	#Create one label for each family
	for i, family in enumerate(families):
		label_id = "p4_label_" + str(i)
		label_text = family
		label_color = str(_featColors[family]).replace(" ", ",")
		label_xpos = str(0.02)
		label_ypos = str(0.95 - (0.035)*i)
		label_size = str(16)
		label_command = "2dlabels create " + label_id + " text '\u2588 " + label_text + "' color " + label_color + " size " + label_size + " xpos " + label_xpos + " ypos " + label_ypos
		runCommand(label_command)

def calc_p4map(molecules, families=('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic', 'LumpedHydrophobe'), mergeMetric=1, mergeTol=1.5, dirMergeMode=1, minRepeats=1, showVectors=True):
	rdkit_mols = []
	rdkit_maps = []
	for mol in molecules:
		_del_chimeraHs(mol)
		rdkit_mol, rdkit_map = _chimera_to_rdkit(mol)
		rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
		rdkit_mols.append(rdkit_mol)
		rdkit_maps.append(rdkit_map)

	fdef = BuildFeatureFactory(FEATURES_FILE)
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
				#elif f.GetFamily() == 'Aromatic':
				#	f.featDirs = FeatDirUtils.GetAromaticFeatVects(m.GetConformer(-1), f.GetAtomIds(), f.GetPos(-1), scale=(mergeTol))				
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

def chimera_p4(molecules_sel, mergeTol=1.5, minRepeats=1, showVectors=True, families=('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic', 'LumpedHydrophobe'), showLegend=True, _gui=None):
	chimera.openModels.remove(chimera.openModels.list(id=100))
	registerAttribute(chimera.Bond, "order")
	
	if _gui:
		molecules = molecules_sel
	else:
		molecules = molecules_sel.molecules()

	if not len(molecules) > 0:
		raise chimera.UserError("At least 1 molecule is needed to do a pharmacophore")

	p4map = calc_p4map(molecules, families=families, mergeTol=mergeTol*mergeTol, minRepeats=minRepeats, showVectors=showVectors)

	for feat in p4map._feats:
		if feat.GetFamily() != 'Donor':
			p4_elem = p4_element(shape="sphere", size=(mergeTol/2), origin=chimera.Point(feat.GetPos()[0],feat.GetPos()[1],feat.GetPos()[2]), color=_featColors[feat.GetFamily()])
			p4_elem.draw()
		elif feat.featDirs:
			ps, fType = feat.featDirs			
			for tail, head in ps:
				if fType == 'linear':
					p4_elem = p4_element(shape="arrow", origin=tail, end=head, color=_featColors[feat.GetFamily()])
				elif fType =='cone':
					p4_elem = p4_element(shape="cone", origin=head, end=tail, color=_featColors[feat.GetFamily()], size=(1.33))
				p4_elem.draw()
	
	if showLegend:	
		draw_p4legend(families)
	else:
		runCommand("2dlabels delete *")

	if not _gui:
		msg = "Chimera pharmacophore done"
		chimera.statusline.show_message(msg)
	
	return True

#Function based on plume subalign
def align_o3a(reference, probe, sanitize=True, nConformers=0, **kwargs):
	_del_chimeraHs(reference)
	_del_chimeraHs(probe)
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

		o3as = GetO3AForProbeConfs(rdk_probe, rdk_reference, numThreads=0)
	
		highest_conf_score = 0.0
		for conf_id, o3a in enumerate(o3as):
			score_conf = o3a.Score()
			if score_conf > highest_conf_score:
				highest_conf_score = score_conf
				o3a_result = o3a
				highest_conf_id = conf_id
	else:
		o3a_result = GetO3A(rdk_probe, rdk_reference, probe_params, reference_params)
		highest_conf_id = 0

	o3a_result.Align()
	new_probe_pos = _return_atom_positions(rdk_probe, probe, probe_map, rdkit_confId=highest_conf_id)

	return o3a_result.Score(), new_probe_pos  #return the alignment score and the new atom positions

#New function
def open3align(molecules_sel, transform=True, nConformers=0, reference=None, _gui=None):
	#Delete possible previous pharmacophores
	chimera.openModels.remove(chimera.openModels.list(id=100))
	runCommand("2dlabels delete *")
	
	registerAttribute(chimera.Bond, "order")
	
	if _gui:
		molecules = molecules_sel
		references = [reference] if reference else molecules
	else:
		molecules = molecules_sel.molecules()
		if reference and reference.molecules() not in molecules:
			raise chimera.UserError("Reference must belong to selected molecules")
		references = [reference.molecules()] if reference else molecules

	if not len(molecules) > 1:
		raise chimera.UserError("At least 2 molecules are needed to do an alignment")
	
	#Calculating the best scored alignment
	max_score = 0.0
	for i, reference in enumerate(references):
		msg = "Processing molecule {} of {}".format(i+1, len(references))
		if _gui:
			_gui.status(msg, blankAfter=0)
		else:
			chimera.statusline.show_message(msg)
		align_score = 0.0
		new_atom_pos = {}
		for probe in molecules:
			score, new_atom_pos[probe] = align_o3a(reference, probe, nConformers=nConformers)
			align_score += score
		
		if align_score > max_score:
			max_score = align_score
			max_new_atom_pos = new_atom_pos
		
	if transform:
		for mol in max_new_atom_pos.keys():
			_apply_atom_positions(mol, max_new_atom_pos[mol])

	if not _gui:
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