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

from rdkit import Chem, RDConfig
from rdkit.Chem import SanitizeMol
from rdkit.Chem import AllChem
from rdkit.Chem.FeatMaps import FeatMaps, FeatMapUtils

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
		""".format(self.color, x1, y1, z1, x2, y2, z2, 0.1, 0.2, 0.8)
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

def calc_p4map(molecules, families=('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic'), mergeMetric=1, mergeTol=1.5, dirMergeMode=0, minRepeats=5):
	rdkit_mols = []
	rdkit_maps = []
	for mol in molecules:
		rdkit_mol, rdkit_map = _chimera_to_rdkit(mol)
		rdkit_mols.append(rdkit_mol)
		rdkit_maps.append(rdkit_map)

	fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef'))
	fmParams = {}
	for k in fdef.GetFeatureFamilies():
		fparams = FeatMaps.FeatMapParams()
		fmParams[k] = fparams

	keep = families
	rawFeats=[]
	for m in rdkit_mols:
		for f in fdef.GetFeaturesForMol(m):
			rawFeats.append(f)
	# filter that list down to only include the ones we're intereted in 
	featList = [f for f in rawFeats if f.GetFamily() in keep]

	fmap = FeatMaps.FeatMap(feats = featList,weights=[1]*len(featList),params=fmParams)
	matrix = FeatMapUtils.GetFeatFeatDistMatrix(fmap, mergeMetric=mergeMetric, mergeTol=mergeTol, dirMergeMode=dirMergeMode, compatFunc=FeatMapUtils.familiesMatch)
	p4map = FeatMaps.FeatMap(params=fmParams)
	for i, vector in enumerate(matrix):
		feat_indexs = [vector.index(x) for x in vector if x<100000]
		feat_indexs.append(i)
		if len(feat_indexs) >= minRepeats:
			for feat_index in feat_indexs:
				p4map.AddFeature(fmap._feats[feat_index], weight=1)
	Merge = True
	while Merge == True:
		Merge = FeatMapUtils.MergeFeatPoints(p4map, mergeMetric=mergeMetric, mergeTol=mergeTol)

	return p4map

def chimera_p4(molecules_sel, mergeTol=1.5):
	molecules = molecules_sel.molecules()

	p4map = calc_p4map(molecules)

	for feat in p4map._feats:
		size_sphere = mergeTol/4
		if feat.GetFamily() == 'Acceptor':
			p4_elem = p4_element(shape="sphere", size=size_sphere, origin=chimera.Point(feat.GetPos()[0],feat.GetPos()[1],feat.GetPos()[2]), color='red')
		elif feat.GetFamily() == 'Donor':
			p4_elem = p4_element(shape="sphere", size=size_sphere, origin=chimera.Point(feat.GetPos()[0],feat.GetPos()[1],feat.GetPos()[2]), color='blue')
		elif feat.GetFamily() == 'PosIonizable':
			p4_elem = p4_element(shape="sphere", size=size_sphere, origin=chimera.Point(feat.GetPos()[0],feat.GetPos()[1],feat.GetPos()[2]), color='green')
		elif feat.GetFamily() == 'NegIonizable':
			p4_elem = p4_element(shape="sphere", size=size_sphere, origin=chimera.Point(feat.GetPos()[0],feat.GetPos()[1],feat.GetPos()[2]), color='gray')
		elif feat.GetFamily() == 'Aromatic':
			p4_elem = p4_element(shape="sphere", size=size_sphere, origin=chimera.Point(feat.GetPos()[0],feat.GetPos()[1],feat.GetPos()[2]), color='yellow')
		
		p4_elem.draw()

	msg = "Chimera pharmacophore is working"
	chimera.statusline.show_message(msg)
	
	return True


	#p4_elem3 = p4_element(shape="arrow", origin=(0,0,0), end=(1,1,1), color='green')