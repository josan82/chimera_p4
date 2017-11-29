#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chimera
import cStringIO as StringIO
from SplitMolecule.split import molecule_from_atoms
from Matrix import chimera_xform

import numpy
from rdkit import Chem, RDConfig
from rdkit.Chem import MolFromPDBBlock, FastFindRings, SanitizeMol
from rdkit.Chem.rdMolAlign import GetAlignmentTransform, GetO3A, AlignMol, GetO3AForProbeConfs
from rdkit.Chem.AllChem import EmbedMolecule, ETKDG, MMFFOptimizeMolecule, MMFFGetMoleculeProperties, EmbedMultipleConfs, AddHs, RemoveHs

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

#Function without testing
def _rdkit_to_chimera(molecule, confId=0, name_res="RES"):
	from chimera import openModels, Molecule, Element, Coord, Point

	chimera_mol = Molecule() #creates a blank Chimera molecule
	res = chimera_mol.newResidue(name_res, " ", 1, " ")

	#Creating atoms
	for atom in molecule.GetAtoms():
		Chimera_atom = chimera_mol.newAtom(atom.GetSymbol(), Element(atom.GetSymbol()))
		Chimera_atom.serialNumber = atom.GetIdx()

		atom_position = Point()
		atom_position[0] = list(molecule.GetConformer(id=confId).GetAtomPosition(atom.GetIdx()))[0]
		atom_position[1] = list(molecule.GetConformer(id=confId).GetAtomPosition(atom.GetIdx()))[1]
		atom_position[2] = list(molecule.GetConformer(id=confId).GetAtomPosition(atom.GetIdx()))[2]

		Chimera_atom.setCoord(Coord(atom_position))
		
		res.addAtom(Chimera_atom)

	#Creating bonds
	for bond in molecule.GetBonds():
		Chimera_bond = chimera_mol.newBond(_find_by_sn(chimera_mol, bond.GetBeginAtom().GetIdx()), _find_by_sn(chimera_mol, bond.GetEndAtom().GetIdx()))

	chimera_mol.computeIdatmTypes()

	return chimera_mol

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

'''
	rdk_reference, ref_map = _chimera_to_rdkit(molecules[0])
	for rdkit_atom in rdk_reference.GetAtoms():
		chimera_atom = list(ref_map.keys())[list(ref_map.values()).index(rdkit_atom.GetIdx())]
		if rdkit_atom.GetIdx() == 0:
			chimera_atom.setCoord(Coord(0.414, 30.0, 4.103))
	msg = "Coords {}".format(ref_map)
	chimera.statusline.show_message(msg, blankAfter=5)
	for atom in Chimera_mol.atoms:
		g.add_node(atom.serialNumber, atom=str(atom.element), 
			pos=list(atom.coord()))
	
	#Creating atoms
	for atom in ChemGraph.nodes():
		RDKit_id[atom] = em.AddAtom(Chem.Atom(ChemGraph.node[atom]['atom']))
		em.GetAtomWithIdx(RDKit_id[atom]).SetIntProp('Id_cg', atom)
		if 'pos' in ChemGraph.node[atom]:
			
			em.GetConformer().SetAtomPosition(RDKit_id[atom], 
				ChemGraph.node[atom]['pos'])
	#Creating atoms (nodes)
	for atom in RDKit_mol.GetAtoms():
		if(RDKit_mol.GetConformers()): #Exists a conformer with 3D pos
			g.add_node(atom.GetIntProp('Id_cg'), atom=atom.GetSymbol(), 
				pos=list(RDKit_mol.GetConformer().GetAtomPosition(atom.GetIdx())))
		else:
			g.add_node(atom.GetIntProp('Id_cg'), atom=atom.GetSymbol())
			
	#Creating bonds
	for bond in RDKit_mol.GetBonds():
			g.add_edge(bond.GetBeginAtom().GetIntProp('Id_cg'), 
				bond.GetEndAtom().GetIntProp('Id_cg'), 
				type=bond.GetIntProp('Order_cg'))
	em = Molecule() #creates a blank Chimera molecule
	res = em.newResidue(name_res, " ", 1, " ")
	#Creating atoms
	for atom in ChemGraph.nodes():
		Chimera_atom = em.newAtom(ChemGraph.node[atom]['atom'], 
			Element(ChemGraph.node[atom]['atom']))
		Chimera_atom.serialNumber = atom
		if 'pos' in ChemGraph.node[atom]:
			Chimera_atom.setCoord(Coord(ChemGraph.node[atom]['pos'][0], 
				ChemGraph.node[atom]['pos'][1],
				ChemGraph.node[atom]['pos'][2]))
		res.addAtom(Chimera_atom)
	#Creating bonds
	for edge in ChemGraph.edges():
		Chimera_bond = em.newBond(find_by_sn(em, edge[0]), 
			find_by_sn(em, edge[1]))
		Chimera_bond.order = ChemGraph.edge[edge[0]][edge[1]]['type']
	em.computeIdatmTypes()
	return em
'''