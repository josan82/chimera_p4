#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chimera
import cStringIO as StringIO
from SplitMolecule.split import molecule_from_atoms
from Matrix import chimera_xform

def chimera_p4(molecules_sel):
	msg = "Chimera pharmacophore is working"
	chimera.statusline.show_message(msg)
	
	return True
