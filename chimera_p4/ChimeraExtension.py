#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
import chimera.extension
from Midas.midas_text import addCommand

def cmd_p4(cmdName, args):
	from Midas.midas_text import doExtensionFunc
	from chimera_p4 import chimera_p4
	doExtensionFunc(chimera_p4, args, specInfo=[("molSpec", "molecules_sel", None)])

def cmd_o3align(cmdName, args):
    from Midas.midas_text import doExtensionFunc
    from chimera_p4 import open3align
    doExtensionFunc(open3align, args, specInfo=[("molSpec", "molecules_sel", None)])

addCommand("p4", cmd_p4)
addCommand("open3align", cmd_o3align)

"""
#Code for GUI implementation
class p4Extension(chimera.extension.EMO):

    def name(self):
        return 'Plume Pharmacophore'

    def description(self):
        return "Calculate pharmacophoric model of a set of ligands"

    def categories(self):
        return ['InsiliChem']

    def icon(self):
        return

    def activate(self):
        self.module('gui').showUI()

chimera.extension.manager.registerExtension(p4Extension(__file__))
"""
