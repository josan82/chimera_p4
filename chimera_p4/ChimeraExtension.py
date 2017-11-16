#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
import chimera.extension

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
