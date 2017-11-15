#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Midas.midas_text import addCommand

def cmd_p4(cmdName, args):
    from Midas.midas_text import doExtensionFunc
    from chimera_p4 import chimera_p4
    doExtensionFunc(chimera_p4, args, specInfo=[("molSpec", "molecules_sel", None)])

addCommand("p4", cmd_p4)
