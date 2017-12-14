#!/usr/bin/env python
# encoding: utf-8

# Get used to importing this in your Py27 projects!
from __future__ import print_function, division
# Python stdlib
import Tkinter as tk
# Chimera stuff
import chimera
from chimera.baseDialog import ModelessDialog
from chimera.widgets import MoleculeScrolledListBox

# Additional 3rd parties

# Own
from core import Controller, Model, open3align
#from prefs import prefs, _defaults
from libplume.ui import PlumeBaseDialog

ui = None
def showUI():
	if chimera.nogui:
		tk.Tk().withdraw()
	global ui
	if not ui: 
		ui = p4Dialog()
	model = Model()
	controller = Controller(gui=ui, model=model)
	ui.enter()

ENTRY_STYLE = {
	'background': 'white',
	'borderwidth': 1,
	'highlightthickness': 0,
	'insertwidth': 1,
}
BUTTON_STYLE = {
	'borderwidth': 1,
	'highlightthickness': 0,
}

class p4Dialog(PlumeBaseDialog):
	buttons = ('Run', 'Close')
	default = None
	help = 'https://www.insilichem.com'

	def __init__(self, *args, **kwargs):
		# GUI init
		self.title = 'Pharmacophore'
		self.controller = None

		#Variables
		self._mergeTol = tk.IntVar()
		self._minRepeats = tk.IntVar()
		self._p4Id = tk.IntVar()

		# Fire up
		super(p4Dialog, self).__init__(resizable=False, *args, **kwargs)
		#if not chimera.nogui:  # avoid useless errors during development
		#	chimera.extension.manager.registerInstance(self)

		# Set Defaults
		self._set_defaults()

	def _set_defaults(self):
		pass

	def fill_in_ui(self, parent):
		#First frame for selecting the molecules to align/calculate pharmacophore
		self.ui_input_frame = tk.LabelFrame(self.canvas, text='Select molecules to align/calculate pharmacophore')
		self.ui_input_frame.rowconfigure(0, weight=1)
		self.ui_input_frame.columnconfigure(1, weight=1)
		self.ui_molecules = MoleculeScrolledListBox(self.ui_input_frame, listbox_selectmode="extended")
		self.ui_molecules.grid(row=0, columnspan=3, padx=5, pady=5, sticky='news')
		self.ui_input_frame.pack(expand=True, fill='both', padx=5, pady=5)

		#Second frame to perform alignments
		self.ui_o3align_frame = tk.LabelFrame(self.canvas, text="Perfom an alignment with open3align")
		self.ui_input_frame.rowconfigure(1, weight=1)
		self.ui_input_frame.columnconfigure(1, weight=1)
		self.ui_o3align_btn = tk.Button(self.ui_o3align_frame, text='Align!', command=self._cmd_o3align_btn)
		self.ui_input_frame.pack(expand=True, fill='both', padx=5, pady=5)

		"""
		#Third frame to perform pharmacophores
		ui_config_frame = tk.LabelFrame(self.canvas, text='Configuration parameters')
		tk.Label(ui_config_frame, text='Merge Tolerance').grid(row=1, column=0, padx=3, pady=3, sticky='e')
		self.ui_merge_tol = tk.Entry(ui_config_frame, textvariable=self._mergeTol, width=6).grid(row=1, column=1, padx=3, pady=3)
		tk.Label(ui_config_frame, text='Minimum appearences for feature').grid(row=2, column=0, padx=3, pady=3, sticky='e')
		self.ui_min_repeats = tk.Entry(ui_config_frame, textvariable=self._minRepeats, width=6).grid(row=2, column=1, padx=3, pady=3)
		tk.Label(ui_config_frame, text='Pharmacophore Id').grid(row=3, column=0, padx=3, pady=3, sticky='e')
		self.ui_p4_Id = tk.Entry(ui_config_frame, textvariable=self._p4Id, width=6).grid(row=3, column=1, padx=3, pady=3)
		ui_config_frame.pack(expand=True, fill='both', padx=5, pady=5)
		"""

	def _cmd_o3align_btn(self):
		molecules = self.ui_molecules.getvalue()
		try:
			open3align(molecules)
		except Exception as e:
			self.status('Could not align the molecules!', color='red', blankAfter=4)


	def Run(self):
		"""
		Default! Triggered action if you click on a Run button
		"""
		self.Close()

	def Close(self):
		"""
		Default! Triggered action if you click on the Close button
		"""
		global ui
		ui = None
		ModelessDialog.Close(self)
		self.destroy()

	# Below this line, implement all your custom methods for the GUI.
	def load_controller(self):
		pass
