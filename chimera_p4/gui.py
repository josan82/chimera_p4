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
from core import Controller, Model, open3align, chimera_p4
#from prefs import prefs, _defaults
from base_gui.ui import PlumeBaseDialog

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
	buttons = ('Close')
	default = None
	help = 'https://github.com/josan82/chimera_p4'

	def __init__(self, *args, **kwargs):
		# GUI init
		self.title = 'Pharmacophore'
		self.controller = None

		#Variables
		self._nConformers = tk.IntVar()
		self._mergeTol = tk.DoubleVar()
		self._minRepeats = tk.IntVar()
		#self._p4Id = tk.IntVar()

		# Fire up
		super(p4Dialog, self).__init__(resizable=False, *args, **kwargs)
		#if not chimera.nogui:  # avoid useless errors during development
		#	chimera.extension.manager.registerInstance(self)

		# Set Defaults
		self._set_defaults()

	def _set_defaults(self):
		self._nConformers.set(0)
		self._minRepeats.set(1)
		self._mergeTol.set(1.5)

	def fill_in_ui(self, parent):
		#First frame for selecting the molecules to align/calculate pharmacophore
		self.ui_input_frame = tk.LabelFrame(self.canvas, text='Select molecules to align/calculate pharmacophore')
		self.ui_input_frame.rowconfigure(0, weight=1)
		self.ui_input_frame.columnconfigure(1, weight=1)
		self.ui_molecules = MoleculeScrolledListBox(self.ui_input_frame, listbox_selectmode="extended")
		self.ui_molecules.grid(row=0, columnspan=3, padx=5, pady=5, sticky='news')
		self.ui_input_frame.pack(expand=True, fill='both', padx=5, pady=5)

		#Second frame to perform alignments
		self.ui_o3align_frame = tk.LabelFrame(self.canvas, text="Perform an alignment with open3align")
		tk.Label(self.ui_o3align_frame, text='Number of conformers:').grid(row=0, column=0, padx=5, pady=5, sticky='w')
		self.ui_nconformers = tk.Entry(self.ui_o3align_frame, textvariable=self._nConformers, bg='white', width=6).grid(row=0, column=1, padx=5, pady=5, sticky='w')
		self.ui_o3align_frame.rowconfigure(0, weight=1)
		self.ui_o3align_frame.columnconfigure(1, weight=1)
		self.ui_o3align_btn = tk.Button(self.ui_o3align_frame, text='Align', command=self._cmd_o3align_btn)
		self.ui_o3align_btn.grid(row=0, column=2, padx=5, pady=5)
		self.ui_o3align_frame.pack(expand=True, fill='both', padx=5, pady=5)

		
		#Third frame to perform pharmacophores
		self.ui_p4_frame = tk.LabelFrame(self.canvas, text='Perform a pharmacophore')
		#tk.Label(ui_config_frame, text='Merge Tolerance').grid(row=1, column=0, padx=3, pady=3, sticky='e')
		#self.ui_merge_tol = tk.Entry(ui_config_frame, textvariable=self._mergeTol, width=6).grid(row=1, column=1, padx=3, pady=3)
		tk.Label(self.ui_p4_frame, text='Minimum of repetitions per feature:').grid(row=0, column=0, padx=5, pady=5, sticky='w')
		self.ui_min_repeats = tk.Entry(self.ui_p4_frame, textvariable=self._minRepeats, bg='white', width=3).grid(row=0, column=1, padx=5, pady=5, sticky='w')
		#tk.Label(ui_config_frame, text='Pharmacophore Id').grid(row=3, column=0, padx=3, pady=3, sticky='e')
		#self.ui_p4_Id = tk.Entry(ui_config_frame, textvariable=self._p4Id, width=6).grid(row=3, column=1, padx=3, pady=3)
		self.ui_p4_frame.rowconfigure(0, weight=1)
		self.ui_p4_frame.columnconfigure(1, weight=1)
		self.ui_p4_btn = tk.Button(self.ui_p4_frame, text='Make pharmacophore', command=self._cmd_p4_btn)
		self.ui_p4_btn.grid(row=0, column=2, padx=5, pady=5)
		self.ui_p4_options_btn = tk.Button(self.ui_p4_frame, text="Advanced Options", command=lambda: self.Open_window('ui_input_opt_window', self._fill_ui_input_opt_window))
		self.ui_p4_options_btn.grid(row=1, column=0, padx=5, pady=5, sticky="w")
		self.ui_p4_frame.pack(expand=True, fill='both', padx=5, pady=5)
		
	def _cmd_o3align_btn(self):
		molecules = self.ui_molecules.getvalue()
		nConformers = self._nConformers.get()
		try:
			max_score = open3align(molecules, nConformers=nConformers, _gui=self)
			msg = "Alignment done! Score: {}".format(max_score)
			self.status(msg, color='green', blankAfter=0)
		except Exception as e:
			if len(molecules) < 2:
				self.status('You have to select at least 2 molecules!', color='red', blankAfter=4)
			else:
				self.status('Could not align the molecules!', color='red', blankAfter=4)
	
	def _cmd_p4_btn(self):
		molecules = self.ui_molecules.getvalue()
		minRepeats = self._minRepeats.get()
		mergeTol = self._mergeTol.get()

		try:
			self.status("Working...", blankAfter=0)
			chimera_p4(molecules, mergeTol=mergeTol, minRepeats=minRepeats, _gui=self)
			self.status("Pharmacophore done!", color='green', blankAfter=4)
		except Exception as e:
			if len(molecules) < 1:
				self.status('You have to select at least 1 molecule!', color='red', blankAfter=4)
			else:
				self.status('Could not perform the pharmacophore!', color='red', blankAfter=4)

	def _accept_adv_btn(self):
		self._mergeTol.set(self.ui_input_opt_window.ui_mergeTol.get())

	def _fill_ui_input_opt_window(self):
		# Create TopLevel window
		self.ui_input_opt_window = tk.Toplevel()
		self.Center(self.ui_input_opt_window)
		self.ui_input_opt_window.title("Advanced Options")
		self.ui_input_opt_window.resizable(False, False)

		#First frame for selecting the feature families that the user wants to calculate
		self.ui_features_frame = tk.LabelFrame(self.ui_input_opt_window, text='Select features to calculate')
		self.ui_features_frame.rowconfigure(0, weight=1)
		self.ui_features_frame.columnconfigure(1, weight=1)
		#self.ui_molecules = MoleculeScrolledListBox(self.ui_input_frame, listbox_selectmode="extended")
		#self.ui_molecules.grid(row=0, columnspan=3, padx=5, pady=5, sticky='news')
		self.ui_features_frame.pack(expand=True, fill='both', padx=5, pady=5)

		#Second frame to configure Merge Tolerance parameter
		self.ui_mergeTol_frame = tk.LabelFrame(self.ui_input_opt_window, text="Merge tolerance")
		text_mergeTol = ('The merge tolerance parameter defines the maximum distance \nin which two features of the same family will be considered close \nenough to be merged by the pharmacophore generator.\n')
		tk.Label(self.ui_mergeTol_frame, text=text_mergeTol).grid(row=0, columnspan=2, padx=5, pady=5, sticky='w')
		tk.Label(self.ui_mergeTol_frame, text='Merge tolerance value:').grid(row=1, column=0, padx=5, pady=5, sticky='e')
		self.ui_input_opt_window.ui_mergeTol = tk.Entry(self.ui_mergeTol_frame, bg='white', width=6).grid(row=1, column=1, padx=5, pady=5, sticky='w')
		self.ui_input_opt_window.ui_mergeTol.set(1.5)
		self.ui_mergeTol_frame.rowconfigure(0, weight=1)
		self.ui_mergeTol_frame.columnconfigure(1, weight=1)
		self.ui_mergeTol_frame.pack(expand=True, fill='both', padx=5, pady=5)

		#Third frame to configure other options
		self.ui_other_frame = tk.LabelFrame(self.ui_input_opt_window, text='Other')
		self.ui_other_frame.rowconfigure(0, weight=1)
		self.ui_other_frame.columnconfigure(1, weight=1)
		self.ui_other_frame.pack(expand=True, fill='both', padx=5, pady=5)

		#Cancel and accept buttons
		self.ui_input_opt_window.accept_btn = tk.Button(self.ui_input_opt_window, text="Accept", command=self._accept_adv_btn)
		self.ui_input_opt_window.accept_btn.pack(side='right', padx=5, pady=5) 
		#self.ui_input_opt_window.accept_btn.grid(row=3, column=2, padx=5, pady=5, sticky="e")


	def Open_window(self, window, fill_function):
		"""
		Get sure the window is not opened
		a second time
		Parameters:
		window: window to open
		fill_function: fillin function for window
		"""
		try:
			var_window = window
			var_window.state()
			if window == self.ui_stages_window:
				self.set_stage_variables()
				self.ui_stage_minimiz_tolerance_Entry.configure(state='disabled')
				self.ui_stage_minimiz_maxsteps_Entry.configure(state ='disabled')
				self.ui_stage_barostat_steps_Entry.configure(state='disabled')
				self.ui_stage_pressure_Entry.configure(state='disabled')
			var_window.deiconify()
		except (AttributeError, tk.TclError):
			return fill_function()

	def Center(self, window):
		"""
		Update "requested size" from geometry manager
		"""
		window.update_idletasks()
		x = (window.winfo_screenwidth() -
			 window.winfo_reqwidth()) / 2
		y = (window.winfo_screenheight() -
			 window.winfo_reqheight()) / 2
		window.geometry("+%d+%d" % (x, y))
		window.deiconify()

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
