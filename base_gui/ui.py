#!/usr/bin/env python
# encoding: utf-8


from __future__ import print_function, division
# Python stdlib
import Tkinter as tk
import Tix
import webbrowser as web
import json
from urllib2 import urlopen, HTTPError, URLError
from distutils.version import LooseVersion
from threading import Thread
# Chimera
import Pmw
import chimera
from chimera.baseDialog import ModelessDialog
from chimera.widgets import MoleculeScrolledListBox

# Fix strange bug that can happen
# with newly created conda environments
Tix._default_root = tk._default_root

CHIMERA_BG = chimera.tkgui.app.cget('bg')
STYLES = {
    tk.Entry: {
        'background': 'white',
        'borderwidth': 1,
        'highlightthickness': 0,
        'insertwidth': 1,
    },
    Pmw.EntryField: {
        'entry_background': 'white',
        'entry_borderwidth': 1,
        'entry_highlightthickness': 0,
        'entry_insertwidth': 1,
    },
    tk.Button: {
        'borderwidth': 1,
        'highlightthickness': 0,
    },
    tk.Checkbutton: {
        'highlightbackground': CHIMERA_BG,
        'activebackground': CHIMERA_BG,
    },
    tk.Radiobutton: {
        'highlightbackground': CHIMERA_BG,
        'activebackground': CHIMERA_BG,
    },
    tk.OptionMenu: {
        'borderwidth': 1,
    },
    Pmw.OptionMenu: {
        'menubutton_borderwidth': 1,
        'menu_relief': 'flat',
        'menu_activeborderwidth': 0,
        'menu_activebackground': '#DDD',
        'menu_borderwidth': 1,
        'menu_background': 'white',
        'hull_borderwidth': 0,
    },
    Pmw.ComboBox: {
        'entry_borderwidth': 1,
        'entry_highlightthickness': 0,
        'entry_background': 'white',
        'arrowbutton_borderwidth': 1,
        'arrowbutton_relief': 'flat',
        'arrowbutton_highlightthickness': 0,
        'listbox_borderwidth': 1,
        'listbox_background': 'white',
        'listbox_relief': 'ridge',
        'listbox_highlightthickness': 0,
        'scrolledlist_hull_borderwidth': 0
    },
    Pmw.ScrolledListBox: {
        'listbox_borderwidth': 1,
        'listbox_background': 'white',
        'listbox_relief': 'ridge',
        'listbox_highlightthickness': 0,
        'listbox_selectbackground': '#DDD',
        'listbox_selectborderwidth': 0
    },
    MoleculeScrolledListBox: {
        'listbox_borderwidth': 1,
        'listbox_background': 'white',
        'listbox_highlightthickness': 0,
    },
    tk.Scale: {
        'borderwidth': 1,
        'highlightthickness': 0,
        'sliderrelief': 'flat'
    },
    tk.Button: {
        'borderwidth': 1,
        'highlightthickness': 0,
    },
    tk.Checkbutton: {
        'highlightbackground': CHIMERA_BG,
        'activebackground': CHIMERA_BG,
    }
}


class PlumeBaseDialog(ModelessDialog):

    default = None
    provideStatus = True
    statusPosition = 'left'
    help = 'https://www.insilichem.com'
    VERSION = None
    VERSION_URL = None
    overMaster = True

    def __init__(self, master=None, with_logo=True, check_version=True,
                 callback=None, *args, **kwargs):
        if master is None:
            master = chimera.tkgui.app
        self.with_logo = with_logo
        self.callback = callback
        self.ui_labels = {}
        # Fire up
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)
        super(PlumeBaseDialog, self).__init__(master=master, *args, **kwargs)
        # Fix styles
        self._fix_styles(*self.buttonWidgets.values())
        self._hidden_files_fix()
        if check_version and None not in (self.VERSION, self.VERSION_URL):
            thread = Thread(target=self.check_version)
            thread.start()
            while thread.isAlive():
                chimera.tkgui.app.update()

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def _fix_styles(self, *widgets):
        for widget in widgets:
            try:
                widget.configure(**STYLES[widget.__class__])
            except Exception as e:
                print('Error fixing styles:', type(e), str(e))

    def _team_logo(self, parent):
        parent = tk.Frame(parent)
        # InsiliChem copyright
        bg = chimera.tkgui.app.cget('bg')
        img_data = LOGO_BASE64
        img = tk.PhotoImage(data=img_data)
        logo = tk.Button(parent, image=img, background=bg, borderwidth=0,
                         activebackground=bg, highlightcolor=bg, cursor="hand2",
                         command=lambda *a: web.open_new(r"http://www.insilichem.com/"))
        logo.image = img
        text = tk.Text(parent, background=bg, borderwidth=0, height=4, width=30)
        hrefs = _HyperlinkManager(text)
        big = text.tag_config("big", font="-size 18", foreground="#367159")
        text.insert(tk.INSERT, "InsiliChem", "big")
        text.insert(tk.INSERT, "\nDeveloped by ")
        text.insert(tk.INSERT, "@jaimergp", hrefs.add(lambda *a: web.open_new(r"https://github.com/jaimergp")))
        text.insert(tk.INSERT, "\nat Maréchal Group, UAB, Spain")
        text.configure(state='disabled')
        logo.grid(row=5, column=0, sticky='we', padx=5, pady=3)
        text.grid(row=5, column=1, sticky='we', padx=5, pady=3)
        return parent

    def fillInUI(self, parent):
        # Create main window
        parent.pack(expand=True, fill='both')
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

        self.fill_in_ui(self.canvas)  # Override this method!

        if self.with_logo:
            self.ui_logo = self._team_logo(parent)
            self.ui_logo.pack(padx=5, pady=5)
        to_fix = [getattr(self, attr) for attr in dir(self) if attr.startswith('ui_')]
        self._fix_styles(*to_fix)

    def fill_in_ui(self, canvas):
        raise NotImplementedError

    def Close(self):
        chimera.extension.manager.deregisterInstance(self)
        if callable(self.callback):
            self.callback()
        ModelessDialog.Close(self)
    Quit = Close

    def auto_grid(self, parent, grid, resize_columns=(1,), label_sep=':', **options):
        """
        Auto grid an ordered matrix of Tkinter widgets.

        Parameters
        ----------
        parent : tk.Widget
            The widget that will host the widgets on the grid
        grid : list of list of tk.Widget
            A row x columns matrix of widgets. It is built on lists.
            Each list in the toplevel list represents a row. Each row
            contains widgets, tuples or strings, in column order.
            If it's a widget, it will be grid at the row i (index of first level
            list) and column j (index of second level list).
            If a tuple of widgets is found instead of a naked widget,
            they will be packed in a frame, and grid'ed as a single cell.
            If it's a string, a Label will be created with that text, and grid'ed.

            For example:
            >>> grid = [['A custom label', widget_0_1, widget_0_2], # first row
            >>>         [widget_1_0, widget_1_1, widget_1_2],       # second row
            >>>         [widget_2_0, widget_2_1, (widgets @ 2_2)]]  # third row

        """
        for column in resize_columns:
            parent.columnconfigure(column, weight=int(100/len(resize_columns)))
        _kwargs = {'padx': 2, 'pady': 2, 'ipadx': 2, 'ipady': 2}
        _kwargs.update(options)
        for i, row in enumerate(grid):
            for j, item in enumerate(row):
                kwargs = _kwargs.copy()
                sticky = 'ew'
                if isinstance(item, tuple):
                    frame = tk.Frame(parent)
                    self.auto_pack(frame, item, side='left', padx=2, pady=2, expand=True, fill='x',
                                   label_sep=label_sep)
                    item = frame
                elif isinstance(item, basestring):
                    sticky = 'e'
                    label = self.ui_labels[item] = tk.Label(parent, text=item + label_sep if item else '')
                    item = label
                elif isinstance(item, tk.Checkbutton):
                    sticky = 'w'
                if 'sticky' not in kwargs:
                    kwargs['sticky'] = sticky
                item.grid(in_=parent, row=i, column=j, **kwargs)
                self._fix_styles(item)

    def auto_pack(self, parent, widgets, label_sep=':', **kwargs):
        for widget in widgets:
            options = kwargs.copy()
            if isinstance(widget, basestring):
                label = self.ui_labels[widget] = tk.Label(parent, text=widget + label_sep if widget else '')
                widget = label
            if isinstance(widget, (tk.Button, tk.Label)):
                options['expand'] = False
            widget.pack(in_=parent, **options)
            self._fix_styles(widget)

    def check_version(self):
        """
        Loads a webpage that returns JSON content with 'tag_name' entry,
        as in https://api.github.com/repos/{owner}/{repo}/releases/latest,
        and compares the latest tag with the current version.
        """
        try:
            response = urlopen(self.VERSION_URL)
        except (HTTPError, URLError) as e:
            print('! Could not obtain version info from', self.VERSION_URL)
            print('  Reason:', str(e))
        else:
            content = json.loads(response.read())
            try:
                latest_version = content['tag_name']
            except KeyError:
                print('! Package contains no release info at', self.VERSION_URL)
                return
            if LooseVersion(latest_version[1:]) > LooseVersion(self.VERSION):
                msg = 'New version {} available! You are using version v{}'
                self.status(msg.format(latest_version, self.VERSION), color='blue',
                            blankAfter=5)

    @staticmethod
    def _hidden_files_fix():
        """
        Call a dummy dialog with an impossible option to initialize the file
        dialog without really getting a dialog window; then set the magic
        variables accordingly
        """
        call = chimera.tkgui.app._root().call
        try:
            try:
                call('tk_getOpenFile', '-foobarbaz')
            except tk.TclError:
                pass
            call('set', '::tk::dialog::file::showHiddenBtn', '1')
            call('set', '::tk::dialog::file::showHiddenVar', '0')
        except:
            pass


class _HyperlinkManager:
    """
    from http://effbot.org/zone/tkinter-text-hyperlink.htm
    """

    def __init__(self, text):

        self.text = text

        self.text.tag_config("hyper", foreground="#367159", underline=1)

        self.text.tag_bind("hyper", "<Enter>", self._enter)
        self.text.tag_bind("hyper", "<Leave>", self._leave)
        self.text.tag_bind("hyper", "<Button-1>", self._click)

        self.reset()

    def reset(self):
        self.links = {}

    def add(self, action):
        # add an action to the manager.  returns tags to use in
        # associated text widget
        tag = "hyper-%d" % len(self.links)
        self.links[tag] = action
        return "hyper", tag

    def _enter(self, event):
        self.text.config(cursor="hand2")

    def _leave(self, event):
        self.text.config(cursor="")

    def _click(self, event):
        for tag in self.text.tag_names(tk.CURRENT):
            if tag[:6] == "hyper-":
                self.links[tag]()
                return


LOGO_BASE64 = r''.join("""
R0lGODlhZABrAOeeADhwWDBzWTlxWTF0WjpyWjtzWzx0XD11XT52Xj93X0B4YEF5YUJ6
YkN7Yz98aEt6Y0R8ZEB9aUx7ZEV9ZU18ZUZ+Zk59Zkt+bE9+Z0x/bVB/aE2AblGAaU6Bb
1KBak+CcFOCa1CDcVSDbFGEclWEbVKFc1OHdFSIdVuGdVWJdlyHdlaKd12Id1eLeF6JeF
iMeV+KeWCLemGMe2KNfGOOfWSPfmWQf2aSgGeTgW6RgWiUgm+SgmmVg3CTg3GUhHKVhXO
WhnSXh3WYiHaZiXebinici3mdjICcjXqejYGdjnufjoKej3ygj4OfkH2hkISgkYGhl4Wh
koKimIaik4OjmYeklIWlm4illYmmloennYqnmIionouomYmpn4ypmoqqoI2qm4uroY6rn
JGqopKro5OspJStpZWuppiuoJavp5ewqJixqZqyqpuzq5y0rJ21rZ62rqS1rp+3r6W2r6
a3sKe4sai5sqm6tKq7tau8tqy9t62+uK6/ubW/uq/BurbAu7DCu7jBvLLDvLnCvbPEvbr
DvrvEv7zFwL3Hwr7Iw7/JxMDKxcHLxsLMx8PNyMrMycTOycbPysvOys3Py87QzM/RztDS
z9HT0NLU0dPV0tTW09XX1NbY1dja1v///////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
///////yH5BAEKAP8ALAAAAABkAGsAQAj+ADsJHEiwoMGDCA1ienTITpgiNEJAMAAhU4Y
PK2wU8SKnECNKmxKKHEmyJME9BgCkaUQAgMuXMGPKBMBAU8uZOGfGoeIyTCaTQDtFUpBz
5iMylIIqJRiJyqSURWMGAjonqssCU5cOvNRokBssRX78GPIEDR9FSbUOzCPAKoAiBX1YV
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
QADs""".splitlines()[1:])
