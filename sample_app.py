#!/usr/bin/env python

import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

import os

class Application(tk.Frame):
	def __init__(self,master=None):
		tk.Frame.__init__(self,master)

		self.grid(sticky=tk.N+tk.S+tk.E+tk.W)
		# self.createWidgets()
		self.createNotebook()

		self.s = ttk.Style()


	def createWidgets(self):
		top = self.winfo_toplevel()
		top.rowconfigure(0,weight=1)
		top.columnconfigure(0,weight=1)
		self.rowconfigure(0,weight=1)
		self.columnconfigure(0,weight=1)

		self.quitButton = tk.Button(self, text='Quit',command=self.quit)
		self.quitButton.grid(row=0,column=0,sticky=tk.N+tk.S+tk.E+tk.W)


	def createNotebook(self):
		self.note = ttk.Notebook(self)
		self.note.grid()

		self.tab1 = ttk.Frame(self.note)
		self.note.add(self.tab1,text='Channels')
		self.tab2 = ttk.Frame(self.note)
		self.note.add(self.tab2,text='Files',state='disabled')
		self.tab3 = ttk.Frame(self.note)
		self.note.add(self.tab3,text='Alignment',state='disabled')

		##### Tab 1 #####
		self.nextB = tk.Button(self.tab1,text='Next',command=lambda:self.note.select(self.tab2),state=tk.DISABLED)
		self.nextB.grid(row=4,column=2)

		self.cs = channelMaster(self.tab1,'Structural channel: ',0)
		self.cs.dButton.bind('<Button-1>',lambda e:self.activate(self.nextB,self.note,self.tab2))
		self.c2 = channelMaster(self.tab1,'Channel 2: ',1)
		self.c3 = channelMaster(self.tab1,'Channel 3: ',2)
		self.c4 = channelMaster(self.tab1,'Channel 4: ',3)

		##### Tab 2 #####
		self.filegenB = tk.Button(self.tab2,text='Generate file list',command=self.filegen)
		self.filegenB.grid(row=0)
		
	def activate(self,button,nb,ntab):
		button['state'] = tk.NORMAL
		button['bg'] = 'green'
		nb.tab(ntab,state='normal')

	def filegen(self):
		Lc = [self.cs,self.c2,self.c3,self.c4]

		for i,c in enumerate(Lc):
			if c.dir != None:
				files = os.listdir(c.dir)
			else:
				files = []

			c.files = []
			for f in files:
				if 'Probabilities' in f and 'h5' in f:
					c.files.append(f)

		self.fnums = {}
		self.Lnums = []
		for f in self.cs.files:
			self.fnums[f.split('_')[1]] = tk.IntVar()
			self.Lnums.append(int(f.split('_')[1]))

		print(self.Lnums)

		self.display_files()

	def display_files(self):

		#Header 
		ttk.Label(self.tab2,text='Use?').grid(row=1,column=0)
		ttk.Label(self.tab2,text='Number').grid(row=1,column=1)
		ttk.Label(self.tab2,text='Structural').grid(row=1,column=2)
		ttk.Label(self.tab3)

		for i,n in enumerate(self.Lnums):
			print(i,n)
			key = str(n)
			if len(key) == 1:
				key = '0'+ key
			cb = tk.Checkbutton(self.tab2, text=n, variable=self.fnums[key])
			cb.select()
			cb.grid(row=i+2,column=0)




class channelMaster:
	def __init__(self,parent,text,n):
		ttk.Label(parent,text=text).grid(row=n,column=1)
		self.dir = None
		self.dButton = ttk.Button(parent,text='Select Folder',command=self.open_dir)
		self.dButton.grid(row=n,column=2)

	def open_dir(self):
		self.dir = filedialog.askdirectory()
		self.dButton['text'] = self.dir

app = Application()
app.master.title('Sample Application')
app.mainloop()