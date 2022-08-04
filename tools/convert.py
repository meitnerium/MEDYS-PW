wnm=2100
cmm1=1E7/wnm
ua=cmm1/2.19475E5
E0wcm2=1E14
E0ua=(E0wcm2/3.51E16)**(0.5)
print("w(nm)=", wnm, "w(cm-1) = " , cmm1, "w(u.a.) = ", ua)
print("E(w/cm2)=", E0wcm2, "E0(u.a.) = " , E0ua)
from tkinter.filedialog import askopenfilename,  asksaveasfilename, askdirectory
from tkinter import messagebox, filedialog, StringVar, OptionMenu, Button, Label, Scrollbar, Tk, \
    Menu, Entry, END, Toplevel, Listbox, VERTICAL, HORIZONTAL, Text, INSERT, StringVar
from tkinter.ttk import Separator
from configparser import ConfigParser
import os
from os import path, listdir
from time import sleep
from numpy import zeros, float, array
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
#from scipy.integrate import simps
import numpy as np


def close(self, _event=None):
    window.destroy()

def callback(*args):
    print( "variable changed!" )


class MainWin(object):
    def __init__(self,window):
        window.wm_title("Pyperspec")
        menubar = Menu(window)


        filemenu = Menu(menubar, tearoff=0)
        #filemenu.add_command(label="Open", command=openconf)
        #filemenu.add_command(label="Config", command=Config)
        #filemenu.add_command(label="Start scan", command=startscan)
        #menubar.add_cascade(label="File", menu=filemenu)
        #menubar.add_cascade(label="Analyze", menu=analysemenu)

        #helpmenu = Menu(menubar, tearoff=0)
        #menubar.add_cascade(label="Help", menu=helpmenu)
        #helpmenu.add_command(label="About", command=About)
        #laseronbtn = Button(window, text='Laser ON',command=laseron).grid(row=0,column=0)
        #laseronbtn.pack(side=RIGHT)
        #laseroffbtn = Button(window, text='Laser OFF',command=laseroff).grid(row=0,column=1)
        #laseroffbtn.pack(side=LEFT)

        wllabel = Label(window, text="Wavelength").grid(row=1,column=0)
        self.nmText = StringVar()
        self.nmText.set( "800" )
        self.nmText.trace("w", self.callback)

        wlentry = Entry(window, textvariable=self.nmText).grid(row=1,column=1)
        nmlabel = Label(window, text="nm").grid(row=1,column=2)
        self.uaText = StringVar()
        cmm1=1E7/float(self.nmText.get())
        ua=cmm1/2.19475E5
        self.uaText.set( str(ua) )
        wlentry2 = Entry(window, textvariable=self.uaText).grid(row=1,column=3)
        ualabel = Label(window, text="u.a.").grid(row=1,column=4)

        intenslabel = Label(window, text="Intensity").grid(row=2,column=0)
        self.Wcm2Text = StringVar()
        self.Wcm2Text.set( "1E14" )
        self.Wcm2Text.trace("w", self.callbackWcm2)
        intensentry = Entry(window, textvariable=self.Wcm2Text).grid(row=2,column=1)
        wcm2label = Label(window, text=r"$W/{cm}^2$").grid(row=2,column=2)
        self.intensuaText = StringVar()
        #cmm1=1E7/float(self.nmText.get())
        #ua=cmm1/2.19475E5
        E0wcm2=1E7/float(self.Wcm2Text.get())
        E0ua=(E0wcm2/3.51E16)**(0.5)
        self.intensuaText.set( str(E0ua) )
        intensuaentry = Entry(window, textvariable=self.intensuaText).grid(row=2,column=3)
        intensualabel = Label(window, text="u.a.").grid(row=2,column=4)



        window.config(menu=menubar)


    

        def on_closing():
            if messagebox.askokcancel("Quit", "Do you want to quit?"):
                window.destroy()


        window.protocol("WM_DELETE_WINDOW", on_closing)
        window.mainloop()


    def callback(self,*args):
        print( "variable changed!" )
        cmm1=1E7/float(self.nmText.get())
        ua=cmm1/2.19475E5
        self.uaText.set( str(ua) )
    def callbackWcm2(self,*args):
        print( "variable changed!" )
        E0wcm2=1E7/float(self.Wcm2Text.get())
        E0ua=(E0wcm2/3.51E16)**(0.5)
        self.intensuaText.set( str(E0ua) )


if __name__ == '__main__':
    window = Tk()
    MainWin(window)
