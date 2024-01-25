import os
import sys
import errno
import pickle

import json
import random
import copy

import numpy as np
from numpy.matlib import repmat
import math

import psychopy.core
import psychopy.event
import psychopy.visual
import psychopy.gui
import psychopy.parallel
import psychopy.tools.monitorunittools


experiment_monitor = psychopy.monitors.Monitor(
            'Experiment Monitor', width=53,
            distance=90)
            
experiment_monitor.setSizePix([1920, 1080])

color=[0,0,255]
color=np.array(color)/127.5 - 1 
            
experiment_window = psychopy.visual.Window(
            monitor=experiment_monitor, fullscr=True, color=color,
            colorSpace='rgb', units='deg', allowGUI=False)
            
experiment_window.flip()
            
psychopy.visual.Circle(
            experiment_window, lineColor=None, fillColor = [-1,-1,-1], 
            fillColorSpace='rgb', radius=.075, units='deg', interpolate=True
        ).draw()
            
experiment_window.flip()
psychopy.core.wait(.3) 
keys = psychopy.event.waitKeys(keyList=['space'])

experiment_window.close()
sys.exit(0)

#        
#checking = True
#color_info = {'r': '0','g': '0','b': '0'}
#
#while checking:
#    
#    col_info = psychopy.gui.DlgFromDict(
#        color_info, title='rgb?',
#        order=['r', 'g', 'b'],
#        screen=1
#    )
#
#    if color_info['r']=='q':
#        checking==False
#        experiment_window.close()
#        sys.exit(0)
#    else:
#        bg_color = np.array([float(color_info['r']), float(color_info['g']), float(color_info['b'])])/127.5 - 1 
#        print(bg_color)
#        
#        backgroundRect = psychopy.visual.Rect(
#            experiment_window, fillColor=bg_color, units='norm', width=2,
#            height=2)
#
#
#        
#
#
##        backgroundRect.draw()
#        experiment_window.color = bg_color
#        experiment_window.clearBuffer(color=True)
#        textObject.draw()
#        psychopy.core.wait(.3) 
#        experiment_window.flip()
#
#        keys = None
#
#        psychopy.core.wait(.2)  # Prevents accidental key presses
#        keys = psychopy.event.waitKeys(keyList=['space'])
##        backgroundRect.draw()
##        textObject.draw()
##        experiment_window.flip()
#            
#        
##        psychopy.visual.Rect(
##            experiment_window,
##            size=5,
##            units='deg',
##            fillColor=colors,
##            fillColorSpace='rgb').draw()
##        psychopy.visual.Circle(
##            experiment_window,
##            lineColor=None,
##            fillColor = [-1,-1,-1], 
##            fillColorSpace='rgb',
##            radius=.075,
##            units='deg'
##        ).draw()
##        experiment_window.flip()
##        
##        waiting_for_space = True
##        while waiting_for_space:
##            keys = psychopy.event.getKeys(keyList=["space"])
##            if 'space' in keys:
##                waiting_for_space=False
#
#