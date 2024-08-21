#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:47:02 2020

@author: sarazin
"""

import numpy as np
import math

class observable:
    def __init__(self,name,block, LineID,position): #Position 1 for 1 entry, position 2 for 2 entries
        self.name = name
        self.value = np.array([-1e15])
        self.block = block
        self.lineID = LineID
        self.position = int(position)
        
    def set_val(self, value):
        self.value[0] = value
        
    
    def reset_val(self):
        self.value[0] = -1e15
        
        
    def collect_value(self,file_out):
        self.reset_val()
        with open(file_out) as data:
            Lines = data.readlines()
            Block_trigger = False
            for line in Lines :
                if self.block in line :
                   Block_trigger = True
                if Block_trigger and (self.lineID in line) :
		   #print line
                   self.set_val(float([a for a in line.split(' ') if a != ''][self.position]))
                   break
            if self.value[0] == -1e15  :
                print(("ERROR in collecting the value of observable :",  self.name))

    def GammaTotH(self,file_out):
        #print("ENTER IN THE FONCTION")
        self.reset_val()
        Block_lines = []
        with open(file_out, 'r') as file:
            slha_content = file.read()
        Lines= slha_content.splitlines()
        for line in Lines :
            if line.startswith(self.block):
                Block_lines.append(line)
                Split_Blocklines = Block_lines[0].split()
               # print("IN THE CLASS OBSERVABLE : \n")
               # print(Split_Blocklines)
                self.value[0] = float(Split_Blocklines[2])
