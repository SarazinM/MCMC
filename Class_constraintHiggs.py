#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:10:16 2020

@author: sarazin
"""

import numpy as np
import math


class ConstraintHiggs:
    def __init__(self,name, Type, exp_value, sigma, Block, position):
        self.name = name
        self.exp_val = float(exp_value)
        self.likelihood = -1e15
        self.sigma = float(sigma)
        self.value = np.array([-1e15])
        self.block = Block
        self.position = position
        self.type = Type

    
    def set_val(self, value):
        self.value[0] = value
        
    
    def reset_val(self):
        self.value[0] = -1e15

        
        
    def collect_value(self,file_out):
        self.reset_val()
        Block_lines = []
        with open(file_out, 'r') as file:
            slha_content = file.read()
        #print "\n===========IN THE CLASS============\n"
        Lines= slha_content.splitlines()
        Block_trigger = False
        #print "Looking for ",self.block
        for line in Lines :
            if line.startswith(self.block):
                #print "!!HEADER FOUND!!"
                Block_trigger = True
            elif (Block_trigger and not line.startswith('#')):
                if self.position in line:
                    # print "POSITION FOUND"
                     Block_lines.append(line)
                     Split_Blocklines = Block_lines[0].split()
                     self.value[0] = float(Split_Blocklines[0])
                    # print "Value : ",self.value
                   
        if (self.value[0] == -1e15 ) :
            print(("ERROR in collecting the value of observable :",  self.name))
            

                  
    
    def compute_likelihood(self):
        self.likelihood = -1e15
        if (self.type == 'G') :
            if (self.sigma == 0.0) :
                    print(('ERROR : No sigma given for constraint : ', self.name))
            else:
                    self.likelihood = math.exp(-(float(self.value[0]) - self.exp_val)**2/(2*self.sigma **2))
            if self.value == -1e15:
                 self.likelihood = 1

        if self.type == 'No':
             self.likelihood = 1
             #print self.name, ": ",self.value


        if (self.type == 'S'):
            if (self.sigma == 0.0) :
                    print(('ERROR : No sigma given for constraint : ', self.name))
        
            if (self.value[0] < self.exp_val):
                self.likelihood = 1.0
            else:
                self.likelihood = math.exp(-(float(self.value[0]) - self.exp_val)**2/(2*self.sigma **2))
