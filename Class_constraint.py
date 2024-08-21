#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:10:16 2020

@author: sarazin
"""

import numpy as np
import math
from DirectDetection import *

class Constraint:
    def __init__(self,name, Type, exp_value, sigma, Block, LineID, position):
        self.name = name
        self.exp_val = float(exp_value)
        self.likelihood = -1e15
        self.sigma = float(sigma)
        self.value = np.array([-1e15])
        self.block = Block
        self.lineID = LineID
        self.position = int(position)
        self.type = Type
        self.DM_mass = 0.0
    
    def set_val(self, value):
        self.value[0] = value
        
    
    def reset_val(self):
        self.value[0] = -1e15
        self.DM_mass = 0.0
        
        
    def collect_value(self,file_out):
        self.reset_val()
        with open(file_out) as data:
            Lines = data.readlines()
            Block_trigger = False
            for line in Lines :
                if( self.block in line) :
                   Block_trigger = True
                if (Block_trigger and self.lineID in line) :
                   self.set_val(float([a for a in line.split(' ') if a != ''][self.position]))
                   break
            if (self.value[0] == -1e15 ) :
                print(("ERROR in collecting the value of observable :",  self.name))
        # Store Dark Matter Mass
            Block_trigger = False
            for line in Lines :
                if ('Block RELIC DENSITY' in line) :
                   Block_trigger = True
                if (Block_trigger and ('   3  ' in line)) :
                   self.DM_mass = float([a for a in line.split(' ') if a != ''][1])	
                   break
         
                
    
    def compute_likelihood(self):
        self.likelihood = -1e15
        if (self.type == 'G') :
            if (self.sigma == 0.0) :
                    print(('ERROR : No sigma given for constraint : ', self.name))
            else:
                    self.likelihood = math.exp(-(float(self.value[0]) - self.exp_val)**2/(2*self.sigma **2))
        if (self.type == 'S'):
            if (self.sigma == 0.0) :
                    print(('ERROR : No sigma given for constraint : ', self.name))
        
            if (self.value[0] < self.exp_val):
                self.likelihood = 1.0
            else:
                self.likelihood = math.exp(-(float(self.value[0]) - self.exp_val)**2/(2*self.sigma **2))
    
        if (self.type == 'Nu') :
            if (self.sigma == 0.0) :
                    print(('ERROR : No sigma given for constraint : ', self.name))
            else:
                    self.likelihood = math.exp(-(float(abs(self.value[0])) - self.exp_val)**2/(2*self.sigma **2))
    
    
        if (self.type == 'DD_SI') :
            self.exp_val = Spin_Independant(self.DM_mass)
            self.sigma = 0.1*self.exp_val
            if(self.value[0] < self.exp_val):
                self.likelihood = 1.0
            else:
                self.likelihood = math.exp(-(float(self.value[0]) - self.exp_val)**2/(2*self.sigma **2))
        #print("\n DM Mass = ",self.DM_mass)
        #print("\n Exp val = ", self.exp_val)
    
        if (self.type == 'DD_SDN') :
            self.exp_val = 	Spin_Dependant_neutron(self.DM_mass)
            self.sigma = 0.1*self.exp_val
            if( self.value[0] < self.exp_val):
                self.likelihood = 1.0
            else:
                self.likelihood = math.exp(-(float(self.value[0]) - self.exp_val)**2/(2*self.sigma **2))
            print("\n DM Mass = ",self.DM_mass)
            print("\n Exp val = ", self.exp_val)
        
        if (self.type == 'DD_SDP') :
            print("\n Experimental value :", self.exp_val)
            self.exp_val = 	Spin_Dependant_proton(self.DM_mass)
            self.sigma = 0.1*self.exp_val
            print("\n Sigma = ",self.sigma)
            if( self.value[0] < self.exp_val):
                self.likelihood = 1.0
            else:
                self.likelihood = math.exp(-(float(self.value[0]) - self.exp_val)**2/(2*self.sigma **2))
            
        
        
        
        
