
import random as rnd
import numpy as np
import copy

class parameter_log :
# reading/writing of the input file. 
#name -> string : name of the variable
# minim, max -> float: range definitions
# width -> float : width in which we will take another value from the first one. taken between 0 and 1.
	def __init__(self, Name, minim, maxi, Width, sign):
		self.name = Name
		self.min = np.float64(minim) # In log scale
		self.max = np.float64(maxi) # in log scale
		self.width = float(Width)
		self.sign = float(sign)
		self.val = np.array([-1e15])
		self.prev_val = np.array([-1e15])

	def set_val(self,value, sign):
		self.val[0] = sign * 10**value

	def initialise(self):
		self.set_val(rnd.uniform(self.min,self.max), self.sign)
		self.prev_val[0] = self.val[0]#10**(self.val[0])*self.sign


	def accepted(self):
		self.prev_val = copy.copy(self.val)


	def jump(self):
		sig = int(np.sign(self.prev_val[0]))
		#print('In jump function')
		step = (self.max - self.min) * self.width
		#'''
		proposal = rnd.gauss(np.log10(self.prev_val[0]/sig), step )
		#print(self.val[0], proposal, step)
		while proposal < self.min or proposal > self.max :	#ensure that gaussian value is within required range            	
			proposal = rnd.gauss(np.log10(float(self.prev_val[0]/sig)),step)
		self.val[0] = sig * 10**proposal   #+ (np.random.random()-0.5)*5e-7
