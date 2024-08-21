
import random as rnd
import numpy as np
import copy

class parameter_CI :
# reading/writing of the input file. 
#name -> string : name of the variable
# minim, max -> float: range definitions
# width -> float : width in which we will take another value from the first one. taken between 0 and 1.
	def __init__(self, Name, minim, maxi, Width):
		self.name = Name
		self.min = float(minim)
		self.max = float(maxi)
		self.width = float(Width)
		self.val = np.array([-1e15])
		self.prev_val = np.array([-1e15])

	def set_val(self,value):
		self.val[0] = value

	def initialise(self):
		self.set_val(rnd.uniform(self.min,self.max))
		self.prev_val[0] = self.val[0]

	def accepted(self):
		self.prev_val = copy.copy(self.val)


	def jump(self):
		step = (self.max - self.min) * self.width
		proposal = rnd.gauss(self.prev_val[0], step )
		while proposal < self.min or proposal > self.max :	#ensure that gaussian value is within required range            	
			proposal = rnd.gauss(float(self.prev_val[0]),step)
		self.val[0] = proposal
