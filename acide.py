import numpy as np
import random as rd
#import proteine as pt

class Acide(object):
	#Variable globale
	nb_links = 5
	
	#Constructeur
	def __init__(self, sequence = np.random.random_integers(0,1,nb_links)):
		self.sequence  = sequence
		#self.rigid = 0
		#self.shearable = 0
		self.rigid = np.random.randint(2)
		if (self.rigid == 1):
			self.shearable = 0
		else :
			self.shearable = np.random.randint(2)
		self.line = 0
		self.column = 0
		self.is_defective = False
		
	def defective(self) :
		self.sequence = [0,0,0,0,0]
		self.is_defective = True
	
	
	def mutation(self, link):
		if self.sequence[link]==0:
			self.sequence[link]=1
		else:
			self.sequence[link]=0
			
	def _shearable(self) :
		s = 0
		if self.shearable == 1 :
			s = 1
		return s
		
	def _rigid(self) :
		s = 0
		if self.rigid == 1 :
			s = 1
		return s
		
		
		
	"""	
	Propagation de la rigidite : 
	Input : r = 0 (row 0)
	pour tous les autres r :
	Propagation de la rigidite sigma :
	
	sig(r,c) = Theta(x) = Theta (sum{k=-2:2} l(rck) sig(r-1, c+k)- sig0)
	
	ou sig0=2 est le nombre minimal de voisins rigides pour que la rigidite soit transmise.
	Et Theta = 1 si x>=0 et 0 si x<0
	l(rck) donne la connectivite d'un aa au soisin de la ligne inferieure (1 ou 0)
	En d'autres termes, avec 2 voisins ou plus dont sig>0 on a sig(r,c) = 1 sinon 0
	
	
	
	Propagation de la "shearablite" :
	Tous les aa fluides de la ligne 0 sont aussi shearables
	Un aa sera shearable si au moins un de ses 3 plus proches voisins (r-1, c) ou (r-1, c+-1) est shearable
	Il ne peut pas etre shearable si il est rigide
	
	s(r,c) = (1-sig(r,c)) Theta sum{k=-1:1}(s(r-1, c+k) - s0)
	
	Avec s0 = 1
	
	"""


	
		

if __name__ == "__main__":
	aminoacid = Acide()
	print(aminoacid.sequence[2])

