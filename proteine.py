import numpy as np
import random
import acide as ac

class Protein(object):
	#Variables globales
	h = 18 #Columns
	w = 30 #Lines
	
	def __init__(self, genome = np.random.random_integers(0,1,h*w*ac.Acide.nb_links), proteome = []):
		self.genome = genome
		self.proteome = proteome
		for i in range(Protein.h*Protein.w):
			aminoacid = ac.Acide(genome[5 * i : 5 * i + 5])
			aminoacid.line = i // Protein.h #La liste commence en bas a gauche
			aminoacid.column = i // Protein.w
			self.proteome.append(aminoacid)
			
	#states est un tableau de 1 et 0 qui donne la rigidite des aa en premiere ligne. sig = 0 entraine s = 0
	def set_input(self, rigid_list):
		for i in range (Protein.w) :
			self.proteome[i].rigid = rigid_list[i]
			if rigid_list[i] == 0:
				self.proteome[i].shearable = 1
			
	#Mutation aleatoire dans le genome --> mise a jour sequence acide amine
	#def mutation
			
	#Mise a jour de la proteine en fonction des aa voisins  (au dela d'une certaine ligne pour eviter les operations inutiles apres mutations)
	
	#A modifier : Faire en sorte que les indices des voisins soient bien sur la ligne precedente sur les bords
	def update_prot(self,indice):
		for aa in range (indice, Protein.h*Protein.w):
			#print(aa)
			rigid = 0
			shearable = 0
			for neighbor in range (aa - Protein.w - 2, aa - Protein.w + 3):
				i=0
				if self.proteome[aa].sequence[i] != 0 :
					if self.proteome[neighbor].rigid == 1:
						#print("abc")
						rigid += 1
						#print(rigid)
					
					i+=1
				print("rigid : " + str(rigid))
					
			for neighbor in range (aa - Protein.w - 1, aa - Protein.w + 2):
				i=0
				if self.proteome[neighbor].shearable == 1:
					shearable += 1
				i+=1
				print("shearable : " + str(shearable))
				print("")
				
				
				
			if  rigid >=2 :
				#print("abc")
				self.proteome[aa].rigid = 1
			if shearable >=1 :
				self.proteome[aa].shearable = 1
			print(self.proteome[aa].rigid)
			
				

	
	
if __name__ == "__main__":


			
			
	proteine = Protein()
	print(len(proteine.genome))#2700
	print(Protein.h*Protein.w)
	print((proteine.proteome[530].sequence))
	print(len(proteine.proteome)) #540
	
	
	#Test methode set_input
	rigid_input = np.random.random_integers(0,1,30)
	#rigid_input = [1,1,1,1,0]*6
	#print(rigid_input)

	
	proteine.set_input(rigid_input)
	"""
	for i in range (Protein.w) :
		print(proteine.proteome[i].rigid)
	"""
	"""
	for aa in proteine.proteome:
		#print(aa.line)
		print(aa.rigid)
	"""
	"""
	for i in range (0,60):
		print(proteine.proteome[i].shearable)
	"""
	
	print("")
	#Test update
	proteine.update_prot(32)
	"""
	for i in range (0,90):
		print(proteine.proteome[i].rigid)
	"""
	
