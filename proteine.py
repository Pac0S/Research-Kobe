import numpy as np
import random
import acide as ac
import time
from copy import deepcopy
from copy import copy

class Protein(object):
	#Variables globales
	h = 9 #Columns
	w = 15 #Lines
	output = [[1]*(w//3)+[0]*5+[1]*(w - w//3 -5),[0]*(w//3)+[1]*5+[0]*(w - w//3 -5)] #[Rigid sequence, Shearable sequence]

	def __init__(self, genome = np.random.randint(2, size=(15*ac.Acide.nb_links, 9)) , proteome = np.empty([15, 9], dtype = object)):
		self.genome = genome
		self.proteome = proteome
		self.mutations = [0,0,0]#Mutations [Positives, egales, deletaires]
		self.fitness =  0
		
		for i in range (Protein.w) :
			for j in range (Protein.h) :
			
				aminoacid = ac.Acide(genome[5 * i : 5 * i + 5, j])
		
				aminoacid.line = j #La liste commence en bas a gauche : de 0 a 17
				aminoacid.column = i #de 0 a 29
				
				self.proteome[i,j] = aminoacid
				
				#print(i, j)
				#print("Nouvel acide")
				#print(aminoacid.sequence)
		
					
	#states est un tableau de 1 et 0 qui donne la rigidite des aa en premiere ligne. sig = 0 entraine s = 0
	def set_input(self, rigid_list):
		for i in range(Protein.w) : #La taille de la liste doit etre egale a w = 30
			self.proteome[i,0].rigid = rigid_list[i]
			if rigid_list[i] == 0:
				self.proteome[i,0].shearable = 1 #Tous les aa fluides en premiere ligne sont shearables
			if rigid_list[i] == 1 :
				self.proteome[i,0].shearable = 0#Il est interdit d'avoir un aa rigide et shearable
			
	
	def update_fitness(self):
		self.fitness = 0
		for i in range (Protein.w):
			aa = self.proteome[i, Protein.h-1]
			if aa.shearable == Protein.output[1][i] and aa.rigid == Protein.output[0][i]:
				self.fitness += 1

			
	#Mise a jour de la proteine en fonction des aa voisins  (au dela d'une certaine ligne pour eviter les operations inutiles apres mutations)
	def update_prot(self):
		for j in range (1, Protein.h) :
			for i in range (Protein.w) :
				rigid = 0
				shearable = 0
				
				
				link = 0	
				#On compte le nombre de voisins (sur 5) rigides lies a chaque acide amine 
				for neighbor in range (i - 2, i + 3):
					if neighbor >= Protein.w-1 :
						neighbor2 = neighbor - Protein.w - 1
					else :
						neighbor2 =neighbor
					#print(neighbor2)
					#print(neighbor2)
					#print(link)
					if self.proteome[i,j].sequence[link]!=0 :
						if self.proteome[neighbor2, j-1].rigid == 1 :
							rigid += 1
							#print("ligne = ", self.proteome[neighbor2, j-1].line)
							#print(j)
					link += 1
					#print("")
				#print("stop")



				link = 0
				#On compte le nombre (sur 3) de voisins shearables lies a chaque acide amine
				for neighbor in range (i - 1, i + 2):
					if neighbor >= Protein.w - 1 :
						neighbor2 = neighbor - Protein.w - 1
						
					else :
						neighbor2 =neighbor
						
					if self.proteome[i,j].sequence[link+1]!=0 :
						if self.proteome[neighbor2, j-1].shearable == 1 :	
							shearable += 1
					link += 1
				
	
				if  rigid >=2 :
					#print("abc")
					self.proteome[i,j].rigid = 1
				else :
					self.proteome[i,j].rigid = 0
					
					
					
				if shearable >=1 and self.proteome[i,j].rigid == 0:
					self.proteome[i,j].shearable = 1
					#if self.proteome[i, j].line == 1:
						#print(i,j)
				else :
					self.proteome[i,j].shearable = 0
					
				if self.proteome[i,j].rigid == 1 :
					self.proteome[i,j].shearable = 0 #Interdit d'avoir rigid et shearable
		self.update_fitness()
					

			
		
	
	
	#Modifier pour n'accepter que les mutations favorables ou neutres				
	def mut_prot(self):
		#Sauvegarde proteine en cas de mutation deletaire
		prot_cop = copy(self)
						
		#Mutation genome et mise a jour sequence acide amine
		line = np.random.random_integers(self.h-1)#On choisit une ligne et une colonne a muter aleatoirement dans le genome
		column = np.random.random_integers(self.w*ac.Acide.nb_links-1) #Attention les colonnes du genome sont 5* plus nombreuses que celles du  proteome
		if self.genome[column,line]==0 :
			self.genome[column,line]==1
		else :
			self.genome[column,line]==0
		column_prot = column//ac.Acide.nb_links
		index = column%ac.Acide.nb_links
		
		"""
		print("Genome column : ", column)
		print("Genome line : ", line)
		print("Proteome column : ", column_prot)
		print("Index of mutation : " ,index)
		"""
		
		#Mise a jour de l'acide amine
		self.proteome[column_prot, line].mutation(index)
		#Mise a jour des proprietes de l'aa
		self.update_prot()
		
		#print(self.fitness, prot_cop.fitness)
		#On refuse les mutations délétaires
		if self.fitness < prot_cop.fitness :
			if self.genome[column,line]==0 :
				self.genome[column,line]==1
			else :
				self.genome[column,line]==0
			#Mise a jour de l'acide amine
			self.proteome[column_prot, line].mutation(index)
			#Mise a jour des proprietes de l'aa
			self.update_prot()
			
			self.mutations[2]+=1
			#print("mutation -\n")
			
		elif self.fitness > prot_cop.fitness:
			self.mutations[0]+=1
			print (self.fitness, prot_cop.fitness)
			#print("mutation +\n")	
			
		else :
			self.mutations[1]+=1
			#print("mutation =\n")
		
		



	def run_once(self) :
		self.mut_prot()
		self.update_prot()
		print(self.mutations)
		

		
if __name__ == "__main__":


			
	#Creation d'une proteine
	
	proteine = Protein()
	#print(proteine.proteome)
	
	#Test methode set_input
	
	rigid_input = np.random.random_integers(0,1,Protein.w)
	proteine.set_input(rigid_input)
	
	"""
	print(len(proteine.genome))#2700
	print(Protein.h*Protein.w)
	print(len(proteine.proteome)) #540
	"""
	
	print("")

	#Mise a jour de la proteine selon l'input
	proteine.update_prot()
	
	
	#Visualisation des caracteristiques shearable et rigid
	#La ligne 0 de la proteine est a gauche
	tabshear = np.zeros((proteine.w, proteine.h))
	for i in range (proteine.w):
		for j in range (proteine.h):
			print(proteine.proteome[i][j].shearable)
			tabshear[i][j] = proteine.proteome[i][j].shearable
			
	tabrig = np.zeros((proteine.w, proteine.h))
	for i in range (proteine.w):
		for j in range (proteine.h):
			print(proteine.proteome[i][j].shearable)
			tabrig[i][j] = proteine.proteome[i][j].rigid
			
	print(tabshear)
	print(tabrig)
	
	
	somme = 0
	for i in range(proteine.w):
		somme += sum(aci.shearable==1 for aci in proteine.proteome[i])
	print (somme) #Nombre d'aa shearables (15)
	
	somme = 0
	for i in range(proteine.w):
		somme += sum(aci.rigid==1 for aci in proteine.proteome[i])
	print(somme) #Nombre d'aa rigides (50)
	

	

