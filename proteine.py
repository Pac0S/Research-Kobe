import numpy as np
import random
import acide as ac
import time
from copy import copy
from os import system
import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


class Protein(object):
	#Variables globales
	h = 18 #Columns
	w = 30 #Lines
	output = [[1]*12+[0]*5+[1]*13,[0]*12+[1]*5+[0]*13] #[Rigid sequence, Shearable sequence]

	#def __init__(self, genome = np.random.randint(2, size=(30*ac.Acide.nb_links, 18)) , proteome = np.empty([30, 18], dtype = object)):
	def __init__(self, genome = np.random.choice(2, size=(30*ac.Acide.nb_links, 18), p = [0.3,0.7]) , proteome = np.empty([30, 18], dtype = object)):
		self.genome = genome
		self.proteome = proteome
		self.mutations = [0,0,0]#Mutations [Positives, egales, deletaires]
		self.mut_plus = []
		self.mut_plus_count = 0
		self.mut_less = []
		self.mut_less_count = 0
		self.mut_none = []
		self.mut_none_count = 0
		self.timeline = []
		self.time = 0
		self.fitness = 0
		
		
		for i in range (Protein.w) :
			for j in range (Protein.h) :		
				aminoacid = ac.Acide(genome[5 * i : 5 * i + 5, j])			
				self.proteome[i,j] = aminoacid
				aminoacid.line = j
				aminoacid.column = i
	
			
			
			
	def get_genome(self):
		list_genome = []
		for j in range (1, Protein.h):
			for i in range(Protein.w*ac.Acide.nb_links) :
				list_genome.append(self.genome[i,j])
		return list_genome
		
	def get_shearable_list(self):
		list_shearable = []
		for j in range (1, Protein.h):
			for i in range(Protein.w) :
				list_shearable.append(self.proteome[i,j].shearable)
		return list_shearable
		
		
				
	#rigid_list est un tableau de 1 et 0 qui donne la rigidite des aa en premiere ligne. sig = 0 entraine s = 1
	def set_input(self, rigid_list):
		for i in range(Protein.w) : #La taille de la liste doit etre egale a w 
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
				
	def defective_ac(self,column,line) :
		for i in range(5) :
			self.genome[column*5 + i , line] = 0 #On modifie aussi le genome
		self.proteome[column,line].defective()
		

			
	#Mise a jour de la proteine en fonction des aa voisins (au dela d'une certaine ligne pour eviter les operations inutiles apres mutations)
	def update_prot(self):
		
		
		for j in range (1, Protein.h) :
			for i in range (Protein.w) :
				#Debugger
				#if i == 10 and j ==10 :
				#	import ipdb; ipdb.set_trace()
				rigid = 0 #Nombre de voisins connectés rigides (<=5)
				shearable = 0 #Nombre de voisins connectés shearables (<=3)
				
				
				link = 0 #Indice du lien avec le voisin (0-4)
				#On compte le nombre de voisins (sur 5) rigides lies a chaque acide amine 
				for neighbor in range (i - 2, i + 3):
					if neighbor >= Protein.w-1 :
						neighbor2 = neighbor - Protein.w - 1
						#Certains voisins du bord droit sont au bord gauche. Pas besoin pour bord gauche car indice negatif = fin de liste
					elif neighbor < 0 :
						neighbor2 = Protein.w + neighbor
					else :
						neighbor2 =neighbor
					if self.proteome[i,j].sequence[link]!=0 :
						if self.proteome[neighbor2, j-1].rigid == 1 and not self.proteome[neighbor2, j-1].is_defective :
							rigid += 1
					link += 1


				#link = 0
				#On compte le nombre (sur 3) de voisins shearables lies a chaque acide amine
				for neighbor in range (i - 1, i + 2):
					if neighbor >= Protein.w - 1 :
						neighbor2 = neighbor - Protein.w - 1
					elif neighbor < 0 :
						neighbor2 = Protein.w + neighbor
					else :
						neighbor2 =neighbor				
					#if self.proteome[i,j].sequence[link+1]!=0 :
					if self.proteome[neighbor2, j-1].shearable == 1 and not self.proteome[neighbor2, j-1].is_defective : #A defective acid will be unable to transfer its shearability and rigidity
						shearable += 1
					#link += 1
				
				#Mise a jour des proprietes de l'acide
				if rigid >=2 :
					self.proteome[i,j].rigid = 1
				else :
					self.proteome[i,j].rigid = 0
					
					
					
				if shearable >=1 and self.proteome[i,j].rigid == 0:
					self.proteome[i,j].shearable = 1
				else :
					self.proteome[i,j].shearable = 0
					
				if self.proteome[i,j].rigid == 1 :
					self.proteome[i,j].shearable = 0 #Interdit d'avoir rigid et shearable
		self.update_fitness()
		
		"""
		#Write in file the values of rigid and sheable for the aa of the 11 first lines
		myfile = open("data.txt", "w+")
		for j in range(10):
			myfile.write("\n\n\n Line " + str(j) + "\n\n")
			for i in range(Protein.w):
				myfile.write(str(i) + " : " + str(self.proteome[i,j].sequence)+ "\n Rigid : " + str(self.proteome[i,j].rigid) + " ||| Shearable : " + str(self.proteome[i,j].shearable) + "\n\n")
				#print(i, self.proteome[i,1].sequence)
		"""




	
	#Modifier pour n'accepter que les mutations favorables ou neutres				
	def mut_prot(self):
		#Sauvegarde proteine en cas de mutation deletaire
		prot_cop = copy(self)
						
		#Mutation genome et mise a jour sequence acide amine
		line = np.random.random_integers(self.h-1)#On choisit une ligne et une colonne a muter aleatoirement dans le genome
		column = np.random.random_integers(self.w*ac.Acide.nb_links-1) #Attention les colonnes du genome sont 5* plus nombreuses que celles du  proteome
		
		column_prot = column//ac.Acide.nb_links
		index = column%ac.Acide.nb_links
		
		if not self.proteome[column_prot, line].is_defective :
			if self.genome[column,line]==0 :
				self.genome[column,line]==1
			else :
				self.genome[column,line]==0

			"""
			print("Genome column : ", column)
			print("Genome line : ", line)
			print("Proteome column : ", column_prot)
			print("Index of mutation : " ,index)
			"""
			
			#Mise a jour de l'acide amine correspondant
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
			self.mut_less_count+=1
			self.mut_plus.append(self.mut_plus_count)
			self.mut_less.append(self.mut_less_count)
			self.mut_none.append(self.mut_none_count)
			#print("mutation -\n")
			
		elif self.fitness > prot_cop.fitness:
			self.mutations[0]+=1
			self.mut_plus_count+=1
			self.mut_plus.append(self.mut_plus_count)
			self.mut_less.append(self.mut_less_count)
			self.mut_none.append(self.mut_none_count)
			print(self.mutations)
			print('\007')
			#print (self.fitness, prot_cop.fitness)
			#print("mutation +\n")	
			
		else :
			self.mutations[1]+=1
			self.mut_none_count+=1
			self.mut_plus.append(self.mut_plus_count)
			self.mut_less.append(self.mut_less_count)
			self.mut_none.append(self.mut_none_count)
			#print("mutation =\n")
		self.time+=1
		self.timeline.append(self.time)
		
		
		
	def run_once(self):
		while(proteine.fitness < proteine.w):
			proteine.mut_prot()
			
		
	#Ecriture d'un nombre donné de solutions de sequence dans un fichier
	def gen_sols(self, number):
		i = 0
		self.run_once()
		print("Solution found")
		found = 0
		while(i<number):
			solutions = open('solutions_gen_S.txt','a+')
			reader = open('solutions_gen_S.txt','r')
			string_gen = str(self.get_genome())
			linelist = reader.readlines()
			line_written = False
			for line in linelist :
				if string_gen in line:
					found += 1
					line_written = True
					break
			if not line_written :
				solutions.write(string_gen +"\n")
				i+=1
				print(i)
			proteine.mut_prot()
			reader.close()
			solutions.close()
		print(found)
				
	#Ecriture d'un nombre donné de solutions de conformation dans un fichier
	
	def shear_sols(self, number):
		i = 0
		self.run_once()
		print("Solution found")
		found = 0
		while(i<number):
			solutions = open('obstacle_shear_S.txt','a+')
			reader = open('obstacle_shear_S.txt','r')
			string_shear = str(self.get_shearable_list())
			linelist = reader.readlines()
			line_written = False
			for line in linelist :
				if string_shear in line:
					found += 1
					line_written = True
					break
			if not line_written :
				solutions.write(string_shear +"\n")
				i+=1
				print(i)
			proteine.mut_prot()
			reader.close()
			solutions.close()
		print(found)
		
		
		
	def svd_gen(self, eigen):
		shear_list = []
		#Stockage des solutions dans un np array
		with open('solutions_gen_S.txt','r') as sols :
			l = 0
			for lines in sols :
				l+=1
				for i in range(Protein.w*5*(Protein.h-1)):
					value = lines[3*i+1]
					shear_list.append(int(value))
		svd_np = np.array(shear_list).reshape(l, Protein.w* 5 * (Protein.h-1));
		svd_np = np.unique(svd_np, axis = 0)
		print("Size of matrix for svd : " + str(svd_np.shape))
		u, s, v = np.linalg.svd(svd_np, full_matrices=True)
		print("Left singular vectors : " + str(u.shape) + "\nMatrix sigma : " + str(s.shape) + "\nRight singular vectors : " + str(v.shape))
		s_vec = v[eigen,:]
		s_vec = s_vec.reshape(Protein.h-1,Protein.w*ac.Acide.nb_links)
		
		fig, ax = plt.subplots()
		im = ax.imshow(s_vec, cmap = "bwr")
		fig.tight_layout()
		plt.show()
		
		
		print("First singular vectors : " + str(s[0:7]))
		s_round = np.round(s, 0)
		s_unique, freq = np.unique(s_round, return_counts = True)
		print(s_unique)
		plt.yscale('log')
		
		plt.bar(s_unique, height = freq, color= "red")
		plt.xlim(0,2000)
		plt.show()
		
			
		
	def svd_shear(self, eigen):
		shear_list = []
		#Stockage des solutions dans un np array
		with open('obstacle_shear_S.txt','r') as sols :
			l = 0
			for lines in sols :
				l+=1
				for i in range(Protein.w*(Protein.h-1)):
					value = lines[3*i+1]
					shear_list.append(int(value))
		svd_np = np.array(shear_list).reshape(l, Protein.w * (Protein.h-1));
		svd_np = np.unique(svd_np, axis = 0)
		print("Size of matrix for svd : " + str(svd_np.shape))
		u, s, v = np.linalg.svd(svd_np, full_matrices=True)
		print("Left singular vectors : " + str(u.shape) + "\nMatrix sigma : " + str(s.shape) + "\nRight singular vectors : " + str(v.shape))
		s_vec = v[eigen,:]
		s_vec = s_vec.reshape(17,30)
		
		fig, ax = plt.subplots()
		im = ax.imshow(s_vec, cmap = "bwr")
		fig.tight_layout()
		plt.show()
		
		
		print("First singular vectors : " + str(s[0:7]))
		s_round = np.round(s, 0)
		s_unique, freq = np.unique(s_round, return_counts = True)
		print(s_unique)
		plt.yscale('log')
		
		plt.bar(s_unique, height = freq, color= "red")
		plt.xlim(0,220)
		plt.show()
		
		
		


		
if __name__ == "__main__":


			
	#Creation d'une proteine
	
	proteine = Protein()
	proteine.defective_ac(12,9)
	
	#Test methode set_input
	
	#rest = proteine.w - proteine.w//3 -5
	#rigid_input = [1]*(proteine.w//3)+[0]*5+[1]*rest
	rigid_input = [1]*12 + [0]*5 + [1]*13
	proteine.set_input(rigid_input)
	
	
	#Mise a jour de la proteine selon l'input
	proteine.update_prot()
	
	
	
	#Temps d'exécution
	
	t0 = time.time()
	
	
	#proteine.gen_sols(48551)
	proteine.shear_sols(6000)
	
	tt = time.time()-t0
	h = tt//3600
	rh = tt%3600
	m = rh//60
	s = rh%60
	print(str(h)+"h " + str(m) + "m " + (str(m) + "s"))
	
	
	
	t0 = time.time()
	
	proteine.svd_shear(0)
	proteine.svd_shear(1)
	proteine.svd_shear(2)
	#proteine.svd_gen(0)
	
	
	tt = time.time()-t0
	h = tt//3600
	rh = tt%3600
	m = rh//60
	s = rh%60
	print(str(h)+"h " + str(m) + "m " + (str(m) + "s"))
	
	
	#Affichage d'un vecteurs singulier choisi (heatmap)
	
	
	
	
	
	#Affichage evolution du nombre de mutations délétaires, benefiques ou neutres
	"""
	proteine.run_once()
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(1, proteine.time)
	#plt
	plt.scatter(proteine.timeline, proteine.mut_plus, s=1, c='g', marker='o', label = 'beneficial')
	plt.scatter(proteine.timeline, proteine.mut_less, s=1, c='r', marker='o', label = 'deleterious')
	plt.scatter(proteine.timeline, proteine.mut_none, s=1, c='black', marker='o', label = 'neutral')
	
	
	plt.legend(loc="upper left")
	plt.show()
	"""
	
	
	#Visualisation des caracteristiques shearable et rigid
	"""
	#La ligne 0 de la proteine est a gauche
	tabshear = np.zeros((proteine.w, proteine.h))
	for i in range (proteine.w):
		for j in range (proteine.h):
			#print(proteine.proteome[i][j].shearable)
			tabshear[i][j] = proteine.proteome[i][j].shearable
			
	tabrig = np.zeros((proteine.w, proteine.h))
	for i in range (proteine.w):
		for j in range (proteine.h):
			#print(proteine.proteome[i][j].shearable)
			tabrig[i][j] = proteine.proteome[i][j].rigid
			
	#print(tabshear)
	#print(tabrig)
	
	
	somme = 0
	for i in range(proteine.w):
		somme += sum(aci.shearable==1 for aci in proteine.proteome[i])
	#print (somme) #Nombre d'aa shearables (15)
	
	somme = 0
	for i in range(proteine.w):
		somme += sum(aci.rigid==1 for aci in proteine.proteome[i])
	#print(somme) #Nombre d'aa rigides (50)
	"""

	

