import tkinter as tk
import proteine as pt
import numpy as np
import time

class WorldView(tk.Canvas):
	
	def __init__(self, window, w=pt.Protein.w * 20, h=pt.Protein.h * 20, color="black"):
		tk.Canvas.__init__(self, window, width = w, height = h, bg = color) 
		self.window = window
		self.width = w
		self.height = h
		self.color = color
		self.spacelines = 0
		self.list_ovals = np.empty([pt.Protein.w, pt.Protein.h], dtype = object)
		self.pause = True
		
	def switch(self):
		if self.pause ==False:
			self.pause = True
		else:
			self.pause = False
		
	def draw_ovals(self, proteome):
		for i in range (proteome.shape[0]):
			for j in range (proteome.shape[1]):
				aa = proteome[i,j]
				if aa.rigid == 0 and aa.shearable == 0 : #Fluid non shearable : Red
					oval = self.create_oval(aa.column*20, aa.line*20, (aa.column+1)*20, (aa.line+1)*20, fill = "Red")
					#print(aa.column*10,aa.line*10, (aa.column+1)*10, (aa.line+1)*10, "Draw red")
					self.list_ovals[i,j] = oval
					
				elif aa.rigid == 0 and aa.shearable == 1 : #Fluid shearable : Blue
					oval = self.create_oval(aa.column*20, aa.line*20, (aa.column+1)*20, (aa.line+1)*20, fill = "Blue")
					#print(aa.column*10,aa.line*10, (aa.column+1)*10, (aa.line+1)*10,"Draw Blue")
					self.list_ovals[i,j] = oval
					
				elif aa.rigid == 1 and aa.shearable == 0 : #Rigid non shearable : Grey
					oval = self.create_oval(aa.column*20, aa.line*20, (aa.column+1)*20, (aa.line+1)*20, fill = "Yellow")
					#print(aa.column*10,aa.line*10, (aa.column+1)*10, (aa.line+1)*10, "Draw grey")
					self.list_ovals[i,j] = oval
					
				else : #Black, should not be possible
					oval = self.create_oval(aa.column*20, aa.line*20, (aa.column+1)*20, (aa.line+1)*20, fill = "Black")
					#print(aa.column*10,aa.line*10, (aa.column+1)*10, (aa.line+1)*10, "Draw black")
					self.list_ovals[i,j] = oval
	
	def update_ovals(self,proteome): #updates the number of agents in each state too
		for i in range (proteine.w):
			for j in range (proteine.h):
				aa = proteome[i,j]
				if aa.rigid == 0 and aa.shearable == 0 : #Fluid non shearable : Red
					self.itemconfig(self.list_ovals[i,j],fill = "Red")
					#print("Red")

				elif aa.rigid == 0 and aa.shearable == 1 : #Fluid shearable : Blue
					self.itemconfig(self.list_ovals[i,j],fill = "Blue")
					#print("Blue")
				elif aa.rigid == 1 and aa.shearable == 0 : #Rigid non shearable : Grey
					self.itemconfig(self.list_ovals[i,j],fill = "Yellow")
					#print("Grey")
				else : #Black, should not be possible
					self.itemconfig(self.list_ovals[i,j],fill = "Black")
	
	
	
if __name__ == "__main__":

	"""
	def run_my_prot():
		proteine.run()
		mysquare.update_ovals(proteine.proteome)
	"""
	def run():
		start_time = time.time()
		global proteine
		global mysquare
		if mysquare.pause == False:
			proteine.run_once()
			mysquare.update_ovals(proteine.proteome)
			
			space = (1/1000)*1000
			mafenetre.after(1,run)
		execution_time = time.time()-start_time
		#print(execution_time)
		
		
	def switch_play():
		mysquare.switch()
		if mysquare.pause ==False:
			run()
			
			
	def draw_aa():
		global proteine, mysquare
		proteine = pt.Protein()
		rest = pt.Protein.w - pt.Protein.w//3 -5
		rigid_input = [1]*(pt.Protein.w//3)+[0]*5+[1]*rest
		proteine.set_input(rigid_input)
		proteine.update_prot()	
		mysquare.draw_ovals(proteine.proteome)
		
		
	
	####Definition of Canvas####
	mafenetre = tk.Tk()
	#Frame dessin agents
	frame3 = tk.Frame(mafenetre)
	frame3.pack(side = 'top')
	#Frame de l'affichage
	frame1 = tk.Frame(mafenetre)
	frame1.pack(side = 'top')
	#Frame du bouton play/pause
	frame2 = tk.Frame(mafenetre)
	frame2.pack(side = 'top')
	
	
		
	mysquare = WorldView(frame1)
	mysquare.pack(side = "left")
	
	boutonpause = tk.Button(frame2,text = "Play/Pause", command = switch_play)
	boutonpause.pack(side = "left")
	
	bt_aas = tk.Button(frame3,text = "Create AAs", command = draw_aa)
	bt_aas.pack(side = "left")
	
	"""
	#####Initialization of protein with input and prescribes output####
	proteine = pt.Protein()
	
	rigid_input = [1]*10+[0]*5+[1]*15
	output = [1]*10+[0]*5+[1]*15
	
	
	proteine.set_input(rigid_input)
	proteine.update_prot()
	
	#proteine.mut_prot()
	#proteine.update_prot()
	"""
	
	
	"""
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
	"""		
			
	#print(tabshear)
	#print(tabrig)
	
	
	
	
	####Display of Canvas####
	"""
	mysquare.draw_ovals(proteine.proteome)
	run()
	"""
	mafenetre.mainloop()
	
	

