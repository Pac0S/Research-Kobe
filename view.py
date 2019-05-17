import tkinter as tk
import proteine as pt
import numpy as np

class WorldView(tk.Canvas):
	
	def __init__(self, window, w=600, h=360, color="black"):
		tk.Canvas.__init__(self, window, width = w, height = h, bg = color) 
		self.window = window
		self.width = w
		self.height = h
		self.color = color
		self.spacelines = 0
		self.list_ovals =[]
		self.pause = True
		
	def draw_ovals(self, proteome):
		for i in range (proteome.shape[0]):
			for j in range (proteome.shape[1]):
				aa = proteome[i,j]
				if aa.rigid == 0 and aa.shearable == 0 : #Fluid non shearable : Red
					oval = self.create_oval(aa.column*20, aa.line*20, (aa.column+1)*20, (aa.line+1)*20, fill = "Red")
					#print(aa.column*10,aa.line*10, (aa.column+1)*10, (aa.line+1)*10, "Draw red")
					
				elif aa.rigid == 0 and aa.shearable == 1 : #Fluid shearable : Blue
					oval = self.create_oval(aa.column*20, aa.line*20, (aa.column+1)*20, (aa.line+1)*20, fill = "Blue")
					#print(aa.column*10,aa.line*10, (aa.column+1)*10, (aa.line+1)*10,"Draw Blue")
					
				elif aa.rigid == 1 and aa.shearable == 0 : #Rigid non shearable : Grey
					oval = self.create_oval(aa.column*20, aa.line*20, (aa.column+1)*20, (aa.line+1)*20, fill = "Grey")
					#print(aa.column*10,aa.line*10, (aa.column+1)*10, (aa.line+1)*10, "Draw grey")
					
				else : #Black, should not be possible
					oval = self.create_oval(aa.column*20, aa.line*20, (aa.column+1)*20, (aa.line+1)*20, fill = "Black")
					#print(aa.column*10,aa.line*10, (aa.column+1)*10, (aa.line+1)*10, "Draw black")
				
	
	
	
	
	
if __name__ == "__main__":
	
	mafenetre = tk.Tk()
	frame1 = tk.Frame(mafenetre)
	frame1.pack(side = 'top')
		
	mysquare = WorldView(frame1)
	mysquare.pack(side = "left")
	
	
	proteine = pt.Protein()
	
	rigid_input = np.random.random_integers(0,1,30)
	#rigid_input = [0,0,0,0,0]*6
	rigid_input = [1]*10+[0]*5+[1]*15
	#rigid_input = [1]*30
	proteine.set_input(rigid_input)
	proteine.update_prot()
	
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
			
	for i in range (proteine.w):
		print(proteine.proteome[i,1].sequence)
			
	#print(tabshear)
	#print(tabrig)
	
	#mysquare.create_oval(120, 60, 200, 120, fill = "Blue")
	mysquare.draw_ovals(proteine.proteome)
	mafenetre.mainloop()
	
	"""
	for aa in proteine.proteome:
		print(aa.column, aa.line)
		print(i)
		i+=1
	"""
