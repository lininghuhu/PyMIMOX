#!/usr/bin/env python  
# coding=utf-8  
  
from __future__ import division
import operator
import sys, urllib2, zlib
import math
import re
from pymol import cmd, stored  
from smtplib import *  
from Tkinter import * 
import ttk as ttk 
import tkMessageBox  
import string 
import tkFont

	#  Install the plugin
def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
							'PyMIMOX',
							label='PyMIMOX',
							command=init) 

class App():
	# Plugin interface function

	def __init__(self,root):
		self.root = root 
		root.title('PyMIMOX')
		#root.geometry('400x460')
		ft = tkFont.Font(family = 'Fixdsys',size = 20,weight = tkFont.BOLD)
		
		frCuang = Frame(root, width=430, relief=GROOVE, bd=3)
		frCuang.grid( padx=10,  pady=10)
		
		
		# Introduction box
		
		self.mainlabel = Message(frCuang,fg = '#7B68EE', text='PyMIMOX is a conformational B-cell epitope prediction algorithm based on phage display. The algorithm divides the antigen surface into patches. The Generalised Jaccard Similarity between each patch and each mimotope is computed. And the top five matches are listed as candidate epitopes. The best one will be visualized in PyMOL by default.',width=420, bg='#fcfbf5', bd=2, padx=10,  pady=5,relief=GROOVE,  justify=LEFT)  
		self.mainlabel.grid(row=0, columnspan=12, padx=5,  pady=5)  

  
		self.pach = Label(frCuang, text='patch:',font = ft, borderwidth=2)  
		self.pach.grid(row=1,column=4, sticky = W) 

		
		self.pachipt= Text(frCuang, width=4, height=1, bg='#fcfbf5',relief=GROOVE)
		self.pachipt.insert(END,12)
		self.pachipt.grid(row=1,column=5, sticky = W)
		
		
		self.chain = Label(frCuang, text='chain:',font = ft)
		self.chain.grid(row=1,column=8)
		
		
		self.chainbox_value = StringVar(root)
		self.chainbox= ttk.Combobox(frCuang,textvariable=self.chainbox_value, state='readonly', width=3, height=1)
		self.chainbox['values'] = getChains()
		self.chainbox.grid(row=1,column=9, sticky = W)
	
		
		# self.kong = Label(frCuang, borderwidth=2)  
		# self.kong.grid(row=2,column=4, sticky = W) 
		
		
		
		self.Inputpeptide = Label(frCuang, text=' Peptide   ', borderwidth=2, bd=2, relief=GROOVE)  
		self.Inputpeptide.grid(row=3,column=0,rowspan=4, columnspan=2,padx=5,  pady=5) 
		
		
		self.iptBox = Text(frCuang, width=35, height=12,bg='#fcfbf5',bd=2, relief=GROOVE)
		self.iptBox.grid(row=3,column=3,columnspan=7, sticky = W)
		self.sb1 = Scrollbar(frCuang)
		self.sb1.grid(row = 3, column = 10, sticky = W+N+S)
		self.iptBox['yscrollcommand'] = self.sb1.set
		self.sb1['command'] = self.iptBox.yview 
		
		
		self.runButton = Button(frCuang, text='    Run    ', borderwidth=2, command=self.run)  
		self.runButton.grid(row=3, column=11,padx=5,  pady=5)
		
		
		
		self.Outputpeptide = Label(frCuang, text=' Prediction', borderwidth=2,bd=2, relief=GROOVE)  
		self.Outputpeptide.grid(row=8,column=0, columnspan=2,padx=5,  pady=5) 
		
		
		self.optBox = Text(frCuang, width=35, height=12,bg='#fcfbf5',bd=2, relief=GROOVE)
		self.optBox.grid(row=7,column=3,rowspan=4,columnspan=7, sticky = W)
		self.sb2 = Scrollbar(frCuang)
		self.sb2.grid(row = 7, column = 10,rowspan=4, sticky = W+N+S)
		self.optBox['yscrollcommand'] = self.sb2.set
		self.sb2['command'] = self.optBox.yview 
		
		self.clearButton = Button(frCuang, text='   Clear   ', borderwidth=2, command=self.clear)  
		self.clearButton.grid(row=8, column=11,padx=5,  pady=5)
		
		
		self.showButton = Button(frCuang, text='show', borderwidth=2, height=1,command=self.show_msg)  
		self.showButton.grid(row=9, column=11,sticky = W,  pady=2) 
	
		
		self.showbox_value = StringVar(frCuang)
		self.showbox= ttk.Combobox(frCuang,textvariable=self.showbox_value, state='readonly', width=1, height=1)
		self.showbox['values'] = (1, 2, 3, 4, 5)
		self.showbox.current(0)
		self.showbox.grid(row=9,column=11,sticky = E,padx=2,  pady=2)
		
	
	
	
	
	def show_msg(self):
		trans=int(self.showbox_value.get())
		allmsg = self.optBox.get('1.0', END).encode("utf-8").strip()
		if len(allmsg) != 0:
			Allmsg = allmsg.split('\n')
			msg = Allmsg[1::2]
			Msg = msg[trans-1][1:-1]
			MSG = Msg.split(', [')
			#print MSG
			name = 'predict_sql'
			resisql =''
			chain = self.getSelectedChain()
			for resi in MSG:
				resi =resi.split(',')[1]
				
				resisql+= " chain "+ chain
				resisql+=" and resi "+ str(resi[1:])
			cmd.show_as("line",name)
			cmd.select(name,resisql)
			cmd.show_as("spheres",resisql)
			cmd.color("yellow", resisql)
		
	
	
	# Get the user input patch length function
	def getSize(self):
		size = self.pachipt.get('1.0', END).encode("utf-8").strip()
		if len(size) == 0:  
			tkMessageBox.showwarning('Warning', 'The input can not be empty')  
			self.clear()
			return
		elif not re.search('\d+',size): 
			tkMessageBox.showwarning('Warning', 'The input peptide contains illegal characters')  
			self.clear()
			return
		else:
			return size 
	
	
		
		
	# 	Get Selected Chain
	def getSelectedChain(self):
		return self.chainbox_value.get()
		
		
	# Get the user input peptide function
	def getPeptide(self):
		peptide = self.iptBox.get('1.0', END).encode("utf-8").strip()
		match = re.sub('\s+(?!$)','',peptide)
		if len(peptide) == 0:  
			tkMessageBox.showwarning('Warning', 'The input peptide can not be empty')  
			self.clear()
			return
		elif re.search('B|J|O|U|X|Z|\*|\d|[a-z]',peptide): 
			tkMessageBox.showwarning('Warning', 'The input peptide contains illegal characters')  
			self.clear()
			return
		else:
			return peptide 
			
			
			
	# Clear the user's incorrect input information function		
	def clear(self):  
		self.iptBox.delete(1.0, END)  
		self.optBox.delete(1.0, END) 
		
		
		
	
	# Execute the main program	
	def run(self):
		
		peptide = self.getPeptide()  #get user input and check input
		distance = int(self.getSize())  # Get the user input patch length 
		cluster = Cluster(peptide)  # Peptide classification
		#deal pdb
		chain = self.getSelectedChain() # Get Selected Chain
		calist = resiCApick(chain) # Obtain possible residual information
		caInfo = caInfoFilter(calist) # Filter the residue information into [residue,index,chain, x,y,z]
		patchSet = makePatchSet(caInfo,distance) # Generation of patchs centered on residues
		predictScoreSet = {}  
		sequencemesg  = {}
		for i in range(0,int(len(cluster)/2)):
			pepset = AnalysisPeptide(cluster[i+1]) # Analytical peptides
			sequence = ''.join(cluster[i+1])
			for Key in patchSet.keys() :
				if type(Key[0])!=tuple:
					changedPatch = AnalysisPeptide(patchSet[Key])  # Analytical patch
				#print changedPatch
					curScore = similarity(pepset,changedPatch) # Comparison of similarity
					if Key not in predictScoreSet:
						predictScoreSet[Key] = curScore
						sequencemesg[Key] = cluster[sequence]
					#predictSet.append([Key,curScore])
					elif Key in predictScoreSet and curScore > predictScoreSet[Key]:
						predictScoreSet[Key] = curScore
						sequencemesg[Key] = cluster[sequence]
					#predictSet.append([Key,curScore])
		#print predictScoreSet
		predict = sorted(predictScoreSet.iteritems(),key = operator.itemgetter(1),reverse = True) # Sort by scoring
		
		# Get the predicted result
		
		finalResiSet = predict[0:10]
		#print finalResiSet
		
		
		
		
		#  Output forecast information
		
		Finalpreset = []
		for score in finalResiSet:
			mima = []
			mima1 = tuple(score[0])
			Sequenceinfo = sequencemesg[mima1]
			mima2 = tuple(patchSet[mima1])
			mima.append(mima1)
			mima.append(mima2)
			MIMA = tuple(mima)
			predict_resi = patchSet[MIMA]
			predict_resi[0] = list(predict_resi[0])
			predict_resi = sorted(predict_resi,key = lambda x:x[1])
			Finalpreset.append(predict_resi)
		#print Finalpreset
			
			
		indexmesg = delrepeat(Finalpreset)[0:5]
		#print indexmesg
		number = 1
		meg = {}
		for id in indexmesg:
			tempmima1 = tuple(finalResiSet[id][0])
			Sequenceinfo = sequencemesg[tempmima1]
			meg[number] = Finalpreset[id]
			
			self.optBox.insert(END,"Forecast_set")
			self.optBox.insert(END,number)
			self.optBox.insert(END,":")
			self.optBox.insert(END,Sequenceinfo)
			self.optBox.insert(END,'\n')
			self.optBox.insert(END,Finalpreset[id])
			self.optBox.insert(END,'\n')
			number+=1
			
			
		# Predictive information visualization	
		Num = int(self.showbox_value.get())
		name = 'predict_sql'
		resisql =''
		for resi in meg[Num]:
			resi = list(resi)
			chain = resi[2]
			resisql+= " chain "+ chain
			resisql+=" and resi "+ str(resi[1])			
		cmd.select(name,resisql)
		cmd.show_as("spheres",resisql)
		cmd.color("yellow", resisql)
		
		return meg
	

	

		
	def runApp(self, root):
		root.mainloop()	


		
		


			
	# Get the pdb file chain information	
def getChains():
		chains=[]
		for x in cmd.get_names():
			for ch in cmd.get_chains(x):
				chains.append(ch)
		return chains		
		
		
	# 	Peptide classification
def Cluster(peptide):
	cluster = {}
	Peptide = list(set(peptide.split('\n')))
	if len(Peptide)==1:
		if Peptide[0][0]=="C" and Peptide[0][-1]=="C":
			cluster[1] = list(Peptide[0][1:-1])
			cluster[Peptide[0][1:-1]] = Peptide[0]
		else:
			cluster[1]=list(Peptide[0])
			cluster[Peptide[0]] = Peptide[0]
	else:
		for i in range(0,len(Peptide)):
			if Peptide[i][0]=="C" and Peptide[i][-1]=="C":
				cluster[i+1] = list(Peptide[i][1:-1])
				cluster[Peptide[i][1:-1]] = Peptide[i]
			else:
				cluster[i+1] = list(Peptide[i])
				cluster[Peptide[i]] = Peptide[i]
	#print cluster
	return cluster
	
	
# 	Obtain  all possible residual information	
def resipick(selectedChain=False,cutoff=10, doShow=False, verbose=False):
    
	objSel="(all)"

    # Obtain pdb information
	tmpObj = "__tmp"
	cmd.create(tmpObj, objSel + " and polymer")

	fullObj = "full_str"
	cmd.create(fullObj, objSel + " and polymer")


	cmd.set("dot_solvent")
	cmd.get_area(selection=tmpObj, load_b=1)

	#print selectedChain
    # Remove unselected chains
	if selectedChain:
		cmd.remove(tmpObj + " and not chain "+ selectedChain)


    # threshold on what one considers an "exposed" atom (in A**2):
	cmd.remove(tmpObj + " and b < " + str(cutoff))


	stored.tmp_dict = {}
	cmd.iterate(tmpObj, "stored.tmp_dict[(chain,resv)]=1")
	exposed = stored.tmp_dict.keys()
	exposed.sort()

    # create sels
	selResi = "exposed_res"
	cmd.select(selResi, "byres " + objSel + " in " + tmpObj)

    # show exposed_resi's sels
	resiStr = cmd.get_pdbstr(selResi)

	if doShow != False:
		cmd.show_as("spheres", objSel + " and poly")
		cmd.color("yellow", selResi)

	cmd.delete(tmpObj)
	cmd.delete(fullObj)
	cmd.delete(selResi)
	return resiStr
	
	

# 	Obtain  all possible ca information rely on "resipick" function
def resiCApick(seletedChain=False):
    
	calist = []
	i = resipick(seletedChain)
	ca = "ca_left"
	test = "test"
	cmd.read_pdbstr(i, test)
	cmd.select(ca, "name ca " + " in " + test)
	calist = cmd.get_pdbstr('ca')
    
    
	cmd.delete(test)
	cmd.delete(ca)
	return calist
	
	



	
 # Filter the complete CA sequence information 
def caInfoFilter(caStr):
	
	# print caStr

	s = ''
	caStrFormated = []
	for i in caStr:
		if (i != '\n'):
			s += str(i)
		else:

            # keep ATOM and drop TER
            # print s+'\n'

			if (s[0:4] == "ATOM"):
				#print s[0:4]
				caStrFormated.append(s[17:55])
				s = ''
			else:
				s = ''


    # print len(caStrFormated)

	caInfoSet = []
	resiChainSet = {}

	for item in caStrFormated:
		index = int((item[6:9]))
		chain = item[4]

		if resiChainSet.has_key(chain):
			pass
		else:
			resiChainSet[chain]=[]


        
		if index not in resiChainSet[chain]:

			resiChainSet[chain].append(index)
			cur = []
			cur.append(condonChange(str(item[0:3])))  # residue
			cur.append(int((item[5:9])))  # index
			
			#chain
			cur.append(item[4])
			
            # position information
			#posStr = str(item[13:])
			# postion = posStr.split()
			# print postion
			cur.append(item[13:21])  # x
			cur.append(item[21:29])  # y
			cur.append(item[29:37])  # z
			#print cur
            

			caInfoSet.append(cur)

	return caInfoSet
	
	

# Make Patch Set
def makePatchSet(caInfo, distance):
	

	patch = {}
	
	for ca in caInfo:
		set = []
		preRsi = []

        # init patch
		a = tuple(ca[0:3])
		patch[a] = []
		patch[a].append(ca[0])
		flagPos = ca[3:6]
		set.append(a)
		preRsi.append(a)
		for compareCA in caInfo:
			if compareCA != ca:
				tmpPos = compareCA[3:6]
				dista = getDistance(flagPos, tmpPos)
                # campare
				if dista <= distance :
					
					rsi = compareCA[0]
					patch[a].append(rsi)
					set.append(compareCA[0:3])
		b = tuple(patch[a])			
		preRsi.append(b)			
		preRsi = tuple(preRsi)			
		patch[preRsi] = set
		#print set
		#print '\n'
	return patch
	
	

# Get Distance
def getDistance(pos1,pos2):
	distance = 0
	x = float(pos1[0])-float(pos2[0])
	y = float(pos1[1])-float(pos2[1])
	z = float(pos1[2])-float(pos2[2])
	distance = math.sqrt(x ** 2 + y ** 2 + z ** 2)
	return distance

		
# # Analytical peptides or patchs		
def AnalysisPeptide(peptide):
	AAs=list('ACDEFGHIKLMNPQRSTVWY')
	AAC = {}
	for AA in AAs:
		AAC[AA]=0
	for c in peptide:
		if c != '':
				#condon = condonNormalize(c.upper())
			if (AAC.has_key(c)):
				AAC[c] += 1
                # print peptideFuzzySet
			else:
				tkMessageBox.showwarning('Warning', 'The input peptide contains illegal characters')
                
	return AAC	


# Similarity		
def similarity(Pepset1,Pepset2):
	AAs=list('ACDEFGHIKLMNPQRSTVWY')
	simil = 0
	tmp1 =0
	tmp2 =0
	for h in AAs:
		aa=int(Pepset1[h])
		bb=int(Pepset2[h])
		tmp1+=min(aa,bb)
		tmp2+=max(aa,bb)
	simil = tmp1/tmp2

	return simil	

	
	# 	Change condon
def condonChange(resi):
	codon = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
			'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
			'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
			'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
			'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
			}

	return codon[resi]	


def delrepeat(List):
	news_ids1 = []
	news_ids2 = []
	i = 0
	for id in List: 
		if id not in news_ids1:
			news_ids1.append(id)
			news_ids2.append(i) 
			i +=1
		else:
			i +=1
	return news_ids2 	
		


def condonNormalize(condon):
	model ={'A': 'A', 'C': 'C', 'G': 'G', 'H': 'H', 'M': 'M', 'P': 'P', 'Y': 'Y',
			'R': 'B', 'K': 'B', 'B': 'B',
			'D': 'J', 'E': 'J', 'J': 'J',
			'S': 'O', 'T': 'O', 'O': 'O',
			'L': 'U', 'I': 'U', 'V': 'U', 'U': 'U',
			'Q': 'X', 'N': 'X', 'X': 'X',
			'F': 'Z', 'W': 'Z', 'Z': 'Z'
			}
	return model[condon]		
		
		
		
		
		
	# Initialization	
def init():
	root = Tk()
	app = App(root)
	app.runApp(root)
	
		
		
		
		
		
		
