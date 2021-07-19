import pickle
import math
import copy
import time
from dircycle2 import dircyc,almostdircyc
from heuristic2 import contract_heur,contract_heur_bdir
import cplex

class BNSLlvInst:
	def __init__(self,instance):
		self.instance = instance
		self.originalScores = None
		self.scores = None
		self.data = None
		self.V = None
		self.cComps = None
		self.iComps = None
		self.dPars = None
		self.iPars = None
		self.biPars = None
		self.m = None
		self.z = None
		self.ind = None
		self.indInv = None
		self.udE = None
		self.bi = None
		self.x = None
		self.clusterIP = None
		self.ConflictNodes = None
		self.ConflictEdges = None

	def readFromPkl(self):
		filename = '../Instances/data/'+self.instance
		file = open(filename, 'rb')
		[self.data,self.originalScores] = pickle.load(file)
		file.close()
		
	def readFromDag(self):
		filename = '../Instances/data/'+self.instance
		file = open(filename, 'r')
		x = file.readlines()
		file.close
		self.originalScores = {}
		Vsize = int(x[0].strip('\n'))
		lineind = 1
		for k in range(Vsize):
			line = x[lineind]
			z = line.strip('\n').split(' ')
			i = int(z[0])
			nPaSets = int(z[1])
			self.originalScores[(i,),()] = {}
			for p in range(nPaSets):
				line = x[lineind+1+p]
				z = line.strip('\n').split(' ')
				PaSize = int(z[1])
				PaSet = []
				for t in range(PaSize):
					PaSet.append(int(z[2+t]))
				self.originalScores[(i,),()][(tuple(PaSet),)] = float(z[0])
			lineind = lineind+1+nPaSets
		if len(x) > lineind:
			lineind = lineind+1
			csize = int(x[lineind].strip('\n'))
			for k in range(csize):
				lineind = lineind+1
				line = x[lineind]
				z = line.strip('\n').split()
				nNodes = int(z[0])
				score = float(z[1])
				cComp = []
				for k in range(nNodes):
					cComp.append(int(z[2+k]))
				cComp = tuple(cComp)
				# Only consider size-2 c-components now...
				biEdges = (cComp,)
				if cComp not in self.originalScores.keys():
					self.originalScores[cComp,biEdges] = {}
				parsets = []
				for p in range(len(cComp)):
					lineind = lineind+1
					z = x[lineind].strip('\n').split()
					parset = [int(z[i]) for i in range(2,len(z))]
					parsets.append(tuple(parset))
				self.originalScores[cComp,biEdges][tuple(parsets)] = score
		
	def Prune_scores(self,prune_more=False):
		t0P = time.time()
		sum1 = 0
		sum2 = 0
		Dpars = {}
		for D in self.originalScores.keys():
			if len(D[0]) == 1:
				Dpars[D] = list(self.originalScores[D].keys())
				Dparscopy = Dpars[D].copy()
				for ind in range(len(Dpars[D])):
					for i in range(ind+1,len(self.originalScores[D])):
						if set(Dpars[D][ind][0]).difference(set(Dpars[D][i][0])) == set() and self.originalScores[D][Dpars[D][ind]] >= self.originalScores[D][Dpars[D][i]]:
							if Dpars[D][i] in Dparscopy:
								Dparscopy.remove(Dpars[D][i])
				sum1 = sum1+len(Dpars[D])
				Dpars[D] = Dparscopy.copy()
				sum2 = sum2+len(Dpars[D])

			# size 2 c-component
			elif len(D[0]) == 2:
				Dpars[D] = list(self.originalScores[D].keys())
				Dparscopy = Dpars[D].copy()
				for ind in range(len(Dpars[D])):
					# a = D[0][0], b = D[0][1]
					# a<->b vs a<->b
					delInd = False
					for i in range(ind):
						if set(Dpars[D][i][0]).difference(set(Dpars[D][ind][0])) == set() and set(Dpars[D][i][1]).difference(set(Dpars[D][ind][1])) == set() and self.originalScores[D][Dpars[D][ind]] <= self.originalScores[D][Dpars[D][i]]:
							Dparscopy.remove(Dpars[D][ind])
							delInd = True
							break
					if delInd == True:
						continue
					# a<->b vs a,b
					Da = ((D[0][0],),())
					Db = ((D[0][1],),())
					DaPars = list(self.originalScores[Da].keys())
					DbPars = list(self.originalScores[Db].keys())
					maxa = -float('inf')
					maxb = -float('inf')
					for ind1 in range(len(DaPars)):
						if set(DaPars[ind1][0]).difference(set(Dpars[D][ind][0])) == set() and self.originalScores[Da][DaPars[ind1]] > maxa:
							maxa = self.originalScores[Da][DaPars[ind1]]
					for ind2 in range(len(DbPars)):
						if set(DbPars[ind2][0]).difference(set(Dpars[D][ind][1])) == set() and self.originalScores[Db][DbPars[ind2]] > maxb:
							maxb = self.originalScores[Db][DbPars[ind2]]
					if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb:
						Dparscopy.remove(Dpars[D][ind])
						continue
					tol = 1e-10
					# a<->b vs a->b
					Da = ((D[0][0],),())
					Db = ((D[0][1],),())
					DbPars = list(self.originalScores[Db].keys())
					maxa = self.originalScores[Da][((),)]
					maxb = -float('inf')
					for ind2 in range(len(DbPars)):
						if set(DbPars[ind2][0]).difference(set(Dpars[D][ind][1])) == {D[0][0]} and self.originalScores[Db][DbPars[ind2]] > maxb:
							maxb = self.originalScores[Db][DbPars[ind2]]
					if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb+tol:
						Dparscopy.remove(Dpars[D][ind])
						continue
					# a<->b vs a<-b
					Da = ((D[0][0],),())
					Db = ((D[0][1],),())
					DaPars = list(self.originalScores[Da].keys())
					maxa = -float('inf')
					maxb = self.originalScores[Db][((),)]
					for ind1 in range(len(DaPars)):
						if set(DaPars[ind1][0]).difference(set(Dpars[D][ind][0])) == {D[0][1]} and self.originalScores[Da][DaPars[ind1]] > maxa:
							maxa = self.originalScores[Da][DaPars[ind1]]
					if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb+tol:
						Dparscopy.remove(Dpars[D][ind])
						continue
				sum1 = sum1+len(Dpars[D])
				Dpars[D] = Dparscopy.copy()
				sum2 = sum2+len(Dpars[D])
			elif len(D[0]) == 3:
				Dpars[D] = list(self.originalScores[D].keys())
				Dparscopy = Dpars[D].copy()
				for ind in range(len(Dpars[D])):
					# a<->b<->c
					if len(D[1]) == 2:
						Nodeb = 0
						for i in D[0]:
							if i in D[1][0] and i in D[1][1]:
								Nodeb = i
								break
						Nodea = 0
						for i in D[0]:
							if i != Nodeb and i in D[1][0]:
								Nodea = i
								break
						Nodec = 0
						for i in D[0]:
							if i != Nodeb and i in D[1][1]:
								Nodec = i
								break
						if prune_more == True and (Nodea in Dpars[D][ind][D[0].index(Nodec)] or Nodec in Dpars[D][ind][D[0].index(Nodea)]):
							Dparscopy.remove(Dpars[D][ind])
							continue
						Da = ((Nodea,),())
						Db = ((Nodeb,),())
						Dc = ((Nodec,),())
						Dab = (tuple(sorted((Nodea,Nodeb))),(tuple(sorted((Nodea,Nodeb))),))
						Dbc = (tuple(sorted((Nodeb,Nodec))),(tuple(sorted((Nodeb,Nodec))),))
						# a<->b<->c vs a,b,c
						DaPars = list(self.originalScores[Da].keys())
						DbPars = list(self.originalScores[Db].keys())
						DcPars = list(self.originalScores[Dc].keys())
						maxa = -float('inf')
						maxb = -float('inf')
						maxc = -float('inf')
						for ind1 in range(len(DaPars)):
							if set(DaPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and self.originalScores[Da][DaPars[ind1]] > maxa:
								maxa = self.originalScores[Da][DaPars[ind1]]
						for ind2 in range(len(DbPars)):
							if set(DbPars[ind2][0]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and self.originalScores[Db][DbPars[ind2]] > maxb:
								maxb = self.originalScores[Db][DbPars[ind2]]
						for ind3 in range(len(DcPars)):
							if set(DcPars[ind3][0]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dc][DcPars[ind3]] > maxc:
								maxc = self.originalScores[Dc][DcPars[ind3]]
						if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb+maxc:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# a<->b<->c vs a<->b,c
						DabPars = list(self.originalScores[Dab].keys())
						DbcPars = list(self.originalScores[Dbc].keys())
						maxab = -float('inf')
						maxbc = -float('inf')
						abInda = 0
						abIndb = 1
						if Nodea > Nodeb:
							abInda = 1
							abIndb = 0
						for ind1 in range(len(DabPars)):
							if set(DabPars[ind1][abInda]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DabPars[ind1][abIndb]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and self.originalScores[Dab][DabPars[ind1]] > maxab:
								maxab = self.originalScores[Dab][DabPars[ind1]]
						if self.originalScores[D][Dpars[D][ind]] <= maxab+maxc:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# a<->b<->c vs a,b<->c
						bcIndb = 0
						bcIndc = 1
						if Nodeb > Nodec:
							bcIndb = 1
							bcIndc = 0
						for ind1 in range(len(DbcPars)):
							if set(DbcPars[ind1][bcIndb]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DbcPars[ind1][bcIndc]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dbc][DbcPars[ind1]] > maxbc:
								maxbc = self.originalScores[Dbc][DbcPars[ind1]]
						if self.originalScores[D][Dpars[D][ind]] <= maxa+maxbc:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# a<->b<->c vs a<->b<->c
						delInd = False
						for i in range(ind):
							if set(Dpars[D][i][0]).difference(set(Dpars[D][ind][0])) == set() and set(Dpars[D][i][1]).difference(set(Dpars[D][ind][1])) == set() and set(Dpars[D][i][2]).difference(set(Dpars[D][ind][2])) == set() and self.originalScores[D][Dpars[D][ind]] <= self.originalScores[D][Dpars[D][i]]:
								Dparscopy.remove(Dpars[D][ind])
								delInd = True
								break
						if delInd == True:
							continue
					# a<->b<->c<->a
					if len(D[1]) == 3:
						Nodea = D[0][0]
						Nodeb = D[0][1]
						Nodec = D[0][2]
						Dab = ((Nodea,Nodeb),((Nodea,Nodeb),))
						Dbc = ((Nodeb,Nodec),((Nodeb,Nodec),))
						Dac = ((Nodea,Nodec),((Nodea,Nodec),))
						Dabc = ((Nodea,Nodeb,Nodec),((Nodea,Nodeb),(Nodeb,Nodec)))
						if Dabc not in self.originalScores.keys():
							Dabc = ((Nodea,Nodeb,Nodec),((Nodeb,Nodec),(Nodea,Nodeb)))
						Dbca = ((Nodea,Nodeb,Nodec),((Nodeb,Nodec),(Nodea,Nodec)))
						if Dbca not in self.originalScores.keys():
							Dbca = ((Nodea,Nodeb,Nodec),((Nodea,Nodec),(Nodeb,Nodec)))
						Dcab = ((Nodea,Nodeb,Nodec),((Nodea,Nodec),(Nodeb,Nodec)))
						if Dcab not in self.originalScores.keys():
							Dcab = ((Nodea,Nodeb,Nodec),((Nodeb,Nodec),(Nodea,Nodec)))
						DaPars = list(self.originalScores[Da].keys())
						DbPars = list(self.originalScores[Db].keys())
						DcPars = list(self.originalScores[Dc].keys())
						maxa = -float('inf')
						maxb = -float('inf')
						maxc = -float('inf')
						for ind1 in range(len(DaPars)):
							if set(DaPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and self.originalScores[Da][DaPars[ind1]] > maxa:
								maxa = self.originalScores[Da][DaPars[ind1]]
						for ind2 in range(len(DbPars)):
							if set(DbPars[ind2][0]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and self.originalScores[Db][DbPars[ind2]] > maxb:
								maxb = self.originalScores[Db][DbPars[ind2]]
						for ind3 in range(len(DcPars)):
							if set(DcPars[ind3][0]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dc][DcPars[ind3]] > maxc:
								maxc = self.originalScores[Dc][DcPars[ind3]]
						DabPars = list(self.originalScores[Dab].keys())
						DbcPars = list(self.originalScores[Dbc].keys())
						DacPars = list(self.originalScores[Dac].keys())
						maxab = -float('inf')
						maxbc = -float('inf')
						maxac = -float('inf')
						for ind1 in range(len(DabPars)):
							if set(DabPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DabPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and self.originalScores[Dab][DabPars[ind1]] > maxab:
								maxab = self.originalScores[Dab][DabPars[ind1]]
						for ind1 in range(len(DbcPars)):
							if set(DbcPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DbcPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dbc][DbcPars[ind1]] > maxbc:
								maxbc = self.originalScores[Dbc][DbcPars[ind1]]
						for ind1 in range(len(DacPars)):
							if set(DacPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DacPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dac][DacPars[ind1]] > maxac:
								maxac = self.originalScores[Dac][DacPars[ind1]]
						DabcPars = list(self.originalScores[Dabc].keys())
						DbcaPars = list(self.originalScores[Dbca].keys())
						DcabPars = list(self.originalScores[Dcab].keys())
						maxabc = -float('inf')
						maxbca = -float('inf')
						maxcab = -float('inf')
						for ind1 in range(len(DabcPars)):
							if set(DabcPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DabcPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DabcPars[ind1][2]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dabc][DabcPars[ind1]] > maxabc:
								maxabc = self.originalScores[Dabc][DabcPars[ind1]]
						for ind1 in range(len(DbcaPars)):
							if set(DbcaPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DbcaPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DbcaPars[ind1][2]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dbca][DbcaPars[ind1]] > maxbca:
								maxbca = self.originalScores[Dbca][DbcaPars[ind1]]
						for ind1 in range(len(DcabPars)):
							if set(DcabPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DcabPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DcabPars[ind1][2]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dcab][DcabPars[ind1]] > maxcab:
								maxcab = self.originalScores[Dcab][DcabPars[ind1]]
						# vs a,b,c
						if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb+maxc:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# vs a<->b,c
						if self.originalScores[D][Dpars[D][ind]] <= maxab+maxc:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# vs a, b<->c
						if self.originalScores[D][Dpars[D][ind]] <= maxa+maxbc:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# vs b, c<->a
						if self.originalScores[D][Dpars[D][ind]] <= maxac+maxb:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# vs a<->b<->c
						if self.originalScores[D][Dpars[D][ind]] <= maxabc:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# vs b<->c<->a
						if self.originalScores[D][Dpars[D][ind]] <= maxbca:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# vs c<->a<->b
						if self.originalScores[D][Dpars[D][ind]] <= maxcab:
							Dparscopy.remove(Dpars[D][ind])
							continue
						# vs a<->b<->c<->a
						delInd = False
						for i in range(ind):
							if set(Dpars[D][i][0]).difference(set(Dpars[D][ind][0])) == set() and set(Dpars[D][i][1]).difference(set(Dpars[D][ind][1])) == set() and set(Dpars[D][i][2]).difference(set(Dpars[D][ind][2])) == set() and self.originalScores[D][Dpars[D][ind]] <= self.originalScores[D][Dpars[D][i]]:
								Dparscopy.remove(Dpars[D][ind])
								delInd = True
								break
						if delInd == True:
							continue
				sum1 = sum1+len(Dpars[D])
				Dpars[D] = Dparscopy.copy()
				sum2 = sum2+len(Dpars[D])
			else:
				Dpars[D] = list(self.originalScores[D].keys())
				sum1 = sum1+len(Dpars[D])
				sum2 = sum2+len(Dpars[D])
		print(str(sum1)+" vs "+str(sum2)+", pruning time: "+str(time.time()-t0P))
		fileName = '../Results/'+self.instance+'_cplex.log'
		f = open(fileName,"a")
		f.write(str(sum1)+" vs "+str(sum2)+", pruning time: "+str(time.time()-t0P))
		f.close

		self.scores = {}
		for D in self.originalScores.keys():
			for Dpar in Dpars[D]:
				if D not in self.scores.keys():
					self.scores[D] = {}
				self.scores[D][Dpar] = self.originalScores[D][Dpar]

	def Initialize(self,prune=True,dag=False,printsc=False,prune_parentInDistrict=False):
		if prune == True:
			self.Prune_scores(prune_more=prune_parentInDistrict)
		else:
			self.scores = self.originalScores
		self.V = set()
		self.cComps = []
		for D in self.scores.keys():
			self.cComps.append(D)
			self.V = self.V.union(set(D[0]))
		self.iComps = {}
		for i in self.V:
			self.iComps[i] = []
			for D in self.cComps:
				if i in D[0]:
					self.iComps[i].append(self.cComps.index(D))
		self.dPars = {}
		for d in range(len(self.cComps)):
			self.dPars[d] = []
			for par in self.scores[self.cComps[d]].keys():
				self.dPars[d].append(par)
		self.iPars = {}
		for i in self.V:
			self.iPars[i] = []
			for d in self.iComps[i]:
				for W in self.dPars[d]:
					if W[self.cComps[d][0].index(i)] not in self.iPars[i]:
						self.iPars[i].append(W[self.cComps[d][0].index(i)])
		if printsc == True:
			for d in range(len(self.cComps)):
				print(str(self.cComps[d])+':')
				print(self.scores[self.cComps[d]])
				print('\n')
		
		
		self.biPars = {}
		for D in self.cComps:
			for bi in D[1]:
				if bi not in self.biPars.keys():
					self.biPars[bi] = []
				for W in self.dPars[self.cComps.index(D)]:
					biPar = (W[D[0].index(bi[0])],W[D[0].index(bi[1])])
					if biPar not in self.biPars[bi]:
						self.biPars[bi].append(biPar)
		self.m = cplex.Cplex()
		self.z = {}
		for d in range(len(self.cComps)):
			for dp in range(len(self.dPars[d])):
				if dag == False or len(self.cComps[d][0]) <= 1:
					self.z[d,dp] = self.m.variables.add(obj=[self.scores[self.cComps[d]][self.dPars[d][dp]]],lb=[0],ub=[1],types=['B'],names=['z'+str(d)+','+str(dp)])
				else:
					self.z[d,dp] = self.m.variables.add(obj=[self.scores[self.cComps[d]][self.dPars[d][dp]]],lb=[0],ub=[0],types=['B'],names=['z'+str(d)+','+str(dp)])
		
		self.m.objective.set_sense(self.m.objective.sense.maximize)
		
		
		for i in self.V:
			self.m.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = ['z'+str(d)+','+str(dp) for d in self.iComps[i] for dp in range(len(self.dPars[d]))], val = [1]*sum(len(self.dPars[d]) for d in self.iComps[i]))], senses=["E"], rhs=[1])

		self.indInv = []
		for i in self.V:
			for j in range(i+1,len(self.V)):
				self.indInv.append((i,j))

		self.udE = range(len(self.indInv))

		self.bi = {}
		for e in self.udE:
			self.bi[e] = self.m.variables.add(obj=[0],types=['C'],names=['bi'+str(e)])
			zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if self.indInv[e] in self.cComps[d][1]]
			self.m.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set+['bi'+str(e)], val = [1]*len(zindex_set)+[-1])], senses=["E"], rhs=[0])
			
		self.x = {}
		for i in self.V:
			for ip in range(len(self.iPars[i])):
				self.x[i,ip] = self.m.variables.add(obj=[0],types=['C'],names=['x'+str(i)+','+str(ip)])
				zindex_set = ['z'+str(d)+','+str(dp) for d in self.iComps[i] for dp in range(len(self.dPars[d])) if self.dPars[d][dp][self.cComps[d][0].index(i)] == self.iPars[i][ip]]
				self.m.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set+['x'+str(i)+','+str(ip)], val = [1]*len(zindex_set)+[-1])], senses=["E"], rhs=[0])

	def biClusterToIneq(self,C,ii,jj):
		if jj < ii:
			cp = jj
			jj = ii
			ii = cp
		ifLHS = {(d,dp):False for d in range(len(self.cComps)) for dp in range(len(self.dPars[d]))}
		for d in range(len(self.cComps)):
			if len(set(self.cComps[d][0])&set([ii,jj]))!=1:
				vs = [v for v in C if v in self.cComps[d][0]]
				if len(vs) > 0:
					for dp in range(len(self.dPars[d])):
						for v in vs:
							if set(self.dPars[d][dp][self.cComps[d][0].index(v)])&C.union(set([ii,jj])) == set():
								ifLHS[d,dp] = True
								break
				if (ii,jj) in self.cComps[d][1]:
					for dp in range(len(self.dPars[d])):
						if set(self.dPars[d][dp][self.cComps[d][0].index(ii)])&C.union(set([ii,jj])) == set() and set(self.dPars[d][dp][self.cComps[d][0].index(jj)])&C.union(set([ii,jj])) == set():
							ifLHS[d,dp] = True
		return ifLHS

	def ClusterToIneq(self,C):
		ifLHS = {(d,dp):False for d in range(len(self.cComps)) for dp in range(len(self.dPars[d]))}
		for d in range(len(self.cComps)):
			vs = [v for v in C if v in self.cComps[d][0]]
			if len(vs) > 0:
				for dp in range(len(self.dPars[d])):
						for v in vs:
							if set(self.dPars[d][dp][self.cComps[d][0].index(v)])&set(C) == set():
								ifLHS[d,dp] = True
								break
		return ifLHS

	def Solve_with_cb(self,CG=False):
		t0 = time.time()
		ContinueCondt = True
		LPiter = 0
		Objvalue = float('inf')
		nvar = sum(len(self.dPars[d]) for d in range(len(self.cComps)))
		nbi = len(self.udE)
		ncluster = 0
		nbicluster = 0
		out = self.m.set_results_stream(None)
		out = self.m.set_log_stream(None)

		while ContinueCondt == True and LPiter < 100:
			nz_bi = 0
			nz_z = 0
			ContinueCondt = False
			LPiter = LPiter+1
			self.m.parameters.mip.limits.nodes = 0
			model_soln = self.m.solve()
			PrevObjvalue = Objvalue
			Objvalue = self.m.solution.get_objective_value()
			
			for d in range(len(self.cComps)):
				for dp in range(len(self.dPars[d])):
					if self.m.solution.get_values(['z'+str(d)+','+str(dp)])[0] > 0:
						nz_z = nz_z+1
			for e in self.udE:
				if self.m.solution.get_values(['bi'+str(e)])[0] > 0:
					nz_bi = nz_bi+1
			
			print('LP iter '+str(LPiter)+', ObjVal: '+str(Objvalue)+', time: '+str(time.time()-t0)+', frac. of nonzero variables: '+\
				str(nz_z)+'/'+str(nvar)+', frac. of nonzero bidirected edges: '+str(nz_bi)+'/'+str(nbi)+', # cluster: '+str(ncluster)+', # bi-cluster: '+str(nbicluster))
		
			ncluster = 0
			nbicluster = 0
			bi_value = {}
			x_value = {}
			z_value = {}
			for e in self.udE:
				bi_value[e] = self.m.solution.get_values(['bi'+str(e)])[0]
			for i in self.V:
				for ip in range(len(self.iPars[i])):
					x_value[i,ip] = self.m.solution.get_values(['x'+str(i)+','+str(ip)])[0]
			for d in range(len(self.cComps)):
				for dp in range(len(self.dPars[d])):
					z_value[d,dp] = self.m.solution.get_values(['z'+str(d)+','+str(dp)])[0]

			wt = {}
			for i in self.V:
				for ip in range(len(self.iPars[i])):
					if x_value[i,ip] > 1e-6:
						for par in self.iPars[i][ip]:
							if (par,i) not in wt.keys():
								wt[(par,i)] = x_value[i,ip]
							else:
								wt[(par,i)] = wt[(par,i)]+x_value[i,ip]
			telist = []
			for (i,j) in wt.keys():
				if wt[(i,j)] >= 1-1e-6:
					telist.append(i)
					telist.append(j)
			allcyc = []
			
			allcyc = dircyc(len(self.V),int(len(telist)/2),telist)
			for Cluster in allcyc:
				ifLHS = self.ClusterToIneq(Cluster)
				zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
				self.m.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set, val = [1]*len(zindex_set))], senses=["G"], rhs=[1])
				
			if len(allcyc) > 0:
				ncluster = len(allcyc)
				ContinueCondt = True
			else:
				nnode = len(self.V)
				gcnodes = []
				gcweight = []
				gcparents = []
				for d in range(len(self.cComps)):
					for dp in range(len(self.dPars[d])):
						if z_value[d,dp] > 0:
							gcnodes.append(list(self.cComps[d][0]))
							gcweight.append(z_value[d,dp])
							gcparents.append(list(list(pars) for pars in self.dPars[d][dp]))
				for i in range(len(gcparents)):
					for j in range(len(gcparents[i])):
						if len(gcparents[i][j]) == 0:
							gcparents[i][j].append(nnode)
				dircycs = contract_heur(nnode+1, gcnodes, gcparents, gcweight,1.0)
				for Cluster in dircycs:
					ifLHS = self.ClusterToIneq(Cluster)
					zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
					self.m.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set, val = [1]*len(zindex_set))], senses=["G"], rhs=[1])
				ncluster = ncluster+len(dircycs)

				aldircycs = contract_heur_bdir(nnode+1, gcnodes, gcparents, gcweight)
				for Cluster in aldircycs:
					ii = Cluster[0]
					jj = Cluster[1]
					e = self.indInv.index((ii,jj))
					C = set(Cluster[2:])
					ifLHS = self.biClusterToIneq(C,ii,jj)
					zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
					self.m.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set+['bi'+str(e)], val = [1]*len(zindex_set)+[-1])], senses=["G"], rhs=[0])
				nbicluster = nbicluster+len(aldircycs)

				if ncluster+nbicluster > 0 and Objvalue > PrevObjvalue:
					ContinueCondt = True
		self.m.parameters.mip.limits.nodes = 9223372036800000000
		tRoot = time.time()-t0

		class LazyCallback(cplex.callbacks.LazyConstraintCallback):
			def __init__(self, env):
				super().__init__(env)

			def __call__(self):
				NoCluster = False
				bi_value = {}
				x_value = {}
				for e in self.udE:
					bi_value[e] =  self.get_values(['bi'+str(e)])[0]
				for i in self.V:
					for ip in range(len(self.iPars[i])):
						x_value[i,ip] = self.get_values(['x'+str(i)+','+str(ip)])[0]
				
				ActiveEdgeList = []
				
				for i in self.V:
					for ip in range(len(self.iPars[i])):
						if x_value[i,ip] > 0.5:
							for par in self.iPars[i][ip]:
								ActiveEdgeList.append(par)
								ActiveEdgeList.append(i)
				ne = int(len(ActiveEdgeList)/2)

				# Detecting directed cycles and adding cluster inequalities
				cycList = dircyc(len(self.V),ne,ActiveEdgeList)
				for Cluster in cycList:
					ifLHS = self.ClusterToIneq(Cluster)
					zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
					self.add(constraint=cplex.SparsePair(ind = zindex_set, val = [1]*len(zindex_set)), sense="G", rhs=1)
				# Detecting almost directed cycles and adding bi-cluster inequalities
				for e in self.udE:
					if bi_value[e] > 0.5:
						ii = self.indInv[e][0]
						jj = self.indInv[e][1]
						adcycList = almostdircyc(len(self.V),ne,ActiveEdgeList,ii,jj)+almostdircyc(len(self.V),ne,ActiveEdgeList,jj,ii)
						for Cluster in adcycList:
							C = set(Cluster[1:-1])
							ifLHS = self.biClusterToIneq(C,ii,jj)
							zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
							self.add(constraint=cplex.SparsePair(ind = zindex_set+['bi'+str(e)], val = [1]*len(zindex_set)+[-1]), sense="G", rhs=0)
		self.m.parameters.preprocessing.presolve.set(0)
		self.m.parameters.mip.strategy.search.set(1)
		lazyModel = self.m.register_callback(LazyCallback)
		lazyModel.udE = self.udE
		lazyModel.bi = self.bi
		lazyModel.iPars = self.iPars
		lazyModel.x = self.x
		lazyModel.z = self.z
		lazyModel.V = self.V
		lazyModel.dPars = self.dPars
		lazyModel.cComps = self.cComps
		lazyModel.indInv = self.indInv
		lazyModel.ClusterToIneq = self.ClusterToIneq
		lazyModel.biClusterToIneq = self.biClusterToIneq
		self.m.solve()

		fileName = '../Results/'+self.instance+'_cplex.log'
		wrtStr = 'Time at root: '+str(tRoot)+'\n'
		wrtStr = 'Total solution time: '+str(time.time()-t0)+'\n'
		wrtStr = 'Score: '+str(self.m.solution.get_objective_value())+'\n'
		print('Score: '+str(self.m.solution.get_objective_value()))
		print('Bidirected edges: ')
		wrtStr = wrtStr+'Bidirected edges: \n'
		for e in self.udE:
			if self.m.solution.get_values(['bi'+str(e)])[0] > 0.5:
				print(self.indInv[e])
				wrtStr = wrtStr+str(self.indInv[e])+'\n'
		print('Parent sets: ')
		wrtStr = wrtStr+'Parent sets: \n'
		for i in self.V:
			for ip in range(len(self.iPars[i])):
				if self.m.solution.get_values(['x'+str(i)+','+str(ip)])[0] > 0.5:
					print(str(i)+': '+str(self.iPars[i][ip]))
					wrtStr = wrtStr+str(i)+': '+str(self.iPars[i][ip])+'\n'
		print('z solution: ')
		wrtStr = wrtStr+'z solution: '
		for d in range(len(self.cComps)):
			for dp in range(len(self.dPars[d])):
				if self.m.solution.get_values(['z'+str(d)+','+str(dp)])[0] > 0.5:
					print(str(self.cComps[d])+': '+str(self.dPars[d][dp])+', score: '+str(self.scores[self.cComps[d]][self.dPars[d][dp]]))
					wrtStr = wrtStr+str(self.cComps[d])+': '+str(self.dPars[d][dp])+', score: '+str(self.scores[self.cComps[d]][self.dPars[d][dp]])+'\n'
		f = open(fileName,"a")
		f.write(wrtStr)
		f.close



if __name__ == '__main__':
	scoresets = ['score_example']
	for instName in scoresets:
		inst = BNSLlvInst(instName)
		
		inst.readFromPkl()
		inst.Initialize(prune=True,dag=True)
		inst.Solve_with_cb()
		
		inst.readFromPkl()
		inst.Initialize(prune=True,printsc=False)
		inst.Solve_with_cb()
