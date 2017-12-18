#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#Модуль для сравнения крит групп Ищет крит группу в каждом органе по всем точкам карты, ругается,
#если крит группа не одна на карте.
import sys
import os
from math import cos,sin,radians,sqrt,fabs,acos,pi,degrees,floor,log10
import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import xml.etree.ElementTree
from math import cos,sin,radians,sqrt,fabs,acos,pi,degrees,floor
axeRadius = [0.0,0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,
			22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0,42.0,44.0,46.0,48.0,50.0,
		52.0,54.0,56.0,58.0,60.0]
import openpyxl
mpl.rcParams['font.family'] = 'fantasy'
mpl.rcParams['font.fantasy'] = 'Times New Roman', 'Ubuntu','Arial','Tahoma','Calibri'
points=[]
pathToSIdir = "/home/egor/quest/TIC_graph/0_xxx/"
pathToOut =sys.argv[1]# "/home/egor/Programs/stat/VVER_TOI_scenario_3/results/7200 s/out.xml"

pointsForAnalysis = []

def cleanBlock(block):
	for i in range (len(block)-1, -1, -1):
		if(block[i].strip().strip("\t") == ""):
			del block[i]
	return block

def find_element_in_list(element,list_element):
        try:
		index_element=list_element.index(element)
		return index_element
	except ValueError:
		return -1 

def find_substr_in_list(element,list_element):
	index_element = []
	for i, item in enumerate(list_element):
		if item.find(element) >= 0:
			index_element.append(i)
	return index_element
def readLineF(line):
	lst = line.rstrip()
	vals = filter(None,lst.split(" "))
	return [float(x) for x in vals]

def readLineS(line, statement = ";"):
	lst = line.rstrip()
	vals = filter(None,lst.split(statement))
	return [str(x) for x in vals]

def make_sure_path_exists(path):
	try: 
		os.makedirs(path)
	except OSError:
		if not os.path.isdir(path):
			raise

def isEmpty(s):
    return not bool(s or s.strip())

def reaNPFromFile(pathToFile):
	inp = open(pathToFile,"r")
	lines = inp.readlines()
	inp.close()
	lines = [line for line in lines if line.find("#") < 0]
	pointsForAnalysis = []
	for line in lines:
		pointsForAnalysis.append([str(x.strip("\"")) if i == 0 else float(x) for (i,x) in enumerate(filter(None,line.split(",")))])
	return pointsForAnalysis

class GeoGrid:
	def __init__(self,gridFunction):
		self.lonmin = 0.0
		self.lonmax = 0.0
		self.latmin = 0.0
		self.latmax = 0.0
		self.countx = 0.0
		self.county = 0.0
		self.fmin   = 0.0
		self.fmax   = 0.0
		
		self.calcFunctionId = -1
		self.tipObl = -1
		self.isotop = ""
		self.phchForm = ""
		self.organ = ""
		self.age = ""
		self.timeEnd = 0.0
		
		self.data   = []

		if "calcFunctionId" not in gridFunction.attrib:
			self.calcFunctionId = -1
		else:
			self.calcFunctionId = int(gridFunction.get("calcFunctionId").strip())
		
		if "tipObl" not in gridFunction.attrib:
			self.tipObl = -1
		else:
			self.tipObl = int(gridFunction.get("tipObl").strip())
		
		if "isotop" not in gridFunction.attrib:
			self.isotop = ""
		else:
			self.isotop = gridFunction.get("isotop").strip()
			
		if "phchForm" not in gridFunction.attrib:
			self.phchForm = ""
		else:
			self.phchForm = gridFunction.get("phchForm").strip()
			
		if "organ" not in gridFunction.attrib:
			self.organ = ""
		else:
			self.organ = gridFunction.get("organ").strip()
			
		if "age" not in gridFunction.attrib:
			self.age = ""
		else:
			self.age = gridFunction.get("age").strip()
			
		if "timeEnd" not in gridFunction.attrib:
			self.timeEnd = 0.0
		else:
			self.timeEnd = float(gridFunction.get("timeEnd").strip())
				
		lst = gridFunction.text.split("\n")
		lst = cleanBlock(lst)

		self.countx = int(lst[0].split(" ")[1])
		self.county = int(lst[1].split(" ")[1])
		self.lonmin = float(lst[2].split(" ")[1])
		self.latmin = float(lst[3].split(" ")[1])
		self.dlon = float(lst[4].split(" ")[1])
		self.dlat  =  float(lst[5].split(" ")[1])
		self.lonmax = self.lonmin+self.dlon*(self.countx-1)
		self.latmax = self.latmin+self.dlat*(self.county-1)
		self.matrix=np.array([[0.0]*101]*101)
		lst = lst[7:]
		self.data   = [0.0]*self.countx*self.county
		for j in range(len(lst)): #количество строк
			lin = readLineF(lst[j])
			for i in range(len(lin)):
				self.data[i*self.county+(self.county-j-1)] = lin[i]
		for j in range(len(lst)):
			for i in range(len(lin)):
				if self.data[i*self.county+(self.county-j-1)]==0:
					self.matrix[j][i]=None
				else:
					self.matrix[j][i]=self.data[i*self.county+(self.county-j-1)]
		
		self.fmin   = self.findMin()
		self.fmax   = self.findMax()

		#self.printASCIIGRDFile("test.grd")
	
	def findMin(self):
		if(len(self.data) == 0):
			return -1
		ret = self.data[0]
		for v in self.data:
			if(ret > v):
				ret = v
		return ret
		
	def findMax(self):
		if(len(self.data) == 0):
			return -1
		ret = self.data[0]
		for v in self.data:
			if(ret < v):
				ret = v
		return ret			

	def printASCIIGRDFile(self,filename):
		f = open(str(pathToGRD+filename), 'wt')
		f.write("DSAA" + '\r\n')
		f.write(str(self.countx)+"\t"+str(self.county) + '\r\n')
		f.write('%0.12f' % self.lonmin+"\t"+'%0.12f' % self.lonmax + '\r\n')
		f.write('%0.12f' % self.latmin+"\t"+'%0.12f' % self.latmax + '\r\n')
		f.write('%0.12e' % self.findMin()+"\t"+'%0.12e' % self.findMax() + '\r\n')
		for j in range(0,self.county):
			for i in range(0,self.countx):
				f.write('%0.12e' % self.data[i*self.county+j]+'\t')
			f.write('\r\n')
		
		f.close()
		return

	def getValue(self,lon,lat):
		if(lon<self.lonmin or lon>self.lonmax):
			return 0.
		if(lat<self.latmin or lat>self.latmax):
			return 0.
		ci = int(floor((lon-self.lonmin)/self.dlon))
		cj = int(floor((lat-self.latmin)/self.dlat))
		ci1 = 0
		ci2 = 0
		cj1 = 0
		cj2 = 0
		
		if ci == self.countx-1:
			ci2 = ci
			ci1 = ci-1
		else:
			ci1 = ci
			ci2 = ci+1
	    
		if cj == self.county-1:
			cj2 = cj
			cj1 = cj-1
		else:
			cj1 = cj
			cj2 = cj+1
		
		f11 = self.data[ci1*self.county+cj1] #+ -
		f21 = self.data[ci2*self.county+cj1] #+ -

		f12 = self.data[ci1*self.county+cj2] #- +
		f22 = self.data[ci2*self.county+cj2] #- +
		
		
		x1 = self.lonmin+self.dlon*ci1
		x2 = self.lonmin+self.dlon*ci2
		y1 = self.latmin+self.dlat*cj1
		y2 = self.latmin+self.dlat*cj2
		fy1=(f21-f11)/(x2-x1)*lon+(f11-(f21-f11)/(x2-x1)*x1)
		fy2=(f22-f12)/(x2-x1)*lon+(f12-(f22-f12)/(x2-x1)*x1)
		v = (fy2-fy1)/(y2-y1)*lat+(fy1-(fy2-fy1)/(y2-y1)*y1)
		return  v
	def angleOf(self,dX,dY):
		dFi = 0.0
		dR = sqrt(fabs(dX*dX+dY*dY))
		if (fabs(dR) > 0):
			dFi  = acos(dX/dR)
		if (dX <= 0 and dY < 0):
			dFi  = 2.0*pi - dFi
		if (dX >  0 and dY < 0):
			dFi  = 2.0*pi - dFi
		return dFi

#readGridFromASCIIGRDFile
def rg(title):
	#print "title"
	#print title
	#print "title"
	inp = open(title)
	lines = inp.readlines()
	inp.close()
	countx = int(lines[1].split("\t")[0].strip())
	county = int(lines[1].split("\t")[1].strip())
	lonmin = float(lines[2].split("\t")[0].strip())
	lonmax = float(lines[2].split("\t")[1].strip())
	latmin = float(lines[3].split("\t")[0].strip())
	latmax = float(lines[3].split("\t")[1].strip())
	dlon = (lonmax-lonmin)/(countx-1)
	dlat = (latmax-latmin)/(county-1)
	grid = GeoGrid(countx,county,lonmin,latmin,dlon,dlat,lonmax,latmax)
	lines = lines[5:]
	for j,line in enumerate(lines):
		if isEmpty(line):
			continue
		arr = readLineF(line, "\t")
		for i,val in enumerate(arr):
			grid.data[i*county+j] = val
	grid.fmin = grid.findMin()
	grid.fmax = grid.findMax()
	return grid
def printASCIIGRDFile(filename,countx, county, lonmin, lonmax, latmin, latmax, critArr, allk):
		f = open(filename, 'wt')
		f.write("DSAA" + '\r\n')
		f.write(str(countx)+"\t"+str(county) + '\r\n')
		f.write('%0.12f' % (lonmin)+"\t"+'%0.12f' % (lonmax) + '\r\n')
		f.write('%0.12f' % (latmin)+"\t"+'%0.12f' % (latmax) + '\r\n')
		f.write('%0.12e' % 1+"\t"+'%0.12e' % 6 + '\r\n')
		for j in range(0,county):
			for i in range(0,countx):
				f.write('%0.12e' % critArr[i*county+j]+'\t')
			f.write('\r\n')
		
		f.close()
		return
def main():
	mypath=os.path.dirname(os.path.realpath( __file__ ))
	os.chdir(mypath)
	old_settings = np.seterr(all='print')
	print sys.argv[1]
	if os.path.isdir("./critGroup_new/"+sys.argv[1].split('/')[7].split(' ')[0]+"/")==False:
		os.makedirs("./critGroup_new/"+sys.argv[1].split('/')[7].split(' ')[0]+"/")
	#os.makedirs("./critGroup/"+sys.argv[1].split('/')[7].split(' ')[0]+"/")
	os.chdir("./critGroup_new/"+sys.argv[1].split('/')[7].split(' ')[0]+"/")
	lst = ["calcFunctionId","tipObl","isotop","phchForm","organ","age","timeEnd"]
	headers = []
	gridsWithMaxValues = []
	path = sys.argv[1]
	source = open(path, 'rb')
	et = xml.etree.ElementTree.parse(path)
	root = et.getroot()
	grid = root.find('grid')
	positionInXml = -1
	number=range(24, 197, 8)
	listOfGrid=grid.findall('gridFunction')
	#for i, gridFunction in enumerate(grid.findall('gridFunction')):
	#	positionInXml = positionInXml + 1
	#	lst=[]
		#for i in range(1, 7):
		#	#print i+143
		#	lst.append(rg("/home/egor/quest/criticalGroup/GRD/f{}.grd".format(i+143)))
			#print "/home/egor/quest/criticalGroup/GRD/f{}.grd".format(i+143)
		#countx=lst[0].countx
		#county=lst[0].county
		#lonmin=lst[0].lonmin
		#lonmax=lst[0].lonmax
		#latmin=lst[0].latmin
		#latmax=lst[0].latmax
	bb=True
	zero=np.array([0.0,0.0,0.0,0.0,0.0])
	for j in range(24, 197, 8):
		#print j
		critArray=[]
		lst=[]
		now=0
		allcrit=0
		a=True
		b=True
		array=np.array([0.0,0.0,0.0,0.0,0.0])
		for i in range(j, j+5):
			#print i		
			lst.append(GeoGrid(listOfGrid[i]))
			#m=GeoGrid(listOfGrid[i])
		#print i, i, i, i, i
		#print "///////////////////////////"
		for k in range(len(lst[0].data)):
			#print i
			array=np.array([])
			array=np.array([lst[0].data[k], lst[1].data[k],lst[2].data[k],lst[3].data[k], lst[4].data[k]])#, lst[5].data[k]])
			if a==True:
				maxar=np.array(lst[0].data)
				maxp=np.argmax(maxar)#Номер точки с макс значением
				
				#print maxar, maxp, i
				array_max=np.array([lst[0].data[maxp], lst[1].data[maxp],lst[2].data[maxp],lst[3].data[maxp], lst[4].data[maxp]])
				allcrit=np.argmax(array_max)+1
				#print allcrit
				#print allcrit
				a=False
			#print len(array), array
			#print np.argmax(array)+1
			if np.array_equal(array, zero):
				critArray.append(1.70141e+38)
			else:
				critArray.append(np.argmax(array)+1)
				now=np.argmax(array)+1

			if (j==j)and(np.array_equal(array, zero)==False): #or True) and (np.array_equal(array, zero)==False)(b==True):
				if allcrit!=now:
					#print "kek"
					delta=fabs(array[now-1]-array[allcrit-1])
					#print delta, delta/array[allcrit-1]*100.0
					#print delta/array[allcrit-1]*100.0
					#print "Try"+"   "+str(j)+"   "+str(j+4), k
					if (delta/array[allcrit-1]*100.0)>=10.0:
						print "\n"
						print "WoW"+"   "+str(j)+"   "+str(j+4), k
						print "old crit num "+str(allcrit)+" val "+str(array[allcrit-1])+"  new crit  num "+str(now)+" val "+str(array[now-1])
						print "delta   "+str(delta)+" % "+str(delta/array[allcrit-1]*100.0)
						#print array
					#else:
						#print "Zamena  "+"   "+str(j)+"   "+str(j+4)+"   "+ str(k)+"  ///////////////////////////////////////////////////"
						#critArray[len(critArray)-1]=allcrit
						
						
						
					b=False
		if bb==False:	
			fig = plt.figure()
			plt.imshow(m.matrix, interpolation='none')
			ar=np.array(critArray).reshape(101, 101)
			ar=ar.transpose()
			ar[::][::1]=ar[::][::-1]
			plt.imshow(ar, interpolation='none', alpha=0.5)
			plt.title('Simple pcolor plot')
			plt.show()
			
		printASCIIGRDFile("crit{}.grd".format(str(j)),lst[0].countx, lst[0].county, lst[0].lonmin, lst[0].lonmax, lst[0].latmin, lst[0].latmax,critArray , allcrit)
		#if j== 32:
			#break
	
	#group1=rg("/home/egor/quest/criticalGroup/GRD/f{}.grd".format(i))
	#print grid1.data
	#grid1.printASCIIGRDFile("123123.grd")

	print("Ok!")
	return

if __name__ == "__main__":
	sys.exit(main())