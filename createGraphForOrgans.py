#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#Модуль для подготовки таблиц с дозами в населенных пунктах, приведенных в файле np.dat.
#Реализует быструю (но менее точную, чем createVariationalSeriesNP2) интерполяцию и сложение насчитанных модулем createVariationalSeriesGRD сеток, отражающих 
#оценку доз с заданным уровнем доверия. По умолчанию - 95 процентный доверительный интервал
#запуск: 1. запускается модуль createVariationalSeriesFI 2. запускается модуль createVariationalSeriesGRD 3. запускается данный модуль.
#!!!!!важен порядок следования  функционалов!!!!!
#Разработан: 01.06.2017 
#Автор: Киселев А.А.
#Последняя модификация: 05.10.2017 
import sys
import os
import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from math import cos,sin,radians,sqrt,fabs,acos,pi,degrees,floor
axeRadius = [0.0,0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,
			22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0,42.0,44.0,46.0,48.0,50.0,
		52.0,54.0,56.0,58.0,60.0]
name_thyroid =  ["Щитовидная железа (облако, поверхность, ингаляция)", "Щитовидная железа (ингаляция)",
		["f166","f167","f173"], ["f173"],  "Эквивалентная доза за 10 сут",
		["f46","f47","f53"],    ["f53"],   "ОБЭ - взвешенная доза за 30 сут",
		["f102","f103","f109"], ["f109"],  "Ожидаемая эквивалентная доза за 2 сут",
		 "2.0", "Дозовый предел МАГАТЭ для вмешательства с целью недопущения развития гипотиреоза",
		 "Дозовый предел МАГАТЭ для вмешательства с целью недопущения развития гипотиреоза 2 Гр-экв",
		 "thyroid_full", "thyroid_inh", "0.2"]
name_lungs =    ["Легкие (облако, поверхность, ингаляция)", "Легкие (ингаляция)",
		["f150", "f151", "f157"], ["f157"],  "Эквивалентная доза за 10 сут",
		["f38","f39","f45"],    ["f45"],   "ОБЭ - взвешенная доза за 30 сут",
		["f86","f87","f93"], ["f93"],  "Ожидаемая эквивалентная доза за 2 сут",
		 "30.0", "Дозовый предел МАГАТЭ для вмешательства с целью недопущения развития пневмонии",
		 "Дозовый предел МАГАТЭ для вмешательства с целью недопущения развития пневмонии 2 Гр-экв",
		 "lungs_full", "lungs_inh", "0.2"]
name_redMarrow =["Красный костный мозг (облако, поверхность, ингаляция)", "Красный костный мозг (ингаляция)",
		["f142", "f143", "f149"], ["f149"],  "Эквивалентная доза за 10 сут",
		["f30","f31","f37"],    ["f37"],   "ОБЭ - взвешенная доза за 30 сут",
		["f78","f79","f85"], ["f85"],  "Ожидаемая эквивалентная доза за 2 сут",
		 "2.0", "Дозовый предел МАГАТЭ для вмешательства с целью недопущения развития пневмонии",
		 "Дозовый предел МАГАТЭ для вмешательства с целью недопущения развития пневмонии 2 Гр-экв",
		 "redMarrow_full", "redMarrow_inh", "0.2"]

#print name_thyroid
import openpyxl
mpl.rcParams['font.family'] = 'fantasy'
mpl.rcParams['font.fantasy'] = 'Times New Roman', 'Ubuntu','Arial','Tahoma','Calibri'
points=[]
pathToSIdir = "/home/egor/quest/TIC_graph/0_xxx/"

pointsForAnalysis = []

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

def readLineF(line, statement = ";"):
	lst = line.rstrip()
	vals = filter(None,lst.split(statement))
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
	def __init__(self,countx,county,lonmin,latmin,dlon,dlat,lonmax,latmax):
		self.lonmin = lonmin
		self.lonmax = lonmax
		self.latmin = latmin
		self.latmax = latmax
		self.countx = countx
		self.county = county
		self.fmin   = 0.0
		self.fmax   = 0.0		
		self.dlon = dlon
		self.dlat = dlat
		self.data   = [0.0 for o in range(countx*county)]
	
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

	def clean(self):
		for i in range(len(self.data)):
			self.data[i] = 0.0

	def printASCIIGRDFile(self,filename):
		f = open(filename, 'wt')
		f.write("DSAA" + '\r\n')
		f.write(str(self.countx)+"\t"+str(self.county) + '\r\n')
		f.write('%0.12f' % (self.lonmin)+"\t"+'%0.12f' % (self.lonmax) + '\r\n')
		f.write('%0.12f' % (self.latmin)+"\t"+'%0.12f' % (self.latmax) + '\r\n')
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
		dR = sqrt(fabs(dX*dX+dY*dY));
		if (fabs(dR) > 0):
			dFi  = acos(dX/dR);
		if (dX <= 0 and dY < 0):
			dFi  = 2.0*pi - dFi;
		if (dX >  0 and dY < 0):
			dFi  = 2.0*pi - dFi;
		return dFi;

#readGridFromASCIIGRDFile
def rg(title):
	#print "title"
	#print title
	#print "title"
	inp = open(pathToSIdir+"/"+title+".grd","r")
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


#read grids from file //////////////////////////////////////////////////////////////////////////////////////////////////////////////
def sumGridForDoseNew(*args):
	#print type(args), args, len(args)
	args=args[0]
	args=tuple(args)
	#print type(args), args, len(args)
	#print "//////////////////////////////////"
	global points
	value=[]
	p = pathToSIdir
	if len(args) == 0:
		return None
	result = rg(args[0])
	result.clean()
	for arg in args:
		g = rg(arg)
		#print arg
		for i in range(len(result.data)):
			result.data[i] = result.data[i] + g.data[i]
			#print result.data[i]
		del g
	for i in points:
		value.append(result.getValue(i[0], i[1]))
	return value

#read grids from file
def maxGridForDoseNew(*args):
	p = pathToSIdir
	if len(args) == 0:
		return None
	result = rg(args[0])
	result.clean()
	for arg in args:
		g = rg(arg)
		for i in range(len(result.data)):
			if g.data[i] > result.data[i]:
				 result.data[i] = g.data[i]
		del g
	return result


#work with existing grid
def sumGridForDose(*args):
	if len(args) == 0:
		return None
	result = rg("f0")
	result.clean()
	for arg in args:
		for i in range(len(result.data)):
			result.data[i] = result.data[i] + arg.data[i]
		
	for arg in args:
		del arg
	return result
	
#work with existing grid	
def maxGridForDose(*args):
	if len(args) == 0:
		return None
	result = rg("f0")
	result.clean()
	for arg in args:
		for i in range(len(result.data)):
			if arg.data[i] > result.data[i]:
				 result.data[i] = arg.data[i]
	for arg in args:
		del arg
	return result 
def readAxis(title):
	global points
	f=open(title,'r')
	lines = f.readlines()
	f.close()
	points=[]
	#print lines
	for i in lines:
		#print i
		points.append([float(i.split(',')[0].strip()), float(i.split(',')[1].strip())])
	return points
	
def save(name='', fmt='png'):
	pwd = os.getcwd()
	iPath = './pictures/{}'.format(fmt)
	if not os.path.exists(iPath):
		os.mkdir(iPath)
	os.chdir(iPath)
	plt.savefig('{}.{}'.format(name, fmt), fmt='png', bbox_inches='tight', dpi = 200)
	os.chdir(pwd)
def prepToSave(stroka):
	global pointsValue
	fig = plt.figure()
	#print sumGridForDoseNew("f134")
	x=axeRadius
	plt.plot(x, pointsValue[0:40])
	plt.xlabel(u'x (km)')
	plt.plot(x, pointsValue[0:40], label = stroka.decode('utf-8'))
	plt.legend()
	#plt.grid(True)
	save(stroka)
	return
def organ(string, name, lst):
	global points , pointsValue
	x=axeRadius
	maxArr=np.array([])
	points = readAxis('/home/egor/quest/TIC_graph/Axis/maxPoint_f1_95.dat')
	lines=[]
	fig = plt.figure(figsize=(10, 6)) #figsize(horiz, vert)
	gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1]) 
	ax = plt.subplot(gs[0])
	#Разбиение на 2 оси и 2 подписи
	ax.set_ylabel(u'Зв')
	ax2 = ax.twinx()
	ax2.set_ylabel(u'Гр-экв', labelpad=30)
	ax2.tick_params(axis='y', which='major', 
               labelleft='off', labelright='off', 
               left=False, right=False)
	ax.tick_params(axis='y', which='major', 
	labelleft='on', labelright='on', 
               left=True, right=True)
	#Конец 
	if string=="full":
		ax.set_title(lst[0].decode('utf-8'))
	else:
		ax.set_title(lst[1].decode('utf-8'))
	ax.set_xlabel(u'x (km)')
	
	ax_leg = fig.add_subplot(gs[1])
	
	if string=="full":
		pointsValue=sumGridForDoseNew(lst[2]) #"f166","f167","f173"

	else:
		pointsValue=sumGridForDoseNew(lst[3])
	npPoints=np.array(pointsValue)
	maxArr=np.append(maxArr, npPoints.max())
	lab1=lst[4]
	a,= ax.plot(x, pointsValue[0:40], label=lab1.decode('utf-8'))
	lines.append(a)
	
	if string=="full":
		pointsValue=sumGridForDoseNew(lst[5])
	else:
		pointsValue=sumGridForDoseNew(lst[6])
	npPoints=np.array(pointsValue)
	maxArr=np.append(maxArr, npPoints.max())
	lab2=lst[7]
	b,= ax.plot(x, pointsValue[0:40], label=lab2.decode('utf-8'))
	lines.append(b)
	
	if string=="full":
		pointsValue=sumGridForDoseNew(lst[8])
	else:
		pointsValue=sumGridForDoseNew(lst[9])
	npPoints=np.array(pointsValue)
	maxArr=np.append(maxArr, npPoints.max())
	lab2=lst[10]
	c, = ax.plot(x, pointsValue[0:40], label=lab2.decode('utf-8'))
	lines.append(c)
	maxPoint=maxArr.max()
	if maxPoint>=float(lst[11]):
		eff = np.array([float(lst[11])]*len(x))
		lab2=lst[12]
		d, = ax.plot(x, eff,  label=lab2.decode('utf-8'))
		lines.append(d)
		#ax.text(20.0, 3.0, u'Дозовый предел МАГАТЭ для вмешательства с целью недопущения развития гипотиреоза', fontsize=8)
	else:
		eff = np.array([maxPoint+float(lst[16])]*len(x))
		lab2=lst[13]
		d, = ax.plot(x, eff, color='white', label=lab2.decode('utf-8'))
		lines.append(d)
		
		#ax.text(20.0, maxPoint/2.0 , u'Дозовый предел МАГАТЭ для вмешательства с целью недопущения развития гипотиреоза 2 Гр-экв', fontsize=8)
	
	ax.grid(False, color='black', linestyle='-', linewidth=0.2)
	for line in lines:  # just to make the legend plot
		ax_leg.plot([], [], line.get_color(), label=line.get_label())
	ax_leg.legend(loc='upper left', ncol=1, fontsize=9) 
	ax_leg.axis('off')
	if string=="full":
		save("{}_{}_new".format(lst[14], name))
		print "{}_{}_new".format(lst[14], name)
	else:
		save("{}_{}_new".format(lst[15], name))
		print "{}_{}_new".format(lst[15], name)
	return

def main():
	global pathToSIdir, name_thyroid
	k=pathToSIdir
	var=pathToSIdir.split("_")[1].split("/")[0]
	#print var
	
	for var in ['50','95', '995']:
		print var
		m=k.replace('xxx', var)
		pathToSIdir=m
		#lungs("full", var)
		#lungs("inh", var)
		#thyroid("full", var)
		#thyroid("inh", var)
		#redMarrow("full", var)
		#redMarrow("inh", var)
		organ("full", var, name_thyroid)
		organ("inh", var, name_thyroid)
		organ("full", var, name_lungs)
		organ("inh", var, name_lungs)
		organ("full", var, name_redMarrow)
		organ("inh", var, name_redMarrow)
	#plt.xlabel(u'x (km)')
	#pointsValue=sumGridForDoseNew("f134","f135","f141")
	#lab1="f134, f135, f141 Эфф доза от  облака, поверх, ингал, 10 дней"
	#plt.plot(x, pointsValue[0:40])
	#pointsValue=sumGridForDoseNew("f204","f205","f141")
	#lab2="f204, f205, f141 Эфф доза от  облака, поверх, ингал"
	#plt.plot(x, pointsValue[0:40])
	#pointsValue=sumGridForDoseNew("f141")
	#lab3="f141, Эфф доза от ингаляции"
	#plt.plot(x, pointsValue[0:40])
	#plt.legend((lab1.decode('utf-8'), lab2.decode('utf-8'), lab3.decode('utf-8')), frameon=False)
	#plt.legend()
	#plt.grid(True)
	#save("123")

	



	print("Ok!")
	return

if __name__ == "__main__":
	sys.exit(main())