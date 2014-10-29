'''calculates the parameters for barrel distorsion in imagemagick'''

import subprocess
from math import sqrt, cos, sin, pi, atan2

imageconvert = "C:/Programme/ImageMagick-6.8.6-Q16/convert.exe"

logFile = open("c://doc/python2.log", 'a')
def log(line):
	try:
		logFile.write(line + "\n")
		logFile.flush()
	except:
		pass
		
class matrix:
	def __init__(self):
		self.field = [[0]]

	def set(self,y,x,value):
		while len(self.field) <= y:
			self.field.append([0])
		while len(self.field[y]) <= x:
			self.field[y].append(0)
		self.field[y][x] = value

	def get(self,y,x):
		erg = 0
		if y < len(self.field):
			if x < len(self.field[y]):
				erg = self.field[y][x]
		return erg

	def dimX(self):
		erg = 0
		for i in range(len(self.field)):
			if erg < len(self.field[i]):
				erg = len(self.field[i])
		return erg

	def copy(self):
		erg = matrix()
		for i in range(len(self.field)):
			for j in range(len(self.field[i])):
				erg.set(i,j,self.field[i][j])
		return erg

	def dimY(self):
		return len(self.field)

	def t_copy(self):
		# transponierte Matrix
		erg = matrix()
		for i in range(len(self.field)):
			for j in range(len(self.field[i])):
				erg.set(j,i,self.field[i][j])
		return erg

	def write(self):
		for v1 in self.field:
			e = ''
			for v2 in v1:
				e += "%12.5f" % v2
			print (e)
			
	def writeWMaxima(self, name='d'):
		follow1 = False
		erg = name + ':matrix('
		for v1 in self.field:
			if follow1:
				erg += ','
			follow1 = True
			erg += '['
			follow2 = False
			dimX = self.dimX()
			for v2i in range(dimX):
				if follow2:
					erg += ','
				follow2 = True
				v2 = 0.0
				if v2i < len(v1):
					v2 = v1[v2i]
				if abs(v2) < 1e-12:
					erg += '0.0'
				else:
					erg += str(v2)
			erg += ']'
		return erg + ');'

	def mul(self,y,a):
		# multipliziert Zeile y mit a
		if y < len(self.field):
			for i in range(len(self.field[y])):
				self.field[y][i] *= a

	def sub(self,y1,y2):
		# zieht von Zeile y1 die Zeile y2 ab
		if y2 < len(self.field):
			maxl = len(self.field[y2])
			if y1 < len(self.field):
				maxl = max(maxl,len(self.field[y2]))
			for i in range(maxl):
				self.set(y1,i, self.get(y1,i) - self.get(y2,i))

	def change(self,y1,y2):
		# tauscht Zeile y1 mit Zeile y2
		for i in range(self.dimX()):
			h = self.get(y1, i)
			self.set(y1, i, self.get(y2, i))
			self.set(y2, i, h)

	def mul_matrix(self, m):
		# multipliziert zwei Matrizen miteinander
		erg = matrix()
		if self.dimX() == m.dimY():
			for y in range(self.dimY()):
				for x in range(m.dimX()):
					s = 0
					for k in range(self.dimX()):
						s += self.get(y,k) * m.get(k,x)
					erg.set(y,x,s)
		return erg

	def mul_vektor(self, v):
		# multipliziert eine Matriz mit einem Array (= Vektor)
		erg = []
		#print ("mul_vektor", self.dimX(), len(v))
		if self.dimX() == len(v):
			for y in range(self.dimY()):
				s = 0
				for k in range(self.dimX()):
					s += self.get(y,k) * v[k]
				erg.append(s)
		return erg

	def getZeros(self):
		# unbestimmbare Spalten ermitteln
		erg = []
		if self.dimX() == self.dimY():
			for y in range(self.dimY()):
				isZero = True
				for x in range(self.dimY()):
					#print ("zero test", self.get(y,x), self.get(x,y))
					isZero &= abs(self.get(y,x)) < 1e-4
					isZero &= abs(self.get(x,y)) < 1e-4
				if isZero:
					erg.append(y)
		return erg
		
	def appendCol(self, col, neg = False):
		# Spalte anhaengen
		dimX = self.dimX()
		for i in range(len(col)):
			self.set(i,dimX,col[i]);
		return self

	def cholesky(self, e):
		l = matrix()
		for k in range(len(self.field)):
			pivot = self.get(k,k)
			if 0 < pivot: 
				l.set(k,k,sqrt(pivot))
				for i in range(k+1,len(self.field)):
					l1 = self.get(i,k)/l.get(k,k)
					l.set(i,k,l1)
					for j in range(k+1,i+1):
						self.set(i,j,self.get(i,j) - l.get(i,k) * l.get(j,k))
		#logging.info(l.writeWMaxima('chol_l'))
		c = []
		for i in range(len(e)):
			s = e[i]
			for j in range(i):
				s -= l.get(i,j) * c[j]
			if abs(l.get(i,i)) < 1e-12: 
				c.append(0.0)
			else:
				c.append(s / l.get(i,i))
		#logging.info(c)
		x = [0] * len(c)
		for i in range(len(e)-1,-1,-1):
			s = c[i]
			for k in range(i+1, len(e)):
				s += l.get(k,i) * x[k]
			if 1e-12 < abs(l.get(i,i)): 
				x[i] = -s / l.get(i,i)
		#logging.info(x)
		return x

	def gauss(self, e):
		#print ("********************** lineare Gleichung")
		# lineare Gleichung nach Gauss
		while len(e) < len(self.field):
			e.append(0)
		# forward-Schritt
		i = 0
		lenField = len(self.field)
		while i < lenField-1:
			#print ("forward",i)
			#self.write()
			#print (e)
			#logging.info (self.writeWMaxima('zf' + str(i)))
			pivotE = self.get(i,i)
			#print "Pivot", pivotE
			if 1e-5 < abs(pivotE):
				# alles wie es sein soll
				for j in range(i+1, len(self.field)):
					#print "gauss, forward", i,j
					pivotJ = self.get(j,i)
					#print "Pivot2", pivotJ
					if 1e-5 < abs(pivotJ):
						faktor = pivotE/pivotJ
						self.mul(j, faktor)
						e[j] *= faktor
						e[j] -= e[i]
						self.sub(j,i)
				i += 1
			else:
				erg = 0
				# suche eine Zeile zum tauschen
				for j in range(i+1, len(self.field)):
					if i < len(self.field[j]) and 1e-5 < abs(self.field[j][i]):
						erg = j
				if 0 < erg:
					self.change(i,erg)
					h = e[i]
					e[i] = e[erg]
					e[erg] = h
				else:
					# Ausnahme, 0-Zeile
					i += 1
		# test real rank
		rank = len(self.field)-1
		test_ranking = True
		while 0 <= rank and test_ranking:
			for i in range(len(self.field[rank])):
				if 1e-5 < abs(self.field[rank][i]):
					test_ranking = False
			if test_ranking:
				rank -= 1
		#logging.info ('rank: ' + str(rank))
		#if 0 < rank and rank < len(self.field)-1:
		#	subMatrixStart = len(self.field)-rank-1
		#	subMatrix = matrix()
		#	for i1 in range(rank+1):
		#		for i2 in range(rank+1):
		#			subMatrix.set(i1, i2, self.get(i1, i2+subMatrixStart))
		#	#logging.info (subMatrix.writeWMaxima('subm'))
		# backward-Schritt
		lin = False
		for i in range(len(self.field)-1,0,-1):
			#print ("backward",i)
			#self.write()
			#print (e)
			#logging.info (self.writeWMaxima('zb' + str(i)))
			pivotE = self.get(i,i)
			if 1e-5 < abs(pivotE):
				for j in range(0,i):
					#print ("gauss, 2backward", i,j)
					pivotJ = self.get(j,i)
					if 1e-5 < abs(pivotJ):
						faktor = pivotE/pivotJ
						self.mul(j, faktor)
						e[j] *= faktor
						e[j] -= e[i]
						self.sub(j,i)

				pivotK = self.get(i,i)
				if 1e-7 < abs(pivotK):
					e[i] /= pivotK
				else:
					e[i] = 0
					#lin = True
		pivotK = self.get(0,0)
		if 1e-7 < abs(pivotK):
			e[0] /= pivotK
		#else:
		#   lin = True
		#print ("e", e)
		#if lin:
		#   e = None
		return e
		
inkscape = ''
inkscapeID = 100

def lineInkscape(line, color="000000", px = 0.1):
	global inkscape, inkscapeID
	inkscape += '<path style="fill:none;stroke:#' + color + ';stroke-width:' + str(px) + 'px;stroke-linecap:round;stroke-linejoin:miter;stroke-opacity:1" d="M'
	for p in line:
		inkscape += ' ' + str(10 * p[0]) + ',' + str(10 * p[1])
	inkscape += '" id="path ' + ("%06u" % inkscapeID) + '"/>\n'
	inkscapeID += 1

		
def linelimits(value, maxValue = 0.5):
	if value < -maxValue:
		value = -maxValue
	if value > maxValue:
		value = maxValue
	return value
		
def calculateDistorsion(coords, d = 1.0):
	a, b, c = 0.0, 0.0, 0.0
	if 1 == 0:
		# testing with wxmaxima
		log("define (r(x,y),sqrt(x*x + y*y));")
		log("define (t(r,g,x1,y1), r * r(x1,y1) / (cos(g) * x1 + sin(g) * y1));")
		log("define (k(a,b,c,d,r,g,x1,y1), a * t(r,g,x1,y1) * t(r,g,x1,y1) * t(r,g,x1,y1) * t(r,g,x1,y1) + b * t(r,g,x1,y1) * t(r,g,x1,y1) * t(r,g,x1,y1) + c * t(r,g,x1,y1) * t(r,g,x1,y1) + d * t(r,g,x1,y1));")
		log("define (t_r(r,g,x1,y1), diff(t(r,g,x1,y1),r));")
		log("define (t_g(r,g,x1,y1), diff(t(r,g,x1,y1),g));")
		log("define (k_a(a,b,c,d,r,g,x1,y1), diff(k(a,b,c,d,r,g,x1,y1),a));")
		log("define (k_b(a,b,c,d,r,g,x1,y1), diff(k(a,b,c,d,r,g,x1,y1),b));")
		log("define (k_c(a,b,c,d,r,g,x1,y1), diff(k(a,b,c,d,r,g,x1,y1),c));")
		log("define (k_d(a,b,c,d,r,g,x1,y1), diff(k(a,b,c,d,r,g,x1,y1),d));")
		log("define (k_r(a,b,c,d,r,g,x1,y1), diff(k(a,b,c,d,r,g,x1,y1),r));")
		log("define (k_g(a,b,c,d,r,g,x1,y1), diff(k(a,b,c,d,r,g,x1,y1),g));")
	#log()
	normalenPoints = [{'r':2.0, 'g':1.0} for co in coords]
	level = 0	
	cleanAddPoints = 0
	j = 0
	while j < 2500:
		j += 1
		#log("\n\n")
		m = matrix()
		t = []
		# sometimes the linear approximation runs into weird locale minima
		# so we start with two points 
		cleanCoords = []
		maxLength = 0
		if level == 0:
			for i1 in range(len(coords)):
				maxLength = max(maxLength, len(coords[i1]))
				cleanSubCoords = []
				co_start = coords[i1][0]
				co_end = coords[i1][len(coords[i1])-1]
				cleanSubCoords.append({'x':co_start['x'], 'y':co_start['y']})
				cleanSubCoords.append({'x':co_end['x'], 'y':co_end['y']})
				cleanCoords.append(cleanSubCoords)
				if j == 1:
					# set start points
					x1,y1 = co_start['x'], co_start['y']
					x2,y2 = co_end['x'], co_end['y']
					dx, dy = x2 - x1, y2 - y1
					rr = dx*dx + dy*dy
					n = (x1 * dy - y1 * dx) / rr
					mx, my = dy * n, -dx * n
					#log("#> " + str(sqrt(mx*mx + my*my)) + ' ' + str(atan2(my, mx)))
					normalenPoints[i1]['g'] = atan2(my, mx)
					normalenPoints[i1]['r'] = sqrt(mx*mx + my*my)
		else:	
			for i1 in range(len(coords)):
				cleanCoords.append(coords[i1][:])
			
		for i1 in range(len(cleanCoords)):
			r,g = normalenPoints[i1]['r'], normalenPoints[i1]['g']
			points = []
			lpoints = []
			for i2 in range(len(cleanCoords[i1])):
				x1, y1 = cleanCoords[i1][i2]['x'], cleanCoords[i1][i2]['y']
				v = sqrt(x1 * x1 + y1 * y1)
				cr = x1 * cos(g) + y1 * sin(g)
				# if cr is 0 the normal is in a 90 degrees angle to the  measuring point, 
				# this is a lokal maxima, we can take it out of the approximation
				tt = r * v / cr
				dr = v / cr
				dg = r*(sin(g)*x1 - cos(g)*y1)*v / (cr*cr)
				#log(str(x) + ' ' + str(y) + ' | ' + str(x1) + ' ' + str(y1) + ' | ' + str(cr) + ' ' + str(v) + ' ' + str(tt))
				ttt = a * tt * tt * tt * tt + b * tt * tt * tt + c * tt * tt + d * tt
				ttt_a = tt * tt * tt * tt
				ttt_b = tt * tt * tt
				ttt_c = tt * tt
				ttt_d = tt
				ttt_r = 4 * a * tt * tt * tt * dr + 3 * b * tt * tt * dr + 2 * c * tt * dr + d * dr
				ttt_g = 4 * a * tt * tt * tt * dg + 3 * b * tt * tt * dg + 2 * c * tt * dg + d * dg
				
				px, py = x1 * ttt / v, y1 * ttt / v
				points.append([px,py])
				lpoints.append([[px,py],[x1,y1]])
				#print(tt, ttt, v)
				row_i = len(t)
				
				if level == 0 or level == 1:
					m.set(row_i, 2*i1 + 0, ttt_r)
					m.set(row_i, 2*i1 + 1, ttt_g)
					t.append(ttt-v)
					#print (x1,y1,tt, tt)
				if level == 2:
					m.set(row_i, 0, ttt_c)
					m.set(row_i, 2*i1 + 1, ttt_r)
					m.set(row_i, 2*i1 + 2, ttt_g)
					t.append(ttt-v)
				if level == 3:
					m.set(row_i, 0, ttt_b)
					m.set(row_i, 1, ttt_c)
					m.set(row_i, 2*i1 + 2, ttt_r)
					m.set(row_i, 2*i1 + 3, ttt_g)
					t.append(ttt-v)
				if level == 4:
					m.set(row_i, 0, ttt_a)
					m.set(row_i, 1, ttt_b)
					m.set(row_i, 2, ttt_c)
					m.set(row_i, 2*i1 + 3, ttt_r)
					m.set(row_i, 2*i1 + 4, ttt_g)
					t.append(ttt-v)
				if level == 5:
					m.set(row_i, 0, ttt_a)
					m.set(row_i, 1, ttt_b)
					m.set(row_i, 2, ttt_c)
					m.set(row_i, 3, ttt_d)
					m.set(row_i, 2*i1 + 4, ttt_r)
					m.set(row_i, 2*i1 + 5, ttt_g)
					t.append(ttt-v)
				
				if 1 == 0:
					# tests with wxmaxima
					log(str(x1) + ',' + str(y1))
					log("k(" + str(a) + "," + str(b) + "," + str(c) + "," + str(d) + "," + str(r) + "," + str(g) + "," + str(x1) + "," + str(y1) + ");")
					log(str(ttt))
					log("k_r(" + str(a) + "," + str(b) + "," + str(c) + "," + str(d) + "," + str(r) + "," + str(g) + "," + str(x1) + "," + str(y1) + ");")
					log("k_g(" + str(a) + "," + str(b) + "," + str(c) + "," + str(d) + "," + str(r) + "," + str(g) + "," + str(x1) + "," + str(y1) + ");")
					log(str(ttt_r) + ',' + str(ttt_g))
					log("k_a(" + str(a) + "," + str(b) + "," + str(c) + "," + str(d) + "," + str(r) + "," + str(g) + "," + str(x1) + "," + str(y1) + ");")
					log("k_b(" + str(a) + "," + str(b) + "," + str(c) + "," + str(d) + "," + str(r) + "," + str(g) + "," + str(x1) + "," + str(y1) + ");")
					log("k_c(" + str(a) + "," + str(b) + "," + str(c) + "," + str(d) + "," + str(r) + "," + str(g) + "," + str(x1) + "," + str(y1) + ");")
					log("k_d(" + str(a) + "," + str(b) + "," + str(c) + "," + str(d) + "," + str(r) + "," + str(g) + "," + str(x1) + "," + str(y1) + ");")
					log(str(ttt_a) + ',' + str(ttt_b) + ',' + str(ttt_c) + ',' + str(ttt_d))
			if j % 1 == 1:
				colors = ["00ff00", "0000ff", "ff0000","00aaaa", "aaaa00", "00aaaa"]
				color = colors[level % len(colors)]
				lineInkscape(points, color)	
				#for lp in lpoints:
				#	lineInkscape(lp, color, 1)	
		if 1 == 1:
			#m.write()	
			#print ('t', t)
			#logging.info (m.writeWMaxima('d'))
			#log (m.writeWMaxima('d'))
			#m.write()
			ma = m.t_copy().mul_matrix(m)
			#logging.info (ma.writeWMaxima('dd'))
			#log ('ma:transpose(d).d;')
			#print ('ma')
			#ma.write()
			#print ('t', t)
			vb = m.t_copy().mul_vektor(t)
			#log ('t:matrix(' + str(t) + ');')
			#log ('transpose(d).t');
			#log ('vb:\n' + str(vb))
			#print ('vb', vb)
			mt = ma.copy().appendCol(vb)
			ee = ma.cholesky(vb)
			#log ("\n" + mt.writeWMaxima('mt'))
			#log ('linsolve(maplist(first,mt.[a,b,c,-1]),[a,b,c]);\n\n')
			#print ('ee: ' + str(ee))
			for ii in range(len(coords)):
				if level == 0 or level == 1:
					normalenPoints[ii]['r'] += linelimits(0.5 * ee[2*ii + 0], 0.5)
					normalenPoints[ii]['g'] += 0.05 * linelimits(ee[2*ii + 1], 0.1)
				if level == 2:
					normalenPoints[ii]['r'] += linelimits(1.0 * ee[2*ii + 1], 0.1)
					normalenPoints[ii]['g'] += 0.05 * linelimits(ee[2*ii + 2])
				if level == 3:
					normalenPoints[ii]['r'] += linelimits(1.0 * ee[2*ii + 2], 0.1)
					normalenPoints[ii]['g'] += 0.05 * linelimits(ee[2*ii + 3])
				if level == 4:
					normalenPoints[ii]['r'] += linelimits(1.0 * ee[2*ii + 3], 0.1)
					normalenPoints[ii]['g'] += 0.05 * linelimits(ee[2*ii + 4])
				if level == 5:
					normalenPoints[ii]['r'] += 0.25 * linelimits(ee[2*ii + 4])
					normalenPoints[ii]['g'] += 0.01 * linelimits(ee[2*ii + 5])
				if normalenPoints[ii]['r'] < 0:
					normalenPoints[ii]['r'] = -normalenPoints[ii]['r']
					normalenPoints[ii]['g'] += pi
				#if normalenPoints[ii]['r'] < 1e-5 and level == 0:
				#	# this catches local minima quite well
				#	normalenPoints[ii]['r'] = 0.01
				#	normalenPoints[ii]['g'] += pi / 6
				if 2 * pi < normalenPoints[ii]['g']:
					normalenPoints[ii]['g'] -= 2*pi
				
			if level == 2:
				c += 0.1 * ee[0]
			if level == 3:
				b += 0.1 * ee[0]
				c += 0.1 * ee[1]
			if level == 4:
				a += 0.1 * ee[0]
				b += 0.1 * ee[1]
				c += 0.1 * ee[2]
			if level == 5:
				a += 0.5 * ee[0]
				b += 0.5 * ee[1]
				c += 0.5 * ee[2]
				d += 0.5 * ee[3]
				
			meanErrorRaw = 0.0
			for e in ee:
				meanErrorRaw += e*e
			tt = 0
			for tti in t:
				tt += tti*tti
			print ("%5u" % j,'level' , level, 'meanError:', "%21.18f" % meanErrorRaw,tt)
			if meanErrorRaw < 1e-8 and level < 4:
				level += 1
			if level == 4 and meanErrorRaw < 1e-24:
				# premature exit
				j = 500000
			
		#print(normalenPoints)
	return a,b,c,d
	
# [{'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, {'x': , 'y': }, ]
coords_contour1080_1 = [
[{'x': 354, 'y': 955}, {'x': 303, 'y': 879}, {'x': 256, 'y': 802}, {'x': 197, 'y': 686}, {'x': 153, 'y': 573}, {'x': 118, 'y': 467}, {'x': 88, 'y': 344}, {'x': 71, 'y': 228}, {'x': 64, 'y': 124}, {'x': 69, 'y': 5}], 
[{'x': 640, 'y': 8}, {'x': 767, 'y': 8}, {'x': 863, 'y': 14}, {'x': 956, 'y': 26}, {'x': 1096, 'y': 55}, {'x': 1217, 'y': 93}, {'x': 1277, 'y': 118}],
[{'x': 441, 'y': 76}, {'x': 584, 'y': 71}, {'x': 715, 'y': 77}, {'x': 835, 'y': 89}, {'x': 939, 'y': 107}, {'x': 1044, 'y': 132}, {'x': 1155, 'y': 167}, {'x': 1276, 'y': 219} ],
[{'x': 3, 'y': 304}, {'x': 185, 'y': 269}, {'x': 305, 'y': 256}, {'x': 469, 'y': 248}, {'x': 653, 'y': 248}, {'x': 840, 'y': 259}, {'x': 982, 'y': 277}, {'x': 1181, 'y': 316}, {'x': 1277, 'y': 355} ],
[{'x': 63, 'y': 657}, {'x': 285, 'y': 665}, {'x': 492, 'y': 661}, {'x': 704, 'y': 647}, {'x': 960, 'y': 616} ],
[{'x': 243, 'y': 933}, {'x': 413, 'y': 940}, {'x': 584, 'y': 930}, {'x': 741, 'y': 909}, {'x': 1026, 'y': 833} ]
]

coords = coords_contour1080_1


image_x, image_y = 1280, 960

# normalize image
minHalf = min(image_x / 2, image_y / 2)
for i1 in range(len(coords)):
	for i2 in range(len(coords[i1])):
		coords[i1][i2]['x'] -= image_x / 2
		coords[i1][i2]['y'] -= image_y / 2
		coords[i1][i2]['x'] /= minHalf
		coords[i1][i2]['y'] /= minHalf
#print (coords)	

a,b,c,d = calculateDistorsion(coords, 1.1)
	
# for inkscape make lines for the pathes	
for i1 in range(len(coords)):
	li = []
	for i2 in range(len(coords[i1])):
		li.append([coords[i1][i2]['x'], coords[i1][i2]['y']])
	lineInkscape(li)
	
print (a,b,c,d)
log("")
#log(inkscape)
log("")
log ("%3.15f" % a)
log ("%3.15f" % b)
log ("%3.15f" % c)
log ("%3.15f" % d)



logFile.close()