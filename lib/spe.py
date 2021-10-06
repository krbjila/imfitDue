#########################################
# JILA KRb -- bare bones SPE reader
# KM 06/05/2018
#########################################

import numpy as np

# Offset from start of file in bytes
headerDefs = {
	'datatype': 108,
	'xdim': 42,
	'ydim': 656,
	'numFrames': 1446,
	'data': 4100
}

# In bytes
headerSizes = {
	'datatype': 2,
	'xdim': 2,
	'ydim': 2,
	'numFrames': 4
}

datatypeDefs = {
	'0': 4,
	'1': 4,
	'2': 2,
	'3': 2,
	'8': 4
}

class SpeFile(object):
	def __init__(self, filename):
		self.filename = filename
		self.initialize()

	def initialize(self):
		with open(self.filename, 'rb') as f:
			self.xdim = self.readBytes(f, headerDefs['xdim'], 2, True)
			self.ydim = self.readBytes(f, headerDefs['ydim'], 2, True)
			self.numFrames = self.readBytes(f, headerDefs['numFrames'], 4, True)
			self.datatype = self.readBytes(f, headerDefs['datatype'], 2, True)
			self.data = self.getData(f, self.xdim, self.ydim, self.numFrames, datatypeDefs[str(self.datatype)])

	# Little endian
	def readBytes(self, f, position, num_bytes, from_beginning=False):
		if from_beginning:
			from_what = 0
		else:
			from_what = 1
		try:
			f.seek(position, from_what)
			val = 0
			for i in range(0, num_bytes):
				val += ord(f.read(1)) << (8 * i)
			return val
		except:
			print("Error reading from file.")

	def getData(self, f, xdim, ydim, frames, data_length):
		data = np.zeros((frames, ydim, xdim))
		stride = xdim * ydim
		f.seek(headerDefs['data'])

		# Changed KM 2/25/2019
		try:
			img = np.fromfile(f, np.uint32, frames*stride)
			return img.reshape((frames, ydim, xdim))
		except Exception as e:
			print("Error reading data from file: {}".format(e))

		# try:
		# 	for i in range(0, frames):
		# 		for j in range(0, ydim):
		# 			for k in range(0, xdim):
		# 				data[i, j, k] = self.readBytes(f, 0, data_length)
		# 	return data
		# except:
		# 	print("Error reading data from file.")
