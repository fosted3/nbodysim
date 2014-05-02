from PIL import Image
from math import sqrt
import os.path

def clamp(a, x, b):
	if (x < a):
		return a
	if (x > b):
		return b
	return x

def get_tuple(line):
	temp = "".join(line.split(" ")).split(",")
	assert(len(temp) == 3)
	temp[0] = float(temp[0])
	temp[1] = float(temp[1])
	temp[2] = float(temp[2])
	return tuple(temp)

def inc_tuple(tup, inc):
	return (clamp(0, tup[0] + inc, 255), clamp(0, tup[1] + inc, 255), clamp(0, tup[2] + inc, 255))

def gen_data(frame):
	fname = str(frame)
	while (len(fname) < 4):
		fname = "0" + fname
	fname = "data/" + fname + ".txt"
	return fname

def gen_image(frame):
	fname = str(frame)
	while (len(fname) < 4):
		fname = "0" + fname
	fname = "img/" + fname + ".png"
	return fname
	
def render(frame, size, inc):
	try:
		data = open(gen_data(frame), 'r')
	except:
		return False
	im = Image.new("RGB", (size, size))
	pix = im.load()
	i = 0
	for line in data:
		temp = get_tuple(line)
		x = clamp(0, int(temp[0] + size/2), size - 1)
		y = clamp(0, int(temp[1] + size/2), size - 1)
		pix[x, y] = inc_tuple(pix[x, y], inc)
	print("Saving " + gen_image(frame))
	im.save(gen_image(frame))
	return True
	
def render_iso(frame, size, inc):
	if (os.path.isfile(gen_image(frame))):
		return True
	try:
		data = open(gen_data(frame), 'r')
	except:
		return False
	im = Image.new("RGB", (size * 2, size * 2))
	pix = im.load()
	i = 0
	for line in data:
		temp = get_tuple(line)
		x = clamp(0, int((sqrt(3) / 2.0) * (temp[0] - temp[1]) + size), (size * 2) - 1)
		y = clamp(0, int((-0.5) * (temp[0] + temp[1] + 2 * temp[2]) + size), (size * 2) - 1)
		pix[x, y] = inc_tuple(pix[x, y], inc)
	print("Saving " + gen_image(frame))
	im.save(gen_image(frame))
	return True
	
frame = 0
size = 1024
inc = 2
while (render_iso(frame, size, inc)):
	frame += 1
