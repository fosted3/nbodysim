from PIL import Image
from math import sqrt
from numpy import zeros
import os.path
import sys

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

#def inc_tuple(tup, inc):
#	return (clamp(0, tup[0] + inc, 255), clamp(0, tup[1] + inc, 255), clamp(0, tup[2] + inc, 255))

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
	
def render(frame, img_w, img_h, inc, proj, scale):
	if (os.path.isfile(gen_image(frame))):
		return True
	try:
		data = open(gen_data(frame), 'r')
	except:
		return False
	print("Rendering " + gen_data(frame))
	if (proj == "front" or proj == "iso"):
		im = Image.new("RGB", (img_w, img_h))
		arr = zeros((img_w, img_h))
	else:
		print("Unknown projection " + proj)
		raise sys.exit()
	if (proj == "iso"):
		scale *= 2.0
	pix = im.load()
	for line in data:
		temp = get_tuple(line)
		if (proj == "front"):
			x = temp[0] + img_w/2
			y = temp[1] + img_h/2
		elif (proj == "iso"):
			x = ((sqrt(3) / 2.0) * (temp[0] - temp[1]) + img_w) / 2.0
			y = ((-0.5) * (temp[0] + temp[1] + 2 * temp[2]) + img_h) / 2.0
		x += (x - (img_w / 2)) * (scale - 1)
		y += (y - (img_h / 2)) * (scale - 1)
		x = int(clamp(0, x, img_w - 1))
		y = int(clamp(0, y, img_h - 1))
		arr[x, y] += inc
	for x in range(0, img_w):
		for y in range(0, img_h):
			v = int(clamp(0, arr[x, y], 255))
			pix[x, y] = (v, v, v)
	print("Saving " + gen_image(frame))
	im.save(gen_image(frame))
	return True
	
if (len(sys.argv) != 2 and len(sys.argv) != 3):
	print("Usage: python " + sys.argv[0] + " config [frame]")
	raise sys.exit()

frame = 0
img_w = -1
img_h = -1
scale = 1
inc = 255
proj = "iso"
conf = sys.argv[1]
config = open(conf, 'r')
for line in config:
	pair = line.split()
	if (len(pair) == 2):
		if (pair[0] == "img_w"):
			img_w = int(pair[1])
		if (pair[0] == "img_h"):
			img_h = int(pair[1])
		if (pair[0] == "projection"):
			proj = pair[1]
		if (pair[0] == "brightness"):
			inc = float(pair[1])
		if (pair[0] == "scale"):
			scale = float(pair[1])
if (len(sys.argv) == 3):
	render(int(sys.argv[2]), img_w, img_h, inc, proj, scale)
else:
	while (render(frame, img_w, img_h, inc, proj, scale)):
		frame += 1
