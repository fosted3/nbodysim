from PIL import Image

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
	
def render(frame, size):
	try:
		data = open(gen_data(frame), 'r')
	except:
		return False
	im = Image.new("RGB", (size, size))
	pix = im.load()
	for line in data:
		temp = get_tuple(line)
		x = clamp(0, int(temp[0] + size/2), size - 1)
		y = clamp(0, int(temp[1] + size/2), size - 1)
		pix[x, y] = (255, 255, 255)
	print("Saving " + gen_image(frame))
	im.save(gen_image(frame))
	return True
	
frame = 0
size = 1024
while (render(frame, size)):
	frame += 1
