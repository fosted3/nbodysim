import os
import sys
import queue
import threading
from subprocess import call

def clamp(a, x, b):
	if (x < a):
		return a
	if (x > b):
		return b
	return x

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

def worker():
	while True:
		data = q.get()
		call(["python", "src/draw.py", str(data[0]), str(data[1])])
		q.task_done()

if (len(sys.argv) != 2):
	print("Usage: python " + sys.argv[0] + " config")
	raise sys.exit()
config = open(sys.argv[1], 'r')
threads = 2
num_frames = -1
work = []
for line in config:
	pair = line.split()
	if (len(pair) == 2):
		if (pair[0] == "threads"):
			threads = int(pair[1])
		if (pair[0] == "num_frames"):
			num_frames = int(pair[1])
config.close()
assert(num_frames != -1)
q = queue.Queue()
for frame in range(0, num_frames):
	if (os.path.isfile(gen_data(frame)) and not os.path.isfile(gen_image(frame))):
		work.append((sys.argv[1], frame))
for i in range(clamp(0, threads, len(work))):
     t = threading.Thread(target=worker)
     t.daemon = True
     t.start()
while True:
	for item in work:
		q.put(item)
	work = []
	q.join()
	for frame in range(0, num_frames):
		if (os.path.isfile(gen_data(frame)) and not os.path.isfile(gen_image(frame))):
			work.append((sys.argv[1], frame))
	if (len(work) == 0):
		break
