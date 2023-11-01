#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import csv
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

filename = 'benchmark.csv'
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = 506
	colnum = 3
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    

for num in range(rownum):
	content[num] = rows[num]
	
y = [0.0 for i in range(21)]
x = ['55%:45%','56%:44%','57%:43%','58%:42%','59%:41%','60%:40%','61%:39%','62%:38%','63%:37%','64%:36%','65%:35%','66%:34%','67%:33%','68%:32%','69%:31%','70%:30%','71%:29%','72%:28%','73%:27%','74%:26%','75%:25%']

for cycle in range(23): #0-22
	for group in range(1,22):  #1-22
		k = cycle*22+group 
		y[group-1] = (float)(content[k][2])
	print(y)
	plt.plot(x, y, c='red', label="performance")
	plt.scatter(x, y, c='red')

plt.title("benchmark")

plt.savefig("benchmark.png",dpi=600,format='png')
plt.show()
plt.close()
