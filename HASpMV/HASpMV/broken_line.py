#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import csv
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

filename = 'by_cache.csv'
with open(filename) as csvfile:
	reader = csv.reader(csvfile)
	rows = [row for row in reader]    
	rownum = 391
	colnum = 3
	content = [[0.0 for i in range(colnum)] for i in range(rownum+1)]    

for num in range(rownum):
	content[num] = rows[num]
	
y = [0.0 for i in range(16)]
x = ['65%:35%','66%:34%','67%:33%','68%:32%','69%:31%','70%:30%','71%:29%','72%:28%','73%:27%','74%:26%','75%:25%','76%:24%','77%:23%','78%:22%','79%:21%','80%:20%']
a = [0.0 for i in range(18)]
b = [0.0 for i in range(18)]

fig = plt.figure(figsize=(20, 15))
plt.grid(color='r', ls = '-.', lw = 0.2)
yy = range(0,65,5)

for cycle in range(18): #0-22
	max = 0
	for group in range(1,17):  #1-22
		k = cycle*17+group 
		y[group-1] = (float)(content[k][2])
		if y[group-1]-max>=1e-8:
			max = y[group-1]
			max_p =  x[group-1]
	print(y)
	#print(max,max_p)
	a[cycle] = max_p
	b[cycle] = max
	t = content[cycle*17][0]
	plt.text(15,y[15] , t, ha='center', rotation=0, wrap=False)
	plt.plot(x, y, c='red')
	plt.scatter(x, y, c='red',s = 10)

plt.scatter(a, b, c='black',s = 30)

plt.title("cache")
plt.xticks(fontsize = 9)
plt.yticks(yy,fontsize = 10)
plt.xlabel("P:E nnz proportion")
plt.ylabel("Performance (Gflops)")
fig.savefig("benchmark.png",dpi=600,format='png')
plt.show()
plt.close()
