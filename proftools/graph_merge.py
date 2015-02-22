import numpy as np
from matplotlib import pyplot as plt




lines = [line.strip() for line in open("../build/merge_profiling_try_1.txt")]

particle_num = np.zeros((len(lines), 2), dtype='uint32');
merge_num = np.zeros((len(lines), 2));
split_num = np.zeros((len(lines), 2));


index = 0

for data_str in lines:

	data_str = data_str.split()

	# Setting iteration
	particle_num[index][0] = int(data_str[0])
	merge_num[index][0] = int(data_str[0])
	split_num[index][0] = int(data_str[0])

	particle_num[index][1] = int(data_str[1])
	merge_num[index][1] = int(data_str[2])
	split_num[index][1] = int(data_str[3])

	index+=1




plt.plot(particle_num[:,0:1].reshape( (len(lines)) ) , particle_num[:,1:2].reshape( (len(lines)) ), label = "particles")
plt.plot(merge_num[:,0:1].reshape( (len(lines)) ) , merge_num[:,1:2].reshape( (len(lines)) ), label = "merges")
plt.plot(split_num[:,0:1].reshape( (len(lines)) ) , split_num[:,1:2].reshape( (len(lines)) ), label = "splits")
plt.legend()
plt.show()







