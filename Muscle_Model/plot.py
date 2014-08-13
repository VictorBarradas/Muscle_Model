import matplotlib.pyplot as plt
from pylab import show, plot, grid

with open("test.txt", "r") as file:
    data = file.readlines()

sample_list = []
T_1_list = []
T_2_list = []
T_3_list = []
T_4_list = []
T_5_list = []
T_6_list = []

for line in data:
    sample, T_1, T_2, T_3, T_4, T_5, T_6 = line.split()
    T_1_list.append(T_1)
    T_2_list.append(T_2)
    T_3_list.append(T_3)
    T_4_list.append(T_4)
    T_5_list.append(T_5)
    T_6_list.append(T_6)

plt.plot(T_1_list,'r', label='5Hz')
plt.plot(T_2_list,'g', label='10Hz')
plt.plot(T_3_list,'b', label='20Hz')
plt.plot(T_4_list,'b', label='50Hz')
plt.plot(T_5_list,'b', label='60Hz')
plt.plot(T_6_list,'b', label='100Hz')
    
plt.legend(loc='upper right', numpoints=1)
grid('on')
plt.show()	
    