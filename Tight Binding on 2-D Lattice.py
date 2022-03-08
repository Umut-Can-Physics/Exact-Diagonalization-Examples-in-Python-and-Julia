#Libraries
import numpy as np

#Lattice Size
L_x=3
L_y=3

#Create 2-D Lattice. Numbers Represents Each Lattice Sites
lat = np.arange(L_x*L_y).reshape(L_x,L_y)
array = np.append(lat, lat, axis=0)
lattice = array[:-L_x, :]

#Respectively Site Coordinates on Lattice
x_co = np.arange(L_x)
y_co = np.arange(L_y)
arr = []
for j in range(len(x_co)):
    for i in range(len(y_co)):
        arr=np.append(arr, [x_co[i], y_co[j]])
xy = arr.reshape((L_x*L_y,2))

#Find Neighbors Each Sites with Tight-Binding Approximation (But for Hard-Wall Boundary Conditions) 
def HardBC(arr):
    neighbors = {}
    for i in range(len(arr)):
        for j, value in enumerate(arr[i]):
            if i == 0 or i == len(arr) - 1 or j == 0 or j == len(arr[i]) - 1:
                new_neighbors = []
                if i != 0:
                    new_neighbors.append(arr[i - 1][j])  
                if j != len(arr[i]) - 1:
                    new_neighbors.append(arr[i][j + 1]) 
                if i != len(arr) - 1:
                    new_neighbors.append(arr[i + 1][j])  
                if j != 0:
                    new_neighbors.append(arr[i][j - 1])
            else:
                new_neighbors = [
                    arr[i - 1][j],  
                    arr[i][j + 1],  
                    arr[i + 1][j],  
                    arr[i][j - 1]   
                ]
            neighbors[value] = new_neighbors
    return neighbors

#Find Neighbors Each Sites with Tight-Binding Approximation (But for Periodic Boundary Conditions) 
def PerBC(arr):
    neighbors = {}
    for i in range(len(arr)):
        for j, value in enumerate(arr[i]):
            new_neighbors = [
                arr[(i - 1)%L_x][j%L_y],  
                arr[i%L_x][(j + 1)%L_y],  
                arr[(i + 1)%L_x][j%L_y],  
                arr[i%L_x][(j - 1)%L_y]   
            ]
            neighbors[value] = new_neighbors
    return neighbors

#Definiton For Easy Operation Which Show Nearest Neighbors (Tigh-Binding Approximation) Each Sites
HardBCLat = HardBC(lattice)
PerBCLat = PerBC(lattice)

#Libraries
import matplotlib.pyplot as plt

#Lattice Visualization 
fig1, ax1 = plt.subplots()
ax2 = ax1.twiny()
for i in range(L_x):
    for j in range(L_y):
        plt.plot(i, j, 'ro', markersize=7)
ax1.axes.get_xaxis().set_visible(False)        
plt.xticks(x_co)
plt.yticks(y_co)
ax1.invert_yaxis()  
ax1.set_ylabel('Y Coordinates')  
ax2.set_xlabel('X Coordinates')
ax1.grid()
ax2.grid()   

#Find Tight-Binding Neighbors (for Periodic Boundary Conditions) with Lattice Visualization
def PlotPerBC(choose_x, choose_y):
    plt.plot(choose_x, choose_y, 'bo', markersize=7, label='Your Choice')
    plt.plot((choose_x+1)%L_x, choose_y%L_y, 'co', markersize=7, label='Periodic Neighbors')
    plt.plot((choose_x-1)%L_x, choose_y%L_y, 'co', markersize=7)
    plt.plot(choose_x%L_x, (choose_y+1)%L_y, 'co', markersize=7)
    plt.plot(choose_x%L_x, (choose_y-1)%L_y, 'co', markersize=7)
    plt.legend(loc="upper left")
    plt.show()

#Find Tight-Binding Neighbors (for Hard-Wall Boundary Condition) with Lattice Visualization    
def PlotHardBC(choose_x, choose_y):
    plt.plot(choose_x, choose_y, 'bo', markersize=7, label='Your Choice')
    if choose_x==0 or choose_x==L_x-1 or choose_y==0 or choose_y==L_y-1:
        if choose_x!=0:
            plt.plot(choose_x-1, choose_y, 'co', markersize=7)
        if choose_y!=L_y-1:
            plt.plot(choose_x, choose_y+1, 'co', markersize=7)
        if choose_x!=L_x-1:
            plt.plot(choose_x+1, choose_y, 'co', markersize=7)
        if choose_y!=0:
            plt.plot(choose_x, choose_y-1, 'co', markersize=7)
    else:
        plt.plot(choose_x+1, choose_y, 'co', markersize=7)
        plt.plot(choose_x-1, choose_y, 'co', markersize=7)
        plt.plot(choose_x, choose_y+1, 'co', markersize=7)
        plt.plot(choose_x, choose_y-1, 'co', markersize=7)
