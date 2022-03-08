#Libraries
import numpy as np

#Square Lattice Size
L_x=L_y=5

#q Values Have To Be Eaxct Multiple of Lattice Size
q = L_x
p=1

#Magnetic Flux Per Plaquet
alpha=p/q

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
PerBCLat = PerBC(lattice)

#Real Space Hamiltonian with Twisted Angle Phases
def HMat_Theta(alpha, theta_x, theta_y):
    H = np.zeros((L_x*L_y, L_x*L_y), dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
            if m in PerBCLat[n]:

                if np.absolute(xy[m][0]-xy[n][0])==L_x-1:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*alpha*xy[m][1])*np.exp(-1j*theta_x)
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*alpha*xy[m][1])*np.exp(1j*theta_x)

                elif np.absolute(xy[m][1]-xy[n][1])==L_y-1:
                    if xy[m][1] > xy[n][1]:
                        H[m][n] = -np.exp(-1j*theta_y)
                    elif xy[m][1] < xy[n][1]:
                        H[m][n] = -np.exp(1j*theta_y)

                else:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*alpha*xy[m][1])
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*alpha*xy[m][1])
                    else:
                        H[m][n] = -np.exp(0)
        
    return H

#Twisted Angle Space Does'nt Depend on Lattice Space
theta_size=10
dx = 2*np.pi/theta_size
dy = dx
theta_x = np.arange(0, 2*np.pi, dx)
theta_y = np.arange(0, 2*np.pi, dy)

#Calculate Chern Numbes
chern_array = []
for i in np.arange(0, L_x*L_y, q):
    j=i+q
    S=0
    n1=i
    n2=j
    #Algorithm of Calculate 
    for t_x in range(0, len(theta_x)):
        for t_y in range(0, len(theta_y)):
            w1, v1 = np.linalg.eig(HMat_Theta(alpha, theta_x[t_x], theta_y[t_y]))
            idx1 = np.argsort(w1)
            v1_sorted = v1[:,idx1]
            v11 = v1_sorted[:,n1:n2]
            w2, v2 = np.linalg.eig(HMat_Theta(alpha, theta_x[t_x]+dx, theta_y[t_y]))
            idx2 = np.argsort(w2)
            v2_sorted = v2[:,idx2]
            v22 = v2_sorted[:,n1:n2]
            w3, v3 = np.linalg.eig(HMat_Theta(alpha, theta_x[t_x], theta_y[t_y]+dy))
            idx3 = np.argsort(w3)
            v3_sorted = v3[:,idx3]
            v33 = v3_sorted[:,n1:n2]
            w4, v4 = np.linalg.eig(HMat_Theta(alpha, theta_x[t_x]+dx, theta_y[t_y]+dy))
            idx4 = np.argsort(w4)
            v4_sorted = v4[:,idx4]
            v44 = v4_sorted[:,n1:n2]
            U1 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v11)), v22)) #U_x
            U1 = U1 / np.absolute(U1)
            U2 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v22)), v44)) #U_y dx
            U2 = U2 / np.absolute(U2)
            U3 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v33)), v44)) #U_x dy
            U3 = U3 / np.absolute(U3)
            U4 = np.linalg.det(np.matmul(np.conjugate(np.transpose(v11)), v33)) #U_y
            U4 = U4 / np.absolute(U4)
            F = np.log(U1*U2*1/U3*1/U4)
            S = S+F

    C = 1/(2*np.pi*1j)*S
    C.real
    chern_array.append(C.real)
Chern_Numbers = np.array(chern_array)
print('Chern Numbers of q='+str(q)+'\n Seperation of Bands in Real Space: '+str(Chern_Numbers))
