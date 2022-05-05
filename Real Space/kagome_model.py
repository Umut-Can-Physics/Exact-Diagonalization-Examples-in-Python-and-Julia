import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
import math
import matplotlib.patches as patches

nrows = 1
ncols = 1
fig, axs = plt.subplots(nrows, ncols, sharex=True, figsize=(14, 8))
txt_size = 22

pi = math.pi
sin = math.sin
cos = math.cos
sqrt = math.sqrt
exp = math.exp

ax = 1
ay = 0
bx = cos(pi / 3)
by = sin(pi / 3)

#primitive vectors
a_vec = np.array([1,0])
b_vec = np.array([cos(pi / 3),sin(pi / 3)])

x = [0,1,cos(pi/3),0]
y = [0,0,sin(pi/3),0]

#Kagome lattice size
lx = 15
ly = 15

#Empty coordinates to create coordinates array
co_x = []
co_y= []
xp = np.full((len(x), 1), np.nan)
yp = np.full((len(x), 1), np.nan)
for ix in range(lx):
    for iy in range(ly):
        for ip in range(len(x)):
            xp[ip] = x[ip] + ax * ix * 2 + bx * iy * 2
            yp[ip] = y[ip] + ay * ix * 2 + by * iy * 2
            #Coordinates array
            co_x = np.append(co_x, xp.item(ip))
            co_y = np.append(co_y, yp.item(ip))
            #plot coordinate labels
            # for i_x, i_y in zip(co_x, co_y):
            #     plt.text(i_x, i_y,'({}, {})'.format(np.round(i_x,3), np.round(i_y,3)))
        axs.plot(xp, yp, marker=' ', color='crimson', linestyle=' ')
        axs.plot(xp, yp, marker=' ', color='blue', linestyle='-')
        xp1 = [xp[1], xp[1] + ax]
        yp1 = [yp[1], yp[1] + ay]
        axs.plot(xp1, yp1, marker=' ', color='blue', linestyle='--')
        xp2 = [xp[2], xp[2] + bx]
        yp2 = [yp[2], yp[2] + by]
        axs.plot(xp2, yp2, marker=' ', color='blue', linestyle='--')
        xp3 = [xp[2], xp[2] - bx]
        yp3 = [yp[2], yp[2] + by]
        axs.plot(xp3, yp3, marker=' ', color='blue', linestyle='--')

#Make coordinates array
points = np.column_stack((co_x,co_y))
points = np.vstack({tuple(e) for e in points})
#sort according to coordinate
ind = np.lexsort((points[:,1],points[:,0]))
points = points[ind]

#Finding nearest neighboors with scipy 
tree = spatial.cKDTree(points)
#tree.query_ball_point(([0.0,0.0]), 1.1)

#Show primite vectors
style = "Simple, tail_width=3, head_width=10, head_length=15"
kw = dict(arrowstyle=style, color="black")
arrow_origin = np.array(points[0])
arroww=patches.FancyArrowPatch(arrow_origin,np.array([points[1,0],points[1,1]]), connectionstyle="arc3,rad=0", **kw)
arrowww=patches.FancyArrowPatch(arrow_origin,np.array([points[2,0],points[2,1]]), connectionstyle="arc3,rad=0", **kw)
axs.add_patch(arroww)
axs.add_patch(arrowww)

#Plot site labels
# order = np.arange(len(points))
# fig3, ax3 = plt.subplots()
# ax3.set_title('Site Labels')
# ax3.scatter(points[:,0], points[:,1], s=50)
# for i, txt in enumerate(order):
#     ax3.annotate(txt, (points[i,0], points[i,1]), fontsize=15)
    
#Show hopping positions with NN and NNN of specific site 
order = np.arange(len(points))
for i, txt in enumerate(order):
    axs.annotate(txt, (points[i,0], points[i,1]), fontsize=15)
n = 41
style = "Simple, tail_width=0.5, head_width=4, head_length=8"
kw = dict(arrowstyle=style, color="orange")
kww = dict(arrowstyle=style, color="magenta")
for j in tree.query_ball_point([points[n,0],points[n,1]], 1.1):
    arrow=patches.FancyArrowPatch(points[n],points[j], connectionstyle="arc3,rad=.2", **kw, label='NN')
    axs.add_patch(arrow)
a = tree.query_ball_point([points[n,0],points[n,1]], 1.1)
b = tree.query_ball_point([points[n,0],points[n,1]], 1.75)
c = np.intersect1d(a,b)
d = np.setdiff1d(b,c)
for j in d:
    arrow=patches.FancyArrowPatch(points[n],points[j], connectionstyle="arc3,rad=.2", **kww, label='NNN')
    axs.add_patch(arrow)

#Define to hopping phase, we need to pick important points
#Pick Triangle-2 Points
A = np.zeros((lx,2))
for i in range(1,lx):
    A[i] = 2*b_vec+A[i-1]
B = np.empty((0, 2))
for j in range(1,ly):
    B = np.vstack((B, A + j*2*a_vec))
triangle_2 = np.vstack((A,B))
axs.scatter(triangle_2[:,0], triangle_2[:,1], s=100, c='red', label="B")

#Pick Triangle-1 Points
C = np.zeros((lx+1,2))
for i in range(1,lx+1):
    C[i] = b_vec+A[i-1]
C = np.delete(C,0,0)
D = np.empty((0, 2))
for j in range(1,ly):
    D = np.vstack((D, C + j*2*a_vec))
triangle_1 = np.vstack((C,D))
axs.scatter(triangle_1[:,0], triangle_1[:,1], s=100, c='lime', label="A")

#Pick Triangle-3 Points
E = np.zeros((lx+1,2))
for i in range(1,lx+1):
    E[i] = a_vec+A[i-1]
E = np.delete(E,0,0)
F = np.empty((0, 2))
for j in range(1,ly):
    F = np.vstack((F, E + j*2*a_vec))
triangle_3 = np.vstack((E,F))
axs.scatter(triangle_3[:,0], triangle_3[:,1], s=100, c='blue', label="C")

#To stop repeating legends
handles, labels = fig.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys())

#Define Real Space Kagome Hamiltonian Matrix
t1 = -1;L1 = 0.28;t2 = 0.3;L2 = 0.2
#-------
#t1=-1;L1 = 0;t2 = 0;L2 = 0

#Function to search any points in triangles
def pointri(triangle, point):
    _triangle_a = np.array([round(i[0]) for i in triangle])
    _triangle_a = np.where(_triangle_a == round(point[0]))

    _triangle_b = triangle[_triangle_a]
    _triangle_b = [i[1] for i in _triangle_b if i[1] == point[1]]

    if len(_triangle_b) == 1:
        return True

    else:
        return False

def H():
    H = np.zeros((len(points), len(points)), dtype=complex)
    for m in range(len(points)): 
        for n in range(len(points)):
            if m!=n:
                if m in tree.query_ball_point([points[n,0],points[n,1]], 1.1):
                    
                    #A-->B
                    if pointri(triangle_1,points[m]):
                        if pointri(triangle_2,points[n]):
                            H[m][n] = t1+1j*L1
                            #print("[NN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                            
                    #A-->C
                    if pointri(triangle_1,points[m]):
                        if pointri(triangle_3,points[n]):
                            H[m][n] = t1-1j*L1
                            #print("[NN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                            
                    #B-->A
                    if pointri(triangle_2,points[m]):
                        if pointri(triangle_1,points[n]):
                            H[m][n] = t1-1j*L1
                            #print("[NN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                    
                    #B-->C
                    if pointri(triangle_2,points[m]):
                        if pointri(triangle_3,points[n]):
                            H[m][n] = t1+1j*L1
                            #print("[NN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                            
                    #C-->A
                    if pointri(triangle_3,points[m]):
                        if pointri(triangle_1,points[n]):
                            H[m][n] = t1+1j*L1
                            #print("[NN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                    
                    #C-->B
                    if pointri(triangle_3,points[m]):
                        if pointri(triangle_2,points[n]):
                            H[m][n] = t1-1j*L1
                            #print("[NN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                    
                a = tree.query_ball_point([points[n,0],points[n,1]], 1.1)
                b = tree.query_ball_point([points[n,0],points[n,1]], 1.75)
                c = np.intersect1d(a,b)
                d = np.setdiff1d(b,c)
                if m in d: 
                    
                    #A-->B
                    if pointri(triangle_1,points[m]):
                        if pointri(triangle_2,points[n]):
                            H[m][n] = t2-1j*L2
                            #print("[NNN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                    
                    #A-->C
                    if pointri(triangle_1,points[m]):
                        if pointri(triangle_3,points[n]):
                            H[m][n] = t2+1j*L2
                            #print("[NNN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                    
                    #B-->A
                    if pointri(triangle_2,points[m]):
                        if pointri(triangle_1,points[n]):
                            H[m][n] = t2+1j*L2
                            #print("[NNN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                    
                    #B-->C
                    if pointri(triangle_2,points[m]):
                        if pointri(triangle_3,points[n]):
                            H[m][n] = t2-1j*L2
                            #print("[NNN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                    
                    #C-->A
                    if pointri(triangle_3,points[m]):
                        if pointri(triangle_1,points[n]):
                            H[m][n] = t2-1j*L2
                            #print("[NNN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
                    
                    #C-->B
                    if pointri(triangle_3,points[m]):
                        if pointri(triangle_2,points[n]):
                            H[m][n] = t2+1j*L2
                            #print("[NNN] "+str(m)+" --> "+str(n)+" : "+str(H[m][n]))
    return H         
# np.max(np.abs((H()-np.conj(np.transpose(H())))))
# fig4, ax4 = plt.subplots()
# y_eig = np.zeros(len(points))
# for i in range(len(points)):
#     y_eig[i]=0
# x_eig, x_func = np.linalg.eig(H())
# ax4.set_title('Real Space Kagome EigenValues')
# ax4.plot(y_eig, x_eig, 'rx', markersize=7)

x_eig, x_func = np.linalg.eig(H())
fig5, ax5 = plt.subplots()
idx = np.argsort(x_eig)
sorted = x_eig[idx]
x_sort = np.arange(0, len(points), 1)
y_sort = sorted
ax5.set_xlabel("Sıralanmış Öz Değerler")
ax5.set_ylabel(r'$Enerjiler$')
ax5.set_title(r'$N_{site}= $'+str(len(points))+'\n'+r'$t_1: $'+str(t1)+' , '+r'$\lambda_1: $'+str(L1)+' , '+r'$t_2: $'+str(t2)+' , '+r'$\lambda_2: $'+str(L2))
ax5.plot(x_sort, y_sort, 'ro', markersize=0.7)

axs.set_xlabel(r'$x$', fontsize=txt_size)
axs.set_ylabel(r'$y$', fontsize=txt_size)
axs.set_title('Kagome Örgü Modeli')
axs.tick_params(axis='both', which='major', labelsize=txt_size, length=10)
plt.show()
