##############################################################################
### t1=1;L1=0.28;t2=0.3;L2=0.2 iken matris enerjileri 4'ün üstüne çıkıyor! ###
##############################################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
import matplotlib.patches as patches

# constants
pi = np.pi;sin = np.sin;cos = np.cos;sqrt = np.sqrt;exp = np.exp
ax = 1;ay = 0;bx = cos(pi / 3);by = sin(pi / 3)

# primitive vectors (DEĞİŞTİRME!)
a_vec = np.array([ax,ay])
b_vec = np.array([bx,by])
#a_vec=2*a_vec;b_vec=2*b_vec

# Kagome lattice size
lx = 3
ly = 3

# Plot kagome lattice
co_x = []
co_y= []
x = [0,1,cos(pi/3),0];y = [0,0,sin(pi/3),0]
xp = np.full((len(x), 1), np.nan)
yp = np.full((len(x), 1), np.nan)
fig, axs = plt.subplots()
for ix in range(lx):
    for iy in range(ly):
        for ip in range(len(x)):
            xp[ip] = x[ip] + ax * ix * 2 + bx * iy * 2
            yp[ip] = y[ip] + ay * ix * 2 + by * iy * 2
            #Coordinates array
            co_x = np.append(co_x, xp.item(ip))
            co_y = np.append(co_y, yp.item(ip))
        axs.plot(xp, yp, marker=' ', color='brown', linestyle=' ')
        axs.plot(xp, yp, marker=' ', color='brown', linestyle='-')
        xp1 = [xp[1], xp[1] + ax]
        yp1 = [yp[1], yp[1] + ay]
        axs.plot(xp1, yp1, marker=' ', color='brown', linestyle='--')
        xp2 = [xp[2], xp[2] + bx]
        yp2 = [yp[2], yp[2] + by]
        axs.plot(xp2, yp2, marker=' ', color='brown', linestyle='--')
        xp3 = [xp[2], xp[2] - bx]
        yp3 = [yp[2], yp[2] + by]
        axs.plot(xp3, yp3, marker=' ', color='brown', linestyle='--')
axs.set_xlabel(r'$x$')
axs.set_ylabel(r'$y$')
axs.set_title('2D Kagome Lattice \n'+str(lx)+r'$\times$'+str(ly))
axs.tick_params(axis='both', which='major', labelsize=2, length=10)
points = np.column_stack((co_x,co_y))
points = np.vstack({tuple(e) for e in points})
ind = np.lexsort((points[:,1],points[:,0]))
points = points[ind]
order = np.arange(len(points))
for i, txt in enumerate(order):
    axs.annotate(txt, (points[i,0], points[i,1]), fontsize=15)

# Finding nearest neighboors with scipy 
tree = spatial.cKDTree(points)

# Show primite vectors
plt.annotate(r'$\vec{a}_1$', xy=(0, 0),xytext=(2, -1/2), horizontalalignment="center", fontsize=15)
plt.annotate(r'$\vec{a}_2$', xy=(0, 0),xytext=(1/2,sqrt(3)), horizontalalignment="center", fontsize=15)
style = "Simple, tail_width=3, head_width=10, head_length=15"
kw = dict(arrowstyle=style, color="black")
arrow_origin = np.array(points[0])
arroww=patches.FancyArrowPatch(arrow_origin,np.array([2*a_vec[0],2*a_vec[1]]), connectionstyle="arc3,rad=0", **kw, label='Primitive Vectors')
arrowww=patches.FancyArrowPatch(arrow_origin,np.array([2*b_vec[0],2*b_vec[1]]), connectionstyle="arc3,rad=0", **kw)
axs.add_patch(arroww)
axs.add_patch(arrowww)
    
# Show hard-wall hopping directions
n = 2 # site index
style = "Simple, tail_width=0.5, head_width=4, head_length=8"
kw = dict(arrowstyle=style, color="orange")
kww = dict(arrowstyle=style, color="magenta")
for j in tree.query_ball_point([points[n,0],points[n,1]], 1.1):
    arrow=patches.FancyArrowPatch(points[n],points[j], connectionstyle="arc3,rad=.2", **kw, label='Nearest-Neighbors Hopping')
    axs.add_patch(arrow)
a = tree.query_ball_point([points[n,0],points[n,1]], 1.1)
b = tree.query_ball_point([points[n,0],points[n,1]], 1.75)
c = np.intersect1d(a,b)
d = np.setdiff1d(b,c)
for j in d:
    arrow=patches.FancyArrowPatch(points[n],points[j], connectionstyle="arc3,rad=.2", **kww, label='Next-Nearest-Neighbors Hopping')
    axs.add_patch(arrow)
    
###############################################
## ALGORTIHM FOR PERIODIC BOUNDARY CONDITION ##
###############################################
# sanal örgüdeki noktaların koordinatları
move_points=np.zeros((0,2),int)
for i in points:
    move_points=np.append(move_points,
        np.array([
        i+lx*2*a_vec, # ana örgüdeki noktaları sağa ötele
        i-lx*2*a_vec, # ana örgüdeki noktaları sola ötele
        i+ly*2*b_vec, # ana örgüdeki noktaları yukarı ötele
        i-ly*2*b_vec, # ana örgüdeki noktaları aşağı ötele
        i+lx*2*a_vec+ly*2*b_vec, # ana örgüdeki noktaları sağ üst çapraza ötele
        i-lx*2*a_vec+ly*2*b_vec, # ana örgüdeki noktaları sol üst çapraza ötele
        i-lx*2*a_vec-ly*2*b_vec, # ana örgüdeki noktaları sol alt çapraza ötele
        i+lx*2*a_vec-ly*2*b_vec, # ana örgüdeki noktaları sağ alt çapraza ötele
        ])
        ,axis=0)
# plot sanal örgü noktaları
# axs.scatter(move_points[:,0],move_points[:,1], label="Ghost Sites", c='gray')  

# vectorize (hızlı ama işime yaramıyor)
# move_points_2=np.zeros((0,2))
# move_points_2=np.append(move_points_2, np.array(points+lx*2*a_vec), axis=0)
# move_points_2=np.append(move_points_2, np.array(points-lx*2*a_vec), axis=0)
# move_points_2=np.append(move_points_2, np.array(points+ly*2*b_vec), axis=0)
# move_points_2=np.append(move_points_2, np.array(points-ly*2*b_vec), axis=0)
# move_points_2=np.append(move_points_2, np.array(points+lx*2*a_vec+ly*2*b_vec), axis=0)
# move_points_2=np.append(move_points_2, np.array(points-lx*2*a_vec+ly*2*b_vec), axis=0)
# move_points_2=np.append(move_points_2, np.array(points-lx*2*a_vec-ly*2*b_vec), axis=0)
# move_points_2=np.append(move_points_2, np.array(points+lx*2*a_vec-ly*2*b_vec), axis=0)

# select A,B and C atoms
# atom B
A = np.zeros((lx*3,2))
for i in range(len(A)):
    A[i,0]=move_points[6,0];A[i,1]=move_points[6,1]
for i in range(1,lx*3):
    A[i] = 2*b_vec+A[i-1]
B = np.empty((0, 2))
for j in range(1,ly*3):
    B = np.vstack((B, A + j*2*a_vec))
triangle_2 = np.vstack((A,B))
axs.scatter(triangle_2[:,0], triangle_2[:,1], s=50, c='red', alpha=.5, label="B")
# atom A
C = np.zeros((lx*3,2))
for i in range(len(C)):
    C[i,0]=move_points[6,0];C[i,1]=move_points[6,1]
for i in range(1,lx*3):
    C[i] = b_vec+A[i-1]
C = np.delete(C,0,0)
D = np.empty((0, 2))
for j in range(1,ly*3):
    D = np.vstack((D, C + j*2*a_vec))
triangle_1 = np.vstack((C,D))
axs.scatter(triangle_1[:,0], triangle_1[:,1], s=50, c='lime',alpha=.5, label="A")
# atom C
E = np.zeros((lx*3,2))
for i in range(len(C)):
    E[i,0]=move_points[6,0];E[i,1]=move_points[6,1]
for i in range(1,lx*3):
    E[i] = a_vec+A[i-1]
E = np.delete(E,0,0)
F = np.empty((0, 2))
for j in range(1,ly*3):
    F = np.vstack((F, E + j*2*a_vec))
triangle_3 = np.vstack((E,F))
axs.scatter(triangle_3[:,0], triangle_3[:,1], s=50, c='blue',alpha=.5, label="C")

# To stop repeating legends
handles, labels = fig.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys())


# for i in range(0,len(move_points),8):
#     n1=i
#     n2=i+8
#     move_points[n1:n2,:] #taşınmış noktaların her seferinde ilk 8 elememanı sırasıyla sanal
#     #0,1,2,... sitelerini oluşturur.
# move_points[0:8,:] #0'ın sanal site koordinatları
# move_points[8:16,:] #1'in sanal site koordinatları
# move_points[16:32,:] #2'nin sanal site koordinatları
# move_points[32:64,:] #3'ün sanal site koordinatları

# j. sitenin sadece sanal site koordinatları
def per_neighbors(j):
    # (n1,n2) çifti için 
    n1n2 = np.zeros((0,2),int)
    for i in range(0,len(move_points),8):
        n1n2=np.append(n1n2,np.array([[i,i+8]]),axis=0) # bu dizinin birinci/ikinci elemanı n1/n2
    # satırları ise site numarasının n1 ve n2 elemanlarını göstermektedir.
    j=move_points[n1n2[j][0]:n1n2[j][1],:] # j. sitenin sanal site koordinatları
    return j
        
# j. sitenin NN komşulukları tüm koordinatları
abcd=np.append(points,move_points,axis=0)
treee=spatial.cKDTree(abcd)
def coordinates_NN(j):
    co=np.zeros((0,2),int)
    for i in treee.query_ball_point([abcd[j,0],abcd[j,1]], 1.1):
        co=np.append(co,[abcd[i]],axis=0)
    return co

# j. sitenin periyodik NN komşulukları
def find_per_NearestNeighbors(j):
    hopping_NN=np.zeros((0,1),int)
    hopping_NN=np.append(hopping_NN,tree.query_ball_point([points[j,0],points[j,1]], 1.1))
    hopping_NN=np.delete(hopping_NN,np.where(hopping_NN==j))
    a0 = tree.query_ball_point([points[j,0],points[j,1]], 1.1)
    a1 = tree.query_ball_point([points[j,0],points[j,1]], 1.75)
    a2 = np.intersect1d(a0,a1)
    a3 = np.setdiff1d(a1,a2)
    for i in range(len(points)):
        if len(a0)!=5 and len(a3)!=4:
            for k in coordinates_NN(j):
                if any(np.all(np.round(k,2)==elt) for elt in np.round(per_neighbors(i),2)):
                    hopping_NN=np.append(hopping_NN, i)
                    hopping_NN=np.unique(hopping_NN)
    return hopping_NN    

# j. sitenin NNN komşulukları tüm koordinatları
def coordinates_NNN(j):
    co=np.zeros((0,2),int)
    a0 = treee.query_ball_point([points[j,0],points[j,1]], 1.1)
    a1 = treee.query_ball_point([points[j,0],points[j,1]], 1.75)
    a2 = np.intersect1d(a0,a1)
    a3 = np.setdiff1d(a1,a2)
    a3 = np.append(a3, [j], axis=0)
    for i in a3:
        co=np.append(co,[abcd[i]],axis=0)
    return co

# j. sitenin periyodik NNN komşulukları
def find_per_NextNearestNeighbors(j):
    hopping_NNN=np.zeros((0,1),int)
    a0 = tree.query_ball_point([points[j,0],points[j,1]], 1.1)
    a1 = tree.query_ball_point([points[j,0],points[j,1]], 1.75)
    a2 = np.intersect1d(a0,a1)
    a3 = np.setdiff1d(a1,a2)
    hopping_NNN = np.append(hopping_NNN,a3)
    for i in range(len(points)):
        if len(a0)!=5 or len(a3)!=4: 
            for k in coordinates_NNN(j):
                if any(np.all(np.round(k,2)==elt) for elt in np.round(per_neighbors(i),2)):
                    hopping_NNN=np.append(hopping_NNN, i)
                    hopping_NNN=np.unique(hopping_NNN)
    return hopping_NNN

# t1 = 1;L1 = 0.28;t2 = -0.3;L2 = 0.2
t1=1;L1 = 0;t2 = 0;L2 = 0
#Sadece NNN matris değerleri kontrolü için
# t1 = 0;L1 = 0;t2 = -0.3;L2 = 0.2
def H():
    H = np.zeros((len(points), len(points)), dtype=complex)
    for m in range(len(points)): 
        for n in find_per_NearestNeighbors(m):
            #A-->B
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_1,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_2,2)):
                    H[m][n]= -t1+1j*L1
                    # print(str(m)+" --> "+str(n)+" : "+str("pozitif"))
            #A-->C
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_1,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_3,2)):
                    H[m][n] = -t1-1j*L1
                    # print(str(m)+" --> "+str(n)+" : "+str("negatif"))
            #B-->A
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_2,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_1,2)):
                    H[m][n] = -t1-1j*L1
                    # print(str(m)+" --> "+str(n)+" : "+str("negatif"))
            #B-->C
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_2,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_3,2)):
                    H[m][n] = -t1+1j*L1
                    # print(str(m)+" --> "+str(n)+" : "+str("pozitif"))
            #C-->A
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_3,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_1,2)):
                    H[m][n] = -t1+1j*L1
                    # print(str(m)+" --> "+str(n)+" : "+str("pozitif"))
            #C-->B
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_3,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_2,2)):
                    H[m][n] = -t1-1j*L1
                    # print(str(m)+" --> "+str(n)+" : "+str("negatif"))
        for n in find_per_NextNearestNeighbors(m):
            #A-->B
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_1,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_2,2)):
                    H[m][n] = t2-1j*L2
            #A-->C
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_1,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_3,2)):
                    H[m][n] = t2+1j*L2
            #B-->A
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_2,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_1,2)):
                    H[m][n] = t2+1j*L2
            #B-->C
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_2,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_3,2)):
                    H[m][n] = t2-1j*L2
            #C-->A
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_3,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_1,2)):
                    H[m][n] = t2-1j*L2
            #C-->B
            if any(np.all(np.round(points[m],2)==elt) for elt in np.round(triangle_3,2)):
                if any(np.all(np.round(points[n],2)==elt) for elt in np.round(triangle_2,2)):
                    H[m][n] = t2+1j*L2        
    return H         
# wR,uR=np.linalg.eig(H())
# E_real=np.real(np.sort(wR))

x_eig, x_func = np.linalg.eig(H())
fig5, ax5 = plt.subplots()
idx = np.argsort(x_eig)
sorted = x_eig[idx]
x_sort_real = np.arange(0, len(points), 1)
y_sort_real = sorted
ax5.set_xlabel(r"$n$")
ax5.set_ylabel(r'$E$')
ax5.set_title(r'$N_{site}= $'+str(len(points))+'\n'+r'$t_1: $'+str(t1)+' , '+r'$\lambda_1: $'+str(L1)+' , '+r'$t_2: $'+str(t2)+' , '+r'$\lambda_2: $'+str(L2))
ax5.plot(x_sort_real, y_sort_real, 'r-', markersize=0.7)