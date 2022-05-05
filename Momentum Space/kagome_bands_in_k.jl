using TightBinding
#Primitive vectors
a1 = [sqrt(3)/2,-1/2]
a2 = [sqrt(3)/2,1/2]
#set lattice
la = set_Lattice(2,[a1,a2])
#add atoms
add_atoms!(la,[0,0])
add_atoms!(la,[0,1/2])
add_atoms!(la,[1/2,0])

#construct hoppings
t = -1.0
add_hoppings!(la,-t,1,2,[1/2,0])
add_hoppings!(la,-t,1,2,[-1/2,0])
add_hoppings!(la,-t,1,3,[0,-1/2])
add_hoppings!(la,-t,1,3,[0,1/2])
add_hoppings!(la,-t,3,2,[1/2,-1/2])
add_hoppings!(la,-t,3,2,[-1/2,1/2])

using Plots
#show the lattice structure
#kagome_lattice = plot_lattice_2d(la)
plot_lattice_2d(la)
#savefig(kagome_lattice,"kagome_lattice.png")

using Plots
# Density of states
nk = 1000 #numer ob meshes. nk^d meshes are used. d is a dimension.
plot_DOS(la, nk)

nk = 1000 #numer ob meshes. nk^d meshes are used. d is a dimension.
hist = get_DOS(la, nk)
println(hist.weights) #DOS data
println(hist.edges[1]) #energy mesh
using Plots
plot(hist.edges[1][2:end] .- hist.edges[1].step.hi/2,hist.weights)

#show the band structure
klines = set_Klines()

kmin = [0,0]
kmax = [2π/sqrt(3),0]
add_Kpoints!(klines,kmin,kmax,"G","K")

kmin = [2π/sqrt(3),0]
kmax = [2π/sqrt(3),2π/3]
add_Kpoints!(klines,kmin,kmax,"K","M")

kmin = [2π/sqrt(3),2π/3]
kmax = [0,0]
add_Kpoints!(klines,kmin,kmax,"M","G")

calc_band_plot(klines,la)

savefig(calc_band_plot(klines,la),"in_k.png")

using Plots
#We have already constructed atoms and hoppings.
#We add the line to plot
klines = set_Klines()
kmin = [-π]
kmax = [π]
add_Kpoints!(klines,kmin,kmax,"-pi","pi")

#We consider the periodic boundary condition along the primitive vector
direction = 1
#Periodic boundary condition
calc_band_plot_finite(klines,la,direction,periodic=true)

#We introduce the surface perpendicular to the premitive vector
direction = 1
#Open boundary condition
calc_band_plot_finite(klines,la,direction,periodic=false)
