using QuantumOptics
using LinearAlgebra
using SparseArrays

max_occupation_number = 2
occ_array = []
for i in 0:max_occupation_number
    push!(occ_array, i)
end

Nx = Ny = 2
N = Nx*Ny 
b = NLevelBasis(N)

states = bosonstates(b, [occ_array...]) #en az 0 en fazla 2 bozon
for i in 1:length(states)
    println("$i: ", states[i])
end

b_mb_b = ManyBodyBasis(b, states)

using OffsetArrays
lat = range(1,Nx*Ny)
latt = reshape(lat, (Nx,Ny))
lattice = OffsetArray(latt, 0:Nx-1, 0:Ny-1)
x_co = range(0, Nx-1) 
y_co = range(0, Ny-1)
arr = []
for i in 1:length(x_co)
    for j in 1:length(y_co)
       arr = [arr;x_co[i];y_co[j]]
    end
end
xy = reshape(arr, (2,Nx*Ny)) |> transpose
arr = OffsetArray(arr,0:Nx*Ny*2-1)
#xy[latis,1:x 2:y]
lattice

x = []
function PerBC()
    for i in 0:Nx-1
        for j in 0:Ny-1
            cc = [lattice[mod(j,Ny),mod(i-1,Nx)],lattice[mod(j+1,Ny),mod(i,Nx)],lattice[mod(j,Ny),mod(i+1,Nx)],lattice[mod(j-1,Ny),mod(i,Nx)]]
            push!(x,cc)
        end
    end
    return x
end
PerBC() #index 1'den başlıyor

function HHH(m, n, alpha)
    if abs(xy[m,1]-xy[n,1])==Nx-1
        if xy[m,1] > xy[n,1]
            A = -exp(-1im*2*pi*alpha*xy[m,2])
            #println(m," -> ",n," Kenarda x+ ve faz: ", A)
        elseif xy[m,1] < xy[n,1]
            A = -exp(1im*2*pi*alpha*xy[m,2])
            #println(m," -> ",n," Kenarda x- ve faz: ", A)
        end
    elseif m==n
        A = 0
        #println(m, " -> ", n, " Aynı noktalar ve faz: ", A)
    else
        if xy[m,1] > xy[n,1]
            A = -exp(1im*2*pi*alpha*xy[m,2])
            #println(m," -> ",n," x- ve faz: ", A)
        elseif xy[m,1] < xy[n,1]
            A = -exp(-1im*2*pi*alpha*xy[m,2])
            #println(m," -> ",n," x+ ve faz: ", A)
        else
            A = -exp(0)
            #println(m," -> ",n," x sabit ve faz: ", A)
        end
    end
end 

alpha = 1/4
# 2 / (4*4*1/4) = 1 / 2
S = SparseOperator(b_mb_b)
for m in 1:N
    for n in 1:N
        if m in PerBC()[n]
            S = S + HHH(m, n, alpha)*transition(b_mb_b, m, n)
        end
    end
end
S
eigenenergies(dense(S))
  
using Plots
eig = eigenenergies(dense(S))
x=1:length(eig)
gr()
plot(x, eig, seriestype = :scatter, xlabel="Öz-Değer Sırası", ylabel="Enerjiler", markersize = 5, labels="U=0", legend = :outertopleft)

A = eigenstates(dense(S))

C = Array{Float64}(undef, length(states), 2)
for i in 1:length(states)
    exp = round(expect(number(b_mb_b), A[2][i])) #expected values
    C[i] = exp #expected values (first column)
    C[i,2] = A[1][i] #eigen-values (second column)
end
C
CC = sortslices(C, dims=1 ,by = x -> x[1]) #sorted according to particle numbers
show(stdout, "text/plain", CC)

using DataFrames
df = DataFrame(CC, :auto)
rename!(df,:x1 => :Parçacık_Sayısı)
rename!(df,:x2 => :Enerji)

j = 11
A[1][j] #2 parçacıklı baz vektörü enerjisi
A[2][j] #Bu baz vektöre karşılık gelen özvektörde vakum ve tek parçacıklı baz vektörleri
#açılım katsayıları sıfırdır.
#Yukarıdaki hesaplanan özvektör içerisindeki baz vektörlerin ağırlıkları
for i in 1:length(states)    
    println(abs(A[2][j].data[i]))
end
