using QuantumOptics
using LinearAlgebra
using SparseArrays

Nx = Ny = 4
N = Nx*Ny

max_occupation_number = 2
occ_array = []
for i in 0:max_occupation_number
    push!(occ_array, i)
end

b = NLevelBasis(N)
states = bosonstates(b, [occ_array...])
b_mb = ManyBodyBasis(b, states)

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
for i in 1:N
    println("$i: ", PerBC()[i])
end

function HP(m, n, alpha)
    if abs(xy[m,1]-xy[n,1])==Nx-1
        if xy[m,1] > xy[n,1]
            P = exp(-1im*2*pi*alpha*xy[m,2])
        elseif xy[m,1] < xy[n,1]
            P = exp(1im*2*pi*alpha*xy[m,2])
        end
    elseif m==n
        P = 0
    else
        if xy[m,1] > xy[n,1]
            P = exp(1im*2*pi*alpha*xy[m,2])
        elseif xy[m,1] < xy[n,1]
            P = exp(-1im*2*pi*alpha*xy[m,2])
        else
            P = exp(0)
        end
    end
end 

alpha = 1/4
# 2/(4*4*1/4)
t = 1
U = 2

#kinetic term
KT = SparseOperator(b_mb)
#interaction term
IT = SparseOperator(b_mb)
for m in 1:N
    IT = IT + U/2 * number(b_mb, m) * (number(b_mb, m) - identityoperator(b_mb))
    for n in 1:N
        if m in PerBC()[n]
            KT = KT - t * HP(m, n, alpha) * transition(b_mb, m, n)
        end
    end
end
HH = KT + IT
eigenenergies(dense(HH))
  
using Plots
eig = eigenenergies(dense(HH))
x=1:length(eig)
gr()
plot(x, eig, seriestype = :scatter, xlabel="Öz-Değer Sırası", ylabel="Enerjiler", markersize = 5, label="U=2", legend = :outertopleft)

A = eigenstates(dense(HH))
  
C = Array{Float64}(undef, length(states), 2)
for i in 1:length(states)
    exp = round(expect(number(b_mb), A[2][i]))
    C[i] = exp 
    C[i,2] = A[1][i]
end
C
CC = sortslices(C, dims=1 ,by = x -> x[1])

using DataFrames
df = DataFrame(CC, :auto)
rename!(df,:x1 => :Parçacık_Sayısı)
rename!(df,:x2 => :Enerji)

b_hard = NLevelBasis(N)
states_hard = fermionstates(b_hard, [2])
b_mb_hard = ManyBodyBasis(b_hard, states_hard)
states_hard

UU=2
#kinetic term
KT_hard = SparseOperator(b_mb_hard)
#interaction term
IT_hard = SparseOperator(b_mb_hard)
for m in 1:N
    IT_hard = IT_hard + UU/2 * number(b_mb_hard, m) * (number(b_mb_hard, m) - identityoperator(b_mb_hard))
    for n in 1:N
        if m in PerBC()[n]
            KT_hard = KT_hard - t * HP(m, n, alpha) * transition(b_mb_hard, m, n)
        end
    end
end
HH_hard = KT_hard + IT_hard
eigenenergies(dense(HH_hard))
show(stdout, "text/plain", CC)
#Bu enerji seviyeleri, U çok büyükken, yukarıdaki Hamiltonyende (parçacık=2 ve oc. num=2) iki parçacıklı durumların enerjilerine denk düşüyor.
