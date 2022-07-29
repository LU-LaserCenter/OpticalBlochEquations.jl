using LinearAlgebra
using QuantumOptics
using Plots
using WignerSymbols
using Printf 
using QuantumOptics.steadystate
using DelimitedFiles
using SparseArrays
using Dates

include("param.jl")
include("functions.jl")
include("polarizacija.jl")

n2Fm=Dict{Int64,Tuple}()
Fm2n=Dict{Tuple,Int64}()
myIndex=1

#g un e vertību indeksi
function fill_n2Fm_ats_g(par::param)
    n2Fm_ats_g=Dict{Int64,Tuple}()
    myIndex=1

    for Fg=par.grFmax:-1:par.grFmin
        for mg=Fg:-1:-Fg
            get!(n2Fm_ats_g,myIndex,("ground",(Fg,mg)))  
            myIndex+=1
        end
    end

    return n2Fm_ats_g
end

function fill_n2Fm_ats_e(par::param)
    n2Fm_ats_e=Dict{Int64,Tuple}()
    Fm2n_ats_e=Dict{Tuple,Int64}()
    myIndex=1

    for Fe=par.exFmax:-1:par.exFmin
        for me=Fe:-1:-Fe
            get!(n2Fm_ats_e,myIndex,("excited", (Fe,me)))
            myIndex+=1
        end
    end

    return n2Fm_ats_e
end

#gi gj apvienotā dictionary
function fill_gDict(n2Fm_ats_g)
    gDict=Dict()
    myIndex=1

    for i in keys(n2Fm_ats_g) 
        for j in keys(n2Fm_ats_g)
            push!(gDict,myIndex=>(i,j))
            myIndex+=1
        end
    end

    return gDict
end

#ei ej apvienotā dictionary 
function fill_eDict(n2Fm_ats_e)
    eDict=Dict()
    myIndex=1 

    for i in keys(n2Fm_ats_e) 
        for j in keys(n2Fm_ats_e)
            push!(eDict,myIndex=>(i,j))
            myIndex+=1
        end
    end

    return eDict
end


par=param(cezijsD1)
laz=lazers()
n2Fm_ats_g=fill_n2Fm_ats_g(par)
n2Fm_ats_e=fill_n2Fm_ats_e(par)

gDict=fill_gDict(n2Fm_ats_g)
eDict=fill_eDict(n2Fm_ats_e)

function getGi(key)
    gDict[key][1]
end
function getGj(key)
    gDict[key][2]
end
function getEi(key)
    eDict[key][1]
end
function getEj(key)
    eDict[key][2]
end

#Reverse the dictionaries
gDictInv=Dict(value => key for (key, value) in gDict)
eDictInv=Dict(value => key for (key, value) in eDict)


#ierosmes, novērošanas un zondēšanas ģeometrijas definēšana (pol, θ, ϕ)
e_vec_i=ElectricVector(1, π/2, 0).cyclic
e_vec_n=ElectricVector(1, π/2, π/2).cyclic
e_vect_z=ElectricVector(0, π/2, π/4).cyclic
dip_star = lin_fillDipoleMatrix_star(par,e_vec_i) #izveido dipolu matricu
dip=lin_fillDipoleMatrix(par,e_vec_i)
GAMMA=fill_GAMMA(par,laz) #izveido spontānās relaksācijas matricu

#induced relaxation
λ=zeros(Complex{Float64},par.dim_g*par.dim_g+par.dim_e*par.dim_e)
for index in range(1,length=length(gDict))
    gi=getGi(index)
    gj=getGj(index)
    if(gi==gj) 
        λ[index]=laz.γ/par.dim_g
    end
end

I_matr=[]
A_matr=[]
pop_matr=[]

A_Doplera_matr=[]
I_Doplera_matr=[]
Bmin=-40
Bmax=40
Bstep=5
Brange=range(Bmin,Bmax,step=Bstep)
Brange=[-40,-35,-30,-25,-20,-15,-10,-5,-4,-3,-2,-1,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0]

"""Dopler shift"""
step_lit=2 
sigma=sqrt(laz.kB*laz.T/laz.masa)*laz.ω_svitr/laz.c
dscan=2*sigma*step_lit
dstep_nr=150*step_lit #number of steps
dstep_length=dscan/dstep_nr

#dscan=20
#dstep_length=2.0

for B₀ in Brange
        """Creates matrix of dimensions ((number of ground states)²+(number of excited states)²)x((number of ground states)²+(number of excited states)²)"""
        matr=spzeros(Complex{Float64}, par.dim_g*par.dim_g+par.dim_e*par.dim_e,par.dim_g*par.dim_g+par.dim_e*par.dim_e)
        eigvals_gr,eigvects_gr,eigvals_ex,eigvects_ex=magn_field(par,laz,B₀)

        Fm2E_g=Evals_g(par.Izfmbasis_gr,par.Jzfmbasis_gr,par.F²fmbasis_gr,par.dim_g,eigvals_gr,eigvects_gr)
        Fm2E_e=Evals_e(par.Izfmbasis_ex,par.Jzfmbasis_ex,par.F²fmbasis_ex,par.dim_e,eigvals_ex,eigvects_ex)

        """Dopler""" 
        A_Doplera=0
        I_Doplera=0

        println(B₀)
        println("Doplera cikls: ", Dates.format(now(), "HH:MM:SS"))
        println("Doplera cikls: ")
        dshift=0
        @time for dshift in range(-dscan, step=dstep_length, stop=dscan) 
            detune=laz.ω_svitr-dshift
            
            matr=spzeros(Complex{Float64}, par.dim_g*par.dim_g+par.dim_e*par.dim_e,par.dim_g*par.dim_g+par.dim_e*par.dim_e)
            ksi=fill_ksi(par,laz,Fm2E_g,Fm2E_e,par.dim_g, par.dim_e,detune)
            ksi_cc=fill_ksi_cc(par,laz,Fm2E_g,Fm2E_e,par.dim_e,par.dim_g, detune)

            """Calling functions ρDotgg and ρDotee places matrix blocks dρ(gigj)/dt and dρ(eiej)/dt respectively on the diagonal of matrix "matr" thus creating the full coeffitient matrix for the rate equations"""
            ρDotgg(laz,matr,par.dim_g,par.dim_e,ksi,ksi_cc,dip, dip_star,Fm2E_g)
            ρDotee(laz,matr,par.dim_g,par.dim_e,ksi,ksi_cc,dip, dip_star,Fm2E_e)


            rho_0=zeros(Complex{Float64}, par.dim_g*par.dim_g+par.dim_e*par.dim_e)
            rho_0[1]+=1.0+0.0im
            
            """LU decomposition"""
            lu_matr_sparse=lu(matr)

            sol=zeros(par.dim_g*par.dim_g+par.dim_e*par.dim_e)
            """Solving the rate equations to get the density matrix ρ"""
            eigen_vectors = lu_matr_sparse \ (-λ)


            """Fluorescence"""
            rho_e=permutedims(reshape(eigen_vectors[par.dim_g*par.dim_g+1:(par.dim_g*par.dim_g)+par.dim_e*par.dim_e], par.dim_e, par.dim_e))
            rho_e_dimge=spzeros(Complex{Float64}, par.dim_g+par.dim_e, par.dim_g+par.dim_e)
            rho_e_dimge[(par.dim_g+1):(par.dim_g+par.dim_e),(par.dim_g+1):(par.dim_g+par.dim_e)]=rho_e

            d_star_nov_dimge=spzeros(Complex{Float64},par.dim_g+par.dim_e,par.dim_g+par.dim_e)
            d_star_nov=nov_fillDipoleMatrix_star(par,e_vec_n) 
            d_star_nov_dimge[1:par.dim_g,(par.dim_g+1):(par.dim_g+par.dim_e)]=d_star_nov

            d_nov_dimge=spzeros(Complex{Float64},par.dim_g+par.dim_e,par.dim_g+par.dim_e)
            d_nov=nov_fillDipoleMatrix(par,e_vec_n) 
            d_nov_dimge[(par.dim_g+1):(par.dim_g+par.dim_e),1:par.dim_g]=d_nov

            I=0
            I=d_star_nov_dimge*rho_e_dimge*d_nov_dimge
            I_Doplera+=tr(I)*dstep_length*(exp(-dshift^2/(2*sigma^2)))/(sqrt(2*π)*sigma) 
            push!(I_matr,tr(I))


            # """Absorbtion""" 
            rho_g=permutedims(reshape(eigen_vectors[1:par.dim_g*par.dim_g], par.dim_g, par.dim_g))
            rho_g_dimge=spzeros(Complex{Float64}, par.dim_g+par.dim_e, par.dim_g+par.dim_e)
            rho_g_dimge[1:par.dim_g,1:par.dim_g]=rho_g

            GammaR2=((laz.Γ+laz.Δω+laz.γ)/2)^2 
            #d_star_z=nov_fillDipoleMatrix_star(par.J1,par.J2,par.nucI,par.dim_g,par.dim_e,e_vect_z)
            d_star_z=nov_fillDipoleMatrix_star(par,e_vect_z)
            #d_z=nov_fillDipoleMatrix(par.J1,par.J2,par.nucI,par.dim_g,par.dim_e,e_vect_z)
            d_z=nov_fillDipoleMatrix(par,e_vect_z)
            A=0
            for gj in keys(n2Fm_ats_g) 
                for gk in keys(n2Fm_ats_g) 
                    for ei in keys(n2Fm_ats_e) 
                        A+=(d_z[eDictInv[5,ei],gDictInv[5,gj]]*rho_g[gDictInv[5,gj],gDictInv[5,gk]]*d_star_z[gDictInv[5,gk],eDictInv[5,ei]])/(GammaR2+(detune-omega_ge(par,gj,ei,Fm2E_g,Fm2E_e))^2)
                    end
                end
            end

            push!(A_matr,A)

            A_Doplera+=A*dstep_length*(exp(-(dshift^2)/(2*sigma^2)))/(sqrt(2*π)*sigma) 
            push!(pop_matr, tr(rho_e)) #nosaka, vai blīvuma matricas pēda ir 1 ### FHG 2022-07-18
        
        end
        push!(I_Doplera_matr,I_Doplera)
        push!(A_Doplera_matr,A_Doplera)
        

        #push!(pop_matr, tr(rho_e)) #nosaka, vai blīvuma matricas pēda ir 1

end

B₀=Brange

#Fluorescences signāls bez Doplera
X=[]
Z=[]
index=1
for B in B₀
    push!(X,B)
    y=real(I_matr[index])
    push!(Z,y)
    global index=index+1
end
plot(X,Z)

open("I.txt", "w") do io
    writedlm(io, [X Z])
end 
savefig("Cs-I.png")

##Absorbcijas signāls bez Doplera
 X=[]
 Z=[]
 index=1
 for B in B₀
     push!(X,B)
     y=real(A_matr[index])
     push!(Z,y)
     global index=index+1
 end
 plot(X,Z)

 open("A.txt", "w") do io
     writedlm(io, [X Z])
 end 
 savefig("Cs-A.png")

##Absorbcijas signāls ar Dopleru
 X=[]
 Z=[]
 index=1
 for B in B₀
     push!(X,B)
     y=real(A_Doplera_matr[index])
     push!(Z,y)
     global index=index+1
 end
 plot(X,Z)

 open("A_Doplera.txt", "w") do io
     writedlm(io, [X Z])
 end 

 savefig("Cs-A-Doplera.png")

 ##Absorbcijas signāls ar Dopleru
 X=[]
 Z=[]
 index=1
 for B in B₀
     push!(X,B)
     y=real(I_Doplera_matr[index])
     push!(Z,y)
     global index=index+1
 end
 plot(X,Z)

 open("I_Doplera.txt", "w") do io
     writedlm(io, [X Z])
 end 

 savefig("Cs-I-Doplera.png")

# #Populāciju pēda
 X=[]
 Z=[]
 index=1
 for B in B₀
     push!(X,B)
     y=real(pop_matr[index])
     push!(Z,y)
     global index=index+1
 end
 plot(X,Z)

 savefig("Cs-pop.png")