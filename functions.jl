using SparseArrays
include("param.jl")


""" Index - magnetic sublevel dictionaries """
function mF1(g)
    trunc(Int,n2Fm_ats_g[g][2][2])
end

function mF2(e)
    trunc(Int,n2Fm_ats_e[e][2][2])
end


function F1(g)
    trunc(Int,n2Fm_ats_g[g][2][1])
end

function F2(e)
    trunc(Int,n2Fm_ats_e[e][2][1])
end

"""
magn_field(par,B)

DESCRIPTION
Zeeman splitting H0B = H0 + HB
Computes magnetic sublevel energy splitting in magnetic field.

INPUT
par: atom parameters
B: magnetic field

OUTPUT
results_gr[1]: ground state eigenvalues
results_gr[2]: ground state eigenvectors
results_ex[1]: excited state eigenvalues
results_ex[2]: excited state eigenvectors

"""

function magn_field(par::param,laz::lazers,B)
    H=par.A_gr*par.JIfmbasis_gr + laz.muB_MHz*B*(par.gJ_gr*par.Jzfmbasis_gr+par.gI*par.Izfmbasis_gr)   
    results_gr=eigenstates(dense((H+dagger(H))/2))                                              

    i=par.nucI
    j=par.J2
    H = par.A_ex*par.JIfmbasis_ex+laz.muB_MHz*B*(par.gJ_ex*par.Jzfmbasis_ex+par.gI*par.Izfmbasis_ex)
    results_ex=eigenstates(dense((H+dagger(H))/2))
    
    return results_gr[1],results_gr[2],results_ex[1],results_ex[2]
end

"""
DESCRIPTION
Transition matrix dipole elements - computed from Wigner 3j and 6j symbols

INPUT
F,m,J,I: atom quantum numbers
q: light polarization
g,e: magnetic sublevels corresponding to F, m quantum numbers. Called in rate equations

OUTPUT
Transition matrix dipole elements.

"""

function lin_three_j(f1,m1,f2,m2,q)
    if -m2+m1+q == 0
        value=(-1)^(f2-m2)*wigner3j(f2,1,f1,-m2,q,m1)
    else 
        value=0
    end
    return value
end

function lin_six_j(f1,f2,j1,j2,I)
    value=(-1)^(j2+I+f1+1)*sqrt((2*f1+1)*(2*f2+1))*wigner6j(j2,f2,I,f1,j1,1)
    return value
end

function lin_three_j_star(f1,m1,f2,m2,q)
    value=(-1)^(f1-m1)*wigner3j(f1,1,f2,-m1,q,m2)
    return value
end

function lin_six_j_star(f1,f2,j1,j2,I)
    value=(-1)^(j1+I+f2+1)*sqrt((2*f1+1)*(2*f2+1))*wigner6j(j1,f1,I,f2,j2,1)
    return value
end

function lin_dip(j1,j2,I,e,g,q)
    f1=F1(g)
    f2=F2(e)
    m1=mF1(g)
    m2=mF2(e)
    value = lin_six_j(f1,f2,j1,j2,I)*lin_three_j(f1,m1,f2,m2,q)
    return value
end

function lin_dip_star(j1,j2,I,e,g,q)
    f1=F1(g)
    f2=F2(e)
    m1=mF1(g)
    m2=mF2(e)
    value = lin_three_j_star(f1,m1,f2,m2,q)*lin_six_j_star(f1,f2,j1,j2,I)
    return value
end


""" 
Evals_g(Iz_gr,Jz_gr,F²fmbasis_gr,dim_g,eigvals_gr,eigvects_gr); Evals_e(Iz_ex,Jz_ex,F²fmbasis_ex,dim_e,eigvals_ex,eigvects_ex)

DESCRIPTION
Dictionary that connects calculated Zeeman splitting energy of a magnetic sublevel to F and m quantum numbers of energy levels.

INPUT
eigvals_gr: magnetic sublevel energies

OUTPUT
Dictionary for Fm quantum numbers of the level to the energy of the level.

"""

""" 
omega_g(par,Fm2E_g,gi,gj); omega_e(par,Fm2E_e,ei,ej)

DESCRIPTION
Computes energy difference between magnetic sublevels in ground state; in excited state.

INPUT
par: atom parameters
Fm2E: Fm to energy dictionary
gi,gj; ei,ej: index of magnetic sublevels         

OUTPUT
Energy difference between magnetic sublevels

"""

""" Energy difference between GROUND STATE levels """

function Evals_g(Iz_gr,Jz_gr,F²fmbasis_gr,dim_g,eigvals_gr,eigvects_gr)
    F²_gr=F²fmbasis_gr
    Fm2E_g=Dict()


    for i in 1:dim_g
        
        Fz_gr=Iz_gr+Jz_gr

        mFg=dagger(eigvects_gr[i])*Fz_gr*eigvects_gr[i]
        
        if round(real(mFg)) ==-0.0
            mFg=0.0
        end

        F²=dagger(eigvects_gr[i])*F²_gr*eigvects_gr[i]
        Fg = (-1 + sqrt(1+4*F²))/2

        get!(Fm2E_g,(round(real(Fg)),round(real(mFg))),eigvals_gr[i])
        
    end
    
    return Fm2E_g

end

function omega_g(par,Fm2E_g,gi,gj) 
    Fm2E_g[n2Fm_ats_g[gi][2]]-Fm2E_g[n2Fm_ats_g[gj][2]]
end


""" Energy difference between EXCITED STATE levels """
function Evals_e(Iz_ex,Jz_ex,F²fmbasis_ex,dim_e,eigvals_ex,eigvects_ex)
    F²_ex=F²fmbasis_ex
    Fm2E_e=Dict()


    for i in 1:dim_e
        
        Fz_ex=Iz_ex+Jz_ex
        mFe=dagger(eigvects_ex[i])*Fz_ex*eigvects_ex[i]
        if round(real(mFe)) ==-0.0
            mFe=0.0
        end
        F²=dagger(eigvects_ex[i])*F²_ex*eigvects_ex[i]
        Fe = (-1 + sqrt(1+4*F²))/2

        get!(Fm2E_e,(round(real(Fe)),round(real(mFe))),eigvals_ex[i])

    end
    
    return Fm2E_e

end

function omega_e(par,Fm2E_e,ei,ej)
    Fm2E_e[n2Fm_ats_e[ei][2]]-Fm2E_e[n2Fm_ats_e[ej][2]]
end


""" Energy difference between GROUND STATE and EXCITED STATE levels """
function omega_ge(par::param,g,e,Fm2E_g,Fm2E_e) 
    par.Efs+Fm2E_e[n2Fm_ats_e[e][2]]-Fm2E_g[n2Fm_ats_g[g][2]]
end


"""
ksi_f(par::param,laz::lazers,g,e,Fm2E_g,Fm2E_e,detune) and ksi_cc_f(par::param,laz::lazers,e,g,Fm2E_g,Fm2E_e,detune)

DESCRIPTION
Creates Ξ matrix elements and fills them in a matrix.

INPUT
g: ground state level index, iterated through in rate equations
e: excited state level index, iterated through in rate equations
Fm2E_g: Fm to E ground state level dictionary 
Fm2E_e: Fm to E excited state level dictionary
detune: detuning (shift in frequency) due to Dopler effect

OUTPUT
Ξ matrix element

"""

function ksi_f(par::param,laz::lazers,g,e,Fm2E_g,Fm2E_e,detune) 
    laz.Ωᵣ^2 /(((laz.Γ+laz.Δω)/2) + laz.γ + 1.0im*(detune-omega_ge(par,g,e,Fm2E_g,Fm2E_e)))
end

function ksi_cc_f(par::param,laz::lazers,e,g,Fm2E_g,Fm2E_e,detune)
    laz.Ωᵣ^2 /(((laz.Γ+laz.Δω)/2) + laz.γ - 1.0im*(detune-omega_ge(par,g,e,Fm2E_g,Fm2E_e)))
end

"""
fill_ksi(par::param,laz::lazers,Fm2E_g,Fm2E_e,dim_g,dim_e,detune); fill_ksi_cc(par::param,laz::lazers,Fm2E_g,Fm2E_e,dim_e,dim_g,detune)

DESCRIPTION
Fills Ξ elements in a matrix.

INPUT
Fm2E_g: Fm to ground state energy dictionary 
Fm2E_e: Fm to excited state energy dictionary
dim_g: number of ground state levels
dim_e: number of excited state levels
detune: detuning (shift in frequency) due to Dopler effectv

OUTPUT
Ξ matrix; Ξ matrix complex conjugate

"""

function fill_ksi(par::param,laz::lazers,Fm2E_g,Fm2E_e,dim_g,dim_e,detune)
    ksi=spzeros(Complex{Float64},dim_g,dim_e)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            ksi[g,e]=ksi_f(par,laz,g,e,Fm2E_g,Fm2E_e,detune) 
        end
    end
    return ksi
end

function fill_ksi_cc(par::param,laz::lazers,Fm2E_g,Fm2E_e,dim_e,dim_g,detune)
    ksi_cc=spzeros(Complex{Float64},dim_e,dim_g)
    for e in 1:1:dim_e
        for g in 1:1:dim_g
            ksi_cc[e,g]=ksi_cc_f(par,laz,e,g,Fm2E_g,Fm2E_e,detune)
        end
    end
    return ksi_cc
end


"""
GAMMA_f(Γ,g1,g2,e1,e2,j1,j2,I) 

DESCRIPTION
Calculates Γ matrix elements. Takes into account that two transitons should be excited with the same photon (mF2(e1)-mF1(g1) == mF2(e2)-mF1(g2)).

INPUT
Γ: spontaneous relaxation constant
g1,g2: ground state level index, iterated through in rate equations 
e1,e2: excited state level index, iterated through in rate equations
j1: quantum number of ground state total electron angular momentum
j2: quantum number of excited state total electron angular momentum
I: quantum number of total nuclear angular momentum

OUTPUT
Γ matrix element
"""

function GAMMA_f(Γ,g1,g2,e1,e2,j1,j2,I) 
    sum = 0
    if mF2(e1)-mF1(g1) == mF2(e2)-mF1(g2) 
        for q in -1:+1:+1
            sum+=lin_dip(j1,j2,I,e1,g1,q)*lin_dip_star(j1,j2,I,e2,g2,-q)*(-1)^q
        end    
    end 
    sum*=(2*j2+1)*Γ
    return sum
end

"""
DESCRIPTION
Fills Γ matrix elements in a 4 dimensional matrix, where each element is described by 4 indexes of ground and excited state levels of two magnetic sublevel tranitions

INPUT
par: atom parameters
laz: laser parameters

OUTPUT
Γ matrix
"""
function fill_GAMMA(par::param,laz::lazers)
    GAMMA=zeros(Complex{Float64},par.dim_g,par.dim_g,par.dim_e,par.dim_e)
    for g1 in 1:1:par.dim_g
        for g2 in 1:1:par.dim_g
            for e1 in 1:1:par.dim_e 
                for e2 in 1:1:par.dim_e
                    GAMMA[g1,g2,e1,e2]=GAMMA_f(laz.Γ,g1,g2,e1,e2,par.J1,par.J2,par.nucI) 
                end
            end
        end
    end
    return GAMMA

end


"""
lin_fillDipoleMatrix(par::param,e_vec) and lin_fillDipoleMatrix_star(par::param,e_vec)

DESCRIPTION
Calculates dipole matrix elements using excited light vector components and fills them in a matrix. The dimensions of the matrix are the number of excited state levels x the number of ground state levels

INPUT
par: atom parameters
e_vec: excitation light E vector 

OUTPUT
Transition dipole matrix and its complex conjugate matrix.

"""

function lin_fillDipoleMatrix(par::param,e_vec)
    j1=par.J1
    j2=par.J2
    I=par.nucI
    dim_g=par.dim_g
    dim_e=par.dim_e

    dip_p=spzeros(Complex{Float64},dim_e,dim_g)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_p[e,g]=e_vec[1]*lin_dip(j1,j2,I,e,g,+1) 
        end
    end
    dip_0=spzeros(Complex{Float64},dim_e,dim_g)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_0[e,g]=e_vec[2]*lin_dip(j1,j2,I,e,g,0)
        end
    end
    dip_m=spzeros(Complex{Float64},dim_e,dim_g)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_m[e,g]=e_vec[3]*lin_dip(j1,j2,I,e,g,-1)
        end
    end
    dip = dip_p + dip_0 + dip_m
    return dip 
end

function lin_fillDipoleMatrix_star(par::param,e_vec)
    j1=par.J1
    j2=par.J2
    I=par.nucI
    dim_g=par.dim_g
    dim_e=par.dim_e

    dip_p_star=spzeros(Complex{Float64},dim_g,dim_e) 
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_p_star[g,e]=conj(e_vec[3])*(-1)*lin_dip_star(j1,j2,I,e,g,1) 
        end
    end
    dip_0_star=spzeros(Complex{Float64},dim_g,dim_e)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_0_star[g,e]=conj(e_vec[2])*(-1)^0*lin_dip_star(j1,j2,I,e,g,0)
        end
    end
    dip_m_star=spzeros(Complex{Float64},dim_g,dim_e)
    for g in 1:1:dim_g
        for e in 1:1:dim_e
            dip_m_star[g,e]=conj(e_vec[1])*(-1)*lin_dip_star(j1,j2,I,e,g,-1)
        end
    end
    dip_star = dip_p_star + dip_0_star + dip_m_star
    return dip_star 
end


"""
ρDotgg(laz::lazers,matr,dim_g,dim_e,ksi,ksi_cc, dip, dip_star,Fm2E_g)

DESCRIPTION
Creates rate equation ground state coeficient matrix C:
dρ/dt=Cρ
Each row of the matrix represents one dρ(gigj)/dt equation, with each column representing coeffitient for each ρ element. This gives the matrix a dimension of (number of ground states)²x(number of ground states)². 

OUTPUT
Ground state coeficient matrix C.

"""

function ρDotgg(laz::lazers,matr,dim_g,dim_e,ksi,ksi_cc, dip, dip_star,Fm2E_g)
    for index in range(1,length=length(gDict))
    
        gi=getGi(index)
        gj=getGj(index)
    
        for ek in range(1,length=dim_e)
            for em in range(1,length=dim_e)
                if dip_star[gi,ek]*dip[em,gj] != 0
                    matr[index,length(gDict)+eDictInv[ek,em]]+=(ksi[gi,em]+ksi_cc[ek,gj])*dip_star[gi,ek]*dip[em,gj] 
                end
            end

            for gm in range(1,length=dim_g)
                if dip_star[gi,ek]*dip[ek,gm] != 0
                    matr[index,gDictInv[gm,gj]] += -ksi_cc[ek,gj]*dip_star[gi,ek]*dip[ek,gm]
                end

                if dip_star[gm,ek]*dip[ek,gj] !=0
                    matr[index,gDictInv[gi,gm]] += -ksi[gi,ek]*dip_star[gm,ek]*dip[ek,gj]
                end
            end

            for el in range(1,length=dim_e)
                matr[index,length(gDict)+eDictInv[ek,el]] += GAMMA[gi,gj,ek,el]
            end
        end


        matr[index,gDictInv[gi,gj]]  += -1.0im*omega_g(par,Fm2E_g,gi,gj)
        matr[index,gDictInv[gi,gj]]  += -laz.γ

    end
end    

"""
ρDotee(laz::lazers,matr,dim_g,dim_e,ksi,ksi_cc,dip, dip_star,Fm2E_e)

DESCRIPTION
Creates rate equation excitedd state coeficient matrix C:
dρ/dt=Cρ
Each row of the matrix represents one dρ(eiej)/dt equation, with each column representing coeffitient for each ρ element. This gives the matrix a dimension of (number of excited states)²x(number of excited states)². 

OUTPUT
Excited state coeficient matrix C.

"""

function ρDotee(laz::lazers,matr,dim_g,dim_e,ksi,ksi_cc,dip, dip_star,Fm2E_e)
    for index in range(1,length=length(eDict))
        ei=getEi(index)
        ej=getEj(index)
        
    
        for gk in range(1,length=dim_g)
            for gm in range(1,length=dim_g)
                if dip_star[gm,ej]*dip[ei,gk] !=0
                    matr[length(gDict)+index,gDictInv[gk,gm]] += (ksi[gk,ej]+ksi_cc[ei,gm])*dip_star[gm,ej]*dip[ei,gk] 
                end
            end

            for em in range(1,length=dim_e)
                if dip_star[gk,em]*dip[ei,gk] !=0
                    matr[length(gDict)+index,length(gDict)+eDictInv[em,ej]] += -ksi[gk,ej]*dip_star[gk,em]*dip[ei,gk]
                end

                if dip_star[gk,ej]*dip[em,gk] !=0
                    matr[length(gDict)+index,length(gDict)+eDictInv[ei,em]] += -ksi_cc[ei,gk]*dip_star[gk,ej]*dip[em,gk]
                end
            end
        end

        matr[length(gDict)+index,length(gDict)+eDictInv[ei,ej]] += -1.0im*omega_e(par,Fm2E_e,ei,ej)
        matr[length(gDict)+index,length(gDict)+eDictInv[ei,ej]] += -(laz.Γ+laz.γ)
    end
end 

"""
nov_fillDipoleMatrix(j1,j2,I,dim_g,dim_e,e_vec) and nov_fillDipoleMatrix_star(j1,j2,I,dim_g,dim_e,e_vec).

DESCRIPTION
Calculates dipole matrix elements for observing using observation light vector components and fills them in a matrix. The dimensions of the matrix are the number of excited state levels x the number of ground state levels

INPUT
par: atom parameters
e_vec: observation light E vector 

OUTPUT
Transition dipole matrix and its complex conjugate matrix.

"""

function nov_fillDipoleMatrix(par::param,e_vec) 
    j1=par.J1
    j2=par.J2
    I=par.nucI
    dim_g=par.dim_g
    dim_e=par.dim_e

    dip_p=spzeros(Complex{Float64},dim_e,dim_g)
    for g in keys(n2Fm_ats_g) 
        for e in keys(n2Fm_ats_e) 
            dip_p[eDictInv[5,e],gDictInv[5,g]]=e_vec[1]*lin_dip(j1,j2,I,e,g,+1) 
        end
    end
    dip_0=spzeros(Complex{Float64},dim_e,dim_g)
    for g in keys(n2Fm_ats_g)
        for e in keys(n2Fm_ats_e)
            dip_0[eDictInv[5,e],gDictInv[5,g]]=e_vec[2]*lin_dip(j1,j2,I,e,g,0)
        end
    end
    dip_m=spzeros(Complex{Float64},dim_e,dim_g)
    for g in keys(n2Fm_ats_g)
        for e in keys(n2Fm_ats_e)
            dip_m[eDictInv[5,e],gDictInv[5,g]]=e_vec[3]*lin_dip(j1,j2,I,e,g,-1)
        end
    end
    dip = dip_p + dip_0 + dip_m
    return dip 
end

function nov_fillDipoleMatrix_star(par::param,e_vec)
    j1=par.J1
    j2=par.J2
    I=par.nucI
    dim_g=par.dim_g
    dim_e=par.dim_e

    dip_p_star=spzeros(Complex{Float64},dim_g,dim_e)
    for g in keys(n2Fm_ats_g)
        for e in keys(n2Fm_ats_e)
            dip_p_star[gDictInv[5,g],eDictInv[5,e]]=conj(e_vec[3])*(-1)*lin_dip_star(j1,j2,I,e,g,1)
        end
    end
    dip_0_star=spzeros(Complex{Float64},dim_g,dim_e)
    for g in keys(n2Fm_ats_g)
        for e in keys(n2Fm_ats_e)
            dip_0_star[gDictInv[5,g],eDictInv[5,e]]=conj(e_vec[2])*(-1)^0*lin_dip_star(j1,j2,I,e,g,0)
        end
    end
    dip_m_star=spzeros(Complex{Float64},dim_g,dim_e)
    for g in keys(n2Fm_ats_g)
        for e in keys(n2Fm_ats_e)
            dip_m_star[gDictInv[5,g],eDictInv[5,e]]=conj(e_vec[1])*(-1)*lin_dip_star(j1,j2,I,e,g,-1)
        end
    end
    dip_star = dip_p_star + dip_0_star + dip_m_star
    return dip_star 
end