import LinearAlgebra
import Plots
using Plots
using LinearAlgebra

function Wpoly6(xi,xj)
    sigma = 315/(64*pi*h^9)
    q = norm(Differenz(xi,xj))
    if 0<=q<=h
        W = sigma*(h^2-q^2)^3
    else
        W=0
    end
    return W
end

function dWspiky(xi,xj)
    r = Differenz(xi,xj)
    q = norm(r)
    sigma = -45/(q*pi*h^6)
    if 0<=q<=h
        W = (sigma * (h-q)^2 )
        W= W.*r
    else
        W=[0 0]
    end
    return W
end

function d2Wviscosity(xi,xj)
    sigma = 45/(pi*h^6)
    q = norm(Differenz(xi,xj))
    if 0<=q<=h
        W = sigma *(h-q)
    else
        W=0
    end
    return W
end

function Gitter()
    ax = partikelrand_links+((abs(partikelrand_rechts)+(abs(partikelrand_links)))/(2*gittergroeße))
    bx = partikelrand_rechts-((abs(partikelrand_rechts)+(abs(partikelrand_links)))/(2*gittergroeße))
    abx = bx-ax
    xsmall=collect(ax:abx/(gittergroeße-1):bx)
    ay = partikelrand_unten+((abs(partikelrand_oben)+(abs(partikelrand_unten)))/(2*gittergroeße))
    by = partikelrand_oben - ((abs(partikelrand_oben)+(abs(partikelrand_unten)))/(2*gittergroeße))
    aby = by-ay
    ysmall=collect(ay:aby/(gittergroeße-1):by)
    x = zeros(length(xsmall)*length(ysmall))
    y = zeros(length(xsmall)*length(ysmall))
    for i in 1:length(xsmall)
        for j in 1:length(ysmall)
            x[(i-1)*length(xsmall)+j]=xsmall[i]
            y[(i-1)*length(xsmall)+j]=ysmall[j]
        end
    end
    return x,y,[x y],xsmall,ysmall
end

function Wandpartikel(wandwert_links, wandwert_unten, wandwert_rechts)
    wy_unten = zeros(length(ysmall))
    wy_unten = wy_unten .+wandwert_unten
    wand_unten = [xsmall wy_unten]
    wx_links = zeros(length(xsmall))
    wx_links = wx_links.+wandwert_links
    wand_links = [wx_links ysmall]
    wx_rechts = zeros(length(xsmall))
    wx_rechts = wx_rechts.+wandwert_rechts
    wand_rechts = [wx_rechts ysmall]
    return wand_links,wand_unten,wand_rechts
end

function Wandpartikel2()
    wand_x = collect(wandwert_links:0.05:wandwert_rechts)
    wand_y = collect(wandwert_unten:0.05:wandwert_oben)
    hilfswand_unten = zeros(length(wand_x)) .+ wandwert_unten
    hilfswand_links = zeros(length(wand_y)) .+ wandwert_links
    hilfswand_rechts = zeros(length(wand_y)) .+ wandwert_rechts
    wand_unten = [wand_x hilfswand_unten]
    wand_links = [hilfswand_links wand_y]
    wand_rechts = [hilfswand_rechts wand_y]
    return wand_links,wand_unten,wand_rechts
end

function Differenz(xi,xj)
    xi_x = xi[1,1]
    xi_y = xi[2,1]
    xj_x = xj[1,1]
    xj_y = xj[2,1]
    minimum_x = xi_x - xj_x
    minimum_y = xi_y - xj_y
    return [minimum_x minimum_y]
end

function Dichte(pos)
    rho = zeros(size(pos,1))
    for i in 1:partikelanzahl
        sum = 0
        xi = pos[i,:]
        for j in 1:partikelanzahl
            xj = pos[j,:]
            sum = sum + m*Wpoly6(xi,xj)
        end

        "Dichte zu den Wandpartikeln"
        
        for p in 1:size(wand_unten,1)
            xp = wand_unten[p,:]
            sum = sum+m*Wpoly6(xi,xp)
        end
        for k in 1:size(wand_links,1)
            xk = wand_links[k,:]
            sum = sum+m*Wpoly6(xi,xk)
        end
        for l in 1:size(wand_rechts,1)
            xl = wand_rechts[l,:]
            sum = sum+m*Wpoly6(xi,xl)
        end
        rho[i] = sum
    end
    return rho
end

function Druck(rho)
    P = k*(rho.*rho)
    "P = k * (((rho./reference_density).^7)-1)"
    return P
end

function F_pressure(rho,P,pos)
    f_pres = zeros(partikelanzahl,2)
    for i in 1:partikelanzahl
        sum = zeros(1,2)
        xi = pos[i,:]
        for j in 1:partikelanzahl
            if i!=j
                xj = pos[j,:]
                sum = sum + m*((P[i]/(rho[i].*rho[i]))+(P[j]/(rho[j].*rho[j]))) .* dWspiky(xi,xj)
            end
        end
        f_pres[i,1] = sum[1,1]
        f_pres[i,2] = sum[1,2]
    end
    return -f_pres
end

function F_viscosity(rho,pos,v)
     f_vis = zeros(partikelanzahl,2)
     for i in 1:partikelanzahl
         sum = zeros(2,1)
         xi = pos[i,:]
         for j in 1:partikelanzahl
             if i !=j
                 xj = pos[j,:]
                 sum = sum + (m/rho[j])*d2Wviscosity(xi,xj).*v[j,:]
                 "sum = sum + (m/rho[j])*d2Wviscosity(xi,xj).*(v[j,:]-v[i,:])"
             end
         end
         f_vis[i,1] = nu*m.*sum[1,1]
         f_vis[i,2] = nu*m.*sum[2,1]
     end
     return f_vis
end

"Initialisierung ds Gitters"

partikelrand_links = -0.95
partikelrand_rechts = -0.35
partikelrand_unten = -0.95
partikelrand_oben = 0.95
gittergroeße = 20
wandwert_rechts = 1
wandwert_links = -1
wandwert_unten = -0.95
wandwert_oben = 4
x,y,anfangspos,xsmall,ysmall = Gitter()
partikelanzahl = length(x)
wand_links,wand_unten,wand_rechts = Wandpartikel2()

"Initialisierung der Konstanten"

star_mass = 2
m = star_mass/partikelanzahl
h = 0.1
nu = 1.0087
n = 1
k = 0.1
reference_density = 2.861

"Initialisierung der Variablen"

t_end = 0
dt = 0.004
v = ones(partikelanzahl,2)
v[:,1] = v[:,1].*0
v[:,2] = v[:,2].*0
global pos = anfangspos
global f_grav = zeros(partikelanzahl,2)
f_grav[:,2] .= -9.81

"Startbedingungen"

global rho = Dichte(pos)
global P = Druck(rho)
global f_pres = F_pressure(rho,P,pos)
global f_vis = F_viscosity(rho,pos,v)
global a = -nu*v + f_pres + f_grav


"Leap-frog-time-integration"
global v_minus = v - (dt/2)*a

for i in 1:t_end
    global a = -nu*v + f_pres + f_grav
    global v_plus = v_minus + a*dt
    global v = 0.5*(v_minus+v_plus)
    global pos = pos + v_plus*dt
    global v_minus = v_plus

    "Statische Ränder"
    for j in 1:partikelanzahl
        if pos[j,1] .> wandwert_rechts
            pos[j,1] = 2*wandwert_rechts - pos[j,1]
            v[j,1] = -v[j,1]
        end
        if pos[j,1] .< wandwert_links
            pos[j,1] = 2*wandwert_links - pos[j,1]
            v[j,1] = -v[j,1]
        end
        if pos[j,2] .< wandwert_unten
            pos[j,2] = 2*wandwert_unten - pos[j,2]
            "v[j,2] = -v[j,2]"
        end
    end
end


x = pos[:,1]
y = pos[:,2]
scatter(x,y,title="Partikelbewegung")
