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
    ax = linker_rand+((abs(rechter_rand)+(abs(linker_rand)))/(2*gittergroeße))
    bx = rechter_rand-((abs(rechter_rand)+(abs(linker_rand)))/(2*gittergroeße))
    abx = bx-ax
    xsmall=collect(ax:abx/(gittergroeße-1):bx)
    ay = unterer_rand+((abs(oberer_rand)+(abs(unterer_rand)))/(2*gittergroeße))
    by = oberer_rand - ((abs(oberer_rand)+(abs(unterer_rand)))/(2*gittergroeße))
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

function Differenz(xi,xj)
    xi_x = xi[1,1]
    xi_y = xi[2,1]
    xj_x = xj[1,1]
    xj_y = xj[2,1]
    minimum_x = xi_x - xj_x
    minimum_y = xi_y - xj_y
    if xj_x < xi_x
        if (xj_x - linker_rand)+(rechter_rand - xi_x) < abs(minimum_x)
            minimum_x = -((xj_x-linker_rand)+(rechter_rand - xi_x))
        end
    else
        if (rechter_rand-xj_x)+(xi_x-linker_rand) < abs(minimum_x)
            minimum_x = ((rechter_rand-xj_x)+(xi_x-linker_rand))
        end
    end
    if xj_y < xi_y
        if (xj_y - unterer_rand)+(oberer_rand - xi_y) < abs(minimum_y)
            minimum_y = -((xj_y-unterer_rand)+(oberer_rand - xi_y))
        end
    else
        if (oberer_rand-xj_y)+(xi_y-unterer_rand) < abs(minimum_y)
            minimum_y = ((oberer_rand-xj_y)+(xi_y-unterer_rand))
        end
    end
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
        rho[i] = sum
    end
    return rho
end

function Druck(rho)
    P = k*(rho.*rho)
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
             end
         end
         f_vis[i,1] = nu*m.*sum[1,1]
         f_vis[i,2] = nu*m.*sum[2,1]
     end
     return f_vis
end

"Initialisierung ds Gitters"

linker_rand = -2
rechter_rand = 2
unterer_rand = 0
oberer_rand = 1.5
gittergroeße = 20
x,y,anfangspos,xsmall,ysmall = Gitter()
partikelanzahl = length(x)

"Initialisierung der Konstanten"

star_mass = 2
m = star_mass/partikelanzahl
h = 0.1
nu = 1.0087
n = 1
k = 0.1

"Initialisierung der Variablen"

t_end = 0
dt = 0.01
v = ones(partikelanzahl,2)
v = v.*0.1
global pos = anfangspos

"Startbedingungen"

if t_end == 0
    global rho = Dichte(pos)
    global P = Druck(rho)
    global f_pres = F_pressure(rho,P,pos)
    global f_vis = F_viscosity(rho,pos,v)
end

"Zeitintegration"

for i in 1:t_end
    global rho = Dichte(pos)
    global P = Druck(rho)
    global f_pres = F_pressure(rho,P,pos)
    global f_vis = F_viscosity(rho,pos,v)
    global v = v + (dt/1).*(f_vis+f_pres)
    global pos = pos + dt.*v

    "Periodische Ränder"
    for i in 1:partikelanzahl
        if pos[i,1] .> rechter_rand
            pos[i,1] = pos[i,1].-(abs(rechter_rand)+(abs(linker_rand)))
        end
        if pos[i,1] .< linker_rand
            pos[i,1] = pos[i,1].+(abs(rechter_rand)+(abs(linker_rand)))
        end
        if pos[i,2] .> oberer_rand
            pos[i,2] = pos[i,2].-(abs(oberer_rand)+(abs(unterer_rand)))
        end
        if pos[i,2] .< unterer_rand
            pos[i,2] = pos[i,2].+(abs(oberer_rand)+(abs(unterer_rand)))
        end
    end

end

x = pos[:,1]
y = pos[:,2]
scatter(x,y,title="Partikelbewegung")
