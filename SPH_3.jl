import LinearAlgebra
import Plots
using Plots
using LinearAlgebra


function Wcub(r)
    sigma = 40 / (7*pi*h^2)
    q = 1/h * norm(r)
    if 0<=q<=0.5
        return sigma*(6*(q^3-q^2)+1)
    elseif 0.5<q<=1
        return sigma*(2*(1-q)^3)
    else
        return 0
    end
end

function dWcub(r)
    sigma = 40 / (7*pi*h^2)
    p = r
    q = 1/h * norm(r)
    if 0<=q<=0.5
        W= -sigma*(6*(3*q^2-2*q))
    elseif 0.5<q<=1
        W= sigma*(6*(1-q)^2)
    else
        W= 0
    end
    if p>=0
        W=W*(-1)
    end
    return W
end

function d2Wcub(r)
    sigma = 40 / (7*pi*h^2)
    q = 1/h * norm(r)
    if 0<=q<=0.5
        W= -sigma*(-36*q+12)
    elseif 0.5<q<=1
        W= sigma*(12*(1-q))
    else
        W= 0
    end
    return W
end

function Wpoly6(r)
    sigma = 315/(64*pi*h^9)
    q = norm(r)
    if 0<=q<=h
        W = sigma*(h^2-q^2)^3
    else
        W=0
    end
    return W
end

function dWpoly6(r)
    sigma =945/(32*pi*h^9)
    p=r
    q = norm(r)
    if 0<=q<h
        W = sigma*(-r)*(h^2-q^2)^2
    else
        W=0
    end
    if p>=0
        W=W*(-1)
    end
    return W
end

function d2Wpoly6(r)
    sigma = 945/(8*pi*h^9)
    q = norm(r)
    if 0<=q<h
        W = sigma * (q^2 - 0.75(h^2-q^2))*(h^2-q^2)
    else
        W=0
    end
    return W
end

function Wspiky(r)
    sigma = 15/(pi*h^6)
    q = norm(r)
    if 0<=q<h
        W = sigma * (h-q)^3
    else
        W=0
    end
    return W
end

function dWspiky(r)
    sigma = -45/(pi*h^6)
    q = norm(r)
    if 0<=q<=h
        W = sigma * (h-q)^2 * r/q
    else
        W=0
    end
    if r.>=0
        W = W*(-1)
    end
    if q<=0.000000001
        W = 14.32394487824193
    end
    return W
end

function Wviscosity(r)
    sigma = 15/(2*pi*h^3)
    q = norm(r)
    if 0<=q<=h
        W = sigma * ((-q^3)/(2*h^3)+(q^2/h^2)+(h/2*q)-1)
    else
        W=0
    end
    return W
end

function d2Wviscosity(r)
    sigma = 45/(pi*h^6)
    q = norm(r)
    if 0<=q<=h
        W = sigma *(h-q)
    else
        W=0
    end
    return W
end

function Gitter()
    xsmall=collect(-1.9:3.8/(gittergroeße-1):1.9)
    ysmall=collect(0.0375:1.425/(gittergroeße-1):1.4625)
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

function Abstand(xi,xj)
    if xi[1]>=1.9 && xj[1]<=(-1.9)
        xi[1] = xi[1] .- 4
    end
    if xi[1] <=(-1.9) && xj[1]>=1.9
        xi[1] = xi[1] .+4
    end
    if xi[2]>=1.4 && xj[2]<=0.1
        xj[2] = xj[2] .+ 1.5
    end
    if xi[2]<=0.1 && xj[2]>=1.4
        xj[2] = xj[2] .- 1.5
    end
    return norm(xi-xj)
end


function Dichte(pos)
    rho = zeros(size(pos,1))
    for i in 1:partikelanzahl
        sum = 0
        xi = pos[i,:]
        for j in 1:partikelanzahl
            xj = pos[j,:]
            sum = sum + m*Wpoly6(Abstand(xi,xj))
        end
        rho[i] = sum
    end
    return rho
end

function Druck(rho)
    P = k*(rho .- rho_ruhe)
    return P
end

function F_pressure(rho,P,pos)
    f_pres = zeros(partikelanzahl)
    for i in 1:partikelanzahl
        sum = 0
        xi = pos[i,:]
        for j in 1:partikelanzahl
            xj = pos[j,:]
            sum = sum + m*((P[i]+P[j])/(2*rho[j])) * dWspiky(Abstand(xi,xj))
        end
        f_pres[i] = -sum
    end
    return f_pres
end

function F_viscosity(rho,pos,v)
     f_vis = zeros(partikelanzahl)
     for i in 1:partikelanzahl
         sum = 0
         xi = pos[i,:]
         for j in 1:partikelanzahl
             xj = pos[j,:]
             sum = sum + m*((v[j]-v[i])/rho[j])*d2Wviscosity(Abstand(xi,xj))
         end
         f_vis[i] = sum*nu
     end
     return f_vis
end


"Initialisierung der Konstanten"
m = 0.005
h = 0.1
rho_ruhe = 7.833407355304225
nu = 1.0087
n = 1
k = 0.1

"Initialberechnung des Gitters"
gittergroeße= 20
x,y,anfangspos,xsmall,ysmall = Gitter()
partikelanzahl = length(x)

"Initialisierung der Variablen"
t_end = 12
dt = 0.04
v = ones(partikelanzahl,2)
v[:,2] = v[:,2] .*0
global pos = anfangspos

"Main-Methode"

x_werte = zeros(partikelanzahl,t_end+1)
y_werte = zeros(partikelanzahl,t_end+1)
pos = [x[:] y[:]]
rho = Dichte(pos)
for j in 1:t_end
    x_werte[:,j] = x[:]
    y_werte[:,j] = y[:]
    global pos = [x[:] y[:]]
    global rho = Dichte(pos)
    global P = Druck(rho)
    global f_pres = F_pressure(rho,P,pos)
    global f_vis = F_viscosity(rho,pos,v)
    for i in 1:partikelanzahl
        global f_komplett = f_pres + f_vis
        global acc = f_komplett ./ rho
        v[i,1] = v[i,1] + (dt/2)*acc[i]
        x[i] = x[i] + dt*v[i,1]
        y[i] = y[i] + dt*v[i,2]

        "Periodische Ränder"
        if x[i] .> 2
            x[i] = x[i].-4
        end
        if x[i] .< -2
            x[i] = x[i].+4
        end
        if y[i] .> 1.5
            y[i] = y[i].-1.5
        end
        if y[i] .< 0
            y[i] = y[i].+1.5
        end
    end
end
x_werte[:,t_end+1] = x[:]
y_werte[:,t_end+1] = y[:]
global pos = [x[:] y[:]]
global rho = Dichte(pos)
global P = Druck(rho)
global f_pres = F_pressure(rho,P,pos)
global f_vis = F_viscosity(rho,pos,v)

scatter(x_werte[:,t_end+1], y_werte[:,t_end+1],title="Partikelbewegung")
