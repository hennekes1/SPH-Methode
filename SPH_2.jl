import LinearAlgebra
import Plots
using Plots
using LinearAlgebra

function Wcub(r,h)
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

function dWcub(r,h)
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

function d2Wcub(r,h)
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

function Wpoly6(r,h)
    sigma = 315/(64*pi*h^9)
    q = norm(r)
    if 0<=q<=h
        W = sigma*(h^2-q^2)^3
    else
        W=0
    end
    return W
end

function dWpoly6(r,h)
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

function d2Wpoly6(r,h)
    sigma = 945/(8*pi*h^9)
    q = norm(r)
    if 0<=q<h
        W = sigma * (q^2 - 0.75(h^2-q^2))*(h^2-q^2)
    else
        W=0
    end
    return W
end

function Wspiky(r,h)
    sigma = 15/(pi*h^6)
    q = norm(r)
    if 0<=q<h
        W = sigma * (h-q)^3
    else
        W=0
    end
    return W
end

function dWspiky(r,h)
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

function Wviscosity(r,h)
    sigma = 15/(2*pi*h^3)
    q = norm(r)
    if 0<=q<=h
        W = sigma * ((-q^3)/(2*h^3)+(q^2/h^2)+(h/2*q)-1)
    else
        W=0
    end
    return W
end

function d2Wviscosity(r,h)
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
    xsmall=collect(-2:4/(gittergroeße-1):2)
    ysmall=collect(0:1.5/(gittergroeße-1):1.5)
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
    return norm(xi-xj)
end

function Dichte(pos,m,h)
    rho = zeros(size(pos,1))
    for i in 1:length(rho)
        sum = 0
        xi = pos[i,:]
        for j in 1:size(pos,1)
            xj = pos[j,:]
            sum = sum + m*Wpoly6(Abstand(xi,xj),h)
        end
        rho[i] = sum
    end
    return rho
end


function Druck(rho,k,n)
    P = k * rho.^(1+1/n)
    return P
end

function F_Pressure(pos,v,m,h,k,n)
    rho = Dichte(pos,m,h)
    P = Druck(rho,k,n)
    f_pres = zeros(partikelanzahl,1)
    for i in 1:partikelanzahl
        xi = pos[i,:]
        for j in 1:partikelanzahl
            xj = pos[j,:]
            f_pres[i] = f_pres[i]+m*((P[i] / (rho[i].*rho[i]))+(P[j] / (rho[j].*rho[j])))*dWspiky(Abstand(xi,xj),h)
        end
    end
    return f_pres
end

function F_Viscosity(pos,v,m,h,k,n,nu)
    f_vis = zeros(partikelanzahl,1)
    for i in 1:partikelanzahl
        xi = pos[i,:]
        for j in 1:partikelanzahl
            xj = pos[j,:]
            f_vis[i] = f_vis[i] + m*((v[j,1]-v[i,1])/rho[j])*d2Wviscosity(Abstand(xi,xj),h)
        end
    end
    f_vis = f_vis .* (nu./rho)
    return f_vis
end

m = 0.005
h = 0.1
gittergroeße= 20
x,y,anfangspos,xsmall,ysmall = Gitter()
partikelanzahl = length(x)
"Pressure constant"
k = 0.1
"Polytropic index"
n = 1
"Anfangsgeschwindigkeiten der Partikel"
v = ones(partikelanzahl,2)
v = v*2
v[:,2] .= v[:,2] * 0
"Viskositätskonstante"
nu = 1
"external force constant--Wird hier nicht berechnet!"
lambda = 2.01
"Schwerkraftkonstante"
g = [0 ;-9.81]
"Zeitintegration"
dt = 0.04
time_steps = 10
x_werte = zeros(partikelanzahl,time_steps+1)
y_werte = zeros(partikelanzahl,time_steps+1)
pos = [x[:] y[:]]
rho = Dichte(pos,m,h)
for j in 1:time_steps
    x_werte[:,j] = x[:]
    y_werte[:,j] = y[:]
    global pos = [x[:] y[:]]
    global rho = Dichte(pos,m,h)
    global P = Druck(rho,k,n)
    global f_pres = F_Pressure(pos,v,m,h,k,n)
    global f_vis = F_Viscosity(pos,v,m,h,k,n,nu)
    for i in 1:partikelanzahl
        global delta_v = -f_pres./rho
         ".+ f_vis"
        v[i,1] = v[i,1] + dt*delta_v[i]
        "v[i,2] = v[i,2] + dt*delta_v[i]"

        x[i] = x[i] + dt*v[i,1]
        y[i] = y[i] + dt*v[i,2]

        "Periodische Ränder"
        if x[i] .> 2
            x[i] = 2*2-x[i]
            v[i,1] = v[i,1] * (-1)
        end
        if x[i] .< -2
            x[i] = 2*(-2) - x[i]
            v[i,1] = v[i,1] * (-1)
        end
        if y[i] .> 1.5
            y[i] = 2*1.5 - y[i]
            v[i,2] = v[i,2] *(-1)
        end
        if y[i] .< 0
            y[i] = 2*0 - y[i]
            v[i,2] = v[i,2] *(-1)
        end
        if x[i] .> 2
            x[i] = 2*2-x[i]
            v[i,1] = v[i,1] * (-1)
        end
        if x[i] .< -2
            x[i] = 2*(-2) - x[i]
            v[i,1] = v[i,1] * (-1)
        end
        if y[i] .> 1.5
            y[i] = 2*1.5 - y[i]
            v[i,2] = v[i,2] *(-1)
        end
        if y[i] .< 0
            y[i] = 2*0 - y[i]
            v[i,2] = v[i,2] *(-1)
        end
    end
end
x_werte[:,time_steps+1] = x[:]
y_werte[:,time_steps+1] = y[:]
global pos = [x[:] y[:]]
global rho = Dichte(pos,m,h)
global P = Druck(rho,k,n)
global f_pres = F_Pressure(pos,v,m,h,k,n)

scatter(x_werte, y_werte,layout=(3,2))
scatter(x_werte[:,time_steps+1], y_werte[:,time_steps+1])
