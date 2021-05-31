import LinearAlgebra
import Plots
using Plots
using LinearAlgebra

function W(r,h)
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

function gradW(r,h)
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
    "wx = W *x
    wy = W * y"
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
        for j in 1:size(pos,1)
            sum = sum + m*W(pos[i,:]-pos[j,:],h)
        end
        rho[i] = sum
    end
    return rho
end


function Druck(rho,k,n)
    P = k * rho.^(1+1/n)
    return P
end

function Beschleunigung(pos,v,m,h,k,n)
    rho = Dichte(pos,m,h)
    P = Druck(rho,k,n)
    a = zeros(partikelanzahl,2)
    for i in 1:partikelanzahl
        xi = pos[i,:]
        for j in (i+1):partikelanzahl
            xj = pos[j,:]
            "a[i,1] = a[i,1]-m*((P[i] / (rho[i].*rho[i]))+(P[j] / (rho[j].*rho[j])))*gradW(Abstand(xi,xj),h)
            a[i,2] = a[i,2]-m*((P[i] / (rho[i].*rho[i]))+(P[j] / (rho[j].*rho[j])))*gradW(Abstand(xi,xj),h)"
        end
    end

    "Dämpfung und externe Kräfte"
    "a[:,1] = a[:,1] - lambda .* pos[:,1] - nu .* v[:,1]
    a[:,2] = a[:,2] - lambda .* pos[:,2] - nu .* v[:,2]"
    return a
end

m = 0.005
h = 0.1
sigma = 40 / (7*pi*h^2)
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
v[:,2] .= v[:,2] * 0.75
"damping"
nu = 1
"external force constant--Wird hier nicht berechnet!"
lambda = 2.01
dt = 0.04

"Zeitintegration"

time_steps = 50
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
    global a = Beschleunigung(pos,v,m,h,k,n)
    for i in 1:partikelanzahl
        v[i,1] = v[i,1] + dt*a[i,1]
        v[i,2] = v[i,2] + dt*a[i,2]
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
    end
end
x_werte[:,time_steps+1] = x[:]
y_werte[:,time_steps+1] = y[:]
global pos = [x[:] y[:]]
global rho = Dichte(pos,m,h)
global P = Druck(rho,k,n)
global a = Beschleunigung(pos,v,m,h,k,n)

scatter(x_werte, y_werte,layout=(3,2))
scatter(x_werte[:,time_steps+1], y_werte[:,time_steps+1])
