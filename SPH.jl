import LinearAlgebra
import Plots
using Plots
using LinearAlgebra

"function Wgaus(r, h)
    q = norm(r)^2 / h^2
    return (1.0 / (h^3 * (pi^(3/2)))) * exp(-q)
end

function gradWgaus(r,h)
    q = norm(r)^2 / h^2
    return r*(-2.0 / (h^5 * (pi^(3/2)))) * exp(-q)
end"

function A1(x,y)
    return 0.5*(x+y)
end

function A2(x,y)
    return 0.5*(x^2+y^2)-1
end

function A3(x,y)
    return sin(5*x)*cos(3*y)
end

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
        return W=W*(-1)
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

function Dichte()
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

function A1_dis(x,y)
    sum = 0
    for j in 1:size(pos,1)
        xj = pos[j,:]
        xi = [x ; y]
        sum = sum + A1(xj[1],xj[2])*(m/rho[j])*W(xi-xj,h)
    end
    return sum
end

function A2_dis(x,y)
    sum = 0
    for j in 1:size(pos,1)
        xj = pos[j,:]
        xi = [x ; y]
        sum = sum + A2(xj[1],xj[2])*(m/rho[j])*W(xi - xj,h)
    end
    return sum
end

function A3_dis(x,y)
    sum = 0
    for j in 1:size(pos,1)
        xj = pos[j,:]
        xi = [x ; y]
        sum = sum + A3(xj[1],xj[2])*(m/rho[j])*W(xi - xj,h)
    end
    return sum
end

m = 0.005
h = 0.1
sigma = 40 / (7*pi*h^2)
gittergroeße= 20
x,y,pos,xsmall,ysmall = Gitter()
partikelanzahl = length(x)
rho = Dichte()

s = zeros(gittergroeße,2)
for i in 1:size(s,1)
    s[i,1] = xsmall[i]
    s[i,2] = ysmall[i]
end

a1 = zeros(size(s,1))
a2 = zeros(size(s,1))
a3 = zeros(size(s,1))
a1_dis = zeros(size(s,1))
a2_dis = zeros(size(s,1))
a3_dis = zeros(size(s,1))
for i in 1:size(s,1)
    a1[i] = A1(s[i,1],s[i,2])
    a2[i] = A2(s[i,1],s[i,2])
    a3[i] = A3(s[i,1],s[i,2])
    a1_dis[i] = A1_dis(s[i,1],s[i,2])
    a2_dis[i] = A2_dis(s[i,1],s[i,2])
    a3_dis[i] = A3_dis(s[i,1],s[i,2])
end

a = [a1 a2 a3]
a_dis = [a1_dis a2_dis a3_dis]
a_komplett = [a a_dis]

a1_error = a1 - a1_dis
a2_error = a2 - a2_dis
a3_error = a3 - a3_dis

a_error = [a1_error a2_error a3_error]

function Druck()
    P = k * rho.^(1+1/n)
    return P
end

"Pressure constant"
k = 0.1
"Polytropic index"
n = 1
P = Druck()
"Anfangsgeschwindigkeiten der Partikel"
v = rand(partikelanzahl)
"damping"
nu = 1
"external force constant--Wird hier nicht berechnet!"
lambda = 2.01
lamda = 0
dt = 0.04

function Beschleunigung()
    a = zeros(partikelanzahl)
    for i in 1:partikelanzahl
        xi = pos[i,:]
        for j in (i+1):partikelanzahl
            xj = pos[j,:]
            a[i] = a[i]-m*((P[i] / rho[i].*rho[i])+(P[j] / rho[j].*rho[j]))*gradW(Abstand(xi,xj),h)
        end
    end
    return a
end

a = Beschleunigung()

for i in 1:partikelanzahl
    v[i] = v[i] + dt * a[i]
    x[i] = x[i] + dt*v[i]
end

scatter(x, y, title = "Partikelbewegung", label = "Partikel", lw = 3)
