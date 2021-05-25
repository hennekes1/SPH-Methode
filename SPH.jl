import LinearAlgebra
import Plots
using Plots
using LinearAlgebra

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
    q = 1/h * norm(r)
    if 0<=q<=0.5
        return sigma*(6*(q^3-q^2)+1)
    elseif 0.5<q<=1
        return sigma*(2*(1-q)^3)
    else
        return 0
    end
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
        sum = sum + A1(xj[1],xj[2])*(m/rho[j])*W(xi - xj,h)
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

"function Laplace_op()
    for i in 1:size(pos,1)
        sum = 0
        for j in 1:size(pos,1)
            sum = sum + (18/p[j])*A(xi,xj) *2*2te Ablteing W / Abstand(x[i],x[j])
        end
    end
end"

m = 18
h = 0.3
sigma = 40 / (7*pi*h^2)
gittergroeße= 50
x,y,pos,xsmall,ysmall = Gitter()
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

plot(xsmall,a_komplett)
