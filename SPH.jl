using Plots
using LinearAlgebra
using StaticArrays

function Wpoly6(xi,xj,h)
    #sigma = 315/(64*pi*h^9)
    sigma = 4 / (pi*h^8)
    q = norm(Differenz(xi,xj))
    if 0<=q<=h
        W = sigma*(h^2-q^2)^3
    else
        W=0.0
    end
    return W
end

function dWspiky(xi,xj,h)
    r = Differenz(xi,xj)
    q = norm(r)
    #sigma = -45/(q*pi*h^6)
    sigma = - 30 / (pi * h^5)
    if 0<=q<=h
        sigma2 = (sigma * (h-q)^2 )/q
        W = sigma2 * r
    else
        W=SVector(0.0, 0.0)
    end

    return W
end

function d2Wviscosity(xi,xj,h)
    #sigma = 45/(pi*h^6)
    sigma = 40 / (pi * h^5)
    q = norm(Differenz(xi,xj))
    if 0<=q<=h
        W = sigma *(h-q)
    else
        W=0.0
    end
    return W
end

function Gitter(partikelrand_links, partikelrand_rechts, partikelrand_unten, partikelrand_oben,
                gittergroeße)
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


function Differenz(xi,xj)
    xi_x = xi[1,1]
    xi_y = xi[2,1]
    xj_x = xj[1,1]
    xj_y = xj[2,1]
    minimum_x = xi_x - xj_x
    minimum_y = xi_y - xj_y
    return SVector(minimum_x, minimum_y)
end

function Dichte(pos, m, h)
    rho = zeros(size(pos,1))
    for i in 1:length(rho)
        xi = SVector(pos[i, 1], pos[i, 2])
        rho[i] = rho[i] + m*Wpoly6(xi,xi,h)
        for j in (i+1):length(rho)
            xj = SVector(pos[j, 1], pos[j, 2])
            rho_ij = m*Wpoly6(xi,xj,h)
            rho[i] = rho[i] + rho_ij
            rho[j] = rho[j] + rho_ij
        end
    end
    return rho
end

function Druck(rho,rest_density,gas_konstante)
    #P = k*(rho.-rest_density)
    P = (gas_konstante ^2 * (rest_density/7)) .* ((rho./rest_density).*(rho./rest_density).*(rho./rest_density).*(rho./rest_density).*(rho./rest_density).*(rho./rest_density).*(rho./rest_density) .-1)
    return P
end

function F_pressure(rho,P,pos,m,h)
    partikelanzahl = length(rho)
    f_pres = zeros(partikelanzahl,2)
    for i in 1:partikelanzahl
        sum = SVector(0.0, 0.0)
        xi = SVector(pos[i, 1], pos[i, 2])
        for j in 1:partikelanzahl
            if i!=j
                xj = SVector(pos[j, 1], pos[j, 2])
                sum = sum -0.5* (((P[i]+P[j]))*dWspiky(xi,xj,h))./(rho[j])
            end
        end
        f_pres[i,1] = m*sum[1]
        f_pres[i,2] = m*sum[2]
    end
    return f_pres
end

function F_viscosity(rho,pos,v,m,mu,h)
     partikelanzahl = length(rho)
     f_vis = zeros(partikelanzahl,2)
     for i in 1:partikelanzahl
         sum1 = 0.0
         sum2 = 0.0
         xi = SVector(pos[i, 1], pos[i, 2])
         for j in 1:partikelanzahl
             if i !=j
                 xj = SVector(pos[j, 1], pos[j, 2])
                 sum1 = sum1 + d2Wviscosity(xi,xj,h).*(v[j,1].-v[i,1])/rho[j]
                 sum2 = sum2 + d2Wviscosity(xi,xj,h).*(v[j,2].-v[i,2])/rho[j]
             end
         end
         f_vis[i,1] = m.*mu.*sum1
         f_vis[i,2] = m.*mu.*sum2
     end
     return f_vis
end

function F_gravity(partikelanzahl,rho)
    f_grav = zeros(partikelanzahl,2)
    f_grav[:,2] .= -9.81.*rho
    return f_grav
end

"Initialisierung ds Gitters"

function main()
partikelrand_links = 0
partikelrand_rechts = 2
partikelrand_unten = 0
partikelrand_oben = 2


gittergroeße = 25


wandwert_rechts = 4
wandwert_links = 0
wandwert_unten = 0
wandwert_oben = 4
x,y,anfangspos,xsmall,ysmall = Gitter(partikelrand_links, partikelrand_rechts, partikelrand_unten,
                                      partikelrand_oben, gittergroeße)
partikelanzahl = length(x)
#dy = (wandwert_oben - wandwert_unten) / (2*gittergroeße)
#dx = (wandwert_rechts - wandwert_links) / (2*gittergroeße)
dy = (partikelrand_oben - partikelrand_unten) / gittergroeße
dx = (partikelrand_rechts - partikelrand_links) / gittergroeße

"Initialisierung der Konstanten"

star_mass = 2
#m = star_mass/partikelanzahl
#m = 0.0064
m = dx * dy
#h = 0.0457
h = 2 * dx
mu = 0.1
k = 0.5
rest_density = 1
restitution = -0.5
gas_konstante = 10


"Initialisierung der Variablen"

t_end= 6000
dt = 0.0008
t = 0
v = ones(partikelanzahl,2)
v[:,1] = v[:,1].*2.5
v[:,2] = v[:,2].*(-2)
#v[:,1] = v[:,1].*0
#v[:,2] = v[:,2].*0
pos = anfangspos

"Startbedingungen"
#rho = Dichte(pos, m, h)
rho = ones(partikelanzahl)
P = Druck(rho,rest_density,gas_konstante)
f_pres = zeros(partikelanzahl,2)
f_vis = zeros(partikelanzahl,2)
f_grav = zeros(partikelanzahl,2)
f_gesamt = f_pres + f_vis + f_grav
time_plot = 1


"Leap-frog-time-integration"

for i in 1:t_end

    t = t + dt
    v = v + dt * f_gesamt./rho
    pos = pos + dt * v

    "Statische Ränder"
    for j in 1:partikelanzahl
        if pos[j,1] .> wandwert_rechts
        #if pos[j,1] .> (wandwert_rechts-dx/2)
            pos[j,1] = 2*wandwert_rechts - pos[j,1]
            #pos[j,1] = wandwert_rechts - dx/2
            if v[j,1] .> 0
                #v[j,1] = v[j,1] -(1+restitution)*v[j,1]
                v[j,1] = restitution * v[j,1]
            end
        end
        if pos[j,1] .< wandwert_links
        #if pos[j,1] .< (wandwert_links+dx/2)
            pos[j,1] = 2*wandwert_links - pos[j,1]
            #pos[j,1] = wandwert_links +dx/2
            if v[j,1] .< 0
            #    v[j,1] = v[j,1] -(1+restitution)*v[j,1]
                v[j,1] = restitution * v[j,1]
            end
        end
        if pos[j,2] .< wandwert_unten
        #if pos[j,2] .< (wandwert_unten+dy/2)
            pos[j,2] = 2*wandwert_unten - pos[j,2]
            #pos[j,2] = wandwert_unten +dy/2
            if v[j,2] .< 0
                #v[j,2] = v[j,2] -(1+restitution)*v[j,2]
                v[j,2] = restitution * v[j,2]
            end
        end
        if pos[j,2] .> wandwert_oben
        #if pos[j,2] .> wandwert_oben-dy
        #if pos[j,2] .> (wandwert_oben-dy/2)
            pos[j,2] = 2*wandwert_oben - pos[j,2]
            #pos[j,2] = wandwert_oben - dy/2
            if v[j,2] .> 0
                #v[j,2] = v[j,2] -(1+restitution)*v[j,2]
                v[j,2] = restitution * v[j,2]
            end
        end
    end
    rho = Dichte(pos, m, h)
    P = Druck(rho,rest_density,gas_konstante)
    f_pres = F_pressure(rho,P,pos,m,h)
    f_vis = F_viscosity(rho,pos,v,m,mu,h)
    f_grav = F_gravity(partikelanzahl,rho)
    f_gesamt = f_pres + f_vis + f_grav
    #f_gesamt = f_grav
    if i == time_plot
        x = pos[:,1]
        y= pos[:,2]
        display(scatter(xlims=(wandwert_links,wandwert_rechts), ylims=(wandwert_unten,wandwert_oben), x,y))
        println("Zeit ist gerade bei: " , t)
        #savefig("SPH-Methode\\test" *  ".png")
        time_plot = time_plot + 10
        #println("Rho")
    #    println(sum(rho))

    end
#println("Die Zeit ist gerade:")
#println(t)
#println("Rho")
#println(minimum(rho))
#println(maximum(rho))
#println(sum(rho))
#println("Druck")
#println(minimum(P))
#println(maximum(P))
#println(sum(P))
#println("Kräfte")
#println(minimum(f_gesamt[:,1]))
#println(maximum(f_gesamt[:,1]))
#println(sum(f_gesamt[:,1]))
#println(minimum(f_gesamt[:,2]))
#println(maximum(f_gesamt[:,1]))
#println(sum(f_gesamt[:,2]))

end

return rho,P,f_gesamt[:,1],f_gesamt[:,2],v[:,1],v[:,2],pos[:,1],pos[:,2],dx,dy,xsmall,ysmall
end # function main()

rho,P,fx,fy,vx,vy,x,y,dx,dy,xsmall,ysmall = main()

#savefig("SPH-Methode\\")
