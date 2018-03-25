

function CostMDP(Np,T,U,Pt,g)

D, V = eig(Pt)
#println("V = ", V)

rhot = real(V[:,1])
#println("")
#println("rhot = ", rhot)
rhot = rhot/sum(rhot)
#println("rhot = ", rhot)

const lambda1 = -20

if size(g,1) != 1

	delta=0.1
    Popt=zeros(Np*2,Np*2,T-1)
    phi=zeros(Np*2,T)
    phi[:,T]=U[:,T]

        for t=1:T-1 
        for a=1:(2*Np)
            lambda=lambda1
            for it=1:10000
                norm=0
                for c=1:2*Np
                    norm=norm+exp(-(phi[c,T-t+1]-lambda)/g[c,a])*Pt[c,a]
                end

                norm = norm/exp(1)
                lambda=lambda-delta*(norm-1)
            end
            for b=1:(2*Np)
                Popt[b,a,T-t]=exp(-(phi[b,T-t+1]-lambda)/g[b,a])*Pt[b,a]
                Popt[b,a,T-t]=exp(-(phi[b,T-t+1]-lambda)/g[b,a]-1)*Pt[b,a]

            end
        end
        for a=1:(2*Np)
            aux=0
            for b=1:(2*Np)
                aux=aux+(phi[b,T-t+1]+g[b,a]*log(Popt[b,a,T-t]/Pt[b,a]))*Popt[b,a,T-t]
            end
            phi[a,T-t]=aux+U[a,T-t]
        end
    end

    rho=zeros(2*Np,T)
    rho[:,1]=rhot
    for t=1:(T-1)
        for a=1:(2*Np)
            r=0
            for b=1:(2*Np)
                r=r+Popt[a,b,t]*rho[b,t]
            end
            rho[a,t+1]=r
        end
    end
    
elseif g ==1
    
    u=exp(-U) 
    for t=1:(T-1)
        for a=1:(2*Np)
            uc=0
            for b=1:(2*Np)
                uc=uc+u[b,T-t+1]*Pt[b,a]
            end
            u[a,T-t]=uc*exp(-U[a,T-t])
        end
    end
    
    Popt=zeros(Np*2,Np*2,T-1)
    for t=1:T-1
        for a=1:(2*Np)
            up=0
            for c=1:(2*Np)
                up=up+u[c,T-t+1]*Pt[c,a]
            end
            for b=1:(2*Np)
                Popt[b,a,T-t]=u[b,T-t+1]*Pt[b,a]/up
            end
        end
    end
    
    ### computing optimal P - Popt
    rho=zeros(2*Np,T)
    rho[:,1]=rhot 
    for t=1:(T-1)
        for a=1:(2*Np)
            r=0
            for b=1:(2*Np)
                r=r+Popt[a,b,t]*rho[b,t]
            end
            rho[a,t+1]=r
        end
    end
end


#println("Popt = ", Popt)
println("")
println("rho = ", rho)
println("rho_size", size(rho))
println("")
for t=1:8
println("state $t = ", rho[t,:])
end
return rho
end