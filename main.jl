using JuMP
using DataFrames
using Requests
Sbase = 100
include("Input.jl") #Input data for the distribution network
include("opf.jl") #opf solution (LinDistflow)
include("mdp.jl") #mdp file


buses, lines, generators, rPTDF = 
    load_case_data()
    tic()
    time_det = zeros(4) # getsolvetime(m::Model)

numbuses = collect(keys(buses))
numlines = collect(keys(lines))
gen_set = collect(keys(generators))
gen_bus = []
for g in keys(generators)
    push!(gen_bus, generators[g].bus_idx)
end

T = 24 #Time Horizon
en_price = [91.0777 69.516 96.6204 97.8976 73.4845 50 66 77 80 91 73 68 90 86 54 76 82 91 54 69 82 76 76 82] #Random price for T

bus_MDP = [17 20 23 26] #Location for the mdp ensemble
Np=4 #number of up & down states

load_MDP_o_Pd = zeros(length(bus_MDP)) 
load_MDP_o_Qd = zeros(length(bus_MDP)) 

for i=1:length(bus_MDP)
load_MDP_o_Pd[i] = buses[bus_MDP[i]].d_P
load_MDP_o_Qd[i] = buses[bus_MDP[i]].d_Q
end  
frac_MDP = linspace(.1,2,8)

load_MDP_p = zeros(length(bus_MDP),T)  
load_MDP_q = zeros(length(bus_MDP),T)  
cong_price_p = zeros(length(bus_MDP),T)
cong_price_q = zeros(length(bus_MDP),T) 

g_uni=1; 
g_pen=1*ones(Np*2,Np*2) #replace the factor 1 here with 10 for nonuniform case
g_pen[1,2*Np]=1
for a=1:(2*Np-1)
    g_pen[a+1,a]=1
end 

Pt=[0.2  0.5 0.1 0.03 0.02 0.03 0.1 0.02
    0.02 0.2 0.5 0.1 0.03 0.02 0.03 0.1
    0.1 0.02 0.2 0.5 0.1 0.03 0.02 0.03
    0.03 0.1 0.02 0.2 0.5 0.1 0.03 0.02
    0.02 0.03 0.1 0.02 0.2 0.5 0.1 0.03
    0.03 0.02 0.03 0.1 0.02 0.2 0.5 0.1
    0.1 0.03 0.02 0.03 0.1 0.02 0.2 0.5
    0.5 0.1 0.03 0.02 0.03 0.1 0.02 0.2]


U=zeros(Np*2,T)
U_2=zeros(Np*2,T)

load_MDP_pold = ones(length(bus_MDP),T)
load_MDP_qold = ones(length(bus_MDP),T)
priceold = load_MDP_pold
pricenew = load_MDP_qold

rho_final = zeros(Np*2,T,length(bus_MDP))
rho_old = rho_final

N_iter =4;

conv_P=zeros(N_iter,T,length(bus_MDP))
conv_Q=zeros(N_iter,T,length(bus_MDP))

conv_lambdaP=zeros(N_iter,T,length(bus_MDP))
conv_lambdaQ=zeros(N_iter,T,length(bus_MDP))

opf_obj=zeros(N_iter,T)
rho_opt = zeros(T,Np*2)
rho_opt2 = zeros(T,Np*2)
vpbus = zeros(T,33)
#var_P = zeros(T)
v_v = zeros(T)
g_v = zeros(T)
voltage_profile = zeros(T)
losses = zeros(T)

#price =zeros(Np*2)
for iter = 1:N_iter
        #tic()
for i = 1:length(bus_MDP)
    for t=1:T
            price = en_price[t] + cong_price_p[i,t]*load_MDP_o_Pd[i]*frac_MDP'' + cong_price_q[i,t]*load_MDP_o_Qd[i]*frac_MDP''
            println(price)
            U[:,t]=price

    end

            println("****************************")
            println("Cost for iteration $iter = ", U)
            println("****************************")

        rho = CostMDP(Np,T,U,Pt,g_pen)

        load_MDP_p[i,:] = (load_MDP_o_Pd[i]*frac_MDP')*rho  
        load_MDP_q[i,:] = (load_MDP_o_Qd[i]*frac_MDP')*rho
        rho_final[:,:,i] = rho
        rho_opt = rho
end
println("U = ", U)
println("g = ", g_pen)

    load_MDP_pold = load_MDP_p;
    load_MDP_qold = load_MDP_q;

        for t = 1:T 
for i = 1:length(bus_MDP)
        buses[bus_MDP[i]].d_P = load_MDP_p[i,t]#[:,t]
        buses[bus_MDP[i]].d_Q = load_MDP_q[i,t]#[:,t]
end
        yp, yq, obj_val, vv, gg, voltage_p, v_p_b = opf_RRC(t)
        v_v[t] = vv
        g_v[t] = gg
        opf_obj[iter,t] = obj_val
        voltage_profile = voltage_p
        losses[t] = obj_val
        for i =1:33
        vpbus[t,i] = v_p_b[i]
        end
for i = 1:length(bus_MDP)
        cong_price_p[:,t] = yp[bus_MDP[i]]#*0.0010#*0.1#00 
        cong_price_q[:,t] = yq[bus_MDP[i]]#*0.0010#*0.1#00
        println("")
        println("Pd at time $t for iteration $iter at $bus_MDP = ", buses[bus_MDP[i]].d_P*Sbase)
                println("ypg = ", yp[17])

        println("")
        conv_P[iter,t,i] = buses[bus_MDP[i]].d_P#*Sbase
        conv_Q[iter,t,i] = buses[bus_MDP[i]].d_Q#*Sbase
        conv_lambdaP[iter,t,i] = yp[bus_MDP[i]]
        conv_lambdaQ[iter,t,i] = yq[bus_MDP[i]]
end
        end
println("*************************************************")
end
for i=1:length(bus_MDP)
println("Bus $(bus_MDP[i])")
for t=1:T
println("For T = $t")
println("conv_P = ", conv_P[:,t,i])
println("")
println("conv_Q = ", conv_Q[:,t,i])
println("")
println("conv_lambdaP = ", conv_lambdaP[:,t,i])
println("")
println("conv_lambdaQ = ", conv_lambdaQ[:,t,i])
println("")
println("")
end
println("**************************************************")
end
for t=1:T
println("For T = $t")
println("opf_obj = ", opf_obj[:,t])
println("")
end
println("")
println("rho = ", rho_opt)
println("")
for t=1:8
println("state $t = ", rho_opt[t,:])
end
println("")
for t=1:T
println("For T = $t")
println("No. of voltage violated constraints = ", v_v[t])
@printf("%.2f%% of voltage constraints hold  \n", ((1-(v_v[t]/(2*500*n_buses)))*100))
@printf("%.2f%% of voltage constraints not hold  \n", 100-((1-(v_v[t]/(2*500*n_buses)))*100))
println("")
println("No. of generation violated constraints = ", g_v[t])
@printf("%.2f%% of generation constraints hold \n", ((1-(g_v[t]/(4*500*length(gen_bus))))*100))
@printf("%.2f%% of generation constraints not hold \n", 100-((1-(g_v[t]/(4*500*length(gen_bus))))*100))
println("")
end
println("Time taken (seconds): ", toq())
println("End") 
