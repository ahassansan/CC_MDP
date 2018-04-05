using JuMP, JuMPChance
using Ipopt
using Distributions

buses, lines, generators, rPTDF = 
    load_case_data()
    tic()
    time_det = zeros(4) # getsolvetime(m::Model)
n_buses = 33#length(BUSES)

numbuses = collect(keys(buses))
numlines = collect(keys(lines))
gen_set = collect(keys(generators))
gen_bus = []
for g in keys(generators)
    push!(gen_bus, generators[g].bus_idx)
end

η_v = 0.05 # 1 - Confidence for voltage limit cc
η_g = 0.05 # 1 - Confidence for generation limit cc
sd = 0.3
var = sd^2#0.1
var_true_vector = ones(n_buses)*var
error_variances = var_true_vector

function opf_RRC(time_per)

info("Starting Chance Constraint LinDist Model")

bus_set = collect(keys(buses))
abc= zeros(33)
for b in bus_set
abc[b] = error_variances[b]*(buses[b].p_PV)
end
println("????????????????????????????????????")
println("abc = ",abc)
println("????????????????????????????????????")

line_set = collect(keys(lines))
gen_set = collect(keys(generators))
gen_bus = []
for g in keys(generators)
    push!(gen_bus, generators[g].bus_idx)
end
root_bus = 0
for b in keys(buses)
    if buses[b].is_root
        root_bus = b
    end
end   

lines_to = Dict()
for l in keys(lines)
    lines_to[lines[l].to_node] = lines[l]
end

v_root = 1

var_sum_P = sum(abc)
var_P = abc

error_var_sum = sum(abc)

m = ChanceModel(solver=GurobiSolver())

@variable(m, v[bus_set] >= 0) #variable for voltage square
@variable(m, fp[bus_set]) #variable for active power flow
@variable(m, fq[bus_set]) # variable for reactive power flow
@variable(m, gp[bus_set]) #variable for active power generation
@variable(m, gq[bus_set]) #variable for reactive power generation
@variable(m, α[bus_set] >= 0) # variable to distribute load deviations (affine frequency control)
@variable(m, loadp[bus_set]) 
@variable(m, loadq[bus_set]) 

@indepnormal(m, e[z=1:length(bus_set)], mean=zeros(length(bus_set)), var=var_P)

# Expected loss minimization
@objective(m, Min, 10*(sum( lines_to[b].r*(sum(rPTDF[lines_to[b].index, k] *((var_P[k] + (α[k]^2)*var_sum_P)) for k in bus_set) + (fp[b])^2 ) for b in setdiff(bus_set, [root_bus])) +
sum( lines_to[b].r*(sum(rPTDF[lines_to[b].index, k] *((var_P[k] + (α[k]^2)*var_sum_P)) for k in bus_set) + (fq[b])^2 ) for b in setdiff(bus_set, [root_bus]))  ))

# Control variable
@constraint(m, sum(α) == 1)

#Voltage and flow constraints for root node
@constraint(m, v[root_bus] == v_root)
@constraint(m, fp[root_bus] == 0)
@constraint(m, fq[root_bus] == 0)
@constraint(m, con11[b in bus_set], loadp[b] == buses[b].d_P  )
@constraint(m, con22[b in bus_set], loadq[b] == buses[b].d_Q  ) 

for b in setdiff(bus_set, gen_bus)
# All buses without generation
    @constraint(m, α[b] == 0)
    @constraint(m, gp[b] == 0)
    @constraint(m, gq[b] == 0)
end

for b in bus_set
# All buses
    @constraint(m, loadp[b] - gp[b] + sum(fp[k] for k in buses[b].children) - buses[b].p_PV == fp[b])
    @constraint(m, loadq[b] - gq[b] + sum(fq[k] for k in buses[b].children) - buses[b].p_PV*buses[b].tanphi == fq[b])
end

for b in setdiff(bus_set, [root_bus])
# All buses without root
    b_ancestor = buses[b].ancestor[1]
    @constraint(m, v[b] == v[b_ancestor] - 2*(lines_to[b].r * fp[b] + lines_to[b].x * fq[b]))
end

# CHANCE CONSTRAINTS:
for b in gen_bus
# All buses with generation
    @constraint(m, gp[b] + α[b]*sum(e[bb] for bb in bus_set) <= buses[b].generator.g_P_max, with_probability = (1 - η_g))
    @constraint(m, gp[b] + α[b]*sum(e[bb] for bb in bus_set) >= 0, with_probability = (1 - η_g))

    @constraint(m, gq[b] + α[b]*sum(e[bb] * buses[bb].tanphi for bb in bus_set) <=  buses[b].generator.g_Q_max, with_probability = (1 - η_g))
    @constraint(m, gq[b] + α[b]*sum(e[bb] * buses[bb].tanphi for bb in bus_set) >= -buses[b].generator.g_Q_max, with_probability = (1 - η_g))
end

for b in setdiff(bus_set, [root_bus])
# All buses without root
    l = lines_to[b]
    line_arr = [lines[l] for l in line_set]
    @constraint(m, v[b] - 2 * (sum(rPTDF[l.index, b] * l.r * sum(rPTDF[l.index, j] * (e[j] - α[j]*sum(e[z] for z in bus_set)) for j in bus_set) for l in line_arr) +
                    sum(rPTDF[l.index, b] * l.x * sum(rPTDF[l.index, j] * (e[j]*buses[j].tanphi - α[j]*sum(e[z]*buses[j].tanphi for z in bus_set)) for j in bus_set) for l in line_arr))
                    <= buses[b].v_max, with_probability = (1 - η_v))
    @constraint(m, v[b] - 2 * (sum(rPTDF[l.index, b] * l.r * sum(rPTDF[l.index, j] * (e[j] - α[j]*sum(e[z] for z in bus_set)) for j in bus_set) for l in line_arr) +
                    sum(rPTDF[l.index, b] * l.x * sum(rPTDF[l.index, j] * (e[j]*buses[j].tanphi - α[j]*sum(e[z]*buses[j].tanphi for z in bus_set)) for j in bus_set) for l in line_arr))
                    >= buses[b].v_min, with_probability = (1 - η_v))
end

method = :Cuts
# method = :Reformulate
status = solve(m, method=method)
objective = getobjectivevalue(m)
if status == :Optimal
    info("CC solved to optimality with method $(method)")
    info("Objective: $(objective)")
else
    warn("CC not optimal, termination status $(status)")
end

# Prepare Results
bus_results = DataFrame(bus = Any[], d_P = Any[], d_Q = Any[], cosphi=Any[], tanphi = [],  g_P = Any[], g_Q = Any[], g_cost=Any[], v_squared = Any[], α=Any[], p_PV = Any[], q_PV= Any[])
line_results = DataFrame(line = Any[], from = Any[], to = Any[], f_P = Any[], f_Q = Any[], a = Any[])
for b in bus_set
    cost = 0
    if in(b, gen_bus)
        cost = buses[b].generator.cost
    end

    # Numerical smoothing
    α_res = getvalue(α[b])
    α_res < 1e-10 ? α_res = 0 : nothing

    res = [b, buses[b].d_P, buses[b].d_Q, buses[b].cosphi, buses[b].tanphi,  getvalue(gp[b]), getvalue(gq[b]), cost, getvalue(v[b]), α_res, buses[b].p_PV, buses[b].p_PV*buses[b].tanphi]
    push!(bus_results, res)

    # Calc square current on line
    a_res = 0
    fp_res = getvalue(fp[b])
    fq_res = getvalue(fq[b])
    if b != root_bus
        v_res = getvalue(v[buses[b].ancestor[1]])
        a_res = (fp_res^2 + fq_res^2)/v_res^2
        lres = [lines_to[b].index, buses[b].ancestor[1], b, fp_res, fq_res, a_res]
        push!(line_results, lres)
    end
end

sort!(bus_results, cols=:bus)
sort!(line_results, cols=:to)

Y_P = getdual(con11)
#println("yp", Y_P)
Y_Q = getdual(con22)
#println("yq", Y_Q)
println(line_results)
println(bus_results)
volt_prof = getvalue(v)

n_buses = nrow(bus_results)
n_lines = nrow(line_results)
########### CC_TEST###################
samples = 500
var_true = var_P
    OBSERVATIONS = zeros(length(var_true), samples )
    for (i,v) in enumerate(var_true)
        if v > 0
            dist = Normal(0,sqrt(v))
            OBSERVATIONS[i,:] = rand(dist, samples)
        end
    end
    n_buses = 33 
line_arr = [lines[l] for l in 1:n_lines]

dP = bus_results[:d_P]
dQ = bus_results[:d_Q]
gP = bus_results[:g_P]
gQ = bus_results[:g_Q]
α = bus_results[:α]
v_base = bus_results[:v_squared]
p_PV = bus_results[:p_PV]
q_PV = bus_results[:q_PV]

fP_base = line_results[:f_P]
fQ_base = line_results[:f_Q]

v_violation = 0
v_violation_temp = 0
g_violation = 0
f_violation = 0
v_violateion_in_sample = 0
v_violation_per_bus = zeros(33)

# Some helping Data
gen_bus = []
for g in keys(generators)
    push!(gen_bus, generators[g].bus_idx)
end
lines_to = Dict()
for l in keys(lines)
    lines_to[lines[l].to_node] = lines[l]
end

for s in 1:samples
    # if print_on 

    # end
    # Create Load Deviation from true distribution
    if s > length(OBSERVATIONS[1,:])
        warn("Requested sample size $(samples) larger than available sample set.")
        break
    end

    v_viol_flag = false

    load_dev_P = OBSERVATIONS[:,s]
    tanphi = [buses[b].tanphi for b in keys(buses)]
    load_dev_Q = load_dev_P .* tanphi
    load_dev_P_sum = sum(load_dev_P)
    load_dev_Q_sum = sum(load_dev_Q)

    dP_real = dP + load_dev_P
    dQ_real = dQ + load_dev_Q

    gP_real = gP + (α .* load_dev_P_sum)
    gQ_real = gQ + (α .* load_dev_Q_sum)

    # Calculate Flow and Voltage Based on LinDistFlow
    net_load_P = dP_real - gP_real - p_PV
    net_load_Q = dQ_real - gQ_real - q_PV

    fP_real = rPTDF*net_load_P
    fQ_real = rPTDF*net_load_Q
    R = [l.r for l in line_arr]
    X = [l.x for l in line_arr]

    v_real = ones(n_buses) - 2 * rPTDF' * (R.*fP_real + X.*fQ_real)

    # Check for constraint violations
    for i in gen_bus
        gP_real[i] > buses[i].generator.g_P_max ? g_violation += 1 : nothing
        gQ_real[i] > buses[i].generator.g_Q_max ? g_violation += 1 : nothing
        gP_real[i] < 0 ? g_violation += 1 : nothing
        gQ_real[i] < - buses[i].generator.g_Q_max ? g_violation += 1 : nothing
    end

    for i in 1:n_buses
        v_real[i] > buses[i].v_max ? v_violation += 1 : nothing
        v_real[i] < buses[i].v_min ? v_violation += 1 : nothing
    end

    for i in 1:n_buses
        if v_real[i] > buses[i].v_max
        v_violation_temp += 1
        elseif v_real[i] < buses[i].v_min
        v_violation_temp += 1 
        end
        v_violation_per_bus[i] += v_violation_temp
        v_violation_temp = 0
    end

end
percent_viol = zeros(33)
for i =1:33
percent_viol[i] = ((1-(v_violation_per_bus[i]/(2*samples)))*100)
end

println(" ++ Test Results with $samples samples ++ ")
@printf("%.2f%% of voltage constraints hold (%d violations).  \n", ((1-(v_violation/(2*samples*n_buses)))*100), v_violation)
@printf("%.2f%% of generation constraints hold (%d violations) . \n", ((1-(g_violation/(4*samples*length(gen_bus))))*100), g_violation)
println()
println()


return Y_P, Y_Q, getobjectivevalue(m), v_violation, g_violation, volt_prof, percent_viol
end #Solve_CCLinDist
