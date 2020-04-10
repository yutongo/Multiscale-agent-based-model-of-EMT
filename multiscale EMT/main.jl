include("Oval2Analysis_new.jl")
include("oval2_model_new.jl")
include("PosNormalRandom.jl")
using Oval2Analysis_new, Plots; plotlyjs();HDF5
function sode_emt(initial,t0,te)
    final=te-t0
    cat_multiplicative = 0.01ones(18) #multiplicative noise
    cat_multiplicative[15] = 0 #no multiplicative noise for abundance of TGF/miR200 complex
    prob = oval2_problem_new(u0=initial,noise=true,t_start=0.0,t_final = final[1],epsilon_m = cat_multiplicative)
    @time sol = solve(prob,SOSRI(),abstol=1e-4,reltol=1e-4,progress=true,saveat=1)
    return sol[end]
end
# use SODE track the development
#4 steady_states denoting 4 cell types
i1 = [0.104505, 1.01831, 0.00293619, 0.00267363, 0.00932251, 0.0354959, 0.260837, 0.000751749, 0.000153455, 3.13247e-5, 6.39432e-6, 1.30527e-6, 0.0607589, 1.24215, 0.0488046, 201.485, 135.094, 1.54814];
mes = [0.203673, 2.0135, 0.00243876, 0.00232337, 0.165027, 1.58873, 0.00269959, 3.91302e-5, 9.28936e-9, 2.20526e-12, 5.23519e-16, 1.24281e-19, 0.479603, 2.01778, 0.00226636, 104.181, 209.681, 0.641536];
i2 = [0.130952, 1.28448, 0.00269865, 0.00250373, 0.02894, 0.266073, 0.0435972, 0.000261974, 2.48615e-6, 2.35937e-8, 2.23906e-10, 2.12488e-12, 0.205279, 1.50978, 0.0327468, 157.779, 154.228, 1.45355];
epi =  [0.0226754, 0.0365571, 0.071047, 0.0190197, 0.00838248, 0.0286987, 0.286771, 0.000684817, 0.000157506, 3.62262e-5, 8.33194e-6, 1.91633e-6, 0.0595746, 1.23995, 0.0489362, 463.727, 100.62, 1.54878];

tstart = 0; tend = 5000; dt = 50; #time step size is 50
time = tstart:dt:tend
ntstep = length(time)
#initial 200 cells, max 8000 cells
n_cells = 8000; n_0 = 200;
#track times the cell can divide before terminally differentiated (TD) cell
Div = zeros(n_cells)
Div[1:n_0] = rand(big.(2:7),n_0) #each cell can divive 2-7 times brfore TD
Div_0 = Div
Cells = zeros(n_cells,18,ntstep) #store gene expression of each cell at each time
Track = zeros(n_cells,ntstep,20) #Track cell,time,tree
Age = zeros(n_cells) #keep track of age of terminally differentiated cells
tdiv = zeros(n_cells) #division time
time_div = zeros(n_cells) #time from beginning to the time of last division
time_TD = zeros(n_cells)  #keep track of TD cells expression with respect to time


#initialize  50 cells for each cell type
Cells[1:50,:,1] = epi'.*ones(50,18)
Cells[51:100,:,1] = i1'.*ones(50,18)
Cells[101:150,:,1] = i2'.*ones(50,18)
Cells[151:n_0,:,1] = mes'.*ones(50,18)
Track[1:n_0,1,1] = (1:n_0)'


tdiv[1:n_0] = PosNormalRandom(n_0,1,700,200) #first division time for initial cells

#find the time when cell starting to divide
fastest_div_time = minimum(time_div[1:n_0] + tdiv[1:n_0])
flag = 0
i_div = 1 #before the first round of division


# not end and no max step
while (fastest_div_time <= tend) && (i_div[1] <= ntstep)
    old_i_div = i_div # of step
        i_div = find(time->time >= fastest_div_time,time)[1] #index of first div time
    #existing cells before this division
    s = sum(Cells[:,:,old_i_div],2)
    ind_r = find(s->s>0,s) #existing cells before this division
    if (length(ind_r) == n_cells) || (isempty(ind_r))
        break
    end

    # find the division cell
    ind_d = find(x->x<=time[i_div],time_div[1:ind_r[end]] + tdiv[1:ind_r[end]])
    ind_d = intersect(ind_r',ind_d')

    # find TD cell before this division
    ind_div0 = find(x->x<=0,Div)
    ind_TD = intersect(1:ind_r[end],ind_div0')

    for j in 1:1:length(ind_r)
        if isempty(intersect(ind_r[j],ind_TD)) #not TD cell
            ii = old_i_div[1]+1
            for i = ii:i_div #last time to this division
                Cells[ind_r[j],:,i] = sode_emt(Cells[ind_r[j],:,i-1],time[i-1],time[i])
                Track[ind_r[j],i,:] = Track[ind_r[j],i-1,:]
            end
        end
    end

    # cells from this current generation that can divide or not
    for k in 1:1:length(ind_d)
        if ind_r[end]+1>n_cells # if reach max cell
            flag = 1
            break
        elseif Div[ind_d[k]]>0
            y_temp = Cells[ind_d[k],:,i_div]
            Cells[ind_r[end] + k,:,i_div] = y_temp
            a = find(x-> x==0, Track[ind_d[k],i_div,:])[1]
            Track[ind_d[k],i_div,a] = Track[ind_d[k],i_div,1]+0.1
            Track[ind_r[end]+k,i_div,:] = Track[ind_d[k],i_div,:]
            Track[ind_r[end]+k,i_div,a] = Track[ind_d[k],i_div,1]+0.2

            # noise at division
            noise_coef = PosNormalRandom(1,18,0,0.7)
            noise_coef[15] = 0
            for kk = 1:18
                while abs(noise_coef[kk]) >0.5   #  0.5
                    noise_coef[kk] = PosNormalRandom(1,1,0,0.7)[1]
                end
            end
            noise_div = Cells[ind_d[k],:,i_div].*noise_coef'
            for kk = 1:18
                coin = rand(big.(0:1),1)[1]
                if coin == 0
                    noise_div[kk] = (-1)*noise_div[kk]
                end
            end
            Cells[ind_d[k],:,i_div] = Cells[ind_d[k],:,i_div] + noise_div
            Cells[ind_r[end]+k,:,i_div] = Cells[ind_r[end]+k,:,i_div] - noise_div
            # update how many times left where cell can divide brfore TD
            Div[ind_d[k]]= Div[ind_d[k]] - 1
            Div[ind_r[end]+k] = Div[ind_d[k]]
            # update time from begin to most recent division
            time_div[ind_d[k]] = time[i_div]#time_div[k] + tdiv[k]
            time_div[ind_r[end]+k] = time[i_div]#time_div[k] + tdiv[k]

                if Div[ind_d[k]] >0 #assign division time to new daughter cells
                    tdiv[ind_d[k]] = PosNormalRandom(1,1,700,200)[1]
                    tdiv[ind_r[end] + k] = PosNormalRandom(1,1,700,200)[1]
                elseif Div[ind_d[k]] <= 0 #daughet cells are TD, no division happens
                    time_TD[ind_d[k]] = time[i_div]
                    time_TD[ind_r[end] + k] = time[i_div]
                    tdiv[ind_d[k]] = 0
                    tdiv[ind_r[end] + k] = 0
                end


        elseif Div[ind_d[k]]<=0
            continue
        end
    end

    t_div_next = time_div + tdiv
    ID = find(x->x>0,tdiv)
    dummy_div = find(x-> x >0, Div)
    ID = intersect(ID,dummy_div)
    if isempty(ID) #if no division any more, just growth all TD
        break
    else
    fastest_div_time = minimum(t_div_next[ID])
    end
end

# all TD cells
ind_div0 = find(x->x<=0,Div)
s = sum(Cells[:,:,i_div],2)
ind_r = find(s->s>0,s)
ind_TD = intersect(ind_r,ind_div0)
death_rate = PosNormalRandom(1,n_cells,1000,100)
for k in 1:1:length(ind_TD)
    Age[ind_TD[k]] = time[i_div]-time_TD[ind_TD[k]]
    if time[end]>=time_TD[ind_TD[k]]+death_rate[ind_TD[k]] #live before the end
        death_time = find(x->x>=time_TD[ind_TD[k]]+death_rate[ind_TD[k]],time)[1]
            for kk = Int(time_TD[ind_TD[k]]/dt) + 2:death_time
                Cells[ind_TD[k],:,kk] = sode_emt(Cells[ind_TD[k],:,kk-1],time[kk-1],time[kk])
                Track[ind_TD[k],kk,:] = Track[ind_TD[k],kk-1,:]
            end
            for kk = death_time+1:ntstep
                Age[ind_TD[k]] = 0
                Cells[ind_TD[k],:,kk] = zeros(1,18)
                Track[ind_TD[k],kk,:] = zeros(20)
            end
    else #live until the end
        for kk = Int(time_TD[ind_TD[k]]/dt)+2:ntstep
            Cells[ind_TD[k],:,kk] = sode_emt(Cells[ind_TD[k],:,kk-1],time[kk-1],time[kk])
            Track[ind_TD[k],kk,:] = Track[ind_TD[k],kk-1,:]
        end
    end
end

using MAT
using HDF5


using JLD
@save "data.jld" # save all data
#save all variables to mat
file = matopen("emt4_0522.mat", "w")
write(file, "i1", i1)
write(file, "mes", mes)
write(file, "i2", i2)
write(file, "epi", epi)
write(file, "tend", tend)
write(file, "ntstep", ntstep)
write(file, "n_cells", n_cells)
write(file, "n_0", n_0)
write(file, "Div", Div)
write(file, "Cells", Cells)
write(file, "Track", Track)
write(file, "Age", Age)
write(file, "death_rate", death_rate)
write(file, "tdiv", tdiv)
write(file, "time_div", time_div)
write(file, "time_TD", time_TD)
write(file, "stem_marker", stem_marker)
write(file, "fastest_div_time", fastest_div_time)
write(file, "i_div", i_div)
write(file, "Div_0", Div_0)
write(file, "ind_TD", ind_TD)
write(file, "ind_r", ind_r)
close(file)
