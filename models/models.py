import pandas as pd
import numpy as np
import random
import os

##########################################
# Game dynamics              
##########################################
def i_payoff(A, NN, node_i, actual_state):
    '''
    Returns payoff of node i given its Neares Neighboors (NN).
    @ params: A - payoff matrix (numpy array)
              NN - nearest neighboor obtained from adjacency matrix (numpy array)
              node_i - actual node to calculate payoff
              actual_state - current node states
    @ return: u - node payoff (float)
    '''
    u = 0; # reset payoff
    for j, j_nn in enumerate(NN):
        u += A[ int(actual_state[node_i]), int(actual_state[j_nn]) ]    # A[actual_node_state, actual_NN_state]
    u /= NN.size

    return u

def Pij(p_i, p_j):
    ''' 
    Fermi-Dirac distribution defines how player i change strategy given j player given payoff's
    '''
    K = 0.5     # Irrationality of players
    Pij = 1/(1 + np.exp(-(p_j - p_i)/K))

    return Pij

##########################################################################
# Models                                    
######################################################################### 
def sis_replicator(graph, max_time, params, local=False):

    N = len( graph.nodes() )
    G = graph

    id_node = {node:idx for idx,node in enumerate(G.nodes())}   # In case nodes are not ints

    ####################
    # Initial params 
    ####################
    # Payoff matrix (of each two-player game)
    S = -0.5
    T = 1.5

    sigma = params.get('awareness', 0)

    probInfect = params['probInfect']
    maxRecTime = params.get('time2Recover', 7)
    initInfected = params.get('initInfected', None)
    initDefectors = params.get('initDefectors', None)

    p_inf = probInfect

    #####################
    # Initial states   
    #####################
    epid_states = np.zeros([N, max_time])     # 1: infected, 0: susceptible
    game_states = np.zeros([N, max_time])     # 1: defectors, 0: coperator
    
    recTime = np.zeros(N,)
    payoffs = np.zeros([N, max_time])         # Payoff's of each iteration

    # Initial infected nodes
    if initInfected is None:    # If initial conditions are given
        initInfected = random.sample(range(int(N)), int(N/(N-1)))
        epid_states[initInfected,0] = 1

    else:
        epid_states[initInfected,0] = 1
    # In recovery times
    recTime[initInfected] = maxRecTime
    
    # Initial defectors nodes
    if initDefectors is None:   # If initial conditions are given
        initDefectors = random.sample(range(int(N)), int(N/2))
        game_states[initDefectors,0] = 1
    else:
        game_states[initDefectors,0] = 1

    ##########################
    # Start simulation 
    ##########################
    for t in range(1, max_time):

        curInf_nodes = np.nonzero(epid_states[:,t-1])[0]   # Get current infected nodes
        epid_states[curInf_nodes, t] = 1                         # Update infected nodes

        if local is False:
            FracInf = np.sum(epid_states[:,t])/N         # Get fraction of infected nodes
            sigma_infection = sigma*FracInf              # Calculate sigma infection

        recTime = np.maximum(recTime-1,0)                # Update infection state depending on recovery time
        epid_states[ recTime>0, t] = 1
        epid_states[ recTime==0, t] = 0


        ##### Run Game 
        for idx, node_i in enumerate(G.nodes()):
            node_i_nn_ = [neig for neig in G.neighbors(node_i)]          # get neighbors
            node_i_nn  = np.array([id_node[nn] for nn in node_i_nn_ ])
            
            if local is True:
                FracInf = np.sum(epid_states[node_i_nn, t-1])/node_i_nn.size           # Calculate fraction of infected neighbors
                sigma_infection = sigma*FracInf                                        # Infection awearness

            A = np.array([[1,S], [T-sigma_infection, 0-sigma_infection]])              # Game payoff matrix
            
            payoffs[node_i, t-1] = i_payoff(A, node_i_nn, node_i, game_states[:,t-1])  # Calculate payoff

        for idx, node_i in enumerate(G.nodes()):
            node_i_nn_ = [neig for neig in G.neighbors(node_i)]        # get neighbors
            node_i_nn  = np.array([id_node[nn] for nn in node_i_nn_ ])

            node2change = random.randint(1, node_i_nn.size)
            P_change = Pij(payoffs[node_i,t-1], payoffs[node2change,t-1])   # Calculate probability Pij

            th = 0.5
            game_states[node_i, t] = (P_change>th)*game_states[node2change,t-1]+(P_change<=th)*game_states[node_i,t-1]


        ##### Run Infection
        for idx, node_i in enumerate(G.nodes()):
            if epid_states[node_i, t-1] == 1:      # If node is already infected
                continue                           # Skip                                    
            
            cur_node_nn_ = [neig for neig in G.neighbors(node_i)]                     # Get neighbors
            cur_node_nn = np.array([id_node[nn] for nn in cur_node_nn_ ])
            cur_node_inf_nn = np.nonzero(epid_states[cur_node_nn,t-1])[0]             # Get infected neighbors
            cur_game_nn = game_states[cur_node_inf_nn,t-1]                            # Get game states infected neighbors
                
            n_infected_neighbors = np.sum(epid_states[cur_node_nn,t-1])
            # if local is False:
            p_inf_node_i = p_inf*np.exp(-np.sum(game_states[curInf_nodes,t]==0)/N)        # Infection probability due to cooperation state of all network infected nodes
                            
            if local is True:       # If only considering local neighbors
                if n_infected_neighbors > 0:    # If there is infection in my neighbors
                    p_inf_node_i = p_inf*np.exp(-np.sum(cur_game_nn==0)/cur_game_nn.size)  # Infection probability due to cooperation state of infected neighbors
                else:
                    p_inf_node_i = 0.00000001

            if n_infected_neighbors == 0:   # If no neighbor is infected, infection probability ~ 0
                p_inf_node_i = 0.00000001

            # Update new infection
            if random.random() < p_inf_node_i:
                epid_states[node_i, t] = 1
                # Update recovery time                                     
                recTime[node_i] = maxRecTime

    return game_states, epid_states, payoffs

def run_model(model, graph, params, n_iters=10, max_time=300, num_checkpoints=1,  local=False, path_to_save_checkpoints=None):

    df_list = []
    chk_point = 0
    for n_iter in range(n_iters):

        game_states, infection_states, payoffs = model(graph, max_time, params, local=local)

        num_def = np.sum(game_states, axis=0)
        num_inf = np.sum(infection_states, axis=0)        
        n_total = np.array( [len(graph.nodes)]*max_time )
        num_cop = n_total-num_def
        num_sus = n_total-num_inf

        susceptibles = list(num_sus)
        infected     = list(num_inf)
        defectors    = list(num_def)
        cooperators  = list(num_cop)

        df_response         = pd.DataFrame(columns=['sim_id','time','S','I','C','D'])
        df_response['time'] = list(range(max_time))
        df_response['S'] = susceptibles
        df_response['I'] = infected
        df_response['C'] = cooperators
        df_response['D'] = defectors
        df_response['sim_id'] = [n_iter]*len(defectors)

        df_list.append(df_response)
        df_response = df_response.set_index('time')

        if chk_point < num_checkpoints:

            # save disease checkpoints
            np.savetxt(os.path.join(path_to_save_checkpoints,'epid_iter_{}_of_{}_beta_{}_sigma_{}.csv'.format(chk_point,n_iters,params['probInfect'],params['awareness'])),
                        infection_states, delimiter=',')
            # save game checkpoints
            np.savetxt(os.path.join(path_to_save_checkpoints,'game_iter_{}_of_{}_beta_{}_sigma_{}.csv'.format(chk_point,n_iters,params['probInfect'],params['awareness'])),
                        game_states, delimiter=',')

            chk_point+=1

    df_return = pd.concat(df_list)

    return df_return
