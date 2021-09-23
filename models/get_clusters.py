import networkx as nx
import numpy as np
import pandas as pd
import community as community_louvian
import heapq
from operator import itemgetter
from scipy.stats import mode
from tqdm import tqdm
import os


def get_partition(G, n_biggest=3):

    partition = community_louvian.best_partition(G, randon_state=0)     # Get partition/clusters

    n_cluster = max(partition.values())+1               # Get number of clusters

    cluster_nodes = {}                                  # Get nodes per cluster {cluster:[nodes]}
    
    for key, value in partition.items(): 
        if value not in cluster_nodes: 
            cluster_nodes[value] = [key] 
        else: 
            cluster_nodes[value].append(key)

    n_nodes_in_cluster = { cluster:len(nodes) for cluster, nodes in cluster_nodes.items() }       # Get number of nodes per cluster

    top_clusters = dict(heapq.nlargest(n_biggest, n_nodes_in_cluster.items(), key=itemgetter(1))) # Get n biggest cluster

    # small_clusters = dict(heapq.nsmallest(n_biggest-1, n_nodes_in_cluster.items(), key=itemgetter(1))) # Get n smallest cluster
    
    top_cluster_nodes = {cluster:cluster_nodes[cluster] for cluster,_ in top_clusters.items()}    # Get biggest clusters items

    # small_cluster_nodes = {cluster:cluster_nodes[cluster] for cluster,_ in small_clusters.items()}    # Get smallest clusters items
    
    return partition, n_cluster, cluster_nodes, top_clusters, top_cluster_nodes




def top_cluster_dynamics(top_cluster_nodes, disease_checkpoint_path, game_checkpoint_path, iters, beta, sigma):

    for it in range(iters+1):

        disease_checkpoint = np.loadtxt(os.path.join(disease_checkpoint_path,
                             'epid_iter_{}_of_20_beta_{}_sigma_{}.csv'.format(it,beta,sigma), delimiter=','))
        game_checkpoint    = np.loadtxt(os.path.join(game_checkpoint_path,
                             'game_iter_{}_of_20_beta_{}_sigma_{}.csv'.format(it,beta,sigma), delimiter=','))
    
        maxtime = disease_checkpoint.shape[0]          # Get time stamp

        df_list = []

        for cluster, nodes in top_cluster_nodes.items():

            nodes_ = np.array(nodes)

            n_def  = np.sum(game_checkpoint[nodes_], axis=0)
            n_inf  = np.sum(disease_checkpoint[nodes_], axis=0)
            n_tot  = np.array([len(nodes)] * maxtime)
            n_cop  = n_tot - n_def
            n_sus  = n_tot - n_inf

            n_def  = n_def/len(nodes)                  # Normalize
            n_inf  = n_inf/len(nodes)
            n_cop  = n_cop/len(nodes)
            n_sus  = n_sus/len(nodes)

            susceptibles = list(n_sus)
            infected     = list(n_inf)
            defectors    = list(n_def)
            cooperators  = list(n_cop)

            df_cluster_dynamic = pd.DataFrame(columns=['sim_it','time','cluster','n_nodes','S','I','C','D','beta','sigma','R0']) 

            df_cluster_dynamic['sim_it']  = [it] * maxtime
            df_cluster_dynamic['cluster'] = cluster
            df_cluster_dynamic['n_nodes'] = len(nodes)
            df_cluster_dynamic['S']       = susceptibles
            df_cluster_dynamic['I']       = infected
            df_cluster_dynamic['C']       = cooperators
            df_cluster_dynamic['D']       = defectors
            df_cluster_dynamic['beta']    = beta
            df_cluster_dynamic['sigma']   = sigma
            df_cluster_dynamic['R0']      = beta/(1/7)

            df_list.append(df_cluster_dynamic)

    df_return = pd.concat(df_list)

    return df_return


def cluster_dynamics(top_cluster_nodes, type_sim, results_path, iters, beta, sigma):

    maxtime = 150

    ic_list = []
    
    for ic in tqdm(range(0,10), total=10):

        list_game = []
        list_epid = []


        for it in range(iters+1):

            disease_checkpoint = np.loadtxt(os.path.join(results_path, type_sim, 'scale_free', 'checkpoints', 'ic_0{}'.format(ic),
                                'epid_iter_{}_of_20_beta_{}_sigma_{}.csv'.format(it,beta,sigma)), delimiter=',')
            game_checkpoint    = np.loadtxt(os.path.join(results_path, type_sim, 'scale_free', 'checkpoints', 'ic_0{}'.format(ic),
                                'game_iter_{}_of_20_beta_{}_sigma_{}.csv'.format(it,beta,sigma)), delimiter=',')
    
            list_game.append(game_checkpoint)
            list_epid.append(disease_checkpoint)
        
        arr_game = np.dstack(list_game)
        arr_epid = np.dstack(list_epid)

        checks_game_mode = mode(arr_game, axis=2)[0].reshape((5000,150))
        checks_epid_mode = mode(arr_epid, axis=2)[0].reshape((5000,150))

        df_list = []

        for cluster, nodes in top_cluster_nodes.items():

            nodes_ = np.array(nodes)

            n_def  = np.sum(checks_game_mode[nodes_], axis=0)
            n_inf  = np.sum(checks_epid_mode[nodes_], axis=0)
            n_tot  = np.array([len(nodes)] * maxtime)
            n_cop  = n_tot - n_def
            n_sus  = n_tot - n_inf

            n_def  = n_def/len(nodes)                  # Normalize
            n_inf  = n_inf/len(nodes)
            n_cop  = n_cop/len(nodes)
            n_sus  = n_sus/len(nodes)

            susceptibles = list(n_sus)
            infected     = list(n_inf)
            defectors    = list(n_def)
            cooperators  = list(n_cop)

            df_cluster_dynamic = pd.DataFrame(columns=['sim_ic','time','cluster','n_nodes','S','I','C','D','beta','sigma','R0']) 

            df_cluster_dynamic['sim_ic']  = [ic] * maxtime
            df_cluster_dynamic['time']    = list(range(maxtime))
            df_cluster_dynamic['cluster'] = cluster
            df_cluster_dynamic['n_nodes'] = len(nodes)
            df_cluster_dynamic['S']       = susceptibles
            df_cluster_dynamic['I']       = infected
            df_cluster_dynamic['C']       = cooperators
            df_cluster_dynamic['D']       = defectors
            df_cluster_dynamic['beta']    = beta
            df_cluster_dynamic['sigma']   = sigma
            df_cluster_dynamic['R0']      = beta/(1/7)

            df_list.append(df_cluster_dynamic)

        ic_list.append(pd.concat(df_list))


    df_return = pd.concat(ic_list)

    return df_return





