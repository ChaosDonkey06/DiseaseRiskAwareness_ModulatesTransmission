import config as cf
import numpy as np
import pandas as pd
import os
import random 
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from models import sis_replicator
from models.get_clusters import get_partition, cluster_dynamics
assert cf

import argparse 

parser = argparse.ArgumentParser(description='Network simulations.')

parser.add_argument('--network_type', type=str, default='scale_free',
                    help='Network type for storing...')
parser.add_argument('--network_name', type=str, default='scale_free_1000',
                    help='Network type for storing...')
parser.add_argument('--num_nodes', type=int, default=1000,
                    help='Number of nodes for specific network...')
parser.add_argument('--iterations', type=int, default=10,
                    help='Selected iterations for plotting...')
parser.add_argument('--n_clusters', type=int, default=3,
                    help='Selected iterations for plotting...')
parser.add_argument('--beta_select', type=float, default=0.6,
                    help='Selected beta for plotting...')
parser.add_argument('--beta_key', type=str, default='060',
                    help='Selected beta for plotting...')
parser.add_argument('--sigma_select', type=float, default=1.0,
                    help='Selected beta for plotting...')
parser.add_argument('--sigma_key', type=str, default='100',
                    help='Selected beta for plotting...')

args = parser.parse_args()


main_path = os.path.split(os.getcwd())[0] + '/Epidemiology_behavior_dynamics'
config_path = main_path + '/config.csv'
config_data = pd.read_csv(config_path, sep=',', header=None, index_col=0)

networks_path = config_data.loc['networks_dir'][1]
results_path  = config_data.loc['results_dir'][1]
figures_path  = config_data.loc['figures_dir'][1]
num_nodes     = args.num_nodes

# Path to checkpoints
game_checkpoint_path_global = os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', 'global', args.network_type,'checkpoints')
disease_checkpoint_path_global = os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', 'global', args.network_type,'checkpoints')

game_checkpoint_path_local = os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', 'local', args.network_type,'checkpoints')
disease_checkpoint_path_local = os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', 'local', args.network_type,'checkpoints')

# Load network
G = nx.read_gpickle( os.path.join(main_path, networks_path, str(num_nodes), args.network_name) )


# Get communities
partition, n_cluster, cluster_nodes, top_clusters, top_cluster_nodes = get_partition(G, n_biggest=args.n_clusters) #n_biggest=3)


# Get cluster dynamics
df_cluster_dyncs_global = cluster_dynamics(top_cluster_nodes, disease_checkpoint_path_global, game_checkpoint_path_global, beta=args.beta_select, sigma=args.sigma_select)
df_cluster_dyncs_local = cluster_dynamics(top_cluster_nodes, disease_checkpoint_path_local, game_checkpoint_path_local, beta=args.beta_select, sigma=args.sigma_select)


## Plot results

colors_plt = [ 'tab:green', 'tab:red', 'tab:blue'] #, 'tab:purple', 'tab:cyan', 'tab:orange'  ]

n_top_clusters = list(top_cluster_nodes.keys())
n_top_clusters = n_top_clusters[:3]                 # Get biggest 3

fig, ax = plt.subplots(1,2,figsize=(20, 7))
ax = ax.flatten()

for idx, clust in tqdm(enumerate(n_top_clusters), total=len(n_top_clusters)):

    # Get cluster data
    clust_mask_global = df_cluster_dyncs_global['cluster'] == clust
    clust_mask_local = df_cluster_dyncs_local['cluster'] == clust

    # For global
    df_clust_i_glob = pd.DataFrame(df_cluster_dyncs_global[clust_mask_global])
    df_clust_i_glob['type'] = ['Global'] * len(df_clust_i_glob)
    df_clust_i_loc  = pd.DataFrame(df_cluster_dyncs_local[clust_mask_local])
    df_clust_i_loc['type'] = ['Local'] * len(df_clust_i_loc)

    df_res = [df_clust_i_loc, df_clust_i_glob]
    df_res_c = pd.concat(df_res)


    # Plot global
    sns.lineplot( ax = ax[0],
                  data = df_res_c,
                  x = 'time', y = 'I',
                  label = r'Cluster {}: {} individuals'.format(clust, top_clusters[clust]),
                  style='type',
                  color=colors_plt[idx],
                  alpha=0.5)
    ax[0].get_legend().remove()
    #ax[0].lines[0].set_linestyle("--")
    ax[0].set_title(r'Clustered disease dynamics $R_0=4.2$ $\sigma={}$ '.format(str(args.sigma_select)),fontsize=22)
    ax[0].set_xlabel(r'Days',fontsize=21)
    ax[0].xaxis.set_tick_params(labelsize=20)
    ax[0].yaxis.set_tick_params(labelsize=20)
    ax[0].set_xlim([-0.1,151])
    ax[0].set_ylabel(r'Cluster Inf. Fraction $I$',fontsize=21)
    ax[0].set_ylim([-0.1,1.1])

    sns.lineplot( ax = ax[1],
                  data = df_res_c,
                  x = 'time', y = 'C',
                  style='type',
                  color=colors_plt[idx],
                  alpha = 0.5)
    ax[1].get_legend().remove()
    #ax[0].lines[0].set_linestyle("--")
    ax[1].set_title(r'Clustered behavioral dynamics $R_0=4.2$ $\sigma={}$ '.format(str(args.sigma_select)),fontsize=22)
    ax[1].set_xlabel(r'Days',fontsize=21)
    ax[1].xaxis.set_tick_params(labelsize=20)
    ax[1].yaxis.set_tick_params(labelsize=20)
    ax[1].set_xlim([-0.1,151])
    ax[1].set_ylabel(r'Cluster Coop. Fraction $c$',fontsize=21)
    ax[1].set_ylim([-0.1,1.1])
    plt.tight_layout()

if not os.path.isdir( os.path.join(figures_path, 'dynamics', str(num_nodes)) ):
        os.makedirs( os.path.join(figures_path, 'dynamics', str(num_nodes)) )
im_path_save = os.path.join(figures_path, 'dynamics', str(num_nodes))
plt.savefig(os.path.join(im_path_save, 'new_cluster_dynamics_sigma_{}_beta_{}.png'.format(str(args.sigma_select),str(args.beta_select))), 
                             dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)


### save legends
fig, ax = plt.subplots(2,1,figsize=(9, 8))
coms = [1,2,3]
for idx, clust in tqdm(enumerate(n_top_clusters), total=len(n_top_clusters)):

    # Get cluster data
    clust_mask_global = df_cluster_dyncs_global['cluster'] == clust

    # For global
    df_clust_i_glob = pd.DataFrame(df_cluster_dyncs_global[clust_mask_global])

    sns.lineplot( ax = ax[0],
                  data = df_clust_i_glob, 
                  x = 'time', y = 'I',
                  label = r'Community {}: {} ind.'.format(coms[idx], top_clusters[clust]),
                  color = colors_plt[idx])
    ax[0].get_legend().remove()
    ax[0].set_title(r'Disease dynamics')
    #ax[0].set_xlabel('')
    ax[0].set_xlabel(r'Days',fontsize=21)
    #ax[0].set_xticks(fontsize=20)
    ax[0].set_xticklabels('')
    ax[0].set_xlim([-0.1,151])
    ax[0].set_ylabel(r'Inf. Fraction $\bar{I}$',fontsize=21)
    #ax[0].set_yticks(fontsize=20)
    ax[0].set_ylim([-0.1,1.1])

    plt.figlegend(bbox_to_anchor=(0.9,0.4), fontsize=22)

if not os.path.isdir( os.path.join(figures_path, 'dynamics', str(num_nodes), 'labels') ):
        os.makedirs( os.path.join(figures_path, 'dynamics', str(num_nodes), 'labels') )
save_path = os.path.join(figures_path, 'dynamics', str(num_nodes), 'labels')
plt.savefig(os.path.join(save_path, 'new_colorlabel_cluster_dynamics.png'), 
                             dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)
#plt.show()

#### save type
fig, ax = plt.subplots(2,1,figsize=(9, 8))


df_clust_t_loc  = df_cluster_dyncs_local.copy()
df_clust_t_loc['type'] = ['Global'] * len(df_cluster_dyncs_local)
df_clust_t_glob = df_cluster_dyncs_global.copy()
df_clust_t_glob['type'] = ['Local'] * len(df_cluster_dyncs_global)

df_res = [df_clust_t_loc, df_clust_t_glob]
df_res_c = pd.concat(df_res)

sns.lineplot( ax = ax[0],
                data = df_res_c, 
                x = 'time', y = 'I',
                style='type', alpha=0.5)
ax[0].get_legend().remove()
ax[0].set_title(r'Disease dynamics')
#ax[0].set_xlabel('')
ax[0].set_xlabel(r'Days',fontsize=21)
#ax[0].set_xticks(fontsize=20)
ax[0].set_xticklabels('')
ax[0].set_xlim([-0.1,151])
ax[0].set_ylabel(r'Inf. Fraction $\bar{I}$',fontsize=21)
#ax[0].set_yticks(fontsize=20)
ax[0].set_ylim([-0.1,1.1])

plt.figlegend(bbox_to_anchor=(0.7,0.4), fontsize=22)

plt.savefig(os.path.join(figures_path, 'dynamics', 'style_cluster_dynamics.png'), 
                             dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)



## Plot graph
from matplotlib import cm
def plot_G(G, communities, not_comms, plot_title, nodecmap,figures_path=figures_path):

    n_comms = max(communities.values())
    pos     = nx.kamada_kawai_layout(G)

    plt.figure(figsize=(12,12))
    #cmap = cm.get_cmap('hsv')#, n_comms+)
    #plt.title(str(plot_title) + ' {} biggest clusters'.format(n_comms), size=15)

    nx.draw(G, pos,
            nodelist    = list(not_comms.keys()),
            node_size   = 12,
            #cmap        = cmap,
            node_color  = 'black',
            edge_color  = 'gray',
            width       = .2,
            with_labels = False
            )

    nx.draw(G, pos,
            nodelist    = list(communities.keys()),
            node_size   = 12,
            #cmap        = cmap,
            node_color  = nodecmap,
            edge_color  = 'gray',
            width       = .2,
            with_labels = False
            )
    plt.savefig(os.path.join(figures_path,'networks',str(num_nodes),'graph_cluster_dynamics_new_6clus.png'), 
                              dpi=400, transparent = True, bbox_inches = 'tight', pad_inches = 0.1)


n_top_clusters = list(top_cluster_nodes.keys())
nodes = []
cluster = []
n_nodes = []
n_cluster = []
for idx, clus in enumerate(n_top_clusters):
    for key, value in partition.items():
        if clus == value:
            nodes.append(key)
            cluster.append(value)
        else:
            n_nodes.append(key)
            n_cluster.append(value)

node_cluster = dict(zip(nodes,cluster))
other_node_cluster = dict(zip(n_nodes,n_cluster))
colors_plt = [ 'tab:green', 'tab:red', 'tab:blue', 'tab:purple', 'tab:cyan', 'tab:orange' ]
node_cmap = []
for n, c in node_cluster.items():
    if c == n_top_clusters[0]:
        node_cmap.append(colors_plt[0])
    if c == n_top_clusters[1]:
        node_cmap.append(colors_plt[1])
    if c == n_top_clusters[2]:
        node_cmap.append(colors_plt[2])
    if c == n_top_clusters[3]:
        node_cmap.append(colors_plt[3])
    if c == n_top_clusters[4]:
        node_cmap.append(colors_plt[4])
    if c == n_top_clusters[5]:
        node_cmap.append(colors_plt[5])

plot_G(G, node_cluster, other_node_cluster, 'Clustered graph', node_cmap)