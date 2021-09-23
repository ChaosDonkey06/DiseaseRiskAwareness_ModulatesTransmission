import networkx as nx
import numpy as np
import pandas as pd
import os 

main_path = os.path.split(os.getcwd())[0] + '/Epidemiology_behavior_dynamics'
config_path = main_path + '/config.csv'
config_data = pd.read_csv(config_path, sep=',', header=None, index_col=0) 
networks_path = config_data.loc['networks_dir'][1]

import argparse

parser = argparse.ArgumentParser(description='Creation of networks using Networkx.')

parser.add_argument('--num_nodes', type=int, default=1000,
                    help='Define if create for all networks or only one in specific')
parser.add_argument('--all', type=bool, default=False,
                    help='Define if create for all networks or only one in specific')
parser.add_argument('--specific_network', type=str, default='scale_free',
                    help='Define the specific network you want to create: scale_free, small_world, grid')
args = parser.parse_args()

if not os.path.isdir( os.path.join(main_path, networks_path, str(args.num_nodes)) ):
        os.makedirs(os.path.join(main_path, networks_path, str(args.num_nodes)))


if (args.all == False) and (args.specific_network == 'scale_free'):
    # Scale-free graph
    graph_sf = nx.scale_free_graph( args.num_nodes ).to_undirected()
    nx.write_gpickle(graph_sf, os.path.join(main_path, networks_path, str(args.num_nodes), 'scale_free_'+str(args.num_nodes) ) )
    print('Scale-Free graph created!\n')

elif (args.all == False) and (args.specific_network == 'small_world'):
    # Watts-Strogatz graph
    graph_ws = nx.watts_strogatz_graph(args.num_nodes, p=0.5, k=5).to_undirected()
    nx.write_gpickle(graph_ws, os.path.join(main_path, networks_path, str(args.num_nodes), 'small_world_'+str(args.num_nodes)) )
    print('Watts-Strogatz graph created!\n')

elif (args.all == False) and (args.specific_network == 'grid'):
   # Grid graph
    graph_g = nx.grid_graph(dim=[int(np.sqrt(args.num_nodes)),int(np.sqrt(args.num_nodes))]).to_undirected()
    nx.write_gpickle(graph_g, os.path.join(main_path, networks_path, str(args.num_nodes), 'grid_'+str(args.num_nodes)) )
    print('Grid graph created!\n')

elif args.all == True:
    # Scale-free graph
    graph_sf = nx.scale_free_graph( args.num_nodes ).to_undirected()
    nx.write_gpickle(graph_sf, os.path.join(main_path, networks_path, str(args.num_nodes), 'scale_free_'+str(args.num_nodes) ) )
    print('Scale-Free graph created!\n')

    # Watts-Strogatz graph
    graph_ws = nx.watts_strogatz_graph(args.num_nodes, p=0.5, k=5).to_undirected()
    nx.write_gpickle(graph_ws, os.path.join(main_path, networks_path, str(args.num_nodes), 'small_world_'+str(args.num_nodes)) )
    print('Watts-Strogatz graph created!\n')

    # Grid graph
    graph_g = nx.grid_graph(dim=[int(np.sqrt(args.num_nodes)),int(np.sqrt(args.num_nodes))]).to_undirected()
    nx.write_gpickle(graph_g, os.path.join(main_path, networks_path, str(args.num_nodes), 'grid_'+str(args.num_nodes)) )
    print('Grid graph created!\n')

else:
    print('Invalid')