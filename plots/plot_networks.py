import networkx as nx
import pandas as pd
import os 
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description='Networks visualizarion')

parser.add_argument('--network_type', type=str,
                    help='Type of network [scale_free,small_world,grid]')
parser.add_argument('--network_name', type=str,
                    help='Name of the network created in the /networks folder')
parser.add_argument('--num_nodes', type=int, default=1000,
                    help='Number of nodes for specific network...')

args = parser.parse_args() 

main_path = os.path.split(os.getcwd())[0] + '/Epidemiology_behavior_dynamics'
config_path = main_path + '/config.csv'
config_data = pd.read_csv(config_path, sep=',', header=None, index_col=0)

figures_path  = config_data.loc['figures_dir'][1]
networks_path = config_data.loc['networks_dir'][1]
num_nodes     = args.num_nodes

read_network = nx.read_gpickle( os.path.join(main_path, networks_path, str(num_nodes), str(args.network_name)) )

if args.network_type == 'grid':
        plt.figure(figsize=(12,12))
        pos = nx.spectral_layout(read_network)
        nx.draw(G=read_network, pos=pos, 
                node_size=12,
                node_color= 'black', 
                edge_color='gray',
                width=.2,
                edge_cmap=plt.cm.Blues, with_labels=False)
        if not os.path.isdir( os.path.join(figures_path, 'networks', str(num_nodes)) ):
                os.makedirs ( os.path.join(figures_path, 'networks', str(num_nodes)) )
        figures_path_save = os.path.join(figures_path, 'networks', str(num_nodes))
        plt.savefig(os.path.join(figures_path_save, 'grid_networks_viz.png'), dpi=400, transparent = True, bbox_inches = 'tight', pad_inches = 1)
        

elif (args.network_type == 'scale_free') or (args.network_type == 'small_world'):
        plt.figure(figsize=(12,12))
        pos = nx.kamada_kawai_layout(read_network)
        nx.draw(G=read_network, pos=pos, 
                node_size=12,
                node_color= 'black',
                edge_color='gray',
                width=.2,
                edge_cmap=plt.cm.Blues, with_labels=False)
        if not os.path.isdir( os.path.join(figures_path, 'networks', str(num_nodes)) ):
                os.makedirs ( os.path.join(figures_path, 'networks', str(num_nodes)) )
        figures_path_save = os.path.join(figures_path, 'networks', str(num_nodes))
        plt.savefig( os.path.join(figures_path_save, '{}.png'.format(str(args.network_name)) ), dpi=400, transparent = True, bbox_inches = 'tight', pad_inches = 0.5)

else:
        print('Invalid')

print('\t Network visualization created!\n')