import config as cf
import numpy as np
import pandas as pd
import os
import sys
import networkx as nx
assert cf

import argparse 

parser = argparse.ArgumentParser(description='Network simulations with defined initial conditions.')

parser.add_argument('--network_type', type=str, default='scale_free',
                    help='Network type for storing...')
parser.add_argument('--network_name', type=str, default='new_scale_free_1000',
                    help='Network type for storing...')
parser.add_argument('--num_nodes', type=int, default=1000,
                    help='Number of nodes for specific network...')
parser.add_argument('--beta', type=float,
                    help='Specify the infection probability')
parser.add_argument('--sigma', type=float,
                    help='Specify the awareness')
parser.add_argument('--type_sim', default='global',type=str, 
                    help='For running local or global simulation')
parser.add_argument('--num_checkpoints', type=int, default=8,
                    help='Number of checkpoints per initial condition so save...')
parser.add_argument('--n_iters', default=20,type=int, 
                    help='Number of iterations')
parser.add_argument('--max_time', default=150,type=int, 
                    help='Number of days of simulation')

args = parser.parse_args()

main_path = os.path.split(os.getcwd())[0] + '/Epidemiology_behavior_dynamics'
config_path = main_path + '/config.csv'
config_data = pd.read_csv(config_path, sep=',', header=None, index_col=0)

networks_path = config_data.loc['networks_dir'][1]
results_path  = config_data.loc['results_dir'][1]
num_nodes     = args.num_nodes

G = nx.read_gpickle( os.path.join(main_path, networks_path, str(num_nodes), args.network_name) )


df_params = pd.DataFrame(columns=['beta', 'sigma', 'R0'])

# Selected sigmas and betas
select_sigmas = [args.sigma]
select_betas  = [args.beta]
gamma         = 1/7
repro_numbers = list(np.array(select_betas)/gamma)

df_params['beta']      = select_betas
df_params['sigma']     = select_sigmas
df_params['R0']        = repro_numbers

sys.path.append('../')
from models import models 


initConditions = pd.read_csv(os.path.join(main_path,'run','init_conditions','initial_conditions.csv'))


if args.type_sim=='local':
    local = True
elif args.type_sim=='global':
    local = False

print('Running simulations for {} network in {}  scheme\n'.format(args.network_type, args.type_sim))

from tqdm import tqdm

for i in tqdm(range(0,len(initConditions.index)), total = len(initConditions.index)):
    print('Solving for '+str(i)+' initial params')

    for idx, r in df_params.iterrows():

        # Get initial parameters
        initInf_i = initConditions['I'][i]
        initDef_i = np.fromstring(initConditions['D'][i], sep = '|').astype(int)

        # Parameters
        model_params = {}
        model_params['time2Recover']  = 7
        model_params['probInfect']    = r['beta']
        model_params['awareness']     = r['sigma']
        model_params['initInfected']  = initInf_i
        model_params['initDefectors'] = initDef_i

        if not os.path.isdir( os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', args.type_sim, args.network_type) ):
            os.makedirs(os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', args.type_sim, args.network_type))
            
        if not os.path.isdir( os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', args.type_sim, args.network_type, 'checkpoints', 'ic_0{}'.format(i)) ):
            os.makedirs(os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', args.type_sim, args.network_type, 'checkpoints', 'ic_0{}'.format(i)))

        path_to_save_checkpoints = os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', args.type_sim, args.network_type, 'checkpoints', 'ic_0{}'.format(i))
        path_to_save_response    = os.path.join(results_path, str(num_nodes)+'_seed_checkpoints_new', args.type_sim, args.network_type, 'dynamics_beta_{}_sigma_{}'.format(r['beta'], r['sigma']) +'.csv')

        df_response = models.run_model(models.sis_replicator, G, params=model_params, n_iters=args.n_iters, max_time=args.max_time, num_checkpoints=args.num_checkpoints, local=local, path_to_save_checkpoints= path_to_save_checkpoints)
        df_response.to_csv( path_to_save_response )

print('\t DONE!\n')