import config as cf
import pandas as pd
import os
import sys
import networkx as nx
from models import models_beta_test
assert cf

import argparse 

parser = argparse.ArgumentParser(description='Network simulations.')

parser.add_argument('--network_type', type=str, default='scale_free',
                    help='Network type for loadind and storing...')
parser.add_argument('--network_name', type=str, default='scale_free_1000',
                    help='Network name for loading and storing...')
parser.add_argument('--beta_function', type=str, default='exponential',
                    help='Network name for loading and storing...')
parser.add_argument('--num_nodes', type=int, default=1000,
                    help='Number of nodes for specific network...')
parser.add_argument('--type_sim', default='global',type=str, 
                    help='For running local or global simulation')
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
awareness_path = config_data.loc['sigma_search_dir'][1]
infection_prob_path = config_data.loc['beta_search_dir'][1]
num_nodes     = args.num_nodes

sigma_search = pd.read_csv(awareness_path, dtype={'key':str, 'value':float})
beta_search  = pd.read_csv(infection_prob_path, dtype={'key':str, 'value':float})

G = nx.read_gpickle( os.path.join(main_path, networks_path, str(num_nodes), args.network_name) )

df = pd.concat([sigma_search, beta_search], axis=1)

df_param_run = pd.DataFrame(columns=['beta_key', 'sigma_key', 'beta_val', 'sigma_val'])

beta_key  = []
sigma_key = []
beta_val  = []
sigma_val = []
for idx_sigma , r_sigma in sigma_search.iterrows():
    for idx_beta , r_beta in beta_search.iterrows():

        beta_key.append( r_beta['key']   )
        sigma_key.append( r_sigma['key'] )
        beta_val.append( r_beta['value'] )
        sigma_val.append( r_sigma['value'] )

df_param_run['beta_key'] = beta_key  
df_param_run['sigma_key'] = sigma_key 
df_param_run['beta_val'] = beta_val  
df_param_run['sigma_val'] = sigma_val 

if args.type_sim=='local':
    local = True
elif args.type_sim=='global':
    local = False
n_iters = args.n_iters
print('Running simulations for {} network in {}  scheme\n'.format(args.network_type, args.type_sim))

sys.path.append('../')
from models import models

from tqdm import tqdm
for idx, r in tqdm(df_param_run.iterrows()):
    # Parameters
    model_params = {}
    model_params['time2Recover']  = 7
    model_params['probInfect']    = r['beta_val']
    model_params['awareness']     = r['sigma_val']
    model_params['initInfected']  = None
    model_params['initDefectors'] = None

    if not os.path.isdir( os.path.join(results_path, args.beta_function, str(num_nodes), args.type_sim, args.network_type) ):
        os.makedirs(os.path.join(results_path, args.beta_function, str(num_nodes), args.type_sim, args.network_type))
        
    if not os.path.isdir( os.path.join(results_path, args.beta_function, str(num_nodes), args.type_sim, args.network_type, 'checkpoints') ):
        os.makedirs(os.path.join(results_path, args.beta_function, str(num_nodes), args.type_sim, args.network_type, 'checkpoints'))

    path_to_save_checkpoints = os.path.join(results_path,
                                            args.beta_function, 
                                            str(num_nodes), 
                                            args.type_sim, 
                                            args.network_type, 
                                            'checkpoints', 
                                            'beta_{}_sigma_{}'.format(r['beta_key'], r['sigma_key'])+'.csv' )
    path_to_save_response    = os.path.join(results_path,
                                            args.beta_function, 
                                            str(num_nodes), 
                                            args.type_sim, 
                                            args.network_type, 
                                            'dynamics_beta_{}_sigma_{}'.format(r['beta_key'], r['sigma_key']) +'.csv')

    if os.path.exists(path_to_save_response):
        continue
    print( 'Running for beta={}, sigma={} \r'.format(r['beta_val'], r['sigma_val']) )

    df_response = models_beta_test.run_model(models_beta_test.sis_replicator, 
                                            G, 
                                            params=model_params, 
                                            n_iters=n_iters, 
                                            max_time=args.max_time, 
                                            num_checkpoints=0, 
                                            local=local, 
                                            beta_fun=args.beta_function, 
                                            path_to_save_checkpoints= path_to_save_checkpoints)
    df_response.to_csv( path_to_save_response )
print('\t DONE!\n')