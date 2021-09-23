import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from tqdm import tqdm
import seaborn as sns
import numpy as np
import pandas as pd
import argparse 
import os


parser = argparse.ArgumentParser(description='Time dynamics figures.')

parser.add_argument('--network_type', type=str, default='scale_free',
                    help='Network type for storing...')
parser.add_argument('--network_name', type=str, default='new_scale_free_1000',
                    help='Network type for storing...')
parser.add_argument('--num_nodes', type=int, default=1000,
                    help='Number of nodes for specific network...')
parser.add_argument('--type_sim', default='local',type=str, 
                    help='For running local or global simulation')
parser.add_argument('--type_hm', default='R0',type=str, 
                    help='Define the yaxis in heatmaps (R0 or beta)')

args = parser.parse_args()

main_path = os.path.split(os.getcwd())[0] + '/Epidemiology_behavior_dynamics'
config_path = main_path + '/config.csv'
config_data = pd.read_csv(config_path, sep=',', header=None, index_col=0)

networks_path = config_data.loc['networks_dir'][1]
results_path  = config_data.loc['results_dir'][1]
figures_path  = config_data.loc['figures_dir'][1]
awareness_path = config_data.loc['sigma_plot_dir'][1]
infection_prob_path = config_data.loc['beta_plot_dir'][1]
num_nodes     = args.num_nodes

sigma_plot = pd.read_csv(awareness_path, dtype={'key':str, 'value':float})
beta_plot  = pd.read_csv(infection_prob_path, dtype={'key':str, 'value':float})

time2Recover = (1/7)  # gamma

df = pd.concat([sigma_plot, beta_plot], axis=1)

df_params = pd.DataFrame(columns=['beta_key', 'sigma_key', 'beta', 'sigma', 'R0'])

beta_key  = []
sigma_key = []
beta  = []
sigma = []
R0    = []
for idx_sigma , r_sigma in sigma_plot.iterrows():
    for idx_beta , r_beta in beta_plot.iterrows():

        beta_key.append( r_beta['key']   )
        sigma_key.append( r_sigma['key'] )
        beta.append( r_beta['value'] )
        sigma.append( r_sigma['value'] )
        R0.append( r_beta['value'] / time2Recover )  # R0 = beta/gamma

df_params['beta_key'] = beta_key  
df_params['sigma_key'] = sigma_key 
df_params['beta'] = beta  
df_params['sigma'] = sigma 
df_params['R0'] = R0

colors_plt = [ 'tab:red', 'royalblue', 'green', 'tab:purple', 'tab:cyan', 'tab:orange' ]

# Read results
fig, ax = plt.subplots(1,2,figsize=(20, 7))

for idx, r in tqdm(df_params.iterrows()):    
    #Read global results
    path_to_results_local = os.path.join(results_path, str(num_nodes), 'local', args.network_type, 'dynamics_beta_{}_sigma_{}'.format(r['beta_key'], r['sigma_key']) +'.csv')
    res_local = pd.read_csv(path_to_results_local, usecols=['sim_id', 'time', 'S', 'I', 'C','D'])
    res_local_plot = res_local.copy()
    res_local_plot[['S','I','C','D']] = res_local_plot[['S','I','C','D']]/num_nodes
    res_local_plot['type'] = ['local'] * len(res_local_plot)

    #Read local results
    path_to_results_glob = os.path.join(results_path, str(num_nodes), 'global', args.network_type, 'dynamics_beta_{}_sigma_{}'.format(r['beta_key'], r['sigma_key']) +'.csv')
    res_glob = pd.read_csv(path_to_results_glob, usecols=['sim_id', 'time', 'S', 'I', 'C','D'])
    res_glob_plot = res_glob.copy()
    res_glob_plot[['S','I','C','D']] = res_glob_plot[['S','I','C','D']]/num_nodes
    res_glob_plot['type'] = ['global'] * len(res_glob_plot)

    df_res = [res_local_plot,res_glob_plot]
    df_res_c = pd.concat(df_res)

    # ········· Plot global············

    # Plot disease dynamics

    sns.lineplot( ax = ax[0],
                  data = df_res_c, 
                  x = 'time', y = 'I',
                  label = r'$R_0$={:.1f},$\sigma={:.1f}$'.format(r['R0'],r['sigma']),
                  style='type',
                  color = colors_plt[idx] ,
                  alpha=0.5)
    #ax[0].get_legend().remove()
    ax[0].lines[0].set_linestyle("--")
    ax[0].set_title(r'Disease dynamics $R_0 = {}$'.format(r['R0']),fontsize=22)
    #ax[0].set_xlabel('')
    ax[0].set_xlabel(r'Days',fontsize=21)
    #ax[0].set_xticks(fontsize=20)
    #ax[0].tick_params(axis='x', labelsize=20)
    ax[0].xaxis.set_tick_params(labelsize=20)
    ax[0].yaxis.set_tick_params(labelsize=20)

    ax[0].set_xlim([-0.1,151])
    ax[0].set_ylabel(r'Inf. Fraction ($I$)',fontsize=21)
    #ax[0].set_yticks(fontsize=20)
    ax[0].set_ylim([-0.1,1.1])

    # Plot game dynamics

    sns.lineplot( ax = ax[1],
                  data = df_res_c, 
                  x = 'time', y = 'C',
                  style='type',
                  color = colors_plt[idx],
                  alpha=0.5)
    #ax[1].get_legend().remove()
    ax[1].lines[0].set_linestyle("--")
    ax[1].set_title(r'Behavioral dynamics $R_0 = {}$'.format(r['R0']),fontsize=22)
    ax[1].set_xlabel(r'Days',fontsize=21)
    ax[1].xaxis.set_tick_params(labelsize=20)
    ax[1].yaxis.set_tick_params(labelsize=20)
    #ax[1].set_xticks(fontsize=20)
    ax[1].set_xlim([-0.1,151])
    ax[1].set_ylabel(r'Coop. Fraction ($c$)',fontsize=21)
    #ax[1].set_yticks(fontsize=20)
    ax[1].set_ylim([-0.1,1.1])
    plt.tight_layout()

if not os.path.isdir( os.path.join(figures_path, 'dynamics', str(num_nodes)) ):
        os.makedirs( os.path.join(figures_path, 'dynamics', str(num_nodes)) )

path_save = os.path.join(figures_path, 'dynamics', str(num_nodes))

plt.savefig(os.path.join(path_save, '{}_beta_{}_dynamics.png'.format(args.network_type,str(r['sigma_key']))), 
                            dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)
plt.show()


# Figure for legend labels
fig, ax = plt.subplots(2,1,figsize=(9, 8))
for idx, r in tqdm(df_params.iterrows()):
    # Read global results
    path_to_results_local = os.path.join(results_path, str(num_nodes), 'local', args.network_type, 'dynamics_beta_{}_sigma_{}'.format(r['beta_key'], r['sigma_key']) +'.csv')
    res_local = pd.read_csv(path_to_results_local, usecols=['sim_id', 'time', 'S', 'I', 'C','D'])
    res_local_plot = res_local.copy()
    res_local_plot[['S','I','C','D']] = res_local_plot[['S','I','C','D']]/num_nodes
    res_local_plot['type'] = ['global'] * len(res_local_plot)


    sns.lineplot( ax = ax[0],
                  data = res_local_plot, 
                  x = 'time', y = 'I',
                  label = r'$\sigma={:.1f}$'.format(r['sigma']),
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

    plt.figlegend(bbox_to_anchor=(0.7,0.4), fontsize=22)

    if not os.path.isdir( os.path.join(figures_path, 'dynamics', str(num_nodes)) ):
        os.makedirs( os.path.join(figures_path, 'dynamics', str(num_nodes)) )

    path_save = os.path.join(figures_path, 'dynamics', str(num_nodes))

    plt.savefig(os.path.join(path_save, 'colorlabel_scale_free_beta_{}_dynamics.png'.format(str(r['beta_key']))), 
                                dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)


    # Figure for legend sim types
    fig, ax = plt.subplots(2,1,figsize=(9, 8))

    # Read global results
    path_to_results_local = os.path.join(results_path, str(num_nodes), 'local', args.network_type, 'dynamics_beta_{}_sigma_{}'.format('060', '100') +'.csv')
    res_local = pd.read_csv(path_to_results_local, usecols=['sim_id', 'time', 'S', 'I', 'C','D'])
    res_local_plot = res_local.copy()
    res_local_plot[['S','I','C','D']] = res_local_plot[['S','I','C','D']]/num_nodes
    res_local_plot['type'] = ['Local'] * len(res_local_plot)

    # Read local results
    path_to_results_glob = os.path.join(results_path, str(num_nodes), 'global', args.network_type, 'dynamics_beta_{}_sigma_{}'.format('060', '100') +'.csv')
    res_glob = pd.read_csv(path_to_results_glob, usecols=['sim_id', 'time', 'S', 'I', 'C','D'])
    res_glob_plot = res_glob.copy()
    res_glob_plot[['S','I','C','D']] = res_glob_plot[['S','I','C','D']]/num_nodes
    res_glob_plot['type'] = ['Global'] * len(res_glob_plot)

    df_res = [res_local_plot,res_glob_plot]
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

    if not os.path.isdir( os.path.join(figures_path, 'dynamics', str(num_nodes), 'labels') ):
        os.makedirs( os.path.join(figures_path, 'dynamics', str(num_nodes), 'labels') )

    path_save = os.path.join(figures_path, 'dynamics', str(num_nodes), 'labels')
    plt.savefig(os.path.join(path_save, 'sylelabel_beta_{}_dynamics.png'.format(str(r['beta_key']))), 
                                dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)
