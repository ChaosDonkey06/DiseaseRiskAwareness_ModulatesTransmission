import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import argparse 
import os
from tqdm import tqdm


parser = argparse.ArgumentParser(description='Heatmaps figures.')

parser.add_argument('--network_type', type=str, default='scale_free',
                    help='Network type for storing...')
parser.add_argument('--network_name', type=str, default='scale_free_1000',
                    help='Network type for storing...')
parser.add_argument('--beta_function', type=str, default='exponential',
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
awareness_path = config_data.loc['sigma_search_dir'][1]
infection_prob_path = config_data.loc['beta_search_dir'][1]
num_nodes     = args.num_nodes

sigma_search = pd.read_csv(awareness_path, dtype={'key':str, 'value':float})
beta_search  = pd.read_csv(infection_prob_path, dtype={'key':str, 'value':float})

time2Recover = (1/7)  # gamma

df = pd.concat([sigma_search, beta_search], axis=1)

df_param_run = pd.DataFrame(columns=['beta_key', 'sigma_key', 'beta_val', 'sigma_val', 'R0'])

beta_key  = []
sigma_key = []
beta_val  = []
sigma_val = []
R0_val    = []
for idx_sigma , r_sigma in sigma_search.iterrows():
    for idx_beta , r_beta in beta_search.iterrows():
        beta_key.append( r_beta['key']   )
        sigma_key.append( r_sigma['key'] )
        beta_val.append( r_beta['value'] )
        sigma_val.append( r_sigma['value'] )
        R0_val.append( r_beta['value'] / time2Recover )  # R0 = beta/gamma

df_param_run['beta_key'] = beta_key  
df_param_run['sigma_key'] = sigma_key 
df_param_run['beta_val'] = beta_val  
df_param_run['sigma_val'] = sigma_val 
df_param_run['R0'] = R0_val

list_df = []
for idx, r in tqdm(df_param_run.iterrows(), total=df_param_run.shape[0]):

    path_to_results = os.path.join(results_path, 
                                   args.beta_function, 
                                   str(num_nodes), 
                                   args.type_sim, 
                                   args.network_type, 
                                   'dynamics_beta_{}_sigma_{}'.format(r['beta_key'], r['sigma_key']) +'.csv')
    res = pd.read_csv(path_to_results, usecols=['sim_id', 'time', 'S', 'I', 'C','D'])

    # Calculate mean over iterations.
    res = res.groupby('time').mean()/num_nodes
    res = res.reset_index()[['time', 'S', 'I', 'C','D']]
    res['beta']  = r['beta_val']
    res['sigma'] = r['sigma_val']
    res['R0']    = r['R0']

    res['beta_key']  = r['beta_key']
    res['sigma_key'] = r['sigma_key']

    list_df.append(res)

df_response = pd.concat(list_df)



def return_pivoted_df(df_to_pivot, param, var='I'):
    df_heat_map = df_to_pivot.copy()
    df_heat_map = df_heat_map.pivot(param, 'sigma', var)
    df_heat_map = df_heat_map.iloc[::-1]

    return df_heat_map


def create_heatmaps(res_df=df_response, yaxis = args.type_hm, type_sim = args.type_sim, net_type = args.network_type, figs_path = figures_path):
    ## Create a heatmap for the infected and cooperation fraction given the results

    df_response_lastweek = df_response.copy()
    df_response_lastweek = df_response_lastweek.query("time >= 142")
    # steady state
    df_response_lastweek = df_response_lastweek.groupby([yaxis, 'sigma']).mean()[['S', 'I', 'C','D']].reset_index()

    heatmap_epid = return_pivoted_df(df_response_lastweek, yaxis, 'I')
    heatmap_game = return_pivoted_df(df_response_lastweek, yaxis, 'C')

    title_epid_hm = r'Inf. fraction $\bar{}$ - {}'.format('I',type_sim)
    title_game_hm = r'Coop. fraction $\bar{}$ - {}'.format('c',type_sim)
    sigma_label   = r'Awareness $\sigma$'
    path_s        = 'heatmaps/{}'.format(type_sim)

    ylabels = [r'$R_{0}$', r'$\beta_{max}$']
    if yaxis == 'R0':
        ylabel = ylabels[0]
    elif yaxis == 'beta':
        ylabel = ylabels[1]

    ## Epidemic heatmap
    fig, ax = plt.subplots(1,1,figsize=(6, 5))
    sns.heatmap(ax= ax, data=heatmap_epid, cmap='gist_heat_r',  cbar=True, vmin=0.0, vmax=1.0)
    cax = plt.gcf().axes[-1]
    cax.tick_params(labelsize=16) 
    ax.set_title(title_epid_hm, fontsize=21) 
    ax.set_xlabel(sigma_label, fontsize=22)
    ax.set_ylabel(ylabel, fontsize=21)

    xticks = heatmap_epid.columns
    keptxticksidx = np.linspace(0,len(xticks),6)
    xtickslabels = list(xticks[ np.maximum(keptxticksidx.astype(int)-1,0) ])
    xtickslabels = ['{:.1f}'.format(l) for l in xtickslabels]
    ax.set_xticks(keptxticksidx)
    ax.set_xticklabels(xtickslabels, fontsize=20, rotation=0)

    yticks = heatmap_epid.index
    keptyticksidx = np.linspace(0,len(yticks),6)
    ytickslabels = list(yticks[ np.maximum(keptyticksidx.astype(int)-1,0) ])
    ytickslabels = ['{:.1f}'.format(l) for l in ytickslabels]
    ax.set_yticks(keptyticksidx)
    ax.set_yticklabels(ytickslabels, fontsize=20)
    plt.tight_layout()

    # Save epidemic stability heatmap 

    if not os.path.isdir( figures_path ):
        os.makedirs( figures_path )
        
    if not os.path.isdir( os.path.join(figures_path, 'test', args.beta_function, path_s, net_type ) ):
        os.makedirs( os.path.join(figures_path, 'test', args.beta_function, path_s, net_type ) )

    path_save = os.path.join(figures_path, 'test', args.beta_function, path_s, net_type )
        
    plt.savefig(os.path.join(path_save, 'epid_heatmap_{}_{}_yaxis_{}.png'.format(net_type,type_sim,yaxis)), 
                        dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)

    ## Game heatmap
    fig, ax = plt.subplots(1,1,figsize=(6, 5))

    sns.heatmap(ax= ax, data=heatmap_game, cmap='RdYlGn',  cbar=True, vmin=0.0, vmax=1.0)
    cax = plt.gcf().axes[-1]
    cax.tick_params(labelsize=16) 
    ax.set_title(title_game_hm, fontsize=21) 
    ax.set_xlabel(sigma_label, fontsize=22) 
    ax.set_ylabel(ylabel, fontsize=21)

    xticks = heatmap_game.columns
    keptxticksidx = np.linspace(0,len(xticks),6)
    xtickslabels = list(xticks[ np.maximum(keptxticksidx.astype(int)-1,0) ])
    xtickslabels = ['{:.1f}'.format(l) for l in xtickslabels]
    ax.set_xticks(keptxticksidx)
    ax.set_xticklabels(xtickslabels, fontsize=20, rotation=0)

    yticks = heatmap_game.index
    keptyticksidx = np.linspace(0,len(yticks),6)
    ytickslabels = list(yticks[ np.maximum(keptyticksidx.astype(int)-1,0) ])
    ytickslabels = ['{:.1f}'.format(l) for l in ytickslabels]
    ax.set_yticks(keptyticksidx)
    ax.set_yticklabels(ytickslabels, fontsize=20)

    plt.tight_layout()
    # Save game dynamic heatmap 
    plt.savefig(os.path.join(path_save, 'game_heatmap_{}_{}_yaxis_{}.png'.format(net_type,type_sim,yaxis)), 
                             dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)


### RUN ###

create_heatmaps()
print('Heatmaps created for {} {}'.format(args.type_sim,args.network_type))

### Save lables
import matplotlib as mpl

fig, ax = plt.subplots(figsize=(11.5,2))
fig.subplots_adjust(bottom=0.5)

norm = mpl.colors.Normalize(vmin=0, vmax=1)

fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='gist_heat_r'),
                                   cax=ax, orientation='horizontal')
ax.tick_params(labelsize=22)
ax.set_xlabel(r'Infected fraction $\bar{I}$', fontsize=22)

if not os.path.isdir( os.path.join(figures_path, 'heatmaps', str(num_nodes), 'labels') ):
    os.makedirs ( os.path.join(figures_path, 'heatmaps', str(num_nodes), 'labels') )

path_save = os.path.join(figures_path, 'heatmaps', str(num_nodes), 'labels')
plt.savefig(os.path.join(path_save, 'epid_label.png'), 
                             dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)

fig, ax = plt.subplots(figsize=(11.5,2))
fig.subplots_adjust(bottom=0.5)

norm = mpl.colors.Normalize(vmin=0, vmax=1)

fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='RdYlGn'),
                                   cax=ax, orientation='horizontal')
ax.tick_params(labelsize=22)
ax.set_xlabel(r'Cooperating fraction $\bar{c}$', fontsize=22)


plt.savefig(os.path.join(path_save, 'game_label.png'), 
                             dpi=400, transparent = False, bbox_inches = 'tight', pad_inches = 0.1)
