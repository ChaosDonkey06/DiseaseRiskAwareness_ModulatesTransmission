# How disease risk awareness modulates transmission: coupling infectious disease models with behavioral dynamics
Code for the preprint https://t.co/cwj59tO9XP?amp=1
𝗧𝗶𝘁𝗹𝗲: How disease risk awareness modulates transmission: coupling infectious disease models with behavioral dynamics.
𝗔𝘂𝘁𝗵𝗼𝗿𝘀: Jaime Cascante-Vega, Samuel Torres-Florez, Juan Cordovez, Mauricio Santos-Vega

## Dependencies
    python 3.7
    pandas
    numpy
    scipy
    networkx
    python-louvian 

## Setup
1. Make sure to have dependencies installed
2. Clone our repository `https://github.com/biomac-lab/Epidemiology_behavior_dynamics.git`
3. Create networks
4. Simulate (may take a lot of time!)
5. Generate figures

## Usage

### Creation of networks using *networkx*

In this investigation, we account for three different types of networks where we ran simulations on: Scale-free, Watts-Strogatz small world and grid topologies. To create these networks that are used in simulation, you need to define the number of nodes (*n*) with `--num_nodes <n>`. If not, it is set to *n=1000* automatically. Hence, if you want to create all the networks specify `--all True`. To do this, run:

    python models/create_networkxs.py --num_nodes 1000 --all True 

Conversely, if you want to create a specific network (say `scale_free`, `small_world` or `grid` ) with a given number of nodes (*n*), set `--specific_network <network_name>`. To create a scale-free network with 5000 nodes, run:

    python models/create_networkxs.py --num_nodes 5000 --specific_network scale_free

### Run simulations

In order to run a simulation, a combination of the infection probability (*beta*) and awareness (*sigma*) must be given. This combination are specified on the files `beta_search.csv` and `sigma_search.csv` located on the `run/param_search` folder. Each file will contain a row with the parameter value and a key for saving the result (as shown below). We test the model on an interval between 0-1, this means the files `<>_search.csv` contains values in the range of 0 and 1 with space intervals of 0.02 as decribed below:

    key,value
    000,0.00,
    002,0.02,
    .
    .
    098,0.98,
    100,1.00,

The default setting for `sigma_search.csv` are:

    key,value
    100,1.00
    070,0.70
    050,0.50

And for `beta_search.csv` are:

    key,value
    060,0.60
    070,0.70

Besides, the other entries to the model are the `network_type`, the `network_name` (output name after creating the network(s)), the `type_sim` which may be specified as `global` or `local`, the number of iterations `n_iters` (20 by default), and the length of simulation given by `max_time` which is set to 150 days. If you want to change the number of iterations (*iters*), add `--n_iters <iters>` to the command. To change the length of simulation (*days*) add `--max_time <days>`. The execution is then: `python run/run_sims.py --network_type <> --network_name <> --num_nodes <> --type_sim <>`.

##### For running simulations over a scale-free network with 1000 nodes in both information transmission scheme (as shown in the paper)

    python run/run_sims.py --network_type scale_free --network_name scale_free_1000 --num_nodes 1000 --type_sim local
    python run/run_sims.py --network_type scale_free --network_name scale_free_1000 --num_nodes 1000 --type_sim global

##### For running simulations over a scale-free network with 5000 nodes in both infomation transmission scheme (as shown in the paper)

    python run/run_sims.py --network_type scale_free --network_name scale_free_5000 --num_nodes 5000 --type_sim local
    python run/run_sims.py --network_type scale_free --network_name scale_free_5000 --num_nodes 5000 --type_sim global

##### For running simulations over a small-world and grid network with 1000 nodes in both infomation transmission scheme (as shown in the paper)

    python run/run_sims.py --network_type small_work --network_name small_work_1000 --num_nodes 1000 --type_sim local
    python run/run_sims.py --network_type small_work --network_name small_work_1000 --num_nodes 1000 --type_sim global
    python run/run_sims.py --network_type grid --network_name grid_1000 --num_nodes 1000 --type_sim local
    python run/run_sims.py --network_type grid --network_name grid_1000 --num_nodes 1000 --type_sim global

##### For running simulations over a ODE (as shown in the paper)

    python run/run_ode.py

#### Simulation for cluster analysis

In order to analyse how the community structures (i.e. clusters or hubs) affected infection and behavior, several initial conditions are tested. To generate this initial conditions for a network with 5000 individuals, run:

    python run/init_conditions/create_init_conditions.py --num_nodes 5000

Each iteration of the simulation records the state of each invididual node and is saved in a `.txt` file. You will need to specify `network_type`, `network_name` and the `type_sim` as mentioned above. The infection probability `beta` and awareness `sigma` are also required. The number of checkpoitns to be saved is 8, however you can change it by addind `--num_checkpoints <num>` to the command. If you want to change the number of iterrations (*iters*), add `--n_iters <iters>` to the command. To change the length of simulation (*days*) add `--max_time <days>`. The execution is then: `python run/init_conditions/run_checkpoints.py --network_type <> --network_name <> --num_nodes <> --type_sim <> --beta <> --sigma <>`. Remember that the number of nodes needs to be specified in the `config.csv` file.

##### For for saving the checkpoint of the simulations executed over a scale-free network with 5000 nodes in both infomation transmission scheme (as shown in the paper)

    python run/run_checkpoints.py --network_type scale_free --network_name scale_free_5000 --num_nodes 5000 --type_sim local --beta 6.0 --sigma 1.0
    python run/run_checkpoints.py --network_type scale_free --network_name scale_free_5000 --num_nodes 5000 --type_sim global --beta 6.0 --sigma 1.0

### Figures generation

#### Networks visualization

In order to visualize the networks you will need to specify the `network_type` (i.e.  `scale_free`, `small_world ` or `grid`) and the `network_name` as it was saved in the `/networks` folder. For visualizing a (already) created scale-free network with 5000 nodes, you would run:

    python plots/plot_networks.py --network_type scale_free --network_name scale_free_5000

#### Heatmaps visualization

Plotting a heatmap means you have runned simulations over a range of values for `beta` and `sigma`. Hence, `beta_search.csv` and `sigma_search.csv` had to be modified previously. The values that will be considered by the function for plotting are the same already specified on the `/param_search` files. The `network_type`, `network_name`, `num_nodes`, `type_sim` and `type_hm` parameters need to be specified. The later defines weather the heatmap if ploted using *beta* of *R0* in their *y* axis. An example is displayed as follows:

    python plots/plot_heatmaps.py --network_type scale_free --network_name --scale_free_5000 --num_nodes 5000 --type_sim local --type_hm R0
    python plots/plot_heatmaps.py --network_type scale_free --network_name --scale_free_5000 --num_nodes 5000 --type_sim global --type_hm R0

#### Disease and behavior dynamics over networks

In order to plot the temporal dynamics of the simulation, there are few things to consider. Each figure is created for a specific *beta* and a set of *sigmas*. These values need to be specified in the `beta_plot.csv` and `sigma_plot.csv` files found in the `/plots/params_plot` folder. Thee defauld values are the ones used in the paper figures. We recommend to plot only 3 values of *sigma* per *beta* for best visualization, and only call the function for a set of parameters (leaving only one value of *beta* in `beta_plot.csv`). The `network_type`, `network_name`, `num_nodes`, and `type_hm` need to be specified. For visualizing the plot presented in our paper, run:

    python plots/plot_dynamics.py --network_type scale_free --network_name --scale_free_5000 --num_nodes 5000 --type_hm R0
    python plots/plot_dynamics.py --network_type scale_free --network_name --scale_free_5000 --num_nodes 5000 --type_hm R0

#### Clustered dynamics over a scale-free network

We performe cluster analysis through ploting each cluster/community dynamics over time, we do this for the 3 biggest clusters. To achive this, one must have runned simulation for each time checkpoint (explained earlier). Here, we need to specify the `network_type`, `network_name`, `num_nodes`, `iterations`, `beta_select`, `beta_key`, `sigma_select`, and `sigma_key`. Note that here you need to specify the *beta* and *sigma* values and keys. Then, for visualizing clusters in the scale-free network with 5000 individuals we should run:

    python plots/plot_cluster_dynamics.py --network_type scale_free --network_name scale_free_5000 --num_nodes 5000 --iterations 10 --beta_select 0.6 --beta_key 060 --sigma_select 1.0 --sigma_key 100
