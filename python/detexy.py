import pysam
import sys
import pandas as pd
import numpy as np
import os
import yaml
import warnings
import networkx as nx
import matplotlib
# matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

class detexy():

    def __init__(self):

        self.breakpoint_distance = 150000
        self.reciprocal_distance = 200
        self.detexy_plots_dir = './'
        self.chain_properties = {
            'single_svtype': ['TRA', 'INV'],
        }

    def configure(self, config, class_name='detexy'):

        """
        python3 configure(config, class_name)

        Input Arguments:

            - config: path to config file
            - class_name: section of config containing class-specific properties

        Overwrites any default properties with those passed in the config.

        """

        with open(config, 'r') as f:
            options = yaml.safe_load(f)
            options = options.get(class_name, {})
            PYSAM_ZERO_BASED = options.get('PYSAM_ZERO_BASED', 1)
        
        self.PYSAM_ZERO_BASED = PYSAM_ZERO_BASED
        for option in options:
            if isinstance(options[option], dict):
                default = getattr(self, option)
                nested = options[option]
                for opt in nested.keys():
                    default_arg = default.get(opt, 'MISSING')
                    if default_arg is 'MISSING':
                        warnings.warn('Argument {arg} specified in config is not an allowed property.'.format(arg=opt))
                        continue
                    if not isinstance(nested[opt], type(default_arg)) and default_arg is not None:
                        warnings.warn('Argument {arg} specified in config is not the same type as the default: {default} = {d}, default type: {t}, passed type: {tp}'.format(arg=opt, default=opt, d=default_arg, t=type(default_arg), tp=type(nested[opt])))
                        continue
                    default[opt] = nested[opt]
                setattr(self, option, default)
            else:
                default_arg = getattr(self, option, 'MISSING')
                if default_arg is 'MISSING':
                        warnings.warn('Argument {arg} specified in config is not an allowed property.'.format(arg=opt))
                        continue
                if not isinstance(options[option], type(default_arg)) and default_arg is not None:
                    warnings.warn('Argument {arg} specified in config is not the same type as the default: {default} = {d}, default type: {t}, passed type: {tp}'.format(arg=opt, default=opt, d=default_arg, t=type(default_arg), tp=type(options[option])))
                    continue
                setattr(self, option, options[option])

    def run_detexy(self, tumor_tab, passed_tumor_df=False, plot_results=True):

        """

        python run_detexy(tumor_tab, passed_tumor_df=False, return_graphs=False, plot_results=True)

        Input Arguments:

            - tumor_tab: path to Pype tumor tab file, or pandas DataFrame.
            - passed_tumor_df: Pass True if tumor_tab passed as pandas DataFrame.
            - plot_results: Pass True if chromoplectic chains should be plotted and saved to disk.


        Returns:

            - chromoplectic_results: Pype tab file with the following columns added:
                event_id, chain_id, reciprocal, overlapping

        """

        if passed_tumor_df is False:
            df_tumor = pd.DataFrame.from_csv(tumor_tab, sep="\t", index_col=False)
        else:
            df_tumor = tumor_tab

        df_tumor = self.clean_chromosomes(df_tumor)
        df_tumor = self.add_event_id(df_tumor)
        df = df_tumor[['Chromosome1', 'Chromosome2', 'Position1', 'Position2', 'event_id', 'SVTYPE']].copy().drop_duplicates()
        
        list_breakpoints = self.unnest_translocations(df)
        edges = self.get_nearby_events(df, list_breakpoints)

        if edges is None:
            df_tumor['chain_id'] = np.nan
            df_tumor['single_svtype'] = np.nan
            df_tumor['reciprocal'] = np.nan
            df_tumor['overlapping'] = np.nan
            warnings.warn("No chromoplectic events detected. Consider increasing the breakpoint distance.")
            return df_tumor

        full_graph, subgraphs = self.create_graph(nodes=df, edges=edges)
        chromoplectic, cgraphs = self.identify_chromoplectic_translocations(df, subgraphs)

        chromoplectic = self.validate_chains(chromoplectic)

        chromoplectic_results = pd.merge(df_tumor, chromoplectic, how='left', on=['Chromosome1', 'Chromosome2', 'Position1', 'Position2'])

        if plot_results is True:
            if isinstance(tumor_tab, str):
                output_file = tumor_tab.strip('.tab') + '.detexy.plots.png'
            else:
                output_file = self.detexy_plots_dir + 'detexy.plots.png'
            self.plot_detexy_results(output_file, chromoplectic_results, cgraphs)

        return chromoplectic_results

    def clean_chromosomes(self, df_tumor):

        """

        python clean_chromosomes(df_tumor)

        Force all chromosome values to strings.

        Input Arguments:

            - df_tumor: Pype tumor tab DataFrame.

        Output Arguments:

            - df_tumor: Pype tumor tab DataFrame.

        """

        df_tumor.Chromosome1 = [str(int(x)) if x not in (['X', 'Y']) else str(x) for x in df_tumor.Chromosome1]
        df_tumor.Chromosome2 = [str(int(x)) if x not in (['X', 'Y']) else str(x) for x in df_tumor.Chromosome2]
        return df_tumor

    def add_event_id(self, df_tumor):

        """

        python add_event_id(df_tumor)

        Add a unique identifier to the Pype tumor tab DataFrame.

        Input Arguments:

            - df_tumor: Pype tumor tab DataFrame.

        Output Arguments:

            - df_tumor: Pype tumor tab DataFrame, with new colomn 'event_id'.

        """

        df_tumor['event_id'] = range(0, len(df_tumor))
        return df_tumor

    def validate_chains(self, chromoplectic_results):

        """

        python validate_chains(chromoplectic_results)

        Remove any chains that do not meet the definition of chromoplexy, despite being nearby to one another.

        Flag any chains that constitue reciprocal events.

        Flag any chains that contain overlapping intra-chromosomal rearrangements.

        """

        ineligible_chains = []
        for chain_id, data in chromoplectic_results.groupby(['chain_id']):
            sv_types = data.SVTYPE.unique()
            if len(sv_types) > 1: continue
            sv_type = sv_types[0]
            if sv_type not in self.chain_properties['single_svtype']:
                ineligible_chains.append(chain_id)
        # Flag ineligible chains
        chromoplectic_results['single_svtype'] = '0'
        chromoplectic_results.ix[chromoplectic_results.chain_id.isin(ineligible_chains), 'single_svtype'] = 'FLAGGED'

        chromoplectic_results = self.flag_reciprocal(chromoplectic_results)
        chromoplectic_results = self.flag_overlapping(chromoplectic_results)
        chromoplectic_results = chromoplectic_results.drop(['SVTYPE', 'event_id'], 1)
        return chromoplectic_results

    def flag_overlapping(self, chromoplectic):

        """
        python flag_overlapping(chromoplectic)

        Flag any chains that contain overlapping intra-chromosomal rearrangements.

        """
        chromoplectic['overlapping'] = False

        intra_chromosomal = chromoplectic[chromoplectic.Chromosome1 == chromoplectic.Chromosome2]
        overlapping_chains = []
        for chain_id, data in intra_chromosomal.groupby(['chain_id']):
            overlapping_rearrangements = data.groupby(['Chromosome1']).apply(lambda x: self._test_overlaps(x))
            if any(overlapping_rearrangements):
                overlapping_chains.append(chain_id)

        chromoplectic.ix[(chromoplectic.chain_id.isin(overlapping_chains)), 'overlapping'] = True
        
        return chromoplectic

    def _test_overlaps(self, rearrangements):

        OVERLAP = False
        for idx, row in rearrangements.iterrows():
            pos1_overlaps = rearrangements[(rearrangements.Position1 > row.Position1) & (rearrangements.Position1 < row.Position2)]
            pos2_overlaps = rearrangements[(rearrangements.Position2 > row.Position1) & (rearrangements.Position2 < row.Position2)]
            if len(pos1_overlaps) >= 1 or len(pos2_overlaps) >= 1:
                OVERLAP = True
                break
        return OVERLAP


    def flag_reciprocal(self, chromoplectic):

        """
        python flag_reciprocal(chromoplectic)

        Flag any chains that constitue reciprocal events.

        Reciprocal chromoplectic chains are those containing two nodes, with breakpoints within n bases of on another.

        n is defined by the reciprocal_distance specified in the config.

        Input Arguments:

            - chromoplectic: Pype output of rearrangements participating in chromoplectic chains.

        Output Arguments:

            - chromoplectic: Pype output of rearrangements participating in chromoplectic chains,
                with new columns 'reciprocal', 'overlapping'

        """

        chromoplectic['reciprocal'] = '0'
        nodes = chromoplectic.groupby(['chain_id']).size().reset_index().rename(columns={0:'nodes'})
        nodes = nodes[nodes.nodes == 2]
        if len(nodes) > 0:
            reciprocal = chromoplectic[chromoplectic.chain_id.isin(nodes.chain_id)]
            chromosomes = reciprocal.groupby(['chain_id']).agg({'Chromosome1': 'nunique', 'Chromosome2': 'nunique'}).reset_index()
            chromosomes = chromosomes[(chromosomes.Chromosome1 == 1)  & (chromosomes.Chromosome2 == 1)]
            if len(chromosomes) < 1: return chromoplectic
            reciprocal = reciprocal[reciprocal.chain_id.isin(chromosomes.chain_id)]
            distance = reciprocal.groupby(['chain_id']).agg({'Position1': lambda x: np.abs(x.min() - x.max()), 'Position2': lambda x: np.abs(x.min() - x.max())}).reset_index()
            within_distance = distance[(distance.Position1 < self.reciprocal_distance ) & (distance.Position2 < self.reciprocal_distance)]
            if len(within_distance) < 1: return chromoplectic
            chromoplectic.ix[chromoplectic.chain_id.isin(within_distance.chain_id), 'reciprocal'] = 'FLAGGED'

        return chromoplectic

    def unnest_translocations(self, df):

        """

        python unnest_translocations(df)

        Create a list of all structural rearrangement breakpoints:

            Split one pype record from a single row:
                Chromosome1 Chromosome2 Position1 Position2

            into two rows:

                Chromosome1 Position1
                Chromosome2 Position2

        Input Arguments:

            - df_tumor: Pype tumor tab DataFrame.

        Output Arguments:

            - df_tumor: Pype tumor tab DataFrame, with new colomn 'event_id'.

        """
        
        d = df[['Chromosome1', 'Position1', 'event_id']].rename(columns={'Chromosome1': 'chromosome', 'Position1': 'position'})
        dd = df[['Chromosome2', 'Position2', 'event_id']].rename(columns={'Chromosome2': 'chromosome', 'Position2': 'position'})
        d = pd.concat([d, dd])
        d.position = d.position.astype(int)
        return d


    def create_graph(self, nodes, edges):

        """
        Create undirected graph of translocations, with edges specifying those which are nearby.

        python create_graph(nodes, edges):

        Input Arguments:

            - nodes: list of all translocations
            - edges: list of translocations that are nearby

        Output Arguments:

            - G: graph containing every translocation in tumor tab
            - graphs: connected component subgraphs of G

        """

        G=nx.Graph()
        nodes.apply(lambda x: G.add_node(x.event_id), 1);
        edges.apply(lambda x: G.add_edge(x.from_node, x.to_node), 1);
        graphs = list(nx.connected_component_subgraphs(G))
        return (G, graphs)


    def identify_chromoplectic_translocations(self, df, graphs):

        chromoplexy = {}
        i = 0
        chromoplectic_graphs = []
        for graph in graphs:
            if len(list(graph.node.keys())) <= 1: continue
            i += 1
            chromoplexy[i] = list(graph.node.keys())
            chromoplectic_graphs.append(graph)
        translocations = self.retrieve_translocations(df, chromoplexy)
        return (translocations, chromoplectic_graphs)

    def retrieve_translocations(self, df, chromoplexy):

        i = 0
        chrom_df = []
        for key in chromoplexy.keys():
            i += 1
            chrom = df[df.event_id.isin(chromoplexy[key])].reset_index(drop=True)
            chrom['chain_id'] = i
            chrom_df.append(chrom)
        chrom_df = pd.concat(chrom_df).reset_index(drop=True)
        return chrom_df

    def get_nearby_events(self, df, unnested_df):

        """
        Iterate through pype output and associate each event with any other rearrangement within n bases of either breakpoint.

        The distance n is specified in the config.

        Returns the event ids of any nearby rearrangements.

        """

        results = []
        for idx, row in df.iterrows():
            result = row.to_frame().transpose().copy()
            result1 = pd.DataFrame()
            result2 = pd.DataFrame()
            
            first = unnested_df[unnested_df.chromosome.isin([row.Chromosome1]) & (np.abs(unnested_df.position - row.Position1) < self.breakpoint_distance) & (unnested_df.event_id != row.event_id)].copy()
            if len(first) > 0:
                first_result = result.copy()
                first_result['junction'] = '1'
                first['junction'] = '1'
                result1 = pd.merge(first_result, first, how='left', on=['junction'])
            
            second = unnested_df[unnested_df.chromosome.isin([row.Chromosome2]) & (np.abs(unnested_df.position - row.Position2) < self.breakpoint_distance) & (unnested_df.event_id != row.event_id) ].copy()
            if len(second) > 0:
                second_result = result.copy()
                second_result['junction'] = '2'
                second['junction'] = '2'
                result2 = pd.merge(second_result, second, how='left', on=['junction'])
            
            results.append(result1)
            results.append(result2)
        if len(results) < 1:
            return None
        r = pd.concat(results).reset_index(drop=True)
        if len(r) < 1:
            return None
        edges = r[['event_id_x', 'event_id_y']].rename(columns={'event_id_x': 'from_node', 'event_id_y': 'to_node'})
        return edges


    def plot_detexy_results(self, output_file, chromoplectic_results, graphs):

        chromoplectic_results['event_label'] = chromoplectic_results.apply(lambda x: x.Chromosome1 + ':' + x.Chromosome2, axis=1)

        colours = [['DEL', '#BC8F8F'], ['TRA', '#6495ED'], ['DUP', '#2E8B57'], ['INV', '#BA55D3']]
        colours = pd.DataFrame(colours, columns=['SVTYPE', 'colour'])

        labels = {}
        for idx, event in chromoplectic_results[['event_id', 'event_label']].iterrows():
            labels[str(event.event_id)] = event.event_label

        ratios = [len(g.nodes()) for g in graphs]
        node_count = np.sum(ratios)
        gs = gridspec.GridSpec(len(graphs), 1, height_ratios=ratios)

        fig = plt.figure()
        for g, axo in zip(graphs, gs):
            nodes = [str(x) for x in g.nodes()]
            nodes_for_index = [int(x) for x in g.nodes()]

            events_in_graph = chromoplectic_results[chromoplectic_results.event_id.isin(nodes_for_index)]
            chain_id = events_in_graph.chain_id.unique()[0]
            reciprocal = events_in_graph.reciprocal.unique()[0]
            events_in_graph = pd.merge(events_in_graph, colours, how='left', on='SVTYPE')
            node_colours = [events_in_graph[events_in_graph.event_id == x].colour.item() for x in nodes_for_index]
            event_type_labels = [events_in_graph[events_in_graph.event_id == x].SVTYPE.item() for x in nodes_for_index]
            graph_labels = {}
            for l in labels.keys():
                if l in nodes:
                    graph_labels[int(l)] = labels[l]
            pos = nx.spring_layout(g , iterations=200)
            ax = plt.subplot(axo)
            nx.draw(g, pos=pos, ax=ax, labels=graph_labels, node_size=1200, node_shape='o', node_color=node_colours)
            ax.set_title('Chain ID: {c}, reciprocal: {rec}'.format(c=chain_id, rec=reciprocal))

        for idx, label in colours[['SVTYPE', 'colour']].drop_duplicates().iterrows():
            ax.bar(0, 0, color=label['colour'], alpha=0.8,
                                    label=label['SVTYPE'], linewidth=0)
            lgd = ax.legend(loc='center right', bbox_to_anchor=(1.45,0.6), title='SVTYPE')

        fig.set_size_inches(7 , (node_count * 0.7))
        fig.savefig(output_file, bbox_inches='tight', bbox_extra_artist=(lgd,), format='png', dpi=300)

