# test latent

#
# import bnlearn
#
# df = bnlearn.import_example()
#
# model = bnlearn.structure_learning.fit(df)
#
# G = bnlearn.plot(model)

from pgmpy.factors.discrete import State
from pgmpy.models.BayesianModel import BayesianModel
from pgmpy.sampling import BayesianModelSampling
from pgmpy.factors.discrete import TabularCPD
#from pgmpy.estimators import K2Score, BdeuScore, BicScore
import pandas as pd
import numpy as np
from scipy.special import comb
from utils import compute_local_BiCScore, simulate_sem, compute_local_BiCScore_gauss, \
    compute_BiC_c_component, score_write_to_file, compute_BiC_c_comps_sets, check_connected_component,\
    find_bi_connected_node, has_cycle_in_c_comp
import networkx as nx
import itertools
import pickle

def build_5_var():
    model = BayesianModel([('D', 'A'), ('C', 'A'), ('C','B'), ('E','B')
                           ])
    cpd_d = TabularCPD('D', 2, [[0.1], [0.9]])
    cpd_e = TabularCPD('E', 2, [[0.15], [0.85]])
    cpd_c = TabularCPD('C', 2, [[0.3], [0.7]])
    cpd_a = TabularCPD('A', 2, [[0.3, 0.25, 0.6, 0.7], [0.7, 0.75, 0.4, 0.3]],
                       ['C', 'D'], [2, 2])

    cpd_b = TabularCPD('B', 2, [[0.6, 0.75, 0.2, 0.7], [0.4, 0.25, 0.8, 0.3]],
                       ['C', 'E'], [2, 2])

    # model = BayesianModel([(4, 1), (3, 1), (3, 2), (5, 2)
    #                        ])
    # cpd_d = TabularCPD(4, 2, [[0.8], [0.2]])
    # cpd_e = TabularCPD(5, 2, [[0.15], [0.85]])
    # cpd_c = TabularCPD(3, 2, [[0.7], [0.3]])
    # cpd_a = TabularCPD(1, 2, [[0.3, 0.25, 0.9, 0.55], [0.7, 0.75, 0.1, 0.45]],
    #                    [3, 4], [2, 2])
    #
    # cpd_b = TabularCPD(2, 2, [[0.6, 0.75, 0.2, 0.7], [0.4, 0.25, 0.8, 0.3]],
    #                    [3, 5], [2, 2])

    model.add_cpds( cpd_a, cpd_b, cpd_c, cpd_d, cpd_e)
    return model

def main_gaus():

    # W = [[0, 0, 0, 0, 0],
    #      [0, 0, 0, 0, 0],
    #      [1, 1, 0, 0, 0],
    #      [1, 0, 0, 0, 0],
    #      [0, 1, 0, 0, 0]]
    G = nx.DiGraph()
    G.add_edges_from([(3, 0), (2, 0), (2, 1), (4, 1)])

    sample_size = 100000
    data = simulate_sem(G, sample_size, 1, 'linear-gauss')
    data = data[:,:,0]

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [2, 3]) + \
                         compute_local_BiCScore_gauss(data, 1, [2, 4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, []) + \
                         compute_local_BiCScore_gauss(data, 2, [])
    print('Full self check  ground truth DAG score', latent_model_score)

    # latent_model_score = BicScore(values).local_score('A', ['D']) + \
    #                      BicScore(values).local_score('B', ['E']) + \
    #                      BicScore(values).local_score('D', []) + \
    #                      BicScore(values).local_score('E', []) + \
    #                      BicScore(values).local_score('C', [])
    # print('dag  without A-B edge  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [3]) + \
                         compute_local_BiCScore_gauss(data, 1, [4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, []) + \
                         compute_local_BiCScore_gauss(data, 2, [])
    print('Full self check dag  without A-B edge  BIC', latent_model_score)

    # latent_model_score = BicScore(values).local_score('A', ['D']) + \
    #                      BicScore(values).local_score('B', []) + \
    #                      BicScore(values).local_score('D', []) + \
    #                      BicScore(values).local_score('E', []) + \
    #                      BicScore(values).local_score('C',[])
    # print(' Full DAG  without A-B edge  and E-B BIC ', latent_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [3]) + \
                         compute_local_BiCScore_gauss(data, 1, []) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, []) + \
                         compute_local_BiCScore_gauss(data, 2, [])
    print('Full self check dag with D-A only BIC', latent_model_score)

    # latent_model_score = BicScore(values).local_score('A', []) + \
    #                      BicScore(values).local_score('B', ['E']) + \
    #                      BicScore(values).local_score('D', []) + \
    #                      BicScore(values).local_score('E', []) + \
    #                      BicScore(values).local_score('C', [])
    # print(' Full DAG  without A-B edge  and D-A BIC ', latent_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, []) + \
                         compute_local_BiCScore_gauss(data, 1, [4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, []) + \
                         compute_local_BiCScore_gauss(data, 2, [])
    print('Full self check dag with E-B only BIC', latent_model_score)

    # latent_model_score = BicScore(values).local_score('A', []) + \
    #                      BicScore(values).local_score('B', []) + \
    #                      BicScore(values).local_score('D', []) + \
    #                      BicScore(values).local_score('E', []) + \
    #                      BicScore(values).local_score('C', [])
    # print(' Full empty graph  WITHOUT latent BIC D', latent_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, []) + \
                         compute_local_BiCScore_gauss(data, 1, []) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, []) + \
                         compute_local_BiCScore_gauss(data, 2, [])
    print('Full self check empty dag  BIC', latent_model_score)

    # 4 variable
    # slice data into observed

    # observed_data = values.loc[:, ['A', 'B', 'D', 'E']]

    # score_function = BicScore(observed_data)

    # latent_model_score = score_function.local_score('A', ['B', 'D', 'E']) + \
    #                  score_function.local_score('B', ['E']) + \
    #                  score_function.local_score('D', []) + \
    #                  score_function.local_score('E', [])
    # print('DAG without C BIC model A ', latent_model_score)
    #
    # latent_model_score = score_function.local_score('A', ['D']) + \
    #                      score_function.local_score('B', ['A', 'D', 'E']) + \
    #                      score_function.local_score('D', []) + \
    #                      score_function.local_score('E', [])
    # print('DAG  without C BIC model B', latent_model_score)

    # latent_model_score =  BicScore(values).local_score('A', ['C','D']) + \
    #                       BicScore(values).local_score('B', ['C','E']) + \
    #                       BicScore(values).local_score('D', []) + \
    #                       BicScore(values).local_score('E', [])
    # print('DAG  without A-B, given C,  edge  BIC', latent_model_score)

    # -------- self compute check : compute_local_BiCScore

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [1, 3, 4]) + \
                         compute_local_BiCScore_gauss(data, 1, [4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, [])
    print('4 self check DAG  best graph 1  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [1, 3, 4]) + \
                         compute_local_BiCScore_gauss(data, 1, [0, 4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, [])
    print('4 self check DAG  best graph 1, +bidirected A-B,  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [3]) + \
                         compute_local_BiCScore_gauss(data, 1, [0, 3, 4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, [])
    print('4 self check DAG  best graph 2  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [1, 3]) + \
                         compute_local_BiCScore_gauss(data, 1, [0, 3, 4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, [])
    print('4 self check DAG  best graph 2, +bidirected A-B,  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [3]) + \
                         compute_local_BiCScore_gauss(data, 1, [4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, [])
    print('4 self check DAG  without A-B, given C,  edge  BIC', latent_model_score)

    # bidrected model score
    # bi_model_score = score_function.local_score('A', ['B', 'D']) +\
    #                     score_function.local_score('B', ['A', 'E']) +\
    #                     score_function.local_score('D', []) +\
    #                     score_function.local_score('E', [])
    # print('bidrected sparse  graph BIC', bi_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [1, 3]) + \
                         compute_local_BiCScore_gauss(data, 1, [0, 4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, [])
    print('4 self sparse bidrected   BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore_gauss(data, 0, [1, 3, 4]) + \
                         compute_local_BiCScore_gauss(data, 1, [0, 3, 4]) + \
                         compute_local_BiCScore_gauss(data, 3, []) + \
                         compute_local_BiCScore_gauss(data, 4, [])
    print('4 self dense bidrected   BIC', latent_model_score)

    # bidrected model score
    # bi_model_score = score_function.local_score('A', ['B', 'D','E']) + \
    #                  score_function.local_score('B', ['A', 'E', 'D']) + \
    #                  score_function.local_score('D', []) + \
    #                  score_function.local_score('E', [])
    # print('bidrected dense  graph BIC', bi_model_score)
    #
    # # bidrected model score
    # bi_model_score = score_function.local_score('A', ['B', 'D']) + \
    #                  score_function.local_score('B', ['A', 'E']) + \
    #                  score_function.local_score('D', []) + \
    #                  score_function.local_score('E', [])
    # print('bidrected dense  graph BIC', bi_model_score)

    # bic_local_score_a_bidirected =  score_function.local_score('A', ['B', 'D'])
    bic_local_score_a_bidirected = compute_local_BiCScore_gauss(data, 0, [1, 3])

    # bic_local_score_a_dense = score_function.local_score('A', ['B', 'D', 'E'])
    bic_local_score_a_dense = compute_local_BiCScore_gauss(data, 0, [1, 3, 4])

    # bic_local_score_b_bidirected = score_function.local_score('B', ['A', 'E'])
    bic_local_score_b_bidirected = compute_local_BiCScore_gauss(data, 1, [0, 4])
    # bic_local_score_b_dense = score_function.local_score('B', ['A', 'D', 'E'])
    bic_local_score_b_dense = compute_local_BiCScore_gauss(data, 1, [0, 3, 4])

    print('bidirect edge with sparse connection:',
          bic_local_score_a_bidirected,
          bic_local_score_b_bidirected,
          bic_local_score_a_bidirected + bic_local_score_b_bidirected)

    print('bidirect edge with dense connection :',
          bic_local_score_a_dense,
          bic_local_score_b_dense,
          bic_local_score_a_dense + bic_local_score_b_dense)
    # print(data)

def main():
    dataset = '5_gauss'
    # dataset = 'data_1000_0.1_3_10_1.txt'
    sample_size = 1000
    file_name = 'synth_5_test_score_' + str(sample_size) + '.pkl'


    num_var = 5
    if dataset == '5_discrete':
        model = build_5_var()
        inference = BayesianModelSampling(model)


        data = inference.forward_sample(size=sample_size, return_type='recarray')

        values = pd.DataFrame(data)
        data = values.to_numpy()
        is_discrete = True

        # -------- self compute check : compute_local_BiCScore
        check_subgraph_scores(values)

    elif dataset == '5_gauss':
        W = [[0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [1, 1, 0, 0, 0],
             [1, 0, 0, 0, 0],
             [0, 1, 0, 0, 0]]
        G = nx.DiGraph(np.array(W))
        # sample_size = 1000
        data = simulate_sem(G, sample_size, 1, 'linear-gauss')
        is_discrete = False
        data = data[:,:,0]

        # -------- self compute check : compute_local_BiCScore
        # check_subgraph_scores(data)

        # visible
        observed_data = data[:, [0, 1, 3, 4]]

    elif dataset == '5_exist':
        with open(file_name, 'rb') as f:
            data, old_scores = pickle.load(f)

        # visible
        observed_data = data

    elif dataset.endswith('.txt'):
        with open(dataset, 'rb') as f:
            data = np.loadtxt(f, skiprows=0)

        # visible
        observed_data = data

        # file_name
        file_name = 'score_' + dataset[:-4]

    # k2_score = K2Score(values).score(model)
    # bic_score =  BicScore(values).score(model)
    # print('Full ground truth BIC', bic_score)

    # check current learner results
    # latent_model_score = BicScore(values).local_score('A', ['C', 'D']) + \
    #                      BicScore(values).local_score('B', ['C', 'E']) + \
    #                      BicScore(values).local_score('D', []) + \
    #                      BicScore(values).local_score('E', []) + \
    #                      BicScore(values).local_score('C', [])
    # print('Self compute ground truth DAG   BIC', latent_model_score)
    c_size = 5
    single_parent_size = 3
    other_c_parent_size = 2

    generate_scores_bidirect(observed_data,
                             single_c_parent_size = single_parent_size,
                             other_c_parent_size = other_c_parent_size,
                             c_size = c_size,
                             file_name = file_name)


def generate_scores_bidirect(data, single_c_parent_size, other_c_parent_size, c_size, file_name):
    # generate data for all bidirected scores
    # TODO change c-comp parent size
    num_sample, num_var = data.shape

    # # generate all possible subgraphs, up to size k
    # graphs = gererate_all_bidirected_graphs(num_var, True)

    # compute scores for each graph
    scores = {}
    # c_comps = {}
    # for iVar in range(num_var):
    all_c_components = gererate_all_c_components(num_var, c_size)

    for iComp_size in all_c_components:
        c_comps_of_size = all_c_components[iComp_size]
        parent_size = single_c_parent_size if iComp_size < 2 else other_c_parent_size
        for each_comp in c_comps_of_size: # in a list

            edges = gererate_all_bidirected_edges(each_comp)
            for each_edge_set in edges:

                if iComp_size > 1:
                    # scores[each_comp]['bi_edges'] = [each_comp]
                    # dic_key = (each_comp, each_edge_set)
                    save_edge_set = each_edge_set
                else:
                    # scores[each_comp]['bi_edges'] = []
                    # dic_key = (each_comp, ())
                    save_edge_set = ()

                dic_key = (each_comp, save_edge_set)
                scores[dic_key] = {}

                parent_set_list = []
                for each_var_in_comp in each_comp:
                    # generate parental size
                    parent_set_list.append( gererate_all_parent_sets(each_var_in_comp,
                                                                     each_comp,
                                                                     num_var,
                                                                     parent_size,
                                                                     save_edge_set) )

                # generate parent set combination
                all_parent_list = list(itertools.product(*parent_set_list))

                # compute scores
                for each_parent_set_config in all_parent_list:

                    # scores[dic_key][each_parent_set_config] = {}

                    # remove certain condition:
                    # 1) a->b and a<-b
                    if not has_cycle_in_c_comp(each_parent_set_config, each_comp):
                    # list all potential bidirected edges
                        scores[dic_key][each_parent_set_config] = compute_BiC_c_comps_sets(data,
                                                                                            each_comp,
                                                                                            each_parent_set_config,
                                                                                           each_edge_set)

    # save files
    print('saving file..\n')
    file_name = './' + file_name # /synth_5_test_score.pkl'
    score_write_to_file([data, scores], file_name)

    return

def generate_scores_bidirect_m3hc(data, single_c_parent_size, other_c_parent_size, c_size, file_name,instance_name):
    # generate data for all bidirected scores
    # TODO change c-comp parent size
    num_sample, num_var = data.shape

    # # generate all possible subgraphs, up to size k
    # graphs = gererate_all_bidirected_graphs(num_var, True)

    # compute scores for each graph
    scores = {}
    # c_comps = {}
    # for iVar in range(num_var):
    all_c_components = gererate_all_c_components(num_var, c_size)

    for iComp_size in all_c_components:
        c_comps_of_size = all_c_components[iComp_size]
        parent_size = single_c_parent_size if iComp_size < 2 else other_c_parent_size
        for each_comp in c_comps_of_size: # in a list

            edges = gererate_all_bidirected_edges(each_comp)
            for each_edge_set in edges:

                if iComp_size > 1:
                    # scores[each_comp]['bi_edges'] = [each_comp]
                    # dic_key = (each_comp, each_edge_set)
                    save_edge_set = each_edge_set
                else:
                    # scores[each_comp]['bi_edges'] = []
                    # dic_key = (each_comp, ())
                    save_edge_set = ()

                dic_key = (each_comp, save_edge_set)
                scores[dic_key] = {}

                parent_set_list = []
                for each_var_in_comp in each_comp:
                    # generate parental size
                    parent_set_list.append( gererate_all_parent_sets(each_var_in_comp,
                                                                     each_comp,
                                                                     num_var,
                                                                     parent_size,
                                                                     save_edge_set) )

                # generate parent set combination
                all_parent_list = list(itertools.product(*parent_set_list))

                # compute scores
                for each_parent_set_config in all_parent_list:

                    # scores[dic_key][each_parent_set_config] = {}

                    # remove certain condition:
                    # 1) a->b and a<-b
                    if not has_cycle_in_c_comp(each_parent_set_config, each_comp):
                    # list all potential bidirected edges
                        scores[dic_key][each_parent_set_config] = compute_BiC_c_comps_sets(data,
                                                                                            each_comp,
                                                                                            each_parent_set_config,
                                                                                           each_edge_set)
    matrix = open('../Instances/m3hcAG/m3hc_'+instance_name).read()
    matrix = [item.split() for item in matrix.split('\n')[:-1]]
    dim = len(matrix)
    bi = []
    di = []
    for i in range(dim):
        for j in range(dim):
            if i<j and matrix[i][j] == '2' and matrix[j][i] == '2':
                bi.append((i,j))
            if matrix[i][j] == '2' and matrix[j][i] == '3':
                di.append((i,j))
    Nodes = range(dim)
    cNodes = []
    cCompEdges = []
    cParents = []
    inComp = [0]*dim
    for biedge in bi:
        if inComp[biedge[0]] == 0 and inComp[biedge[1]] == 0:
            inComp[biedge[0]] = 1
            inComp[biedge[1]] = 1
            cNodes.append([biedge[0],biedge[1]])
            cCompEdges.append([(biedge[0],biedge[1])])
        elif inComp[biedge[0]] == 1:
            if inComp[biedge[1]] == 0:
                inComp[biedge[1]] = 1
                isComp = None
                for k in range(len(cNodes)):
                    if biedge[0] in cNodes[k]:
                        isComp = k
                        break
                if biedge[1] not in cNodes[isComp]:
                    cNodes[isComp].append(biedge[1])
                cCompEdges[isComp].append((biedge[0],biedge[1]))
            else:
                isComp = None
                for k in range(len(cNodes)):
                    if biedge[0] in cNodes[k]:
                        isComp = k
                        break
                cCompEdges[isComp].append((biedge[0],biedge[1]))
                jsComp = None
                for k in range(len(cNodes)):
                    if biedge[1] in cNodes[k]:
                        jsComp = k
                        break
                if isComp != jsComp:
                    for node in cNodes[jsComp]:
                        if node not in cNodes[isComp]:
                            cNodes[isComp].append(node)
                    for biedge in cCompEdges[jsComp]:
                        if biedge not in cCompEdges[isComp]:
                            cCompEdges[isComp].append(biedge)
                    cNodes.remove(cNodes[jsComp])
                    cCompEdges.remove(cCompEdges[jsComp])
        elif inComp[biedge[1]] == 1:
            inComp[biedge[0]] = 1
            for k in range(len(cNodes)):
                if biedge[1] in cNodes[k]:
                    isComp = k
                    break
            if biedge[0] not in cNodes[isComp]:
                cNodes[isComp].append(biedge[0])
            cCompEdges[isComp].append((biedge[0],biedge[1]))
    for i in range(dim):
        if inComp[i] == 0:
            cNodes.append([i])
            cCompEdges.append([])
    for k in range(len(cNodes)):
        cNodes[k] = tuple(sorted(cNodes[k]))
        cCompEdges[k] = tuple(cCompEdges[k])
        cParents.append([])
        for i in cNodes[k]:
            ciParents = []
            for j in range(dim):
                if (j,i) in di:
                    ciParents.append(j)
            cParents[k].append(tuple(ciParents.copy()))
        cParents[k] = tuple(cParents[k])
        dic_key = (cNodes[k],cCompEdges[k])
        each_parent_set_config = cParents[k]
        each_comp = cNodes[k]
        each_edge_set = cCompEdges[k]
        if dic_key not in scores.keys():
            scores[dic_key] = {}
        scores[dic_key][each_parent_set_config] = compute_BiC_c_comps_sets(data,each_comp,each_parent_set_config,each_edge_set)
        #print(str(dic_key)+","+str(each_parent_set_config)+","+str(scores[dic_key][each_parent_set_config]))
        

    # save files
    print('saving file..\n')
    file_name = './' + file_name # /synth_5_test_score.pkl'
    score_write_to_file([data, scores], file_name)

    return

def gererate_all_bidirected_edges(each_comp):
    # generate all possible bi-directed edges within a c component
    # edge_num = len(each_comp)   # max number of edges
    edge_num = len(each_comp) - 1
    bi_edges = []
    edge_num_upper_bound = comb(len(each_comp), 2, exact=True) # scipy.misc.comb

    if edge_num < 1:
        flatten_list = [-1]  # no bidrected edges

    else:
        all_edges = list(itertools.combinations( sorted(set(each_comp)), 2))
        # for iCsize in range(1, edge_num+1): # size 0 to c_size
        for iCsize in range(edge_num, edge_num_upper_bound+1):  # size edge_num to max
            # parent_candidates =  list(set(list(range(num_var))) - set(var_in_comp))
            # parent_candidates.remove(iVar)
            sets = list(itertools.combinations( set(all_edges), iCsize))
            # c_comp[iCsize] = sets

            # check connected components, make sure nodes are connected to each other in c-comp
            sets = check_connected_component(sets, each_comp)
            bi_edges.append(sets)

        flatten_list = [item for sublist in bi_edges for item in sublist]
    return flatten_list

def gererate_all_parent_sets(iVar,
                             var_in_comp,
                             num_var,
                             c_size,
                             comp_edges):
    # generate all possible c components
    # comp_edges: edges in a comp

    # bi_directed_edge: bi directed edges from iVar
    bi_directed_node = find_bi_connected_node(iVar, comp_edges)
    candidate_parent_in_comp =  set(list(var_in_comp)) - set([iVar])- \
                                set(list(bi_directed_node))
    c_comp = []
    for iCsize in range(c_size+1): # size 0 to c_size

        # add candiate parental sets
        parent_candidates =  list(sorted(
            (set(list(range(num_var))) - set(var_in_comp)).union(candidate_parent_in_comp)))
        # parent_candidates.remove(iVar)
        sets = list(itertools.combinations( sorted(set(parent_candidates)), iCsize))
        # c_comp[iCsize] = sets
        c_comp.append(sets)

    flatten_list = [item for sublist in c_comp for item in sublist]
    return flatten_list


def gererate_all_c_components(num_var, c_size):
    # generate all possible c components
    c_comp = {}
    for iCsize in range(1, c_size+1):
        sets = list(itertools.combinations( sorted(set(range(num_var))), iCsize))
        c_comp[iCsize] = sets
    return c_comp


def gererate_all_bidirected_graphs(num_var, bi_directed = True):
    # generate all possible subgraphs with bi-directed edges
    # TODO

    W = [[0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [1, 1, 0, 0, 0],
         [1, 0, 0, 0, 0],
         [0, 1, 0, 0, 0]]

    graph = {}
    graph[0] = W
    return graph


def check_subgraph_scores(values):
    data = values.to_numpy()

    latent_model_score = compute_local_BiCScore(data, 0, [2, 3]) + \
                         compute_local_BiCScore(data, 1, [2, 4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, []) + \
                         compute_local_BiCScore(data, 2, [])
    print('Full self check  ground truth DAG score', latent_model_score)


    # latent_model_score = BicScore(values).local_score('A', ['D']) + \
    #                      BicScore(values).local_score('B', ['E']) + \
    #                      BicScore(values).local_score('D', []) + \
    #                      BicScore(values).local_score('E', []) + \
    #                      BicScore(values).local_score('C', [])
    # print('dag  without A-B edge  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, [ 3]) + \
                         compute_local_BiCScore(data, 1, [4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, []) + \
                         compute_local_BiCScore(data, 2, [])
    print('Full self check dag  without A-B edge  BIC', latent_model_score)

    # latent_model_score = BicScore(values).local_score('A', ['D']) + \
    #                      BicScore(values).local_score('B', []) + \
    #                      BicScore(values).local_score('D', []) + \
    #                      BicScore(values).local_score('E', []) + \
    #                      BicScore(values).local_score('C',[])
    # print(' Full DAG  without A-B edge  and E-B BIC ', latent_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, [3]) + \
                         compute_local_BiCScore(data, 1, []) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, []) + \
                         compute_local_BiCScore(data, 2, [])
    print('Full self check dag with D-A only BIC', latent_model_score)

    # latent_model_score = BicScore(values).local_score('A', []) + \
    #                      BicScore(values).local_score('B', ['E']) + \
    #                      BicScore(values).local_score('D', []) + \
    #                      BicScore(values).local_score('E', []) + \
    #                      BicScore(values).local_score('C', [])
    # print(' Full DAG  without A-B edge  and D-A BIC ', latent_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, []) + \
                         compute_local_BiCScore(data, 1, [4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, []) + \
                         compute_local_BiCScore(data, 2, [])
    print('Full self check dag with E-B only BIC', latent_model_score)

    # latent_model_score = BicScore(values).local_score('A', []) + \
    #                      BicScore(values).local_score('B', []) + \
    #                      BicScore(values).local_score('D', []) + \
    #                      BicScore(values).local_score('E', []) + \
    #                      BicScore(values).local_score('C', [])
    # print(' Full empty graph  WITHOUT latent BIC D', latent_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, []) + \
                         compute_local_BiCScore(data, 1, []) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, []) + \
                         compute_local_BiCScore(data, 2, [])
    print('Full self check empty dag  BIC', latent_model_score)


    # 4 variable
    # slice data into observed

    # observed_data = values.loc[:, ['A','B','D','E']]

    # score_function = BicScore(observed_data)


    # latent_model_score = score_function.local_score('A', ['B', 'D', 'E']) + \
    #                  score_function.local_score('B', ['E']) + \
    #                  score_function.local_score('D', []) + \
    #                  score_function.local_score('E', [])
    # print('DAG without C BIC model A ', latent_model_score)
    #
    # latent_model_score = score_function.local_score('A', ['D']) + \
    #                      score_function.local_score('B', ['A', 'D', 'E']) + \
    #                      score_function.local_score('D', []) + \
    #                      score_function.local_score('E', [])
    # print('DAG  without C BIC model B', latent_model_score)

    # latent_model_score =  BicScore(values).local_score('A', ['C','D']) + \
    #                       BicScore(values).local_score('B', ['C','E']) + \
    #                       BicScore(values).local_score('D', []) + \
    #                       BicScore(values).local_score('E', [])
    # print('DAG  without A-B, given C,  edge  BIC', latent_model_score)

    #-------- self compute check : compute_local_BiCScore

    latent_model_score = compute_local_BiCScore(data, 0, [1, 3, 4]) + \
                         compute_local_BiCScore(data, 1, [4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, [])
    print('4 self check DAG  best graph 1  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, [1, 3, 4]) + \
                         compute_local_BiCScore(data, 1, [0, 4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, [])
    print('4 self check DAG  best graph 1, +bidirected A-B,  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, [3]) + \
                         compute_local_BiCScore(data, 1, [0, 3, 4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, [])
    print('4 self check DAG  best graph 2  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, [1, 3]) + \
                         compute_local_BiCScore(data, 1, [0, 3, 4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, [])
    print('4 self check DAG  best graph 2, +bidirected A-B,  BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, [ 3]) + \
                         compute_local_BiCScore(data, 1, [4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, [])
    print('4 self check DAG  without A-B, given C,  edge  BIC', latent_model_score)

    # bidrected model score
    # bi_model_score = score_function.local_score('A', ['B', 'D']) +\
    #                     score_function.local_score('B', ['A', 'E']) +\
    #                     score_function.local_score('D', []) +\
    #                     score_function.local_score('E', [])
    # print('bidrected sparse  graph BIC', bi_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, [1, 3]) + \
                         compute_local_BiCScore(data, 1, [0, 4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, [])
    print('4 self sparse bidrected   BIC', latent_model_score)

    latent_model_score = compute_local_BiCScore(data, 0, [1, 3, 4]) + \
                         compute_local_BiCScore(data, 1, [0, 3, 4]) + \
                         compute_local_BiCScore(data, 3, []) + \
                         compute_local_BiCScore(data, 4, [])
    print('4 self dense bidrected   BIC', latent_model_score)



    # bidrected model score
    # bi_model_score = score_function.local_score('A', ['B', 'D','E']) + \
    #                  score_function.local_score('B', ['A', 'E', 'D']) + \
    #                  score_function.local_score('D', []) + \
    #                  score_function.local_score('E', [])
    # print('bidrected dense  graph BIC', bi_model_score)
    #
    # # bidrected model score
    # bi_model_score = score_function.local_score('A', ['B', 'D']) + \
    #                  score_function.local_score('B', ['A', 'E']) + \
    #                  score_function.local_score('D', []) + \
    #                  score_function.local_score('E', [])
    # print('bidrected dense  graph BIC', bi_model_score)

    # bic_local_score_a_bidirected =  score_function.local_score('A', ['B', 'D'])
    bic_local_score_a_bidirected =  compute_local_BiCScore(data, 0, [1, 3])

    # bic_local_score_a_dense = score_function.local_score('A', ['B', 'D', 'E'])
    bic_local_score_a_dense = compute_local_BiCScore(data, 0, [1, 3, 4])

    # bic_local_score_b_bidirected = score_function.local_score('B', ['A', 'E'])
    bic_local_score_b_bidirected = compute_local_BiCScore(data, 1, [0, 4])
    # bic_local_score_b_dense = score_function.local_score('B', ['A', 'D', 'E'])
    bic_local_score_b_dense = compute_local_BiCScore(data, 1, [0, 3, 4])

    print('bidirect edge with sparse connection:',
          bic_local_score_a_bidirected,
          bic_local_score_b_bidirected,
          bic_local_score_a_bidirected + bic_local_score_b_bidirected)

    print('bidirect edge with dense connection :',
          bic_local_score_a_dense,
          bic_local_score_b_dense,
          bic_local_score_a_dense + bic_local_score_b_dense)
    # print(data)

if __name__ == "__main__":

    main()
    #
    # main_gaus()