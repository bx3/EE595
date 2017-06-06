import json  # for saving and reading
from collections import defaultdict
from itertools import islice
from sklearn.metrics import confusion_matrix
import numpy as np
from sklearn.cluster import KMeans
from matplotlib.cm import spectral
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import itertools as it
from numpy.linalg import norm

# for now just detect and remove the last cartidge return 
def pre_process_line(line):
    if line[-1] is '\n':
        return line[:-1]
    return line    
    
# read
def write_to_file(data_file, data):
    with open(data_file, 'w') as file:
         file.write(json.dumps(data))

#write            
def read_from_file(data_file):
    with open(data_file) as file:    
        return json.load(file)
            
def get_4_mer_stat(code_hist, gen):
    for four_mer in gen:        
        code_hist[four_mer] += 1  
        
def count_list(seq, Q_hist):
    for dna in list(seq):
        Q_hist[dna] += 1    
        
def window(seq, n=4):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    result = ''.join(result)
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + elem
        yield result
        
# take two int list
def evaluate_result(label_test_l, test_pred_l):
    cm = confusion_matrix(label_test_l, test_pred_l)
    accuracy = (cm[0,0]+cm[1,1])/np.sum(cm)
    precision = cm[0,0]/(cm[0,0]+cm[1,0])
    recall = cm[0,0]/(cm[0,0]+cm[0,1])
    print('accuracy: ',accuracy)
    print('precision: ',precision)
    print('recall: ',recall )
    print(cm)    

def compute_lumping_error(class_list, cluster_list):    
    lumping_err = 0
    list_num = len(class_list)
    for i in range(list_num):
        the_class = class_list[i]
        the_cluster = cluster_list[i]
        
        for j in range(i, list_num):
            other_class = class_list[j]
            other_cluster = cluster_list[j]            
            if the_class!=other_class and the_cluster==other_cluster:
                lumping_err += 1
    return lumping_err       

def compute_splitting_error(class_list, cluster_list):    
    splitting_err = 0
    list_num = len(class_list)
    for i in range(list_num):
        the_class = class_list[i]
        the_cluster = cluster_list[i]
        
        for j in range(i, list_num):
            other_class = class_list[j]
            other_cluster = cluster_list[j]            
            if the_class==other_class and the_cluster!=other_cluster:
                splitting_err += 1
    return splitting_err       

def _joint_probabilities_constant_sigma(D, sigma):
    P = np.exp(-D**2/2 * sigma**2)
    P /= np.sum(P, axis=1)
    return P    
    
    
def classify_with_Kmeans(F_reduced, L, init='k-means++', is_plot=True):
    # Kmeans classification
    nc = 2
    data_num = len(L)

    kmeans =KMeans(nc,init=init).fit(F_reduced)
    
    evaluate_performance(L, kmeans.labels_)
    if is_plot:
        #for i in label draw
        fig, ax = plt.subplots()
        c=[spectral(float(i) /nc) for i in kmeans.labels_]
        ax.scatter(F_reduced[:,0],F_reduced[:,1], c=c)
        txt = [str(i) for i in L] #[label_text_dict[i] for i in label]
        for i in range(data_num):
            ax.annotate(txt[i], (F_reduced[i,0],F_reduced[i,1]),size=10)

        plt.show()
    return kmeans.labels_
    
# PCA
def PCA_project(F,L,nc=2):
    x = F
    mean_x = np.mean(x, axis=0)
    std_x = np.std(x, axis=0)
    x_standard = (x-mean_x)/std_x
    pca = PCA(n_components=nc).fit(x_standard)
    projected_x = pca.transform(x_standard)    

    #for i in label draw
    fig, ax = plt.subplots()

    c=[spectral(float(i) /2) for i in L]
    ax.scatter(projected_x[:,0],projected_x[:,1], c=c)

    plt.show()

def convert_Y_to_Q_by_t_dist(Y):
    n = Y.shape[0]    
    Q_map = np.zeros((n,n))
    den = np.sum(np.array([1/(norm(y - x)**2+1) for x, y in it.combinations(Y, 2)])*2)    
    for i in range(n):
        for j in range(n):
            if i !=j:
                num = 1/(1+norm(Y[i,:]-Y[j,:])**2)                
                Q_map[i,j] = num/den
    return Q_map                    
    
def get_N_neighbor(Q_map, N):
    n = Q_map.shape[0]
    Q_trun_map = np.zeros_like(Q_map)
    neighbor_index = []
    for i in range(n):
        index = (-Q_map[i,:]).argsort()[:N]        
        Q_trun_map[i, index] = Q_map[i,index]
        neighbor_index.append(index)
    return Q_trun_map, neighbor_index
    
    
def get_prob_from_prob_bin(x, y, pmf, x_limits, y_limits):
    x_ind = np.nonzero(1- (x >= x_limits))[0]
    y_ind = np.nonzero(1- (y >= y_limits))[0]      
    return pmf[x_ind[0]-1,y_ind[0]-1]    
    
def evaluate_performance(L,A):
    A_f = np.array(A)==0
    data_num = len(L)
    lumping_error = compute_lumping_error(L, A)
    f_lumping_error = compute_lumping_error(L, A_f)
    lumping_error = min(lumping_error, f_lumping_error)            

    splitting_error = compute_splitting_error(L, A)
    f_splitting_error = compute_splitting_error(L, A_f)
    splitting_error = min(splitting_error, f_splitting_error)    
    
    accuracy = np.sum(L==A)/data_num
    f_accuracy = np.sum(L==A_f)/data_num
    accuracy = max(f_accuracy, accuracy)    
    return accuracy, lumping_error, splitting_error