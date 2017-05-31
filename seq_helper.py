import json  # for saving and reading
from collections import defaultdict
from itertools import islice
from sklearn.metrics import confusion_matrix
import numpy as np

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