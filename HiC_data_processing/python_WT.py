import numpy as np
import os
import sys
import math
import random
import time
import itertools as it
import multiprocessing as mp
import gzip
import glob
import pickle
import cooler
import iced

from scipy.spatial import distance
from scipy.linalg import eig
from scipy import stats
import sklearn as sk
import sklearn.cluster as skclust
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score

os.chdir("/trans-C/code")

from data_processing import *
from rand_walk import *
from base import *

#import files
chrom_sizes_file = "/home/shared_folder/genome_files/chrom_hg38.sizes"
fpath = "/home/shared_folder/transC"
res = 1000000
seed_file = "/home/shared_folder/seeds.bed"

from base import *
process_inputs(chrom_sizes_file, fpath, res)
bins = list(bin_map.keys())

seed = read_seed(seed_file)

WT = np.load(f"{fapth}/output_WT/im_matrix_WT.npy")

n_seed = 1
n = 0
lista = []
rep = 1000

while n < rep:
  rand_clique = random.sample(bins, n_seed)
  weight_r = find_clique_weight(rand_clique, WT)
  if weight_r >= weight_WT:
    n += 1
    print(n)
    try:
      pi = do_random_walk(WT, 0.5, rand_clique, path)
      dictionary = dict(zip(bins, pi))
      sorted_dict = sorted(dictionary.items(), key=lambda x:x[1])
      sorted_dict = dict(sorted_dict)
      lista_keys = list(sorted_dict.keys())
      lista_keys.reverse()
      lista.append(lista_keys)
    except:
      n -= 1
      print("eigenvector calculation error")

table = np.array(lista)
DF = pd.DataFrame(table)
DF.to_csv(f"{fpath}/output_WT/table_random_rankings_TTN.csv")




