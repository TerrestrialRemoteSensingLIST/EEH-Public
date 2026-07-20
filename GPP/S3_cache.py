# Create 

import numpy as np
import os
import re

def cache_S3_with_pattern(name,mount_folder,pattern='_(\\d{8}).*_'):
    pid = str(os.getpid())
    map_stuff = dict()
    all_stuff = np.loadtxt(mount_folder+'utils/'+name+'_S3_cached.txt', dtype=bytes, ndmin=1).astype(str)
    for thing in all_stuff:
      m=re.search(pattern,thing)
      if not m is None:
        key = m.group(1)
        map_stuff[key] = thing
    return map_stuff

def cache_S3_simple(name, mount_folder, real_folder):
    pid = str(os.getpid())
    # Charger le fichier en array NumPy
    all_stuff = np.loadtxt(mount_folder + 'utils/' + name + '_S3_cached.txt',
                           dtype=bytes, ndmin=1).astype(str)

    # Convertir en liste et préfixer chaque ligne par real_folder
    all_stuff_prefixed = [real_folder + line for line in all_stuff.tolist()]

    return all_stuff_prefixed
