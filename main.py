# -*- coding: utf-8 -*-
"""
@author: Tommaso Giacometti
"""
from protein_class import Protein
import random
import time
import configparser
import argparse
import utils
import matplotlib.pyplot as plt
import plots


# Parser to get from terminal the configuration file
parser = argparse.ArgumentParser()
# the filename is optional, the default is the config.txt
parser.add_argument('configuration_file', help='file from which takes the configuration', default = 'config.txt', nargs='?')

args = parser.parse_args() # get the args written in the terminal
filename = args.configuration_file # assign filename

# setting the parameters from configuration file
configuration = configparser.ConfigParser()
configuration.read(filename)

config = utils.Configuration(configuration) # class to get the save the configuration from the file

#Random seed setting
random.seed(config.seed)
print(f'The random seed used is {config.seed}')

prot = Protein(config) # Protein class 

plots.view(protein=prot, tit='Initial configuration') # plot the structure of the protein

start = time.time() 

print('--------------------')
print('Evolution started...')
prot.evolution() # evolve the protein with folds foldings
print('Evolution ended')
print('---------------')

# various plots:
plots.view(protein=prot, save=False, tit='Final configuration')
plots.view_min_en(protein=prot)
plots.view_max_comp(protein=prot)
plots.plot_energy(protein=prot, avg=10)
plots.plot_compactness(protein=prot, avg=10)

print(f'It took {time.time()-start:.3f} seconds')

plt.show() # to let the let plots on screen at the end

if config.gif:
    plots.create_gif(protein=prot)


from mpl_toolkits.mplot3d import Axes3D

def view_3d(protein, save=False, tit='3D Structure'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Assuming the Protein class has these methods to get coordinates
    x = protein.get_x_coordinates()
    y = protein.get_y_coordinates()
    z = protein.get_z_coordinates()

    # Plot the 3D trajectory of the protein
    ax.plot(x, y, z, label='Protein Structure')
    ax.scatter(x, y, z)  # Mark the points for clarity

    # Labeling
    ax.set_title(tit)
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.legend()

    if save:
        plt.savefig(f'{tit}.png')

    plt.show()

plots.view_3d(protein=prot, tit='3D Structure of Protein')