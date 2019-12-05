#!/usr/bin/python3


import argparse

parser = argparse.ArgumentParser(description='Generate a walking pattern in PyBullet')
parser.add_argument('n_steps', type=int, help='Number of steps to be taken')
parser.add_argument('--render', default=1, help='create the pybullet server in GUI mode.')
parser.add_argument('--test', default=1, help='WIP record a value along the walk and print a graph at the end (default: saggital hip torque)')
parser.add_argument('--save_csv', default=0, help='save the gait as a csv containing the leg joint values trajectories')


args = parser.parse_args()

print("args:", args)
