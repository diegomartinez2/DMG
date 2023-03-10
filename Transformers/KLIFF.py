#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#git clone https://github.com/openkim/kliff
#pip install ./kliff
from kliff import nn
from kliff.calculators import CalculatorTorch
from kliff.descriptors import SymmetryFunction
from kliff.dataset import Dataset
from kliff.models import NeuralNetwork
from kliff.loss import Loss
from kliff.utils import download_dataset

# Descriptor to featurize atomic configurations
descriptor = SymmetryFunction(
    cut_name="cos", cut_dists={"Si-Si": 5.0}, hyperparams="set51", normalize=True
)

# Fully-connected neural network model with 2 hidden layers, each with 10 units
N1 = 10
N2 = 10
model = NeuralNetwork(descriptor)
model.add_layers(
    # first hidden layer
    nn.Linear(descriptor.get_size(), N1),
    nn.Tanh(),
    # second hidden layer
    nn.Linear(N1, N2),
    nn.Tanh(),
    # output layer
    nn.Linear(N2, 1),
)

# Training set (dataset will be downloaded from:
# https://github.com/openkim/kliff/blob/master/examples/Si_training_set.tar.gz)
dataset_path = download_dataset(dataset_name="Si_training_set")
dataset_path = dataset_path.joinpath("varying_alat")
train_set = Dataset(dataset_path)
configs = train_set.get_configs()

# Set up calculator to compute energy and forces for atomic configurations in the
# training set using the neural network model
calc = CalculatorTorch(model, gpu=False)
calc.create(configs)

# Define a loss function and train the model by minimizing the loss
loss = Loss(calc)
result = loss.minimize(method="Adam", num_epochs=10, batch_size=100, lr=0.001)

# Write trained model as a KIM model to be used in other codes such as LAMMPS and ASE
model.write_kim_model()
