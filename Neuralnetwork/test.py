#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import numpy as np
import random
import os
import sys

import make_my_data as md


def placeholder(size_x, size_y):
    x = tf.placeholder(tf.float32, [size_x, None], name="X")
    y = tf.placeholder(tf.float32, [size_y, None], name="Y")
    return x, y

# initialize w, b
def initialize_parameters(layer_dims):
    tf.set_random_seed(1)
    parameters = {}
    n_layer = len(layer_dims)
    for i in range(1, n_layer):
        parameters['w'+str(i)] = tf.get_variable('w'+str(i), [layer_dims[i], layer_dims[i-1]], initializer = tf.contrib.layers.xavier_initializer())
        parameters['b'+str(i)] = tf.get_variable('b'+str(i), [layer_dims[i], 1], initializer = tf.zeros_initializer())
    return parameters

def forward(input_0, parameters):
    n_layer = len(parameters) / 2
    for i in range(1, n_layer+1):
        if i==1:
            IN = input_0
        else:
            IN = y # y_prev
        w = parameters['w'+str(i)]
        b = parameters['b'+str(i)]
        A = tf.add(tf.matmul(w,IN), b)
        if i!=n_layer:
            y = tf.nn.relu(A)
        else:
            return A

def costCAL(A, label):
    logits = tf.transpose(A)
    labels = tf.transpose(label)
    costFun = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=logits, labels=labels))
    return costFun
    

def model(x_train, y_train, x_test, y_test, layers_dims, learning_rate, num_iterations, print_cost):
    tf.set_random_seed(1)
    # seed = 1
    costs = []
    x, y = placeholder(n_point*n_ver, n_point*n_ver)
    parameters = initialize_parameters(layers_dims)
    A = forward(x, parameters)
    costFun = costCAL(A, y)
    optimizer = tf.AdamOptimizer(learning_rate = learning_rate).minimize(costFun)
    init = tf.global_variables_initializer()
    with tf.Session() as sess:
        sess.run(init)
        for i in range(num_iterations):
            _, cost = sess.run([optimizer, costFun], feed_dict = {X: x_train, Y: y_train})
            costs.append(cost)
            if print_cost==True and i%100==0:
                print("cost after %i interations is %f" % (i, cost))
        parameters = sess.run(parameters)
        print("parameters have been trained.")
        correct_prediction = tf.equal(tf.argmax(A), tf.argmax(y))
        accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
        print("Train Accuracy: ", accuracy.eval({X: x_train, Y: y_train}))
        print("Test Accuracy: ", accuracy.eval({X: x_test, Y: y_test}))
    return parameters

# test


X, Y = md.import_data()

(n_step, n_point, n_ver) = X.shape

n_test = int(n_step/5)

X_TRAIN = X[0:n_step-n_test, :]
Y_TRAIN = Y[0:n_step-n_test, :]
X_TEST = X[n_step-n_test:n_step, :]
Y_TEST = Y[n_step-n_test:n_step, :]
        
parameters = model(X_TRAIN, Y_TRAIN, X_TEST, Y_TEST, [[10,3],[10, 10], [10,3]], 0.001, 1000, True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
