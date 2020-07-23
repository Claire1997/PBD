#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import tensorflow as tf 

def import_data():
    data = np.loadtxt('positions_10_0.txt',dtype='float',delimiter=' ')
    data = np.reshape(data,(1003,30))
    print (data)
    X = data[0:1002, :]
    Y = data[1:1003, :]
    X = tf.convert_to_tensor(X, dtype=tf.float32)
    Y = tf.convert_to_tensor(Y, dtype=tf.float32)
    # print (data)
    # print(Y)
    return X,Y

# import_data()
    