import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
import make_my_data as md

X, Y = md.import_data()

(n_step, n_ver) = X.shape
n_test = int(n_step/5)

X_TRAIN = X[0:n_step-n_test, :]
Y_TRAIN = Y[0:n_step-n_test, :]
X_TEST = X[n_step-n_test:n_step, :]
Y_TEST = Y[n_step-n_test:n_step, :]

train_db = tf.data.Dataset.from_tensor_slices((X_TRAIN, Y_TRAIN)).batch(128)
test_db = tf.data.Dataset.from_tensor_slices((X_TEST, Y_TEST)).batch(128)
train_iter = iter(train_db)
sample = next(train_iter)
print('batch:', sample[0].shape, sample[1].shape)

"""
w1 = tf.Variable(tf.random.truncated_normal([n_ver, 256], stddev=0.1))
b1 = tf.Variable(tf.zeros([256]))
w2 = tf.Variable(tf.random.truncated_normal([256, 128], stddev=0.1))
b2 = tf.Variable(tf.zeros([128]))
w3 = tf.Variable(tf.random.truncated_normal([128, n_ver], stddev=0.1))
b3 = tf.Variable(tf.zeros([30]))

lr = 1e-3

for epoch in range(1000):  
    for step, (x, y) in enumerate(train_db):  # for every batch
        # x:[128, 30]
        # y:[128, 30]
        with tf.GradientTape() as tape: 
            h1 = x@w1 + tf.broadcast_to(b1, [x.shape[0], 256])
            h1 = tf.nn.relu(h1)

            h2 = h1@w2 + b2
            h2 = tf.nn.relu(h2)
            
            out = h2@w3 + b3

            # compute loss
            loss = tf.square(y - out)
            # mean: scalar
            loss = tf.reduce_mean(loss)

        # compute gradients
        grads = tape.gradient(loss, [w1, b1, w2, b2, w3, b3])
        # print(grads)
        # w1 = w1 - lr * w1_gradip3 install keras and it worked.

￼￼ 13￼￼ 5

        w1.assign_sub(lr * grads[0])
        b1.assign_sub(lr * grads[1])
        w2.assign_sub(lr * grads[2])
        b2.assign_sub(lr * grads[3])
        w3.assign_sub(lr * grads[4])
        b3.assign_sub(lr * grads[5])

        if step % 100 == 0:
            print(epoch, step, 'loss:', float(loss))
        

    # test/evaluation
    # [w1, b1, w2, b2, w3, b3]
    if epoch % 100 == 0:
        total_loss = 0
        for step, (x, y) in enumerate(test_db):
            h1 = tf.nn.relu(x@w1 + b1)
            h2 = tf.nn.relu(h1@w2 + b2)
            out = h2@w3 + b3
            total_loss += tf.reduce_mean(tf.square(y - out))
        print('epoch: %i test loss: %f' %(epoch, total_loss))
"""

model = Sequential()
model.add(Dense(256, input_shape=(30,)))
model.add(Activation('relu'))
model.add(Dense(128))
model.add(Activation('relu'))
model.add(Dense(30))

model.compile(optimizer='rmsprop', loss='mse')

model.fit(X_TRAIN, Y_TRAIN, epochs=100, batch_size=32)
score = model.evaluate(X_TEST, Y_TEST, batch_size=32)

print(model.weights)