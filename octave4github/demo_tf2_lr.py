## tensorflow 2.4.1

import tensorflow as tf
import numpy as np

# 批次采样        
def create_batch(data,batch_size):
    """
    不重复采样
    """
    (m,n) = data.shape
    idx = np.random.choice(m,batch_size,replace=False)
    xs = data[idx,:]
    return xs
## end of def create_batch(...)

## A = 0.6c + 0.2 + randn*0.05 
data = np.loadtxt('./data/demo_lr_data_3.csv',delimiter=',')

## 定义初始化器
initializer = tf.ones_initializer()

## 定义变量：斜率 k 和截距 b
k = tf.Variable(initializer(shape=(1, ), dtype=tf.float64), name="k")
b = tf.Variable(initializer(shape=(1, ), dtype=tf.float64), name="b")

batch_num  = 10
batch_size = 20    ## 批次尺寸
epochs = 2000      ## 训练次数
epshow = 200

## 定义优化器
optimizer = tf.optimizers.Adam()

## 进入训练轮次
for epoch in range(epochs):
    for j in range(batch_num):
        ## 采样
        data_j = create_batch(data,batch_size)
        c_batch = data_j[:,0]
        A_batch = data_j[:,1]
        
        ## 进入模型训练
        with tf.GradientTape() as tape:
            hat_A = k * c_batch + b
            loss_fun = tf.reduce_mean(tf.math.pow(hat_A - A_batch,2))
            
        train_variables = [k,b]
        grads = tape.gradient(loss_fun,train_variables)
        optimizer.apply_gradients(zip(grads,train_variables))
        
    if epoch > 0 and (epoch+1)%epshow == 0:
        print("epoch %d，k = %f, b = %f, 损失函数 = %f" % (
            epoch+1,k.numpy(),b.numpy(),loss_fun) )
        
"""
运行结果：
epoch 200，k = 0.404950, b = 0.309575, 损失函数 = 0.003631
epoch 400，k = 0.531663, b = 0.240670, 损失函数 = 0.001865
epoch 600，k = 0.593967, b = 0.207933, 损失函数 = 0.002698
epoch 800，k = 0.599608, b = 0.206620, 损失函数 = 0.002394
epoch 1000，k = 0.599595, b = 0.203257, 损失函数 = 0.003920
epoch 1200，k = 0.595982, b = 0.206595, 损失函数 = 0.003487
epoch 1400，k = 0.600214, b = 0.201276, 损失函数 = 0.001113
epoch 1600，k = 0.598944, b = 0.209028, 损失函数 = 0.002180
epoch 1800，k = 0.601349, b = 0.204697, 损失函数 = 0.004140
epoch 2000，k = 0.599302, b = 0.201845, 损失函数 = 0.001862
"""
