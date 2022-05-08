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

class Linearfit(tf.keras.Model):
    def __init__(self):
        super().__init__()
        self.dense = tf.keras.layers.Dense(
            units = 1,
            activation=None,
            kernel_initializer=tf.zeros_initializer(),
            bias_initializer=tf.zeros_initializer()
        )

    def call(self, input):
        output = self.dense(input)
        return output
## end of class Linearfit(...)

## 调入数据
## A = 0.6c + 0.2 + randn*0.05 
data = np.loadtxt('./data/demo_lr_data_3.csv',delimiter=',')

batch_num  = 10
batch_size = 20    ## 批次尺寸
epochs = 2000      ## 训练轮数
epshow = 200

## 定义模型
model = Linearfit()
optimizer = tf.optimizers.Adam()

## 进入训练轮次
for epoch in range(epochs):
    for j in range(batch_num):
        ## 采样
        data_j = create_batch(data,batch_size)
        c_batch = data_j[:,0].reshape([batch_size,1])
        A_batch = data_j[:,1].reshape([batch_size,1])
        
        ## 进入模型训练
        with tf.GradientTape() as tape:
            hat_A = model(c_batch)  ## 不再显式调用 k * c_batch + b
            loss_fun = tf.reduce_mean(tf.math.pow(hat_A - A_batch,2))

        grads = tape.gradient(loss_fun,model.variables)
        optimizer.apply_gradients(zip(grads,model.variables))
        
    if epoch > 0 and (epoch+1)%epshow == 0:
        print("epoch %d, k = %f, b = %f, 损失函数 = %f" % (
            epoch+1,model.variables[0].numpy(),model.variables[1].numpy(),
            loss_fun))
        
"""
运行结果：
epoch 200, k = 0.482823, b = 0.268453, 损失函数 = 0.003506
epoch 400, k = 0.592939, b = 0.209317, 损失函数 = 0.002622
epoch 600, k = 0.600647, b = 0.206782, 损失函数 = 0.004041
epoch 800, k = 0.602118, b = 0.208749, 损失函数 = 0.002075
epoch 1000, k = 0.594694, b = 0.206938, 损失函数 = 0.001784
epoch 1200, k = 0.599350, b = 0.207490, 损失函数 = 0.001724
epoch 1400, k = 0.599712, b = 0.205404, 损失函数 = 0.002063
epoch 1600, k = 0.600660, b = 0.204496, 损失函数 = 0.002282
epoch 1800, k = 0.604246, b = 0.205373, 损失函数 = 0.002938
epoch 2000, k = 0.601915, b = 0.209566, 损失函数 = 0.002335
"""
