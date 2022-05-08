## 演示用 TensorFlow 建立线性回归
## 本例在 TensorFlow 1.14 版下运行
import tensorflow as tf
import numpy as np

## 1. 准备数据。
## 本数据用 A = 0.6c + 0.2 + randn*0.05 模拟浓度与吸光度关系。 
data = np.loadtxt('./data/demo_lr_data_3.csv',delimiter=',')
sample_num = data.shape[0]   ## 样本数
c = data[:,0]                ## 浓度值 mol/L
A = data[:,1]                ## 吸光度值
c = c.reshape(sample_num,1)  ## 转化为列向量
A = A.reshape(sample_num,1)

lr = 0.01  ## 学习率
k0 = 0.0   ## 校正曲线斜率
b0 = 0.0   ## 校正曲线截距

## 2. 构建模型参数
## 构建浓度变量占位符
x = tf.placeholder(dtype=tf.float32, shape=(None,1))     
## 构建吸光度计算值占位符
hat_y = tf.placeholder(dtype=tf.float32, shape=(None,1))     
## 斜率（对应于权重）
k = tf.Variable(initial_value = [[1]], dtype = tf.float32)   
## 截距（对应于偏置）
b = tf.Variable(initial_value = [[1]],dtype=tf.float32)   
## 定义计算，对应于线性方程 A = k*c + b
y = tf.matmul(x,k) + b    
## 以均方差为损失函数
loss_fun = tf.reduce_sum(tf.pow((y - hat_y),2)) / sample_num
## 采用梯度下降法优化网络
optimizer = tf.train.GradientDescentOptimizer(learning_rate = lr)
## 优化操作
train_op = optimizer.minimize(loss = loss_fun)

## 3. 定义会话并初始化全部变量
sess = tf.Session()
sess.run(tf.global_variables_initializer())

## 4. 进入模型训练
for i in range(0,10000):
    ## 训练
    _,k_val,b_val,loss_val = sess.run([train_op,k,b,loss_fun],
                                      feed_dict={x:c,hat_y:A})
    ## 每经过 1000 次训练，显示当前结果。
    if i > 0 and (i+1)%1000 == 0:
        print('训练次数 %d, k = %f, b = %f, loss_fun = %f' % \
              (i+1,k_val,b_val,loss_val))
    ## 返回斜率和截距    
    k0 = k_val
    b0 = b_val
        
## 5. 训练完之后关闭会话
sess.close()

"""
训练次数 1000, k = 0.593996, b = 0.207984, loss_fun = 0.002292
训练次数 2000, k = 0.596560, b = 0.206645, loss_fun = 0.002291
训练次数 3000, k = 0.597651, b = 0.206075, loss_fun = 0.002291
训练次数 4000, k = 0.598111, b = 0.205835, loss_fun = 0.002290
训练次数 5000, k = 0.598306, b = 0.205733, loss_fun = 0.002290
训练次数 6000, k = 0.598390, b = 0.205689, loss_fun = 0.002290
训练次数 7000, k = 0.598429, b = 0.205668, loss_fun = 0.002290
训练次数 8000, k = 0.598429, b = 0.205668, loss_fun = 0.002290
训练次数 9000, k = 0.598429, b = 0.205668, loss_fun = 0.002290
训练次数 10000, k = 0.598429, b = 0.205668, loss_fun = 0.002290
"""
