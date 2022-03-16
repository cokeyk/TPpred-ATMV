import numpy as np

def normalize_square(data, dim='row'):
    if data.ndim == 1:
        # 只有一条数据
        data = np.expand_dims(data, axis=0)
    if dim == 'row':
        # 每一行求平方和
        fea_norm = np.sum(np.square(data), axis=1)
        fea_norm[fea_norm == 0] = 1e-14
        # 每个元素除以该行的根下的平方和，等价于下式
        res = np.matmul(np.diag(1 / np.sqrt(fea_norm)), data)
        # res = np.mat(np.diag(1 / np.sqrt(fea_norm))) * np.mat(data)
        return res
    else:
        # 每一行求平方和
        fea_norm = np.sum(np.square(data), axis=0)
        fea_norm[fea_norm == 0] = 1e-14
        # 每个元素除以该行的根下的平方和，等价于下式
        # res = np.mat(data) * np.mat(np.diag(1 / np.sqrt(fea_norm)))
        res = np.matmul(data, np.diag(1 / np.sqrt(fea_norm)))
        return res

