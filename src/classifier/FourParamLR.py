"""
4 Parameter Logistic Regression
"""

import numpy as np
import pandas as pd
from math import log
import matplotlib.pyplot as plt
from scipy.optimize import minimize, Bounds
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from itertools import islice


class FourParamLogisticRegression:

    """ 4 Paramater Logistic Regression Model as described in:
        https://www.statforbiology.com/nonlinearregression/usefulequations#logistic_curve

    Parameters
    ----------
    penalty : {'l1', 'l2', 'elasticnet'}, default=None
        Used to specify the norm used in the penalization. The 'newton-cg',
        'sag' and 'lbfgs' solvers support only l2 penalties. 'elasticnet' is
        only supported by the 'saga' solver.

    max_iter : int, default=100
        Maximum number of iterations of the optimization algorithm.

    lr : float, default = 1e-5

    solver: str, {'SGD', 'mini-batch''}, default='SGD'
        Type of optimisation algorithm

    Attributes
    ----------

    coefs_ : 1darray of the coefficients a, b, c and d.

    """
    def __init__(self, penalty=None, lr=1e-3, max_iter=1000, solver='SGD'):
        self.penalty = penalty
        self.lr = lr
        self.max_iter = max_iter
        self.solver = solver
        self.batch_size = 512

        # Initial parameters
        self.a, self.b, self.c, self.d = 0.84177913, 0.22684613, -6.49857856, 0.00797804

    @staticmethod
    def sgd_cross_entropy(actual, predicted):
        # calculate binary cross entropy
        if actual == 1:
            return -log(predicted)
        else:
            return -log(1 - predicted)

    def four_param_sigmoid(self,z):
        """ Forward pass of LR function """
        denom_ = (1 + np.exp( - (self.b*(z[0])) - self.c))
        numer_ = (self.a - self.d)
        result = self.d + numer_/denom_
        return result

    def fit_online(self, X, y):
        """Fit the model using a data stream
        Parameters
        ----------
        X : Generator {array-like, sparse matrix} of shape (n_features)
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
        y : Generator, Target vector relative to X.

        Returns
        -------
        self : object
        """
        not_empty=True
        n= 1 if self.solver == 'SGD' else self.batch_size
        while not_empty:
            xhat = list(islice(X, n))
            ytrue = list(islice(y, n))
            if len(xhat)==0:
                not_empty=False
            else:
                delta_a, delta_b, delta_c, delta_d = 0, 0, 0, 0
                for d_point in range(len(xhat)):
                    yhat1 = self.four_param_sigmoid(xhat[d_point])
                    y_fac = (yhat1 - ytrue[d_point]) / (yhat1 * (1 - yhat1)) # nb negative log likihood for minimisation

                    delta_a += (1/self.batch_size) * y_fac * (yhat1 - self.d) / (self.a - self.d)
                    delta_b += (1/self.batch_size) * y_fac * np.array(xhat[d_point][0]) * ((yhat1 - self.d) * ((self.a - yhat1) / (self.a - self.d)))
                    delta_c += (1/self.batch_size) * y_fac * (yhat1 - self.d) * ((self.a - yhat1) / (self.a - self.d))
                    delta_d += (1/self.batch_size) * y_fac * (self.a - yhat1) / (self.a - self.d)

                self.a = min(1 - 1e-6, max(1e-6, self.a - (self.lr * delta_a)))
                self.b = self.b - (self.lr * delta_b)
                self.c = self.c - (self.lr * delta_c)
                self.d = min(1 - 1e-16, max(1e-6, self.d - (self.lr * delta_d)))

        return self

    def fit(self, X, y):
        """Fit the model according to the given training data.
        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
        y : array-like of shape (n_samples,)
            Target vector relative to X.

        Returns
        -------
        self : object
        """

        for _ in range(self.max_iter):
            if self.solver == 'SGD':
                rnd_index = np.random.randint(len(X), size=1)[0]
                xhat = X[rnd_index]
                ytrue = y[rnd_index]
                yhat = self.four_param_sigmoid(xhat)

                y_fac = (yhat - ytrue)/(yhat*(1-yhat)) #nb negative log likihood for minimisation

                delta_a = y_fac*(yhat - self.d)/(self.a - self.d)
                delta_b = y_fac*np.array(xhat[0]) * ((yhat - self.d) * ((self.a - yhat) / (self.a - self.d)))
                delta_c = y_fac*(yhat - self.d)*((self.a - yhat) / (self.a - self.d))
                delta_d = y_fac*(self.a - yhat)/(self.a - self.d)

                self.a = min(1-1e-6,max(1e-6,self.a - (self.lr * delta_a)))
                self.b = self.b - (self.lr * delta_b)
                self.c = self.c - (self.lr * delta_c)
                self.d = min(1-1e-16,max(1e-6,self.d - (self.lr * delta_d)))

            elif self.solver == 'mini-batch':
                rnd_index = np.random.randint(len(X), size=self.batch_size)
                xhat = [X[inp] for inp in rnd_index]
                ytrue = [y[inp] for inp in rnd_index]
                delta_a, delta_b, delta_c, delta_d = 0, 0, 0, 0

                for d_point in range(len(xhat)):
                    yhat1 = self.four_param_sigmoid(xhat[d_point])
                    y_fac = (yhat1 - ytrue[d_point]) / (yhat1 * (1 - yhat1)) # nb negative log likihood for minimisation

                    delta_a += (1/self.batch_size) * y_fac * (yhat1 - self.d) / (self.a - self.d)
                    delta_b += (1/self.batch_size) * y_fac * np.array(xhat[0]) * ((yhat1 - self.d) * ((self.a - yhat1) / (self.a - self.d)))
                    delta_c += (1/self.batch_size) * y_fac * (yhat1 - self.d) * ((self.a - yhat1) / (self.a - self.d))
                    delta_d += (1/self.batch_size) * y_fac * (self.a - yhat1) / (self.a - self.d)

                self.a = min(1 - 1e-6, max(1e-6, self.a - (self.lr * delta_a)))
                self.b = self.b - (self.lr * delta_b)
                self.c = self.c - (self.lr * delta_c)
                self.d = min(1 - 1e-16, max(1e-6, self.d - (self.lr * delta_d)))

        return self

    def predict(self, X_):

        self.probas = [self.four_param_sigmoid(x_) for x_ in X_]
        self.preds = [1 if x > 0.5 else 0 for x in self.probas]

        return np.array(self.preds)

    def fit_batch(self, X, y):
        """Optimise theta offline for the entire training data using an optimizer.
        Parameters
        ----------
        theta : 1D numpy array with shape (4,)
                Initial guess for the model parameters
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
        y : array-like of shape (n_samples,)
            Target vector relative to X.

        Returns
        -------
        theta : optimized parameters
        """
        theta_0 = np.array([self.a, self.b, self.c, self.d])

        def neg_log_likelihood(theta, X, y):
            m = X.shape[0]
            yhat = theta[3] + ((theta[0] - theta[3])/(1 + np.exp(-(theta[1]*X + theta[2]))))
            return -(1 / m) * np.sum(y*np.log(yhat) + (1 - y)*np.log(1 - yhat))

        def optimize_theta(theta, X, y):
            bounds = Bounds([0, -5, -10, 0], [1, 5, 20, 1])
            opt_weights = minimize(neg_log_likelihood, theta, method='L-BFGS-B',bounds=bounds, args=(X, y.flatten()))

            return opt_weights.x

        self.a, self.b, self.c, self.d = optimize_theta(theta_0, X, y)

        return self


# mdl = FourParamLogisticRegression()
# df1 = pd.read_csv('../../data/SS_only.csv', index_col=[0])
# # train_x = [[ss] for ss in df1['Cn0DbHz']]
# train_x1 = np.array(df1['Cn0DbHz'].values)
# train_y = np.array(df1['Label'].values)
# theta_0 = np.array([0.6,.1,1,0.5])
# mdl = FourParamLogisticRegression().fit_batch(train_x1, train_y)
# # mdl.a, mdl.b, mdl.c, mdl.d = results
# print(results)
# # X_train, X_test, y_train, y_test = train_test_split(train_x, train_y, test_size=0.2, random_state=42)
# # mdl = mdl.fit(X_train,y_train)
# # print(mdl.a, mdl.b, mdl.c, mdl.d)
# #
# #plot
# bins = np.linspace(3, 57, 20)
# df1['bin']=pd.cut(df1['Cn0DbHz'],bins)
# chart=pd.crosstab(df1['bin'],df1['Label'])
# chart['proportion'] = chart[1]/(chart[0]+chart[1])
# chart['midpt']=bins[1:]/2+bins[0:-1]/2
# chart['model']=mdl.predict([[ss] for ss in chart['midpt'].values])

# plt.plot(chart['midpt'].values,chart['proportion'].values,'ro')
# plt.plot(chart['midpt'].values,chart['model'].values,'b--')
# plt.show()
# # #
# print(mdl.a, mdl.b, mdl.c, mdl.d)
# # y_out = mdl.predict(X_test)
# print(y_out)
# print(roc_auc_score(y_test, y_out))