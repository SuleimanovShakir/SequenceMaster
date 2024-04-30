import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier

class RandomForestClassifierCustom(BaseEstimator):
    def __init__(self, n_estimators=10, max_depth=None, max_features=None, random_state=None):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []
        self.classes_ = None

    def fit(self, X, y):
        self.classes_ = sorted(np.unique(y))

        for i in range(0, self.n_estimators):
            np.random.seed(self.random_state + i)
            feature_idx = np.array(np.random.choice(X.shape[1], self.max_features, replace = False))
            self.feat_ids_by_tree.append(feature_idx)
            bootstrapped_idx = np.array(np.random.choice(X.shape[0], X.shape[0], replace=True))
            bootstrapped_X = np.array([X[i] for i in bootstrapped_idx])
            bootstrapped_X_with_selected_features = np.array([bootstrapped_X[:, k] for k in self.feat_ids_by_tree[i]]).T
            bootstrapped_y = np.array([y[i] for i in bootstrapped_idx])
            current_tree_model = DecisionTreeClassifier(max_depth=self.max_depth, random_state=self.random_state)
            current_tree = current_tree_model.fit(bootstrapped_X_with_selected_features, bootstrapped_y)
            self.trees.append(current_tree)
        return self

    def predict_proba(self, X):
        '''
        Function to predict probabilities of the class using the random forest.
        '''
        total_pred = []
        for i in range(len(self.trees)):
            y_pred = self.trees[i].predict_proba(np.array([X[:, k] for k in self.feat_ids_by_tree[i]]).T)
            total_pred.append(y_pred)
        return np.mean(total_pred, axis = 0)

    def predict(self, X):
        '''
        Function to predict classes using random forest.
        '''
        probas = self.predict_proba(X)
        predictions = np.argmax(probas, axis=1)
        return predictions
