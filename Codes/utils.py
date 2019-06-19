from keras import backend as K


def sensitivity(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    return true_positives / (possible_positives + K.epsilon())


def specificity(y_true, y_pred):
    true_negatives = K.sum(K.round(K.clip((1 - y_true) * (1 - y_pred), 0, 1)))
    possible_negatives = K.sum(K.round(K.clip(1 - y_true, 0, 1)))
    return true_negatives / (possible_negatives + K.epsilon())


def weighted_cross_entropy(weights):
    def wce(y_true, y_pred):
        binary_crossentropy = K.binary_crossentropy(y_true, y_pred)
        weight_vector = y_true * weights["1"] + (1. - y_true) * weights["0"]
        weighted_bce = weight_vector * binary_crossentropy
        return K.mean(weighted_bce)

    return wce