from flatten_dict import flatten, unflatten
import numpy as np


def params_flatten(params_unflattened):
    '''
    Returns a flattened (=non-nested) dict of parameters, after converting arrays to (nested) dict substructures with
    string-formatted keys.
    '''

    for key in params_unflattened.keys():
        if type(params_unflattened[key]) == np.ndarray:
            dict_num = dict(enumerate(params_unflattened[key]))
            dict_str = {}
            for subkey in dict_num.keys():
                dict_str[str(subkey)] = dict_num[subkey]
            params_unflattened[key] = dict_str

    params_flattened = flatten(params_unflattened, reducer='dot')

    return params_flattened


def params_unflatten(params_flattened):
    '''
    Performs the reverse of params_flatten.
    '''

    params_unflattened = unflatten(params_flattened, splitter='dot')

    for key in params_unflattened.keys():
        if type(params_unflattened[key]) == dict:
            flag = True
            for subkey in params_unflattened[key].keys():
                if not subkey.isdigit():
                    flag = False
            if flag:
                # reconvert keys to numeric format
                dict_num = {}
                for subkey in params_unflattened[key].keys():
                    dict_num[int(subkey)] = params_unflattened[key][subkey]
                # convert dict to np.ndarray
                params_unflattened[key] = np.array([dict_num[s] for s in dict_num.keys()])

    return params_unflattened
