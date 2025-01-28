from functools import wraps
import json

def get_config():
    with open('analysis_pipeline/config.json') as config_file:
        config = json.load(config_file)

    return config

debug_behavior = get_config()["debug_behavior"]
verbose = eval(debug_behavior['verbose'])
time_func = eval(debug_behavior['time_functions'])

if time_func:
    import time

def debug(func):

    @wraps(func)
    def wrapper(*args, **kwargs):
        if time_func:
            start_time = time.time()

        result = func(*args, **kwargs)

        if time_func:
            end_time = time.time()
            print(f"Function {func.__name__}: {(end_time-start_time):.4f} seconds")
        elif verbose:
            print(f"Function {func.__name__}, ", end='')

        if verbose:
            print("Return value:\n", result)
            print("---" * 30)

        return result
    
    return wrapper