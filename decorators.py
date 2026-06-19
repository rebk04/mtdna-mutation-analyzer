import time
import functools

def measure_time(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        print(f"{func.__name__} completed in {end_time - start_time:.6f} seconds.")
        return result
    return wrapper
