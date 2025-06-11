def __getattr__(name):
    def func(*args, **kwargs):
        return None
    return func
