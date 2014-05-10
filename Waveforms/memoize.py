import functools
def memoize(obj):
    """Decorator for memoization (caching)

    This function serves as a decorator for expensive functions whose
    return values are better being cached.  This code comes from
    <https://wiki.python.org/moin/PythonDecoratorLibrary#Memoize>.

    """
    cache = obj.cache = {}
    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer
