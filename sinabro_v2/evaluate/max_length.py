from sinabro.core import integer

def eval_maxlen(records, mutinfo, **kwargs):
    return (len(records) - 1 >= integer(kwargs.get('maxlen', 10), 'maxlen'), None)
