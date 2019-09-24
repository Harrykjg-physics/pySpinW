import inspect, dis

def expecting(offset=3):
    """Return how many values the caller is expecting"""
    f = inspect.currentframe()
    f = f.f_back.f_back
    c = f.f_code
    i = f.f_lasti + offset
    bytecode = c.co_code
    instruction = bytecode[i]
    if instruction == dis.opmap['UNPACK_SEQUENCE']:
        howmany = bytecode[i + 1]
        return howmany
    elif instruction == dis.opmap['POP_TOP']:
        return 0
    return 1