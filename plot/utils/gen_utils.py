import numpy as np
import re

# =================================
# Utility to turn parameters into a command line string
# =================================
def str_to_bool(s_val):
    """If s_val is already a bool, do nothing.
    Otherwise, convert to a bool.
    """
    if isinstance(s_val, bool):
        return s_val

    # int conversion
    if s_val == 1:
        return True
    if s_val == 0:
        return False

    # string conversion
    if s_val.lower() in ['true', 't', 'yes', 'y', '1']:
        return True
    elif s_val.lower() in ['false', 'f', 'no', 'n', '0']:
        return False
    elif ' ' in s_val:
        return [str_to_bool(s_val_val)
                for s_val_val in s_val.split(' ')]
    elif ', ' in s_val:
        return [str_to_bool(s_val_val)
                for s_val_val in s_val.split(', ')]
    elif ',' in s_val:
        return [str_to_bool(s_val_val)
                for s_val_val in s_val.split(',')]
    else:
        raise ValueError(f'Cannot convert {s_val} to a boolean.')


def str_to_list(s_val):
    assert isinstance(s_val, str), \
        f'Cannot convert the non-string {s_val} to a list '\
        f'using str_to_list() (it is of type {type(s_val)}).'
    if s_val[0] == '[' and s_val[-1] == ']':
        # If it is a string that looks like a list
        s_val = s_val[1:-1].replace('\'', '').\
            replace('\"', '').\
            replace(' ', '').split(',')
    if ' ' in s_val:
        s_val = s_val.split(' ')
    elif ', ' in s_val:
        s_val = s_val.split(', ')
    elif ',' in s_val:
        s_val = s_val.split(',')


# =================================
# Utility to get bins and hist values from files
# =================================
def get_hist(filename):
    """Gets histogram bin edges, bin centers (xs), and histogram
    values (ys) from a file produced by the executables in the
    `write/` folder.
    """
    edges, xs , ys = None, None, None
    with open(filename, 'r', encoding='utf-8') as file:
        for line in file:
            if line.startswith('bin_edges'):
                edges = np.array(next(file).split(', '))
            if line.startswith('xs'):
                xs = np.array(next(file).split(', '))
            if line.startswith('ys'):
                ys = np.array(next(file).split(', '))

    return edges,\
        xs.astype(float),\
        ys.astype(float)


def mathematica_to_python(filename):
    # Open the file and read its contents
    with open(filename, 'r') as file:
        data_str = file.read()

    # Remove backticks and excess white spaces from the input string
    cleaned_data = re.sub(r'`', '', data_str)

    # Convert string to a Python list of tuples
    data = eval(cleaned_data)

    # Extract xs and ys
    xs, ys = zip(*data)

    return list(xs), list(ys)

