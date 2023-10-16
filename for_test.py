
def open_files(path:str):
    with open(path) as open_file:
        lines = open_file[::4]
        print(lines)

        lines = open_file[::4]
        print(lines)








#locus = 'LOCUS       NODE_1_length_2558431_cov_75.1851642558431 bp   DNA linear'
#print(locus.split())


def get_last_dict_value(dict: dict):
    keys = list(dict.keys())
    return dict[keys[-1]]

dict = {
    'aaa': ['1', '2', '1'],
    'bbb': '2',
    'ccc': '3',
}

print(get_last_dict_value(dict))
