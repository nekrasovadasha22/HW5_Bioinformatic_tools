
def open_files(path:str):
    with open(path) as open_file:
        lines = open_file[::4]
        print(lines)

        lines = open_file[::4]
        print(lines)








locus = 'LOCUS       NODE_1_length_2558431_cov_75.1851642558431 bp   DNA linear'
print(locus.split())