
def open_files(path:str):
    with open(path) as open_file:
        lines = open_file[::4]
        print(lines)

        lines = open_file[::4]
        print(lines)


a = input()
open_files(a)