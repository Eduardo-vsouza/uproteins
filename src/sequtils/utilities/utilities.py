import os


def find_coords(orf):
    pos = orf.rfind('_', 0, len(orf)-1)
    coords = orf[pos + 1:-1]
    coords = coords.split("-")
    if 'reverse' in orf:
        start = coords[1]
        end = coords[0]
    else:
        start = coords[0]
        end = coords[1]
    return start, end

def findnth(string, substring, n):
    parts = string.split(substring, n + 1)
    if len(parts) <= n + 1:
        return -1
    return len(string) - len(parts[-1]) - len(substring)

def check_dir(folder):
    if not os.path.exists(folder):
        os.system(f'mkdir {folder}')