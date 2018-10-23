import requests

with open("source/pfam/id_list.txt", "w") as fh:
    handle = requests.get("http://pfam.xfam.org/families?output=text")
    for line in handle.iter_lines():
        line = line.decode()
        if not line.startswith("#"):
            row = line.split("\t")
            if len(row[0]) > 1:
                print(row[0], file=fh)
