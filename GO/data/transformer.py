import re

# {ID:[], name:"", namespace:""}
def print_dict(one_dict):
    for one_id in one_dict["ID"]:
        arr = [one_id, one_dict["name"], one_dict["namespace"]]
        print("\t".join(arr))
with open('go-basic.obo') as f:
    tmp_dict = {}
    for line in f:
        if re.search(r'^\[Term\]', line):
            if tmp_dict:
                print_dict(tmp_dict)
                tmp_dict = {} # clear
        elif re.search(r'^id', line):
            Id = re.search(r'^id:\s(.+)', line).group(1)
            tmp_dict["ID"] = [Id]
        elif re.search(r'^alt_id', line):
            alt_id = re.search(r'^alt_id:\s(.+)', line).group(1)
            tmp_dict["ID"].append(alt_id)
        elif re.search(r'^name:', line): # name namespace
            name = re.search(r'^name:\s(.+)', line).group(1)
            tmp_dict["name"] = name
        elif re.search(r'^namespace', line):
            namespace = re.search(r'^namespace:\s(.+)', line).group(1)
            tmp_dict["namespace"] = namespace
