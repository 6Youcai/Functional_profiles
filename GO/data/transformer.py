import re

# {ID:[], name:"", namespace:""}
def print_dict(one_dict):
    if "ID" not in one_dict:
        return 0
    for one_id in one_dict["ID"]:
        arr = [one_id, one_dict["name"], one_dict["namespace"]]
        print("\t".join(arr))

with open('go-basic.obo') as f:
    tmp_dict = {}
    for line in f:
        if re.search(r'^\[Term|Typedef\]', line): # ensure the last go id dict correctly and be the end
            if "ID" in tmp_dict: # the last few dict has not the "id" key
                print_dict(tmp_dict)
                tmp_dict = {} # clear
        elif re.search(r'^id:\sGO:\d+', line): # go id
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
    # useless here, but......
    print_dict(tmp_dict)
