from itertools import product

def generate_library(peptideScaffold, buildingBlockList):
    
    count = 0
    charList = []

    for char in peptideScaffold:
        if char == "*":
            charList.append(buildingBlockList[count])
            count += 1
        else:
            charList.append([char])
    
    for string in product(*charList):
        yield "-".join(string)
