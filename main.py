import peptic_macrocycle_generator as pmg
import pandas as pd


def main():

    macrocycle_scaffolds = pmg.load.macrocycle_scaffolds()

    print("Checkpoint: Loaded Scaffolds")

    building_blocks = pmg.load.building_blocks()

    print("Checkpoint: Loaded Building Blocks")

    combinations = pmg.gen.combinations(macrocycle_scaffolds, building_blocks)

    print("Checkpoint: Generated Combinations")

    macrocycle_library = pmg.gen.library(combinations)

    print("Checkpoint: Generated Macroycle Library")

    properties = pmg.gen.properties(macrocycle_library)

    print("Checkpoint: Generated Properties")

    table = pd.DataFrame(properties)

    print("Checkpoint: Table constructed")

    table.sort_values('Score', ascending=False)

    table.to_csv("./outs/library.csv", index=False)
    


if __name__ == "__main__":
    main()
    exit()
