from pathlib import Path


input_folder = Path("extracted_chains")
output_file_path = "chains.txt"

with open(output_file_path, "w", encoding="utf8") as output:

    for file in input_folder.glob("*txt"):
        if file.name.endswith("_seq.txt"):
            continue

        sequence = ''

        with open(file, "r", encoding="utf8") as input_file:
            for line in input_file:
                sequence += line.split()[-1]

        file_name = file.stem[:4]

        print(f">{file_name}", file=output)
        print(sequence, file=output)
