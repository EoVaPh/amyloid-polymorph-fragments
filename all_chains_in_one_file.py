from pathlib import Path


input_folder = Path("extracted_chains")

output_folder = Path(r"C:\Users\User\Documents\bioinf\smtb")

output_folder.mkdir(exist_ok=True)

output_file = output_folder / "all_chains_from_cifs.txt"


with open(output_file, "w", encoding="utf8") as output:

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
