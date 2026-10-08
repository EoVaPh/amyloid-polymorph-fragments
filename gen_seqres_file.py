from pathlib import Path


seqs_folder = Path('extracted_chains')
seqs_file = open("seqres.txt", "w", encoding="utf8")

for seq_file_path in seqs_folder.glob("*_seq.txt"):
    id = seq_file_path.stem[:4]

    seq_file = open(seq_file_path, 'r', encoding='utf8')
    seq = seq_file.read().strip()
    seq_file.close()

    seqs_file.write('>' + id + '\n')
    seqs_file.write(seq + '\n')

seqs_file.close()
