from Bio import SeqIO

out = open("phylo/viral_TE_clean.fasta", "w")
for rec in SeqIO.parse("rm_correct/TE_long.fasta", "fasta"):
    # rid вида: NC_002654.1.fasta|Chlamys-1_bBra::NC_002654.1:3712-4028(+)
    parts = rec.id.split("|")
    acc = parts[0].replace(".fasta", "").replace(".", "_")
    tail = parts[1]  # "Chlamys-1_bBra::NC_002654.1:3712-4028(+)"
    name = tail.split("::")[0]      # "Chlamys-1_bBra"
    coords = tail.split("::")[1]    # "NC_002654.1:3712-4028(+)"
    pos = coords.split(":")[1]      # "3712-4028(+)"
    new_id = f"viral_{acc}_{name}_{pos}"
    out.write(f">{new_id}\n{rec.seq}\n")
out.close()
print("Unique names written to phylo/viral_TE_clean.fasta")
