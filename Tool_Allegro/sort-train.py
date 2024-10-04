# ratio of test data : if 9:1, rate_seq = 0.1, if 8:2, rate_seq = 0.2
rate_seq = 0.1

# Step 1: Read the train.xyz file
with open("train.xyz", 'r') as file:
    lines = file.readlines()

# Step 2: Parse the file into blocks and count the number of blocks
blocks = []
index = 0

while index < len(lines):
    num_atoms = int(lines[index].strip())
    block = lines[index:index + num_atoms + 2]
    blocks.append(block)
    index += num_atoms + 2

print(f"There are {len(blocks)} blocks in the file.")

# Step 3: Split the blocks into train_data and test_data based on rate_seq
train_data = []
test_data = []

for i in range(len(blocks)):
    if (i + 1) % int(1 / rate_seq) == 0:
        test_data.append(blocks[i])
    else:
        train_data.append(blocks[i])

# Step 4: Write the blocks to train-seq.xyz file
with open("train-seq.xyz", 'w') as file:
    for block in train_data:
        for line in block:
            file.write(line)
    for block in test_data:
        for line in block:
            file.write(line)

print(f"Data has been written to train-seq.xyz")
