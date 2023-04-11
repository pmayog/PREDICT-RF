import os

total= open('total.csv','w')
if os.path.isdir('csv/'):
    for file in os.listdir('csv/'):
        print(file)
        with open(file, 'r') as protein:
            for line in protein:
                if 'residue' not in line:
                    total.write(line)




