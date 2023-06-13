import os,sys

for i in range(len(sys.argv)):
    "sys.argv[%d] = %s" % (i, sys.argv[i])

if len(sys.argv) != 3:
    print('''python PyExport input output''')
    exit(0)

fin = open(sys.argv[1],"r")
fout = open(sys.argv[2], "w")

fout.write("probe" + "," + "gene" + "\n")

for line in fin:
    line = line.rstrip()
    if(line.startswith("Probe")):
        continue
    else:
        s = line.split(",")
        probe = s[0]
        gene = s[1]
        if("///" in gene):
            gene_s = gene.split(" /// ")
            for i in range(0, len(gene_s)):
                fout.write(probe + "," + gene_s[i] + "\n")
        else:
            fout.write(probe + "," + gene + "\n")
    
fin.close()
fout.close()

