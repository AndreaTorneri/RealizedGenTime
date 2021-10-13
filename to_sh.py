import csv

rows = csv.DictReader(open("SimScens.csv"))

qsub = "qsub -A lp_h_vsc33502 ReGTs.pbs "

for row in rows:
    vars = []
    for key in row:
        vars.append(key + "=" + row[key])
    print(qsub + "-v \"" + ",".join(vars) + "\"")
