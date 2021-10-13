import csv

rows = csv.DictReader(open("SimScens_2V.csv"))

qsub = "qsub -A lp_h_vsc33502 ReGTs_2V.pbs "

for row in rows:
    vars = []
    for key in row:
        vars.append(key + "=" + row[key])
    print(qsub + "-v \"" + ",".join(vars) + "\"")
