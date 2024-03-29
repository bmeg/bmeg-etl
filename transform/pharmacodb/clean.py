
def cleanEmptyFields(row):
    out = {}
    for k, v in row.items():
        if v != "":
            out[k] = v
    return out

def cleanDoseViability(row):
    #group the individual dose and viabilities together:
    row["dose"] = {}
    row["viability"] = {}
    for key in row:
        if key[0:4] == "rawD":
            row["dose"][key[4:]] = row[key]
        elif key[0:4] == "rawV":
            row["viability"][key[4:]] = row[key]
    #connect the doses and viabilitiess together
    dosePairs = []
    newDose = []
    newResponse = []
    for dose in row["dose"]:
        try:
            nd = float(row["dose"][dose])
            nv = float(row["viability"][dose])
        except:
            continue
        else:
            dosePairs.append((nd, nv))
    if len(dosePairs) > 1:
        dosePairs.sort(key=lambda pair: pair[0])
        pairs = list(zip(*dosePairs))
        row["dose_um"] = list(pairs[0])
        row["response"] = list(pairs[1])
        
    #find the profile measurements and make sure to grab the right ones
    responseMetrics = []
    for field in row:
        if field[0:12] == "profileInfo_":
            responseMetrics.append(field)
    recompute = []
    for field in responseMetrics:
        if len(row[field]) > 0:
            if field[-10:] == "RECOMPUTED":
                recompute.append(field)
        else:
            id = field[12:15]
            if id == "IC5":
                row["ic50"] = float(row[field])
            elif id == "EC5":
                row["ec50"] = float(row[field])
            elif id == "AAC":
                row["aac"] = float(row[field])
            elif id == "AUC":
                row["auc"] = float(row[field])
            elif id == "EIN" or id == "E_I" or id == "EMA":
                row["einf"] = float(row[field])
            elif id[:2] == "HS" or id == "SLO":
                row["hs"] = float(row[field])
    for field in recompute:
        id = field[12:15]
        if id == "IC5":
            row["ic50"] = float(row[field])
        elif id == "EC5":
            row["ec50"] = float(row[field])
        elif id == "AAC":
            row["aac"] = float(row[field])
        elif id == "AUC":
            row["auc"] = float(row[field])
        elif id == "EIN" or id == "E_I" or id == "EMA":
            row["einf"] = float(row[field])
        elif id[:2] == "HS" or id == "SLO":
         row["hs"] = float(row[field])
    
    # for drugnum in range(0, len(UNIQUEid)):
    #   aUnique = UNIQUEid[drugnum]
    #   if aUnique[0:3] == "NSC" and aUnique[3] != "-":
    #     aUnique = "NSC-" + aUnique[4:]
    #   row["drug" + str(drugnum)] = UNIQUEid[drugnum]
    
    row["id"] = row["project"] + ":pharmacodb:" + row["sampleid"] + ":" + row["treatmentid"]
    row["project_id"] = row["project"]
    row["submitter_id"] = row["experimentID"]
    if 'chemblid' in row:
        row["compounds"] = [{"id": row['chemblid']}]
    row["aliquot"] = [{"id": "pharmacodb:" + row["project"] + ":aliquot:" + row["sampleid"]}]
    return row

def finalclean(row):
    row["compounds"] = []
    for drugNum in range(0,3):
        if "chembl"+str(drugNum) in row:
            row["compounds"].append({"id": row["chembl"+str(drugNum)]})
        elif len(row["drug"+str(drugNum)]) > 1:
            if row["drug"+str(drugNum)][0:6] == "Chembl" or row["drug"+str(drugNum)][0:6] == "CHEMBL":
                row["compounds"].append({"id": "CHEMBL"+row["drug"+str(drugNum)][6:]})
            else:
                print(row["drug"+str(drugNum)] + "\n")
    return row
