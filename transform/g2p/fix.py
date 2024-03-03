

def feature_scan(row):
    out = []
    for feat in row["features"]:
        if "start" in feat and feat["start"] is not None and feat["chromosome"] is not None:
            out.append( "chr%s:%d-%d" % ( feat["chromosome"], int(feat["start"])+1, int(feat["end"]) ) )
    row["locs"] = out
    return row


def fix_record(row):
    if "association" in row:
        row["description"] = row["association"]["description"]
        if "drug_labels" in row["association"]:
            row["compounds"] = [ { "id" : row["association"]["drug_labels"] } ]
            if 'response_type' in row["association"]:
                row["response_type"] = row['association']['response_type']
    row["submitter_id"] = row["text_object"]["id"]

    row["project_id"] = row["source"]
    g = []
    for i in row["genes"]:
        g.append( { "id" : i } )
    row["genes"] = g

    a = []
    for feat in row["features"]:
        if "start" in feat and feat["start"] is not None and feat["ref"] is not None and feat["alt"] is not None:
            name = "chr%s:%d-%d" % ( feat["chromosome"], int(feat["start"])+1, int(feat["end"]))
            for i in row["locs"]:
                o = None
                if isinstance(i, dict):
                    if i["name"] == name:
                        ref = feat["ref"].replace("-", ".")
                        alt = feat["alt"].replace("-", ".")
                        o = {
                            "reference_bases" : ref,
                            "alternate_bases" : alt,
                            "genome" : "GRCh38",
                            "chromosome": i["chromosome"],
                            "start": int(i["start"]),
                            "end" : int(i["end"]),
                            "id" : "%s:chr%s:%d:%s:%s" % ("GRCh38", feat["chromosome"], int(i["start"]), ref, alt )
                        }
                if o is not None:
                    a.append(o)
    if a:
        row["alleles"] = a

    return row


def prepChembl(row):
    if "compounds" in row:
        row["treatmentid"] = row["compounds"][0]["id"]
    return row

def fix_compound_field(row):
    if "chemblid" in row:
        row["compounds"] = [{"id": row["chemblid"]}]
    return row

def flatten(row):
    out = []
    for a in row["alleles"]:
        a["project_id"] = row["project_id"]
        out.append(a)
    return out