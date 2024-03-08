

def date_extract(data):
    if 'Year' in data and 'Month' in data and 'Day' in data:
        return "%s-%s-%s" % (data['Year'], data['Month'], data['Day'])
    elif 'Year' in data and 'Month' in data:
        return "%s-%s" % (data['Year'], data['Month'])
    elif 'Year' in data:
        return data['Year']
    else:
        return ''


def _get_large_strings_helper(array, strings_progress):
    if isinstance(array, str) and len(array) > 18:
        strings_progress += (" "+array)

    if isinstance(array, dict):
        for key, value in array.items():
            _get_large_strings_helper(value, strings_progress)

    if isinstance(array, list):
        for value in array:
            _get_large_strings_helper(value, strings_progress)

    return strings_progress


def get_large_strings(array):
    large_strings = ""
    return _get_large_strings_helper(array, large_strings)[1:]


def update(x):
    x["pmid"] = x["MedlineCitation"]["PMID"]["#content"]

    #pull out abstract information
    if "Abstract" in x["MedlineCitation"]["Article"]:
        abstract = x["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
        full_text = get_large_strings(abstract)
        x["abstract"] = full_text


    # pull out author information
    if "AuthorList" in x["MedlineCitation"]["Article"]:
        author = x["MedlineCitation"]["Article"]["AuthorList"]["Author"]
        authors = []
        if isinstance(author, dict):
            author = [author]
        for i in author:
            if "LastName" in i and "ForeName" in i and "Initials" in i:
                authors.append( "%s, %s %s" % (i["LastName"], i["ForeName"], i["Initials"]) )
        x["authors"] = authors
        x["title"] = x["MedlineCitation"]["Article"]["ArticleTitle"]
        if isinstance(x["title"], dict) and "#content" in x["title"]:
            x["title"] = x["title"]["#content"]
        if "ArticleDate" in x["MedlineCitation"]["Article"]:
            x["date"] = date_extract(x["MedlineCitation"]["Article"]["ArticleDate"])
    return x
