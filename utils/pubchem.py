def download_compound_properties(CIDs : list, properties : list) -> list:
    """
    Downloads compound properties from PubChem.

    Parameters:
        - CIDs (list)       : List of compound identifiers.
        - properties (list) : List of properties to be fetched.

    Returns:
        - list : List of compounds containing the requested properties.

    Dependencies:
        - requests
    """

    import requests

    CIDs       = ','.join(map(str, CIDs))
    properties = ','.join(properties)

    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{ CIDs }/property/{ properties }/JSON'

    response = requests.get(url)
    response = response.json()
    response = response['PropertyTable']['Properties']

    return response
