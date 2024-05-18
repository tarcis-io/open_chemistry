def download_compound_properties(CIDs : list, properties : list) -> list:
    """
    Downloads properties of chemical compounds from PubChem.

    Parameters:
        - CIDs       (list) : List of the compound identifiers.
        - properties (list) : List of the properties to be fetched.

    Returns:
        - list : List of the compounds containing the requested properties.

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


def download_compound_images(CIDs : list, image_size : int, directory : str):
    """
    Downloads images of chemical compounds from PubChem.

    Parameters:
        - CIDs       (list) : List of the compound identifiers.
        - image_size (int)  : Size of the images.
        - directory  (str)  : Directory where the images will be saved.

    Returns
        - None

    Dependencies:
        - os
        - requests
    """

    import os
    import requests

    for CID in CIDs:

        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{ CID }/PNG?image_size={ image_size }x{ image_size }'

        response = requests.get(url)
        response = response.content

        image_file = f'{ CID }_PubChem.png'
        image_file = os.path.join(directory, image_file)

        with open(image_file, 'wb') as file:

            file.write(response)
