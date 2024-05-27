def download_compound_properties(compounds : list, properties : list) -> list:
    """
    Downloads properties of chemical compounds from PubChem.

    Parameters:
        - compounds  (list) : List of the compound dictionaries with the key 'CID'.
        - properties (list) : List of the properties to be fetched.

    Returns:
        - list : List of the compounds containing the requested properties.

    Dependencies:
        - requests
    """

    import requests

    CIDs       = ','.join(map(str, [compound['CID'] for compound in compounds]))
    properties = ','.join(properties)

    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{ CIDs }/property/{ properties }/JSON'

    response = requests.get(url)
    response = response.json()
    response = response['PropertyTable']['Properties']

    return response


def download_compound_images(compounds : list, directory : str, image_size : int):
    """
    Downloads images of chemical compounds from PubChem.

    Parameters:
        - compounds  (list) : List of the compound dictionaries with the key 'CID'.
        - directory  (str)  : Directory where the images will be saved.
        - image_size (int)  : Size of the images.

    Returns:
        - None

    Dependencies:
        - os
        - requests
    """

    import os
    import requests

    for compound in compounds:

        CID = compound['CID']
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{ CID }/PNG?image_size={ image_size }x{ image_size }'

        response = requests.get(url)
        response = response.content

        image_file = f'{ CID }_PubChem.png'
        image_file = os.path.join(directory, image_file)

        with open(image_file, 'wb') as file:

            file.write(response)


def print_compound_properties(compounds : list):
    """
    Prints the properties of the chemical compounds downloaded from PubChem.

    Parameters:
        - compounds (list) : List of the compound dictionaries.

    Returns:
        - None

    Dependencies:
        - json
    """

    import json

    for compound in compounds:

        print(json.dumps(compound, indent = 4))


def plot_compound_images(compounds : list, directory : str, title : str):
    """
    Plots the images of the chemical compounds downloaded from PubChem.

    Parameters:
        - compounds (list) : List of the compound dictionaries with the keys 'CID' and 'Title'.
        - directory (str)  : Directory where the images are saved.
        - title     (str)  : Title of the image.

    Returns:
        - None

    Dependencies:
        - math
        - matplotlib
        - os
    """

    import math
    import matplotlib.pyplot
    import os

    fig_cols = 5
    fig_rows = math.ceil(len(compounds) / fig_cols)

    fig, axs = matplotlib.pyplot.subplots(fig_rows, fig_cols, figsize = (16, 16), constrained_layout = True)
    fig.suptitle(title, fontsize = 16)

    axs = axs.flatten()
    for ax in axs: ax.axis('off')

    for index, compound in enumerate(compounds):

        CID   = compound['CID']
        Title = compound['Title']

        image_file = f'{ CID }_PubChem.png'
        image_file = os.path.join(directory, image_file)

        ax = axs[index]
        ax.set_title(Title)
        ax.imshow(matplotlib.pyplot.imread(image_file))
