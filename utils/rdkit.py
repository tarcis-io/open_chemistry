def create_compound_images(compounds : list, image_size : int, directory : str):
    """
    Creates images of chemical compounds using the RDKit library.

    Parameters:
        - compounds  (list) : List of compound dictionaries with keys 'CID' and 'CanonicalSMILES'.
        - image_size (int)  : Size of images.
        - directory  (str)  : Directory where the images will be saved.

    Returns:
        - None

    Dependencies:
        - os
        - rdkit
    """

    import os
    from rdkit import Chem

    for compound in compounds:

        CID             = compound['CID']
        CanonicalSMILES = compound['CanonicalSMILES']

        image_file = f'{ CID }_RDKit.png'
        image_file = os.path.join(directory, image_file)

        molecule = Chem.MolFromSmiles(CanonicalSMILES)
        Chem.Draw.MolToFile(molecule, image_file, size = (image_size, image_size))


def plot_compound_images(title : str, directory : str, compounds : list):
    """
    Plots the images of chemical compounds created by the RDKit library.

    Parameters:
        - title     (str)  : The title of the image.
        - directory (str)  : Directory where the images are saved.
        - compounds (list) : List of compound dictionaries with keys 'CID' and 'Title'.

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

    fig, axs = matplotlib.pyplot.subplots(math.ceil(len(compounds) / 5), 5, figsize = (16, 16), constrained_layout = True)
    fig.suptitle(title, fontsize = 16)

    axs = axs.flatten()
    for ax in axs: ax.axis('off')

    for index, compound in enumerate(compounds):

        CID   = compound['CID']
        Title = compound['Title']

        image_file = f'{ CID }_RDKit.png'
        image_file = os.path.join(directory, image_file)

        ax = axs[index]
        ax.set_title(Title)
        ax.imshow(matplotlib.pyplot.imread(image_file))
